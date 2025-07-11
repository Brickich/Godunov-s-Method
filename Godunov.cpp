#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdexcept>


#define GODUNOV_API  extern "C" __declspec(dllexport)



class State
{
public:
    double x , r, u, p, g, c, e, E, H;
    void Update()
    {
        c = sqrt(g * p / r);
        e = p / r / (g - 1);
        E = e + 0.5 * u * u;
        H = E + p / r;
    }
    State() :x(0), r(0), u(0), p(0), g(1.4) {Update();};
    State(double x, double r, double u, double p, double g) :x(x),  r(r), u(u), p(p), g(g)
    {
        Update();
    }
};

class U
{
public:
    double x, r, ru, rE, g;
    U(): x(0), r(0), ru(0), rE(0), g(1.4){};
    U(const State &state) 
    {
       x = state.x, r = state.r, ru = state.r * state.u, rE = state.r * state.E, g = state.g;
    }
};

class F
{
public:
    double x, ru, ruu_p, ruH, g;
    F():x(0), ru(0), ruu_p(0), ruH(0), g(1.4) {};
    F(const State &state)
    {
        x = state.x , ru = state.r * state.u, ruu_p = state.r * state.u * state.u + state.p, ruH = state.r * state.u * state.H, g = state.g;
    }
};

State UtoState(U& u){
    State state;
    state.x = u.x;
    state.r = u.r;
    state.u = u.ru / u.r ;
    state.p = (u.g -1) * (u.rE - 0.5 * state.r * state.u * state.u) ;
    state.g = u.g;
    state.Update();
    return state;
}

void BuildGrid(std::vector<State>& states , const State& left , const State& right, double dx, double left_border){
    for(int i=0 ; i<states.size() ; ++i) {
        double x = left_border + i*dx;

        states[i] = (i < states.size() /2) ? State(x , left.r, left.u , left.p , left.g) : State(x , right.r, right.u , right.p , right.g);
    }
}

double ComputeTimeStep(std::vector<State>& states , double CFL , double dx){
    double u_max = 0;
    for(const auto& state : states){
        double speed = fabs(state.u) + state.c;
        u_max = (u_max > speed) ? u_max : speed;
    }
    return CFL * dx / u_max;
}

void BoundaryCondition(std::vector<State>& states){ 
    double left_x = states[0].x;
    double right_x = states[states.size()-1].x;
    
    states[0] = states[1];
    states[0].x = left_x;

    states[states.size()-1] = states[states.size()-2];
    states[states.size()-1].x = right_x;
    
}

double ComputeF(const State &state, double p_current, bool left_side)
{
    double r = state.r, p = state.p, g = state.g, c = state.c;
    double F;
    if (p_current > p)
    {
        F = (p_current - p) * sqrt(2.0 / (r * ((g + 1.0) * p_current + (g - 1.0) * p)));
    }
    else
    {
        F = 2.0 / (g - 1.0) * c * (pow(p_current / p, (g - 1.0) / (2.0 * g)) -1.0);
    }
    return F ;
}

double ComputeDF(const State &state, double p_current, bool left_side, double tolerance)
{
    double DF = (ComputeF(state, p_current + tolerance, left_side) - ComputeF(state, p_current - tolerance, left_side)) / (2.0 * tolerance);
    return DF;
}


void GetContactParameters(const State& left, const State& right , const double tolerance , double& u_contact , double& p_contact ){
    double Fl, DFl;
    double Fr, DFr;
    double ul = left.u;
    double ur = right.u;
    double crit;
    double p_current = 0.5 * (left.p + right.p);
    for (int iter =0 ; iter < 100; ++iter)
    {
        if (p_current < 0){
            throw std::runtime_error("P < 0");
        }   
        Fl = ComputeF(left, p_current, true), DFl = ComputeDF(left, p_current, true, tolerance);
        Fr = ComputeF(right, p_current, false), DFr = ComputeDF(right, p_current, false, tolerance);
        p_contact = p_current - (Fl + Fr + ur - ul) / (DFl + DFr);
        crit = fabs(p_contact - p_current);
        if (crit < tolerance) break;
        p_current = p_contact;
    }
    u_contact = 0.5 * (Fr - Fl + ul + ur);
}

State RiemanSolver(const State& left, const State& right, double tolerance){
    double r, u , p , g , c;
    double r_new , u_new , p_new;
    double u_contact , p_contact;
    double x = 0.5 * (left.x + right.x);
    if ( fabs(left.r - right.r) < tolerance && fabs(left.p - right.p) < tolerance)
    {
        return State(x, 0.5 * (left.r + right.r), 0.5 * (left.u + right.u), 0.5 * (left.p + right.p), 0.5 * (left.g + right.g));
    }

    GetContactParameters(left, right, tolerance, u_contact, p_contact);
    double D , D_ , c_;
    if (u_contact <=0){//RIGHT
        r = right.r , u = right.u, p=right.p , g = right.g , c=right.c;
        if (p_contact > p) {
            double alpha = sqrt(r * ( (g+1.0) / 2.0 * p_contact + (g-1.0) /2.0 *p));
            D = u + alpha/ r;
            if (0 >= D) {
                r_new = r;
                u_new = u;
                p_new = p;
            }
            else{
                r_new = (r* alpha) / (alpha + r* (u - u_contact));
				u_new = u_contact;
				p_new = p_contact;
            }
        }
        else{
            D = u + c;
            c_ = c - (g-1.0) / 2.0 * (u - u_contact);
            D_ = u_contact + c_;
            if ( 0 >= D){
                r_new = r;
				u_new =u;
				p_new =p;
            }
            else if (D_ <=0){
                u_new = (g-1.0) / (g+1.0) * u - 2.0 / (g+1.0) * c;
                r_new = r * pow( u_new / c , 2.0 / (g-1.0));
                p_new = r_new * u_new * u_new / g;
            }
            else{
                r_new = r * pow(p_contact / p , 1.0 / g);
                u_new = u_contact;
                p_new = p_contact;
            }
        }
    }
    else{ //LEFT
        r = left.r , u = left.u, p=left.p , g = left.g , c=left.c;
        if (p_contact > p) {
            double alpha = sqrt(r * ( (g+1.0) / 2.0 * p_contact + (g-1.0) /2.0 *p));
            D = u - alpha/ r;
            if (0 <= D) {
                r_new = r;
                u_new = u;
                p_new = p;
            }
            else{
                r_new = (r* alpha) / (alpha - r* (u - u_contact));
				u_new = u_contact;
				p_new = p_contact;
            }
        }
        else{
            D = u - c;
            c_ = c + (g-1.0) / 2.0 * (u - u_contact);
            D_ = u_contact - c_;
            if ( 0 <= D){
                r_new = r;
				u_new =u;
				p_new =p;
            }
            else if (0 <= D_){
                u_new = (g-1.0) / (g+1.0) * u  + 2.0 / (g+1.0) * c;
                r_new = r * pow(u_new / c , 2.0 / (g-1.0));
                p_new = r_new * u_new * u_new / g;
            }
            else{
                r_new = r * pow(p_contact / p , 1.0 / g);
                u_new = u_contact;
                p_new = p_contact;
            }
        }
    }

    return State(x , r_new , u_new , p_new , g);
}


void GodunovStep(std::vector<State>& states, double tolerance , double dt, double dx){
    std::vector<F> fluxes(states.size()-1);
    std::vector<State> new_states = states;
    for (size_t i=0 ; i<states.size()-1 ; ++i){
        const State RiemanSolution = RiemanSolver(states[i] , states[i+1] , tolerance);
        fluxes[i] = F(RiemanSolution);
    }

    for(size_t i=1 ; i< states.size() -1; ++i){
        U u_old = U(states[i]);
        u_old.r -=dt/dx * (fluxes[i].ru - fluxes[i-1].ru);
        u_old.ru -=dt/dx * (fluxes[i].ruu_p - fluxes[i-1].ruu_p);
        u_old.rE -= dt/dx * (fluxes[i].ruH - fluxes[i-1].ruH);
        new_states[i] = UtoState(u_old);
        new_states[i].Update();
    }
    states = new_states;
}

State left_state(0.0 , 1.0 , 0.0 , 1.0 , 1.4);
State right_state(0.0 , 0.125 , 0.0 , 0.1, 1.4);
int N = 1000;
std::vector<State> states(N);
double left_border = 0.0,
        right_border = 1.0,
        dx = (right_border - left_border) / N, 
        t_end = 0.2, 
        CFL = 0.5, 
        tolerance = 1e-8;

GODUNOV_API void  InitialStates(double r_left, double u_left, double p_left, double g_left,
                                         double r_right, double u_right, double p_right, double g_right)
{
    left_state = State(0.0 , r_left, u_left, p_left , g_left);
    right_state = State(0.0, r_right, u_right, p_right, g_right);
};

GODUNOV_API void  InitialGridCondition(int N_, double left_border_, double right_border_, double t_end_, double CFL_, double tolerance_)
{
    N = N_;
    left_border = left_border_;
    right_border = right_border_;
    dx = (right_border - left_border) / N;
    t_end = t_end_;
    CFL = CFL_;
    tolerance = tolerance_;
    states.resize(N);
    BuildGrid(states , left_state , right_state , dx , left_border);
}

GODUNOV_API void  Run(){
    double t =0;
    double dt;
    while (t < t_end)
    {
        dt = ComputeTimeStep(states, CFL , dx);
        if (t + dt > t_end) dt = t_end -t;
        GodunovStep(states , tolerance , dt , dx);

        BoundaryCondition(states);
        t+=dt;
    }
}

GODUNOV_API void  GetSolution(double* x , double* r, double* u , double* p , double* g){
    for(int i=0 ; i<states.size() ; ++i){
        x[i] = states[i].x;
        r[i] = states[i].r;
        u[i] = states[i].u;
        p[i] = states[i].p;
        g[i] = states[i].g;
    }
}

