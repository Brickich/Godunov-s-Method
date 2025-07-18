
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#define EXACT_SOLUTION_API  extern "C" __declspec(dllexport)

class State
{
public:
    double x, r, u, p, g, c, e, E, H;
    void Update()
    {
        c = sqrt(g * p / r);
        e = p / r / (g - 1);
        E = e + 0.5 * u * u;
        H = E + p / r;
    }
    State() :x(0), r(0), u(0), p(0), g(1.4) { Update(); };
    State(double x, double r, double u, double p, double g) :x(x), r(r), u(u), p(p), g(g)
    {
        Update();
    }
};

State left(0.0, 1.0, 0.0, 1.0, 1.4);
State right(0.0, 0.125, 0.0, 0.1, 1.4);
int N = 1000;
double left_border = -0.5;
double right_border = 0.5;
double dx = (right_border - left_border) / N;
double t_end = 0.2;
double tolerance = 1e-8;
double u_contact, p_contact;
std::vector<State> states(N);
std::vector<State> solutions(N);


void BuildGrid(std::vector<State>& states, const State& left, const State& right, double dx, double left_border) {
    for (int i = 0; i < states.size(); ++i) {
        double x = left_border + i * dx;
        states[i] = (i < states.size() / 2) ? State(x, left.r, left.u, left.p, left.g) : State(x, right.r, right.u, right.p, right.g);
    }
}

double ComputeF(const State& state, double p_current, bool left_side)
{
    double r = state.r, p = state.p, g = state.g, c = state.c;
    double F;
    if (p_current > p)
    {
        F = (p_current - p) * sqrt(2.0 / (r * ((g + 1.0) * p_current + (g - 1.0) * p)));
    }
    else
    {
        F = 2.0 / (g - 1.0) * c * (pow(p_current / p, (g - 1.0) / (2.0 * g)) - 1.0);
    }
    return F;
}

double ComputeDF(const State& state, double p_current, bool left_side, double tolerance)
{
    double DF = (ComputeF(state, p_current + tolerance, left_side) - ComputeF(state, p_current - tolerance, left_side)) / (2.0 * tolerance);
    return DF;
}

void GetContactParameters(const State& left, const State& right, const double tolerance, double& u_contact, double& p_contact) {
    double Fl, DFl;
    double Fr, DFr;
    double ul = left.u;
    double ur = right.u;
    double crit;
    double p_current = 0.5 * (left.p + right.p);
    int iter = 0;
    for (iter; iter < 100; ++iter)
    {
        if (p_current < 0) {
            throw std::runtime_error("P < 0");
        }
        Fl = ComputeF(left, p_current, true), DFl = ComputeDF(left, p_current, true, tolerance);
        Fr = ComputeF(right, p_current, false), DFr = ComputeDF(right, p_current, false, tolerance);
        p_contact = p_current - (Fl + Fr + ur - ul) / (DFl + DFr);
        crit = fabs(p_contact - p_current);
        if (crit < tolerance) break;
        p_current = p_contact;
    }
    std::cout << "Iterations " << iter << "\n";
    u_contact = 0.5 * (Fr - Fl + ul + ur);
}


State RiemanSolver(const State& left, const State& right, double tolerance, double s, double& u_contact, double& p_contact, double x) {
    double r, u, p, g, c;
    double r_new, u_new, p_new;
    double D, D_, c_;
    if (u_contact <= s)
    { // RIGHT
        r = right.r, u = right.u, p = right.p, g = right.g, c = right.c;
        if (p_contact > p)
        {
            double alpha = sqrt(r * ((g + 1.0) / 2.0 * p_contact + (g - 1.0) / 2.0 * p));
            D = u + alpha / r;
            if (s >= D)
            {
                r_new = r;
                u_new = u;
                p_new = p;
            }
            else
            {
                r_new = (r * alpha) / (alpha + r * (u - u_contact));
                u_new = u_contact;
                p_new = p_contact;
            }
        }
        else
        {
            D = u + c;
            c_ = c - (g - 1.0) / 2.0 * (u - u_contact);
            D_ = u_contact + c_;
            if (s >= D)
            {
                r_new = r;
                u_new = u;
                p_new = p;
            }
            else
            {
                if (D_ <= s)
                {
                    u_new = 2.0 / (g + 1.0) * (-c + 0.5 * (g - 1.0) * u + s);
                    double c_new = 2.0 / (g + 1.0) * (c - 0.5 * (g - 1.0) * (u - s));
                    r_new = r * pow(c_new / c, 2.0 / (g - 1));
                    p_new = p * pow(c_new / c, 2.0 * g / (g - 1));
                    // u_new = (g-1.0) / (g+1.0) * u - 2.0 / (g+1.0) * c;
                    // r_new = r * pow( u_new / c , 2.0 / (g-1.0));
                    // p_new = r_new * u_new * u_new / g;
                }
                else
                {
                    r_new = r * pow(p_contact / p, 1.0 / g);
                    u_new = u_contact;
                    p_new = p_contact;
                }
            }
        }
    }
    else { //LEFT
        r = left.r, u = left.u, p = left.p, g = left.g, c = left.c;
        if (p_contact > p) {
            double alpha = sqrt(r * ((g + 1.0) / 2.0 * p_contact + (g - 1.0) / 2.0 * p));
            D = u - alpha / r;
            if (s <= D) {
                r_new = r;
                u_new = u;
                p_new = p;
            }
            else {
                r_new = (r * alpha) / (alpha - r * (u - u_contact));
                u_new = u_contact;
                p_new = p_contact;
            }
        }
        else {
            D = u - c;
            c_ = c + (g - 1.0) / 2.0 * (u - u_contact);
            D_ = u_contact - c_;
            if (s <= D) {
                r_new = r;
                u_new = u;
                p_new = p;
            }
            else if (s <= D_) {
                u_new = 2.0 / (g + 1.0) * (c + 0.5 * (g - 1.0) * u + s);
                double c_new = 2.0 / (g + 1.0) * (c + 0.5 * (g - 1.0) * (u - s));
                r_new = r * pow(c_new / c, 2.0 / (g - 1));
                p_new = p * pow(c_new / c, 2.0 * g / (g - 1));
                // u_new = (g-1.0) / (g+1.0) * u  + 2.0 / (g+1.0) * c;
                // r_new = r * pow(u_new / c , 2.0 / (g-1.0));
                // p_new = r_new * u_new * u_new / g;
            }
            else {
                r_new = r * pow(p_contact / p, 1.0 / g);
                u_new = u_contact;
                p_new = p_contact;
            }
        }
    }

    return State(x, r_new, u_new, p_new, g);
}

EXACT_SOLUTION_API void InitialStates(double r_left, double u_left, double p_left, double g_left,
    double r_right, double u_right, double p_right, double g_right)
{
    left = State(0.0, r_left, u_left, p_left, g_left);
    right = State(0.0, r_right, u_right, p_right, g_right);
};

EXACT_SOLUTION_API void InitialGridCondition(int N_, double left_border_, double right_border_, double t_end_, double tolerance_) {
    N = N_;
    left_border = left_border_;
    right_border = right_border_;
    dx = (right_border - left_border) / N;
    t_end = t_end_;
    tolerance = tolerance_;
    states.resize(N);
    solutions.resize(N);
    BuildGrid(states, left, right, dx, left_border);
}



EXACT_SOLUTION_API void Run() {
    GetContactParameters(left, right, tolerance, u_contact, p_contact);
    for (size_t i = 0; i < states.size(); ++i) {
        double s = states[i].x / t_end;
        State solution = RiemanSolver(left, right, tolerance, s, u_contact, p_contact, states[i].x);
        solutions[i] = solution;
    }
}

EXACT_SOLUTION_API void GetSolution(double* x, double* r, double* u, double* p, double* g) {
    for (size_t i = 0; i < solutions.size(); ++i)
    {
        x[i] = solutions[i].x;
        r[i] = solutions[i].r;
        u[i] = solutions[i].u;
        p[i] = solutions[i].p;
        g[i] = solutions[i].g;
    }

}

int main(int argc, char const *argv[])
{
     
    return 0;
}
