#include <iostream>
#include <Windows.h>
#include <fstream>
#include <math.h>


int main(int argc, char const *argv[])
{
    HINSTANCE hDll = LoadLibrary("../dll/Godunov2.dll");
    if (!hDll){
        std::cerr<<" Can't load Dll\n";
        return 1;
    }
    else{
        std::cout<<"Loaded a library\n";
    }

    typedef void (*InitialStates)(double , double , double , double ,
                                     double , double, double , double);
    InitialStates init = (InitialStates)GetProcAddress(hDll , "InitialStates");

    if (!init){
        std::cerr<<"Function not found\n";
        return 1;
    }

    typedef void (*InitialGridCondition) (int N_, double left_border_, double right_border_, double t_end_, double CFL_, double tolerance_);

    InitialGridCondition initGrid = (InitialGridCondition)GetProcAddress(hDll , "InitialGridCondition");

    if(!initGrid){
        std::cerr<<"Function not found\n";
        return 1;
    }

    typedef void (*Run) ();
    Run run = (Run)GetProcAddress(hDll, "Run");
    if(!run){
        std::cerr<<"Function not found\n";
        return 1;
    }
    
    typedef void(*Solution) (double *x, double *r, double *u, double *p, double *g);
    Solution solution = (Solution)GetProcAddress(hDll , "GetSolution");
    if(!solution){
        std::cerr<<"Function not found\n";
        return 1;
    }

    init(1.0 , 0.0 , 1.0 , 1.4 , 0.125 , 0.0 , 0.1 , 1.4);
    initGrid(100 , -0.5 , 0.5 , 0.2 , 1.0 , 1e-12);
    run();
    double *xs, *rs , *us , *ps , *gs;
    // solution(xs , rs , us , ps, gs);
    std::ofstream result("result.dat");
    for(int i=0 ; i<100 ; ++i){
        result<<xs[i] << "\t" <<rs[i] << "\t" <<us[i] << "\t"<<ps[i] << "\n";   
    }
    return 0;
}
