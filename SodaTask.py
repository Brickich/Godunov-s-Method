from ctypes import c_double, CDLL, c_int
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os


os.makedirs("output" , exist_ok=True)

Godunov = CDLL("./dll/Godunov.dll")
ExactSolution = CDLL("./dll/ExactSolution.dll")


Godunov.InitialStates.argtypes= [
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
]
Godunov.InitialGridCondition.argtypes = [
    c_int,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
]
Godunov.Run.argtypes = []
Godunov.GetSolution.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
]

Godunov.InitialStates.restype = None
Godunov.InitialGridCondition.restype = None
Godunov.Run.restype = None
Godunov.GetSolution.restype = None

ExactSolution.InitialStates.argtypes= [
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
]
ExactSolution.InitialGridCondition.argtypes = [
    c_int,
    c_double,
    c_double,
    c_double,
    c_double,
]
ExactSolution.Run.argtypes = []
ExactSolution.GetSolution.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
]

ExactSolution.InitialStates.restype = None
ExactSolution.InitialGridCondition.restype = None
ExactSolution.Run.restype = None
ExactSolution.GetSolution.restype = None

def RunGodunov(parameters: dict, N , left_border, right_border , t_end , CFL, tolerance):
    r_left, u_left, p_left, g_left = parameters["left"]
    r_right, u_right, p_right, g_right = parameters["right"]
    x = np.zeros(N)
    r = np.zeros(N)
    u = np.zeros(N)
    p = np.zeros(N)
    g = np.zeros(N)

    Godunov.InitialStates(
        r_left, u_left, p_left, g_left, r_right, u_right, p_right, g_right
    )
    Godunov.InitialGridCondition(N, left_border, right_border, t_end, CFL, tolerance)
    Godunov.Run()
    Godunov.GetSolution(x, r, u, p, g)

    return x, r, u, p, g

def RunExactSolution(parameters: dict , N, left_border, right_border , t_end , tolerance):
    r_left, u_left, p_left, g_left = parameters["left"]
    r_right, u_right, p_right, g_right = parameters["right"]
    x = np.zeros(N)
    r = np.zeros(N)
    u = np.zeros(N)
    p = np.zeros(N)
    g = np.zeros(N)
    ExactSolution.InitialStates( r_left, u_left, p_left, g_left, r_right, u_right, p_right, g_right)
    ExactSolution.InitialGridCondition(N , left_border , right_border , t_end, tolerance)
    ExactSolution.Run()
    ExactSolution.GetSolution(x , r, u, p ,g)
    return x , r, u, p , g

Parameters = {
    "cvr-uv": {
        "left": [1.0, 0.0, 1.0, 1.4],
        "right": [0.125, 0.0, 0.1, 1.4],
    },
    "uv-cvr": {
        "left": [0.125, 0.0, 0.1, 1.4],
        "right": [1.0, 0.0, 1.0, 1.4],
    },
    "uv-uv":{
        "left": [1.0, 0.25, 1.0, 1.4],
        "right": [1.0, -0.25, 1.0, 1.4],
    },
    "cvr-cvr" : {
        "left": [1.0, -0.25, 1.0, 1.4],
        "right": [1.0, 0.25, 1.0, 1.4],
    },
    "vacuum" : {
        "left":[1.0, -10.0, 1.0, 1.4],
        "right": [1.0, 10.0, 1.0, 1.4],
    },
}

def PrintData(filepath:str , data:list):
    file = open(filepath, mode="w+")
    file.write(F"x \t r \t u \t p\n")
    for i in range(len(data[0])):
        file.write(f"{data[0][i]} \t {data[1][i]} \t {data[2][i]} \t {data[3][i]}\n")


def RunDiscontinuitySimulation(type:str = "cvr-uv", N = 1000 , 
                               left_border = -0.5 , right_border = 0.5 , 
                               t_end =0.2 , CFL =0.5, tolerance = 1e-10):
    x_Godunov, r_Godunov, u_Godunov, p_Godunov, g_Godunov = RunGodunov(Parameters[type],N=N,
                                                                       left_border=left_border,
                                                                       right_border=right_border,
                                                                       t_end=t_end, 
                                                                       CFL=CFL, tolerance=tolerance)
    x_Exact, r_Exact , u_Exact , p_Exact , g_Exact = RunExactSolution(Parameters[type], N=N,
                                                                       left_border=left_border,
                                                                       right_border=right_border,
                                                                       t_end=t_end, tolerance=tolerance)

    figure, axes = plt.subplots(nrows=1, ncols=3, squeeze=False)
    column_titles = [r"$\rho$", "u", "p"]
    figure.set_figwidth(15)
    axes[0, 0].plot(x_Godunov, r_Godunov, "r-", linewidth=1.0)
    axes[0, 0].plot(x_Exact, r_Exact, "b:", alpha=0.7, linewidth=2.0)
    axes[0, 0].set_title(column_titles[0])

    axes[0, 1].plot(x_Godunov, u_Godunov, "r-", linewidth=1.0)
    axes[0, 1].plot(x_Exact, u_Exact, "b:", alpha=0.7, linewidth=2.0)
    axes[0, 1].set_title(column_titles[1])

    axes[0, 2].plot(x_Godunov, p_Godunov, "r-", linewidth=1.0)
    axes[0, 2].plot(x_Exact, p_Exact, "b:", alpha=0.7, linewidth=2.0)
    axes[0, 2].set_title(column_titles[2])
    plt.savefig("Soda.png" , dpi = 300 , bbox_inches = "tight")

    PrintData(f"./output/Godunov_{type}.dat", [x_Godunov , r_Godunov , u_Godunov , p_Godunov])
    PrintData(f"./output/Exact_{type}.dat", [x_Exact , r_Exact , u_Exact, p_Exact])
    Image.open("Soda.png").show()
    # plt.show()


RunDiscontinuitySimulation("vacuum" , N=100 , 
                           left_border=-0.5 , right_border=0.5, 
                           t_end=0.15, CFL=1.0 , tolerance=1e-15)


