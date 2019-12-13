Numerical simulation of the Schramm-Loewner Evolution (SLE)


This is a Python implementation of the numerical method presented in "An asymptotic radius of convergence for the Loewner equation and simulation of SLE traces via splitting" by James Foster, Terry Lyons and Vlad Margarint.

The script sle.py generates chordal SLE traces using the Ninomiya-Victoir scheme with a variable step size.

The code requries the libraries NumPy, Matplotlib, cmath and math.

Three example plots of SLE_k on [0,1] are included (for k = 8/3, k = 4 and k = 6).