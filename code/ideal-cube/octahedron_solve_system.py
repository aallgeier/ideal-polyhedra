import numpy as np
from sympy import symbols
from sympy.solvers.solveset import linsolve

def get_x():
        
    found = False
    while not found:
        x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11 = symbols('x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11')

        sol = linsolve([x0 - np.random.rand()*np.pi,
                x1 - np.random.rand()*np.pi,
                x2 - np.random.rand()*np.pi,
                x3 - np.random.rand()*np.pi,
                x0+ x1+ x4 - 2*np.pi,
                x0+ x2+ x5 - 2*np.pi,
                x1+ x3+ x8 - 2*np.pi,
                x10+ x2+ x3 - 2*np.pi,
                x4+ x6+ x7 - 2*np.pi,
                x5+ x6+ x9 - 2*np.pi,
                x11+ x7+ x8 - 2*np.pi,
                x10+ x11+ x9 - 2*np.pi,
                x11 - np.random.rand()*np.pi], (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11))

        if len(sol.args) > 0:
            x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11 = sol.args[0]


            count = 0
            if x0+ x3+ x4+ x8 > 2*np.pi:
                count += 1
            if x0+ x10+ x3+ x5 > 2*np.pi:
                count += 1
            if x0+ x1+ x6+ x7 > 2*np.pi:
                count += 1
            if x0+ x2+ x6+ x9 > 2*np.pi:
                count += 1
            if x0+ x11+ x3+ x6 > 2*np.pi:
                count += 1
            if x1+ x2+ x4+ x5 > 2*np.pi:
                count += 1
            if x1+ x2+ x7+ x9 > 2*np.pi:
                count += 1
            if x1+ x11+ x3+ x7 > 2*np.pi:
                count += 1
            if x1+ x10+ x2+ x8 > 2*np.pi:
                count += 1
            if x11+ x2+ x3+ x9 > 2*np.pi:
                count += 1
            if x0+ x2+ x4+ x7+ x9 > 2*np.pi:
                count += 1
            if x0+ x11+ x3+ x4+ x7 > 2*np.pi:
                count += 1
            if x0+ x10+ x2+ x4+ x8 > 2*np.pi:
                count += 1
            if x0+ x1+ x5+ x7+ x9 > 2*np.pi:
                count += 1
            if x0+ x11+ x3+ x5+ x9 > 2*np.pi:
                count += 1
            if x0+ x1+ x10+ x5+ x8 > 2*np.pi:
                count += 1
            if x0+ x3+ x6+ x7+ x8 > 2*np.pi:
                count += 1
            if x0+ x10+ x3+ x6+ x9 > 2*np.pi:
                count += 1
            if x0+ x1+ x11+ x6+ x8 > 2*np.pi:
                count += 1
            if x0+ x10+ x11+ x2+ x6 > 2*np.pi:
                count += 1
            if x1+ x10+ x3+ x4+ x5 > 2*np.pi:
                count += 1
            if x1+ x2+ x4+ x6+ x9 > 2*np.pi:
                count += 1
            if x1+ x11+ x3+ x4+ x6 > 2*np.pi:
                count += 1
            if x1+ x2+ x5+ x6+ x7 > 2*np.pi:
                count += 1
            if x1+ x10+ x3+ x7+ x9 > 2*np.pi:
                count += 1
            if x1+ x10+ x11+ x2+ x7 > 2*np.pi:
                count += 1
            if x1+ x11+ x2+ x8+ x9 > 2*np.pi:
                count += 1
            if x2+ x3+ x4+ x5+ x8 > 2*np.pi:
                count += 1
            if x11+ x2+ x3+ x5+ x6 > 2*np.pi:
                count += 1
            if x2+ x3+ x7+ x8+ x9 > 2*np.pi:
                count += 1
            if x0+ x10+ x3+ x4+ x7+ x9 > 2*np.pi:
                count += 1
            if x0+ x10+ x11+ x2+ x4+ x7 > 2*np.pi:
                count += 1
            if x0+ x11+ x2+ x4+ x8+ x9 > 2*np.pi:
                count += 1
            if x0+ x3+ x5+ x7+ x8+ x9 > 2*np.pi:
                count += 1
            if x0+ x1+ x11+ x5+ x8+ x9 > 2*np.pi:
                count += 1
            if x0+ x1+ x10+ x11+ x5+ x7 > 2*np.pi:
                count += 1
            if x0+ x10+ x2+ x6+ x7+ x8 > 2*np.pi:
                count += 1
            if x0+ x1+ x10+ x6+ x8+ x9 > 2*np.pi:
                count += 1
            if x1+ x11+ x3+ x4+ x5+ x9 > 2*np.pi:
                count += 1
            if x1+ x10+ x3+ x4+ x6+ x9 > 2*np.pi:
                count += 1
            if x1+ x10+ x11+ x2+ x4+ x6 > 2*np.pi:
                count += 1
            if x1+ x10+ x3+ x5+ x6+ x7 > 2*np.pi:
                count += 1
            if x1+ x11+ x2+ x5+ x6+ x8 > 2*np.pi:
                count += 1
            if x11+ x2+ x3+ x4+ x5+ x7 > 2*np.pi:
                count += 1
            if x2+ x3+ x5+ x6+ x7+ x8 > 2*np.pi:
                count += 1
            if x2+ x3+ x4+ x6+ x8+ x9 > 2*np.pi:
                count += 1
            if x4+ x5+ x7+ x9 > 2*np.pi:
                count += 1
            if x10+ x4+ x5+ x8 > 2*np.pi:
                count += 1
            if x11+ x4+ x6+ x8 > 2*np.pi:
                count += 1
            if x10+ x11+ x5+ x6 > 2*np.pi:
                count += 1
            if x10+ x11+ x4+ x5+ x7 > 2*np.pi:
                count += 1
            if x10+ x4+ x6+ x8+ x9 > 2*np.pi:
                count += 1
            if x11+ x4+ x5+ x8+ x9 > 2*np.pi:
                count += 1
            if x10+ x5+ x6+ x7+ x8 > 2*np.pi:
                count += 1
            if x10+ x7+ x8+ x9 > 2*np.pi:
                count += 1
            if 0 < x0 and x0 < np.pi:
                count += 1
            if 0 < x1 and x1 < np.pi:
                count += 1
            if 0 < x2 and x2 < np.pi:
                count += 1
            if 0 < x3 and x3 < np.pi:
                count += 1
            if 0 < x4 and x4 < np.pi:
                count += 1
            if 0 < x5 and x5 < np.pi:
                count += 1
            if 0 < x6 and x6 < np.pi:
                count += 1
            if 0 < x7 and x7 < np.pi:
                count += 1
            if 0 < x8 and x8 < np.pi:
                count += 1
            if 0 < x9 and x9 < np.pi:
                count += 1
            if 0 < x10 and x10 < np.pi:
                count += 1
            if 0 < x11 and x11 < np.pi:
                count += 1
            
            if count == 67:
                found = True

    return np.array([x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11]).astype(np.float64)

        