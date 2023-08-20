import numpy as np
import math
from sympy import Symbol, solve
from sympy import Matrix, S, linsolve, symbols
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from mpl_toolkits.mplot3d import Axes3D
from sympy import Matrix, S, linsolve, symbols
import random
from octahedron_solve_system import get_x

# needed due to possible floating point error
eps = 0.00000001

"""
These codes are for visualizing a convex ideal CUBE in HYPERBOLIC 3-space in the upper half-space model
and the ball model given a set
of external dihedral angles obtained from Rivin’s theorem.

Example angles:
externalAngles = [1.9748981268459183, 2.7384076996659408, 2.2979863709366652, 1.4516735513314263, 1.5698794806677308,
2.0103008093970063, 2.322152710378157, 2.391153116133702, 2.0931040561822227, 1.9507317874044263, 2.5335253849114983,
1.7989281348636652]
"""

externalAngles = get_x()

"""
This function returns the internal dihedral angles given the external dihedral angles.
"""
def get_internal_angles(extAngles):
    intAngles = [np.pi - extAngle for extAngle in extAngles]
    return intAngles

"""
Function representing the left side of equation 6.2 from project paper.
"""
def F(f):
    return np.sin(f)*np.sin(y[5]-y[7]+f)*np.sin(y[8]-y[0]+f) - np.sin(y[0]-f)*np.sin(y[7]-f)*np.sin(y[11]-y[5]+y[7]-f)

"""
Get possible solutions to equation F(f) = 0.
Since the solving function is unstable(nonlinear), we must take this step.
The solution function will take initial value in range(-10, 10).
Although it is not 100% guaranteed that it will always give the f we want, it is reasonable to assume
it will do so most of the time.
"""
def get_possible_fValues():
    possible_f_sols = []
    for i in range(-10, 10):
        possible_f_sols.append(list(fsolve(F, i)))
    return possible_f_sols

"""
The f values obtained from the previous function may not be between
0 and pi. This function makes sure f is in that range.
The solving function itself may give multiple f values. From my tests
they seem to be almost the same values and can be ignored.
"""
def check_f_appropriate():
    possible_f_sols = get_possible_fValues()
    appropriate_f_sols = []
    for i in range(len(possible_f_sols)):
        if 0 < possible_f_sols[i][0] and possible_f_sols[i][0] < np.pi:
            appropriate_f_sols.append(possible_f_sols[i][0])
    return appropriate_f_sols

"""
Give a, b, c, d, e, f corresponding to the f value obtained from solving the equation.
a, b, c, d, e, f, must be between 0 and pi.
"""
def check_abc():
    appropriate_f_sols = check_f_appropriate()
    fGivesabc = []
    for i in range(len(appropriate_f_sols)):
        f = appropriate_f_sols[i]
        a = y[11] - y[5] + y[7] - f
        b = y[8] - y[0] + f
        c = y[7] - f
        d = y[5] - y[7] + f
        e = y[0] - f
        count = 0
        if 0 < a and a < np.pi:
            count += 1
        if 0 < b and b < np.pi:
            count += 1
        if 0 < c and c < np.pi:
            count += 1
        if 0 < d and d < np.pi:
            count += 1
        if 0 < e and e < np.pi:
            count += 1
        if count == 5:
            fGivesabc.append([a, b, c, d, e, f])
    return fGivesabc

#Get the internal dihedral angles
y =  get_internal_angles(externalAngles)
#Get only the first one because the distinct solution sets are practically the same
abc =  check_abc()[0]
a = abc[0]
b = abc[1]
c = abc[2]
d = abc[3]
e = abc[4]
f = abc[5]

#############################################################################
# VERTICES ON COMPLEX PLANE IN UPPER HALF_SPACE
#############################################################################
"""
This function rotates a vector v = [v1, v2] with an angle of k (anticlockwise is positive).
"""
def rotate_vector(v1, v2, k):
    return np.array([math.cos(k)*v1 - math.sin(k)*v2, math.sin(k)*v1 + math.cos(k)*v2])

"""
Set triangle A, B, F.
"""
def get_A():
    return np.array([0, 0])
def get_B():
    return np.array([5, 0])
#line containing A, B. This will be the x-axis. This lines is parametrized.
def f1(t):
    return np.array([t, 0])
#line containing A, F. This lines is parametrized.
def f3(t):
    A = get_A()
    vecAF = rotate_vector(1, 0, y[2])
    return np.array([A[0] + t*vecAF[0], A[1] + t*vecAF[1]])
#line containing B, F. This lines is parametrized.
def g1(t):
    B = get_B()
    vecBF = rotate_vector(-1, 0, -y[3])
    return np.array([B[0] + t * vecBF[0], B[1] + t * vecBF[1]])
#Get the coordinates of F.
def get_F():
    s = Symbol("s")
    t = Symbol("t")
    sol = list(solve([f3(s)[0]-g1(t)[0], f3(s)[1]-g1(t)[1]], s, t, set=True)[1])[0]
    s = sol[0]
    t = sol[1]
    return [f3(s), s, t]

"""
Get triangle O, F, B.
"""
#line containing F, O. This lines is parametrized.
def h3(t):
    B = get_B()
    F = get_F()[0]
    vecFB = [B[0]-F[0], B[1]-F[1]]
    vecFO = rotate_vector(vecFB[0], vecFB[1], a)
    return np.array([F[0] + t*vecFO[0], F[1] + t*vecFO[1]])
#line containing B,O. This lines is parametrized.
def h1(t):
    A = get_A()
    B = get_B()
    vecBA = [A[0]-B[0], A[1]-B[1]]
    vecBO = rotate_vector(vecBA[0], vecBA[1], -(y[3]+b))
    return np.array([B[0] + t*vecBO[0], B[1] + t*vecBO[1]])
def get_O():
    s = Symbol("s")
    t = Symbol("t")
    sol = list(solve([h1(s)[0]-h3(t)[0], h1(s)[1]-h3(t)[1]], s, t, set=True)[1])[0]
    s = sol[0]
    t = sol[1]
    return [h1(s), s, t]

"""
Get triangle O, F, D.
"""
#line containing F, D. This lines is parametrized.
def g3(t):
    F = get_F()[0]
    O = get_O()[0]
    vecFO = [O[0]-F[0], O[1]-F[1]]
    vecFD = rotate_vector(vecFO[0], vecFO[1], f)
    return np.array([F[0] + t*vecFD[0], F[1] + t*vecFD[1]])
#line containing B, D. This lines is parametrized.
def g2(t):
    A = get_A()
    B = get_B()
    vecBA = [A[0]-B[0], A[1]-B[1]]
    vecBD = rotate_vector(vecBA[0], vecBA[1], -(y[3]+b+c))
    return np.array([B[0]+t*vecBD[0], B[1]+t*vecBD[1]])
def get_D():
    s = Symbol("s")
    t = Symbol("t")
    sol = list(solve([g3(s)[0]-g2(t)[0], g3(s)[1]-g2(t)[1]], s, t, set=True)[1])[0]
    s = sol[0]
    t = sol[1]
    return [g3(s), s, t]

"""
Get h2
"""
#line containing O, D. This lines is parametrized.
def h2(t):
    O = get_O()[0]
    F = get_F()[0]
    vecOF = [F[0]-O[0], F[1]-O[1]]
    vecOD = rotate_vector(vecOF[0], vecOF[1], -(np.pi-e-f))
    return np.array([O[0] + t*vecOD[0], O[1] + t*vecOD[1]])

def get_t_for_h2_O_to_D():
    D_x = get_D()[0][0]
    t = Symbol("t")
    return solve(h2(t)[0]-D_x, t)[0]


"""
Complete the big triangle
"""
#line containing E, C. This lines is parametrized.
def f2(t):
    A = get_A()
    F = get_F()[0]
    D = get_D()[0]
    #Note that E, F, A are on the same line.
    vecFA = [A[0]-F[0], A[1]-F[1]]
    vecEC = rotate_vector(vecFA[0], vecFA[1], y[0])
    return np.array([D[0] + t*vecEC[0], D[1] + t*vecEC[1]])
def get_E():
    s = Symbol("s")
    t = Symbol("t")
    sol = list(solve([f3(s)[0]-f2(t)[0], f3(s)[1]-f2(t)[1]], s, t, set=True)[1])[0]
    s = sol[0]
    t = sol[1]
    return [f3(s), s, t]
def get_C():
    s = Symbol("s")
    t = Symbol("t")
    sol = list(solve([f1(s)[0]-f2(t)[0], f1(s)[1]-f2(t)[1]], s, t, set=True)[1])[0]
    s = sol[0]
    t = sol[1]
    return [f1(s), s, t]


"""
This function calculates the distance between two points.
"""
def distance(P, Q):
    return math.sqrt((P[0]-Q[0])**2 + (P[1]-Q[1])**2)
"""
This function returns the center of the circle that goes through the given three points.
DETAILS:
Two unknons t and s (parameters).
line1 = t*vecNorm1 + midpt1
line2 = s*vecNorm2 + midpt2
#x coordinate equal:
t*vecNorm1[0] + midpt1[0] = s*vecNorm2[0] + midpt2[0]
#y coordinates equal:
t*vecNorm1[1] + midpt1[1] = s*vecNorm2[1] + midpt2[1]
#So the system of equations to solve is:
t*vecNorm1[0] - s*vecNorm2[0] = midpt2[0] - midpt1[0]
t*vecNorm1[1] - s*vecNorm2[1] = midpt2[1] - midpt1[1]
"""
def get_Circle_center(P, Q, R):
    vecPQ = [Q[0]-P[0], Q[1]-P[1]]
    vecQR = [R[0]-Q[0], R[1]-Q[1]]
    #vector normal to vecPQ
    vecNorm1 = rotate_vector(vecPQ[0], vecPQ[1], np.pi/2)
    #vector normal to vecQR
    vecNorm2 = rotate_vector(vecQR[0], vecQR[1], np.pi/2)
    #midpoint  of PQ
    midpt1 = [(P[0] + Q[0]) / 2, (P[1] + Q[1]) / 2]
    #midpoint of QR
    midpt2 = [(Q[0] + R[0]) / 2, (Q[1] + R[1]) / 2]
    #line1 = t*vecNorm1 + midpt1
    #line2 = s*vecNorm2 + midpt2
    t, s = symbols("t, s")
    #Solve system of equations.
    A = Matrix([[vecNorm1[0], -vecNorm2[0]],
                [vecNorm1[1], -vecNorm2[1]]])
    b = Matrix([midpt2[0] - midpt1[0], midpt2[1] - midpt1[1]])
    sols = list(list(linsolve((A, b), [t, s]))[0])
    #put t = sols[0] into line1
    #Get center
    return [sols[0]*vecNorm1[0] + midpt1[0], sols[0]*vecNorm1[1] + midpt1[1]]


"""
This function will return the parametric function of the unique circle
that goes through points P, Q, R.
"""
def circle_function(t, P, Q, R):
    center = get_Circle_center(P, Q, R)
    radius = distance(P, center)
    return radius*np.cos(t)+center[0], radius*np.sin(t)+center[1]

"""
Below is the function that plots the vertices of the cube on the complex plane.
"""
def vertices_on_complex_plane():
    A = get_A()
    B = get_B()
    C = get_C()[0]
    D = get_D()[0]
    E = get_E()[0]
    F = get_F()[0]
    O = get_O()[0]

    #since circle is in parametric form, we want the input to be from 0 to 2pi.
    t_circ = np.arange(0, 2 * np.pi, 0.02)
    #Appropriate input value to show each line segment.
    t_f1 = np.arange(0, get_C()[1], 0.02)
    t_f3 = np.arange(0, get_E()[1], 0.02)
    t_g1 = np.arange(0, get_F()[2], 0.02)
    t_h1 = np.arange(0, get_O()[1], 0.02)
    t_h3 = np.arange(0, get_O()[2], 0.02)
    t_g3 = np.arange(0, get_D()[1], 0.02)
    t_g2 = np.arange(0, get_D()[2], 0.02)
    t_h2 = np.arange(0, get_t_for_h2_O_to_D(), 0.02)
    t_f2 = np.arange(get_E()[2],get_C()[2], 0.02)

    plt.figure(1)
    plt.subplot(211)
    plt.plot(A[0], A[1], 'ro')
    plt.plot(B[0], B[1], 'ro')
    plt.plot(F[0], F[1], 'ro')
    plt.plot(O[0], O[1], 'ro')
    plt.plot(D[0], D[1], 'ro')
    plt.plot(E[0], E[1], 'ro')
    plt.plot(C[0], C[1], 'ro')
    plt.plot(t_f1, [0]*len(t_f1), "k") #f1 is defined to be on the x-axis.
    plt.plot(f2(t_f2)[0], f2(t_f2)[1], "k")
    plt.plot(f3(t_f3)[0], f3(t_f3)[1], "k")
    plt.plot(g1(t_g1)[0], g1(t_g1)[1], 'c')
    plt.plot(g2(t_g2)[0], g2(t_g2)[1], 'c')
    plt.plot(g3(t_g3)[0], g3(t_g3)[1], 'c')
    plt.plot(h1(t_h1)[0], h1(t_h1)[1], "k")
    plt.plot(h2(t_h2)[0], h2(t_h2)[1], "k")
    plt.plot(h3(t_h3)[0], h3(t_h3)[1], "k")
    circle1 = circle_function(t_circ, O, A, B )
    circle2 = circle_function(t_circ, D, C, B)
    circle3 = circle_function(t_circ, F, E, O)
    plt.plot(circle1[0], circle1[1], "b")
    plt.plot(circle2[0], circle2[1], "b")
    plt.plot(circle3[0], circle3[1], "b")
    plt.gca().set_aspect('auto')
    plt.grid(True)
    #plt.legend()
    plt.show()

vertices_on_complex_plane()

#############################################################################
# UPPER HALF SPACE MODEL
#############################################################################
"""
This function will give the parametric function for a spherical triangle.
IMPORTANT DETAILS:
Given 3 points T, U, V onj the complex plane, there is a unique circle that goes through
the three points and thus a unique half-sphere that has this circle as the base. Draw 3
planes that are orthogonal to the complex plane which base lines are line TU, UV, and VT.
The spherical triangle we want is the regioun surrounded by these 3 planes.
parameters:
    *Note: P is a point on the boundary or the interior of the triangle TUV
    mesh - the number of columns in a mesh
    points_on_TU - Number of points to have on line segment TU to make a mesh
    points_on_UP - Number of points to have on line segment UP to make a mesh
"""
def spherical_triangle(T, U, V):
    points_on_TU = 50
    points_on_UP = 50
    # We are building a mesh here. tt and ss will be the mesh
    # interval is the numbers in [0, 1] that t will take
    interval = np.linspace(0, 1, points_on_TU)
    tt = []
    ss = []
    for i in range(len(interval)):
        #By the property of similar triangles, s can take values 0 <= s <= t with 0 <= t <= 1.
        #Here, we are fixing t and obtaining points on the base of the triangle similar to TUV, which
        #side is length t*TU. Then we have UP = s*UV.
        s = list(np.linspace(0, interval[i], points_on_UP))
        t = list(np.array([interval[i]]*points_on_UP))
        tt.append(t)
        ss.append(s)
    #converting to arrays
    tt = np.array([np.array(x) for x in tt])
    ss = np.array([np.array(x) for x in ss])
    #P is an arbitrary point in the interior or the boundary of triangle TUV.
    #vecTP = t*vecTU + s *vecTV, 0<=t<=1, 0<=s<=t
    #Let O(0, 0, 0) and then
    #vecOP = t*vecTU + s *vecTV + vecOT, 0<=t<=1, 0<=s<=t
    vecTU = [U[0]-T[0], U[1]-T[1]]
    vecUV = [V[0] - U[0], V[1] - U[1]]
    #vecOP = [x, y, z].
    #The x, y, values P will take
    x = tt * vecTU[0] + ss * vecUV[0] + T[0]
    y = tt * vecTU[1] + ss * vecUV[1] + T[1]
    #Now, we will get the point on the sphere which has the same x, y coordinates as P.

    #Below is for obtaining the equation of the sphere which base is the circle going through T, U, V.
    #center of the circle containing T, U, V
    circleCenter = get_Circle_center(T, U, V)
    r = distance(T, circleCenter)
    
    #This is for obtaining the point on the sphere that has the same x, y coordinates as P.
    z = (r**2 - (x - circleCenter[0])**2 - (y - circleCenter[1])**2+ eps)**(1/2)
    return [x.astype(float), y.astype(float), z.astype(float)]


"""
Funtion for cut plane (portion of vertical half plane consisting a face of the cube).
P, Q are the line segment from which the orthogonal cut plane extends.
By cut plane we meen not plotting the circle which is made by the sphere
and the plane intersecting.
"""
def orthogonal_CutPlane(P, Q, TrueIfUpperHalf):
    if (TrueIfUpperHalf == True):
        points_on_PQ = 20
        points_in_z_direction  = 20
        vecPQ = np.array([Q[0] - P[0], Q[1] - P[1]])
        t = np.linspace(0, 1, points_on_PQ)
        tt = np.array([np.array([t[i]] * points_in_z_direction) for i in range(len(t))])
        xx = vecPQ[0] * tt + P[0]
        yy = vecPQ[1] * tt + P[1]
        midpt = [(P[0] + Q[0]) / 2, (P[1] + Q[1]) / 2]
        r = distance(P, midpt)
        zz = []
        for i in range(len(xx)):
            zvalue = (r ** 2 - (xx[i][0] - midpt[0]) ** 2 - (yy[i][0] - midpt[1]) ** 2 + eps) ** (1 / 2)
            zrange = np.linspace(float(zvalue), 5, points_in_z_direction)
            zz.append(zrange)
        zz = np.array(zz)
        return [xx.astype(float), yy.astype(float), zz.astype(float)]
    else:
        points_on_PQ = 50
        points_in_z_direction = 50
        vecPQ = np.array([Q[0] - P[0], Q[1] - P[1]])
        t = np.linspace(0, 1, points_on_PQ)
        tt = np.array([np.array([t[i]] * points_in_z_direction) for i in range(len(t))])
        xx = vecPQ[0] * tt + P[0]
        yy = vecPQ[1] * tt + P[1]
        midpt = [(P[0] + Q[0]) / 2, (P[1] + Q[1]) / 2]
        r = distance(P, midpt)
        zz = []
        for i in range(len(xx)):
            zvalue = (r ** 2 - (xx[i][0] - midpt[0]) ** 2 - (yy[i][0] - midpt[1]) ** 2 + eps) ** (1 / 2)
            zrange = np.linspace(float(zvalue), 10000, points_in_z_direction)
            zz.append(zrange)
        zz = np.array(zz)
        return [xx.astype(float), yy.astype(float), zz.astype(float)]

"""
This function plots a convex ideal cube in the upper half-space.
"""
def plot_upper_half():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    A = get_A()
    B = get_B()
    C = get_C()[0]
    D = get_D()[0]
    E = get_E()[0]
    F = get_F()[0]
    O = get_O()[0]
    #since circle is in parametric form, we want the input to be from 0 to 2pi.
    t_circ = np.arange(0, 2 * np.pi, 0.02)
    #Appropriate input value to show each line segment.
    t_f1 = np.arange(0, get_C()[1], 0.02)
    t_f3 = np.arange(0, get_E()[1], 0.02)
    t_g1 = np.arange(0, get_F()[2], 0.02)
    t_h1 = np.arange(0, get_O()[1], 0.02)
    t_h3 = np.arange(0, get_O()[2], 0.02)
    t_g3 = np.arange(0, get_D()[1], 0.02)
    t_g2 = np.arange(0, get_D()[2], 0.02)
    t_h2 = np.arange(0, get_t_for_h2_O_to_D(), 0.02)
    t_f2 = np.arange(get_E()[2],get_C()[2], 0.02)
    ax.plot([A[0]], [A[1]], [0],'ro', marker='o')
    ax.plot([B[0]], [B[1]], [0],'ro', marker='o')
    ax.plot([F[0]], [F[1]], [0],'ro', marker='o')
    ax.plot([O[0]], [O[1]], [0],'ro', marker='o')
    ax.plot([D[0]], [D[1]], [0],'ro', marker='o')
    ax.plot([E[0]], [E[1]], [0],'ro', marker='o')
    ax.plot([C[0]], [C[1]], [0],'ro', marker='o')
    ax.plot(f1(t_f1)[0], [0] * len(t_f1),[0]*len(t_f1))
    ax.plot(f2(t_f2)[0], f2(t_f2)[1], [0]*len(t_f2), "k")
    ax.plot(f3(t_f3)[0], f3(t_f3)[1], [0]*len(t_f3), "k")
    ax.plot(g1(t_g1)[0], g1(t_g1)[1], [0]*len(t_g1), 'c')

    ax.plot(g2(t_g2)[0], g2(t_g2)[1], [0]*len(t_g2), 'c')
    ax.plot(g3(t_g3)[0], g3(t_g3)[1], [0]*len(t_g3), 'c')
    ax.plot(h1(t_h1)[0], h1(t_h1)[1], [0]*len(t_h1), "k")
    ax.plot(h2(t_h2)[0], h2(t_h2)[1], [0]*len(t_h2), "k")
    ax.plot(h3(t_h3)[0], h3(t_h3)[1], [0]*len(t_h3), "k")
    circle1 = circle_function(t_circ, O, A, B )
    circle2 = circle_function(t_circ, D, C, B)
    circle3 = circle_function(t_circ, F, E, O)
    ax.plot(circle1[0], circle1[1], [0] * len(t_circ), "b")
    ax.plot(circle2[0], circle2[1], [0] * len(t_circ), "b")
    ax.plot(circle3[0], circle3[1], [0] * len(t_circ), "b")
    #planes
    planeAB = orthogonal_CutPlane(A, B, True)
    planeBC = orthogonal_CutPlane(B, C, True)
    planeCD = orthogonal_CutPlane(C, D, True)
    planeDE = orthogonal_CutPlane(D, E, True)
    planeEF = orthogonal_CutPlane(E, F, True)
    planeFA = orthogonal_CutPlane(F, A, True)
    
    ax.plot_surface(planeAB[0], planeAB[1], planeAB[2], color='tab:gray', alpha=0.5, shade=False)
    ax.plot_surface(planeBC[0], planeBC[1], planeBC[2], color='tab:gray', alpha=0.5, shade=False)
    
    ax.plot_surface(planeCD[0], planeCD[1], planeCD[2], color='tab:gray', alpha=0.5, shade=False)
    ax.plot_surface(planeDE[0], planeDE[1], planeDE[2], color='tab:gray', alpha=0.5, shade=False)
    
    ax.plot_surface(planeEF[0], planeEF[1], planeEF[2], color='tab:gray', alpha=0.5, shade=False)
    ax.plot_surface(planeFA[0], planeFA[1], planeFA[2], color='tab:gray', alpha=0.5, shade=False)

    #spheres
    testing = spherical_triangle(O, B, A)
    ax.plot_surface(testing[0], testing[1], testing[2], color="b", shade=False)
    
    testing = spherical_triangle(O, F, A)
    ax.plot_surface(testing[0], testing[1], testing[2], color="b", shade=False)
    
    testing = spherical_triangle(F, E, D)
    ax.plot_surface(testing[0], testing[1], testing[2], color="g", shade=False)
    
    testing = spherical_triangle(F, O, D)
    ax.plot_surface(testing[0], testing[1], testing[2], color="g", shade=False)
    
    testing = spherical_triangle(D, C, B)
    ax.plot_surface(testing[0], testing[1], testing[2], color="r", shade=False)
    
    testing = spherical_triangle(D, O, B)
    ax.plot_surface(testing[0], testing[1], testing[2], color="r", shade=False)
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.show()

plot_upper_half()

#############################################################################
# BALL MODEL
#############################################################################
"""
This function will calcule the distance between points in 3D.
"""
def ThreeD_distance(P, Q):
    return ((P[0] - Q[0]) ** 2 + (P[1] - Q[1]) ** 2 + (P[2] - Q[2]) ** 2) ** (1 / 2)

"""
This function reflects a point across the complex plane.
Returns the coordinates of the image in a list
"""
def reflection_complex_plane(P):
    return [P[0], P[1], -P[2]]
"""
This function inverts a point with respect to the sphere with center (0, 0, 2) and radius 2.
Returns the point as a list.
"""
def inversion(P):
    center_of_inversion = np.array([0, 0, 2])
    vec_center_P = np.array(P) - center_of_inversion
    Q = vec_center_P * 4 / ((ThreeD_distance(center_of_inversion, P)) ** 2) + center_of_inversion
    return Q


"""
This function takes a point in the upper half-space model to the
ball model.
"""
def upper_half_to_ball_point(P):
    return inversion(reflection_complex_plane(P))

"""
This function takes a surface in the upper half-space model to the ball model
"""
def upper_half_to_ball_surface(mesh):
    mesh_listx = mesh[0]
    mesh_listy = mesh[1]
    mesh_listz = mesh[2]
    new_mesh_listx = []
    new_mesh_listy = []
    new_mesh_listz = []
    for i in range(len(mesh_listx)):
        listx = []
        listy = []
        listz = []
        for j in range(len(mesh_listx[i])):
            x = mesh_listx[i][j]
            y = mesh_listy[i][j]
            z = mesh_listz[i][j]
            oldpoint = [x, y, z]
            # reflect and invert
            new_point = upper_half_to_ball_point(oldpoint)
            listx.append(new_point[0])
            listy.append(new_point[1])
            listz.append(new_point[2])
        new_mesh_listx.append(np.array(listx))
        new_mesh_listy.append(np.array(listy))
        new_mesh_listz.append(np.array(listz))
    new_mesh_listx = np.array(new_mesh_listx)
    new_mesh_listy = np.array(new_mesh_listy)
    new_mesh_listz = np.array(new_mesh_listz)
    return [new_mesh_listx, new_mesh_listy, new_mesh_listz]

"""
This is a plot that brings the figure to the Poincare ball model.
"""
def ball_model():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    A = get_A()
    B = get_B()
    C = get_C()[0]
    D = get_D()[0]
    E = get_E()[0]
    F = get_F()[0]
    O = get_O()[0]
    
    # turn A, B, C into 3D coordinates since they only have x, y coordinates
    O_3D = [O[0], O[1], 0]
    A_3D = [A[0], A[1], 0]
    B_3D = [B[0], B[1], 0]
    C_3D = [C[0], C[1], 0]
    D_3D = [D[0], D[1], 0]
    E_3D = [E[0], E[1], 0]
    F_3D = [F[0], F[1], 0]
    
    # sphere of inversion
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v) + 1
    ax.plot_wireframe(x, y, z, color="k", alpha=0.3)
    
    # planes
    planeAB = upper_half_to_ball_surface(orthogonal_CutPlane(A, B, False))
    planeBC = upper_half_to_ball_surface(orthogonal_CutPlane(B, C, False))
    planeCD = upper_half_to_ball_surface(orthogonal_CutPlane(C, D, False))
    planeDE = upper_half_to_ball_surface(orthogonal_CutPlane(D, E, False))
    planeEF = upper_half_to_ball_surface(orthogonal_CutPlane(E, F, False))
    planeFA = upper_half_to_ball_surface(orthogonal_CutPlane(F, A, False))
    
    ax.plot_surface(planeAB[0], planeAB[1], planeAB[2], color='y', shade=False)
    ax.plot_surface(planeBC[0], planeBC[1], planeBC[2], color='y', shade=False)
    
    ax.plot_surface(planeCD[0], planeCD[1], planeCD[2], color='m', shade=False)
    ax.plot_surface(planeDE[0], planeDE[1], planeDE[2], color='m', shade=False)
    
    ax.plot_surface(planeEF[0], planeEF[1], planeEF[2], color='c', shade=False)
    ax.plot_surface(planeFA[0], planeFA[1], planeFA[2], color='c', shade=False)
    
    # spheres
    face1 = upper_half_to_ball_surface(spherical_triangle(O, B, A))
    ax.plot_surface(face1[0], face1[1], face1[2], color="b", shade=False)
    face2 = upper_half_to_ball_surface(spherical_triangle(O, F, A))
    ax.plot_surface(face2[0], face2[1], face2[2], color="b", shade=False)
    
    face3 = upper_half_to_ball_surface(spherical_triangle(F, E, D))
    ax.plot_surface(face3[0], face3[1], face3[2], color="g", shade=False)
    face4 = upper_half_to_ball_surface(spherical_triangle(F, O, D))
    ax.plot_surface(face4[0], face4[1], face4[2], color="g", shade=False)
    
    face5 = upper_half_to_ball_surface(spherical_triangle(D, C, B))
    ax.plot_surface(face5[0], face5[1], face5[2], color="r", shade=False)
    face6 = upper_half_to_ball_surface(spherical_triangle(D, O, B))
    ax.plot_surface(face6[0], face6[1], face6[2], color="r", shade=False)
    
    # vertices after inversion
    ball_O = upper_half_to_ball_point(O_3D)
    ball_A = upper_half_to_ball_point(A_3D)
    ball_B = upper_half_to_ball_point(B_3D)
    ball_C = upper_half_to_ball_point(C_3D)
    ball_D = upper_half_to_ball_point(D_3D)
    ball_E = upper_half_to_ball_point(E_3D)
    ball_F = upper_half_to_ball_point(F_3D)
    ax.scatter([ball_A[0]], [ball_A[1]], [ball_A[2]], s=30, c='k', marker='o')
    ax.scatter([ball_B[0]], [ball_B[1]], [ball_B[2]], s=30, c='k', marker='o')
    ax.scatter([ball_F[0]], [ball_F[1]], [ball_F[2]], s=30, c='k', marker='o')
    ax.scatter([ball_O[0]], [ball_O[1]], [ball_O[2]], s=30, c='k', marker='o')
    ax.scatter([ball_D[0]], [ball_D[1]], [ball_D[2]], s=30, c='k', marker='o')
    ax.scatter([ball_E[0]], [ball_E[1]], [ball_E[2]], s=30, c='k', marker='o')
    ax.scatter([ball_C[0]], [ball_C[1]], [ball_C[2]], s=30, c='k', marker='o')
    
    # inversion center
    ax.scatter([0], [0], [2], 'k', s=30, c='k', marker='o')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_aspect('auto')
    plt.show()

ball_model()