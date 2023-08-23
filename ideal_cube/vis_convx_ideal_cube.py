import numpy as np
import math
from sympy import Symbol, solve
from sympy import Matrix, S, linsolve, symbols
from sympy import Matrix, S, linsolve, symbols
from scipy.optimize import fsolve
import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as pyo
from ideal_cube.octahedron_solve_system import get_x
import pandas as pd
# For floating point error
eps = 0.00000001

"""
These codes are for visualizing a convex ideal CUBE in HYPERBOLIC 3-space in the upper half-space model
and the ball model given a set of external dihedral angles obtained from Rivin’s theorem.
"""

print()
print("starting visualizations...")

# Obtain random set of Rivin's angles
#externalAngles = get_x()
externalAngles = [1.91926884, 2.2427268,  2.15274097, 2.1508076,  2.12118966, 2.21117549,
                    2.02383795, 2.1381577,  1.8896509,  2.04817187, 1.97963673, 2.25537671]


# externalAngles = [3.0798724,  0.9052616,  1.24435147, 2.57417689, 2.2980513,  1.95896144,
#  1.88536715, 2.09976686, 2.80374682, 2.43885672, 2.46465695, 1.37967163]

# looks good
#[1.59390645 1.78715237 2.71021501 2.02250796 2.90212648 1.97906385
# 1.89453637 1.48652245 2.47352498 2.40958509 1.55046234 2.32313788]

# looks very good 
# [1.91926884 2.2427268  2.15274097 2.1508076  2.12118966 2.21117549
#  2.02383795 2.1381577  1.8896509  2.04817187 1.97963673 2.25537671]

# doesn't look good
# [1.90270151 2.24812436 2.42396069 1.03102482 2.13235943 1.9565231
#  2.82496023 1.32586565 3.00403613 1.50170198 2.8281998  1.95328353]

# doesn't look good
# [2.50933417 1.46398714 1.52849315 2.11281127 2.30986399 2.24535798
#  1.9517113  2.02161001 2.70638689 2.08611602 2.64188088 1.5551884 ]

print("Found external dihedral angles (radians):")
print(externalAngles)

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
Although it is not 100% guaranteed that it will always give the f we want, it seems reasonable to assume
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

# Get the internal dihedral angles
y =  get_internal_angles(externalAngles)
# Get only the first one because the distinct solution sets are practically the same
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

# line containing A, B. This will be the x-axis. This lines is parametrized.
def f1(t):
    return np.array([t, 0], dtype=object)

# line containing A, F. This lines is parametrized.
def f3(t):
    A = get_A()
    vecAF = rotate_vector(1, 0, y[2])
    return np.array([A[0] + t*vecAF[0], A[1] + t*vecAF[1]])

# line containing B, F. This lines is parametrized.
def g1(t):
    B = get_B()
    vecBF = rotate_vector(-1, 0, -y[3])
    return np.array([B[0] + t * vecBF[0], B[1] + t * vecBF[1]])

# Get the coordinates of F.
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
# line containing F, O. This lines is parametrized.
def h3(t):
    B = get_B()
    F = get_F()[0]
    vecFB = [B[0]-F[0], B[1]-F[1]]
    vecFO = rotate_vector(vecFB[0], vecFB[1], a)
    return np.array([F[0] + t*vecFO[0], F[1] + t*vecFO[1]])

# line containing B,O. This lines is parametrized.
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
# line containing F, D. This lines is parametrized.
def g3(t):
    F = get_F()[0]
    O = get_O()[0]
    vecFO = [O[0]-F[0], O[1]-F[1]]
    vecFD = rotate_vector(vecFO[0], vecFO[1], f)
    return np.array([F[0] + t*vecFD[0], F[1] + t*vecFD[1]])

# line containing B, D. This lines is parametrized.
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
# line containing O, D. This lines is parametrized.
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
# line containing E, C. This lines is parametrized.
def f2(t):
    A = get_A()
    F = get_F()[0]
    D = get_D()[0]
    # Note that E, F, A are on the same line.
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
This function returns the center of the circle that goes through given three points.

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
    # vector normal to vecPQ
    vecNorm1 = rotate_vector(vecPQ[0], vecPQ[1], np.pi/2)
    # vector normal to vecQR
    vecNorm2 = rotate_vector(vecQR[0], vecQR[1], np.pi/2)
    # midpoint  of PQ
    midpt1 = [(P[0] + Q[0]) / 2, (P[1] + Q[1]) / 2]
    # midpoint of QR
    midpt2 = [(Q[0] + R[0]) / 2, (Q[1] + R[1]) / 2]
    # line1 = t*vecNorm1 + midpt1
    # line2 = s*vecNorm2 + midpt2
    t, s = symbols("t, s")
    # Solve system of equations.
    A = Matrix([[vecNorm1[0], -vecNorm2[0]],
                [vecNorm1[1], -vecNorm2[1]]])
    b = Matrix([midpt2[0] - midpt1[0], midpt2[1] - midpt1[1]])
    sols = list(list(linsolve((A, b), [t, s]))[0])
    # put t = sols[0] into line1
    # Get center
    return [sols[0]*vecNorm1[0] + midpt1[0], sols[0]*vecNorm1[1] + midpt1[1]]


"""
This function will return the parametric function of the unique circle
that goes through points P, Q, R.
"""
def circle_function(P, Q, R):
    center = get_Circle_center(P, Q, R)
    radius = distance(P, center)
    return center, radius

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

    #Appropriate input value to show each line segment.
    t_f1 = np.linspace(0, float(get_C()[1]), 2)
    t_f3 = np.linspace(0, float(get_E()[1]), 2)
    t_g1 = np.linspace(0, float(get_F()[2]), 2)
    t_h1 = np.linspace(0, float(get_O()[1]), 2)
    t_h3 = np.linspace(0, float(get_O()[2]), 2)
    t_g3 = np.linspace(0, float(get_D()[1]), 2)
    t_g2 = np.linspace(0, float(get_D()[2]), 2)
    t_h2 = np.linspace(0, float(get_t_for_h2_O_to_D()), 2)
    t_f2 = np.linspace(float(get_E()[2]), float(get_C()[2]), 2)

    points = {
    'A': A,
    'B': B,
    'F': F,
    'O': O,
    'D': D,
    'E': E,
    'C': C
    }

    # Create a scatter trace for each point
    fig0 = go.Figure()

    for name, coords in points.items():
        fig0.add_trace(
            go.Scatter(
                x=[np.array(coords[0], dtype=np.float32)],
                y=[np.array(coords[1], dtype=np.float32)],
                mode='markers',
                marker=dict(color='black'),
                name=name
            )
        )

    # line segments
    fig1 = px.line(pd.DataFrame(dict(x = t_f1, y = np.array([0]*len(t_f1)))), x="x", y="y", color_discrete_sequence=['black']) #f1 is defined to be on the x-axis.
    fig2 = px.line(pd.DataFrame(dict(x = np.array(f2(t_f2)[0], dtype=np.float32), y = np.array(f2(t_f2)[1], dtype=np.float))), x="x", y="y", color_discrete_sequence=['black']) 
    fig3 = px.line(pd.DataFrame(dict(x = np.array(f3(t_f3)[0], dtype=np.float32), y = np.array(f3(t_f3)[1], dtype=np.float))), x="x", y="y", color_discrete_sequence=['black']) 
    fig4 = px.line(pd.DataFrame(dict(x = np.array(g1(t_g1)[0], dtype=np.float32), y = np.array(g1(t_g1)[1], dtype=np.float))), x="x", y="y", color_discrete_sequence=['cyan']) 
    fig5 = px.line(pd.DataFrame(dict(x = np.array(g2(t_g2)[0], dtype=np.float32), y = np.array(g2(t_g2)[1], dtype=np.float))), x="x", y="y", color_discrete_sequence=['cyan']) 
    fig6 = px.line(pd.DataFrame(dict(x = np.array(g3(t_g3)[0], dtype=np.float32), y = np.array(g3(t_g3)[1], dtype=np.float))), x="x", y="y", color_discrete_sequence=['cyan']) 
    fig7 = px.line(pd.DataFrame(dict(x = np.array(h1(t_h1)[0], dtype=np.float32), y = np.array(h1(t_h1)[1], dtype=np.float))), x="x", y="y", color_discrete_sequence=['black']) 
    fig8 = px.line(pd.DataFrame(dict(x = np.array(h2(t_h2)[0], dtype=np.float32), y = np.array(h2(t_h2)[1], dtype=np.float))), x="x", y="y", color_discrete_sequence=['black']) 
    fig9 = px.line(pd.DataFrame(dict(x = np.array(h3(t_h3)[0], dtype=np.float32), y = np.array(h3(t_h3)[1], dtype=np.float))), x="x", y="y", color_discrete_sequence=['black']) 
 
    fig_all = go.Figure(data=fig0.data + fig1.data + fig2.data + fig3.data + fig4.data + fig5.data + fig6.data + fig7.data + fig8.data + fig9.data)
    
    # circles
    circle1_ctr, circle1_rad = circle_function(O, A, B)
    circle2_ctr, circle2_rad = circle_function(D, C, B)
    circle3_ctr, circle3_rad = circle_function(F, E, O)
    fig_all.add_shape(type="circle", xref="x", yref="y",x0=np.array(circle1_ctr[0]-circle1_rad, dtype=np.float32), y0=np.array(circle1_ctr[1]-circle1_rad, dtype=np.float), x1=np.array(circle1_ctr[0]+circle1_rad, dtype=np.float), y1=np.array(circle1_ctr[1]+circle1_rad, dtype=np.float),line_color='blue')
    fig_all.add_shape(type="circle", xref="x", yref="y",x0=np.array(circle2_ctr[0]-circle2_rad, dtype=np.float32), y0=np.array(circle2_ctr[1]-circle2_rad, dtype=np.float), x1=np.array(circle2_ctr[0]+circle2_rad, dtype=np.float), y1=np.array(circle2_ctr[1]+circle2_rad, dtype=np.float),line_color='blue')
    fig_all.add_shape(type="circle", xref="x", yref="y",x0=np.array(circle3_ctr[0]-circle3_rad, dtype=np.float32), y0=np.array(circle3_ctr[1]-circle3_rad, dtype=np.float), x1=np.array(circle3_ctr[0]+circle3_rad, dtype=np.float), y1=np.array(circle3_ctr[1]+circle3_rad, dtype=np.float),line_color='blue')
    
    fig_all.update_traces(marker={'size': 15})
    fig_all.show()

    #pyo.plot(fig_all, filename='plot1.html', auto_open=False)


    
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
    points_on_TU = 25
    points_on_UP = 25
    # We are building a mesh here. tt and ss will be the mesh
    # interval is the numbers in [0, 1] that t will take
    interval = np.linspace(0, 1, points_on_TU)
    tt = []
    ss = []
    for i in range(len(interval)):
        # By the property of similar triangles, s can take values 0 <= s <= t with 0 <= t <= 1.
        # Here, we are fixing t and obtaining points on the base of the triangle similar to TUV, which
        # side is length t*TU. Then we have UP = s*UV.
        s = list(np.linspace(0, interval[i], points_on_UP))
        t = list(np.array([interval[i]]*points_on_UP))
        tt.append(t)
        ss.append(s)
    # converting to arrays
    tt = np.array([np.array(x) for x in tt])
    ss = np.array([np.array(x) for x in ss])
    # P is an arbitrary point in the interior or the boundary of triangle TUV.
    # vecTP = t*vecTU + s *vecTV, 0<=t<=1, 0<=s<=t
    # Let O(0, 0, 0) and then
    # vecOP = t*vecTU + s *vecTV + vecOT, 0<=t<=1, 0<=s<=t
    vecTU = [U[0]-T[0], U[1]-T[1]]
    vecUV = [V[0] - U[0], V[1] - U[1]]
    # vecOP = [x, y, z].
    # The x, y, values P will take
    x = tt * vecTU[0] + ss * vecUV[0] + T[0]
    y = tt * vecTU[1] + ss * vecUV[1] + T[1]
    # Now, we will get the point on the sphere which has the same x, y coordinates as P.

    # Below is for obtaining the equation of the sphere which base is the circle going through T, U, V.
    # center of the circle containing T, U, V
    circleCenter = get_Circle_center(T, U, V)
    r = distance(T, circleCenter)
    
    # This is for obtaining the point on the sphere that has the same x, y coordinates as P.
    z = (r**2 - (x - circleCenter[0])**2 - (y - circleCenter[1])**2+ eps)**(1/2)
    return [x.astype(float), y.astype(float), z.astype(float)]


"""
Funtion for cut plane (portion of vertical half plane consisting a face of the cube).
P, Q are the line segment from which the orthogonal cut plane extends.
By cut plane we meen not plotting the circle which is made by the sphere
and the plane intersecting.
"""
def orthogonal_plane(P, Q, max_z):
    points_on_PQ = 500
    points_in_z_direction  = 500

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
        zrange = np.linspace(float(zvalue), max_z, points_in_z_direction)
        zz.append(zrange)
    zz = np.array(zz)
    return [xx.astype(float), yy.astype(float), zz.astype(float)]


def generate_circle_points(ctr, rad, t_vals):
    return [ctr[0] + rad * np.cos(t_vals), ctr[1] + rad * np.sin(t_vals)]

"""
This function plots a convex ideal cube in the upper half-space.
"""
def plot_upper_half():

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
    t_f1 = np.linspace(0, float(get_C()[1]), 2)
    t_f3 = np.linspace(0, float(get_E()[1]), 2)
    t_g1 = np.linspace(0, float(get_F()[2]), 2)
    t_h1 = np.linspace(0, float(get_O()[1]), 2)
    t_h3 = np.linspace(0, float(get_O()[2]), 2)
    t_g3 = np.linspace(0, float(get_D()[1]), 2)
    t_g2 = np.linspace(0, float(get_D()[2]), 2)
    t_h2 = np.linspace(0, float(get_t_for_h2_O_to_D()), 2)
    t_f2 = np.linspace(float(get_E()[2]), float(get_C()[2]), 2)

    x_vals = [float(A[0]), float(B[0]), float(F[0]), float(O[0]), float(D[0]), float(E[0]), float(C[0])]
    y_vals = [float(A[1]), float(B[1]), float(F[1]), float(O[1]), float(D[1]), float(E[1]), float(C[1])]
    z_vals = [0, 0, 0, 0, 0, 0, 0]  # All points are on the xy-plane

    labels = ['A', 'B', 'F', 'O', 'D', 'E', 'C']

    fig = go.Figure(data=[go.Scatter3d(
        x=x_vals,
        y=y_vals,
        z=z_vals,
        mode='markers+text',
        text=labels,
        textposition="top center",
        showlegend=False,
        marker=dict(size=5, color='black')
    )])

    fig.add_trace(go.Scatter3d(x=np.array(f1(t_f1)[0], dtype=np.float), y=np.array([0] * len(t_f1), dtype=np.float), z=np.array([0]*len(t_f1), dtype=np.float), mode='lines', line=dict(color='black'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(f2(t_f2)[0], dtype=np.float), y=np.array(f2(t_f2)[1], dtype=np.float), z=np.array([0]*len(t_f2), dtype=np.float), mode='lines', line=dict(color='black'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(f3(t_f3)[0], dtype=np.float), y=np.array(f3(t_f3)[1], dtype=np.float), z=np.array([0]*len(t_f3), dtype=np.float), mode='lines', line=dict(color='black'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(g1(t_g1)[0], dtype=np.float), y=np.array(g1(t_g1)[1], dtype=np.float), z=np.array([0]*len(t_g1), dtype=np.float), mode='lines', line=dict(color='cyan'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(g2(t_g2)[0], dtype=np.float), y=np.array(g2(t_g2)[1], dtype=np.float), z=np.array([0]*len(t_g2), dtype=np.float), mode='lines', line=dict(color='cyan'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(g3(t_g3)[0], dtype=np.float), y=np.array(g3(t_g3)[1], dtype=np.float), z=np.array([0]*len(t_g3), dtype=np.float), mode='lines', line=dict(color='cyan'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(h1(t_h1)[0], dtype=np.float), y=np.array(h1(t_h1)[1], dtype=np.float), z=np.array([0]*len(t_h1), dtype=np.float), mode='lines', line=dict(color='black'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(h2(t_h2)[0], dtype=np.float), y=np.array(h2(t_h2)[1], dtype=np.float), z=np.array([0]*len(t_h2), dtype=np.float), mode='lines', line=dict(color='black'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(h3(t_h3)[0], dtype=np.float), y=np.array(h3(t_h3)[1], dtype=np.float), z=np.array([0]*len(t_h3), dtype=np.float), mode='lines', line=dict(color='black'), showlegend=False))
    
    circle1_ctr, circle1_rad = circle_function(O, A, B)
    circle2_ctr, circle2_rad = circle_function(D, C, B)
    circle3_ctr, circle3_rad = circle_function(F, E, O)
    t_circ = np.linspace(0, 2 * np.pi, 100)  # Generating 100 points to represent the circle
    circle1_points = generate_circle_points(circle1_ctr, circle1_rad, t_circ)
    circle2_points = generate_circle_points(circle2_ctr, circle2_rad, t_circ)
    circle3_points = generate_circle_points(circle3_ctr, circle3_rad, t_circ)

    fig.add_trace(go.Scatter3d(x=np.array(circle1_points[0], dtype=np.float), y=np.array(circle1_points[1], dtype=np.float), z=np.array([0] * len(t_circ), dtype=np.float), mode='lines', line=dict(color='blue'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(circle2_points[0], dtype=np.float), y=np.array(circle2_points[1], dtype=np.float), z=np.array([0] * len(t_circ), dtype=np.float), mode='lines', line=dict(color='blue'), showlegend=False))
    fig.add_trace(go.Scatter3d(x=np.array(circle3_points[0], dtype=np.float), y=np.array(circle3_points[1], dtype=np.float), z=np.array([0] * len(t_circ), dtype=np.float), mode='lines', line=dict(color='blue'), showlegend=False))

    zmax_val = -float("inf")
    #spheres
    testing = spherical_triangle(O, B, A)
    uniform_blue = [[0, 'rgb(0, 0, 255)'], [1, 'rgb(0, 0, 255)']]  # A colorscale that's blue throughout
    fig.add_trace(go.Surface(
        z=testing[2],
        x=testing[0],
        y=testing[1],
        colorscale=uniform_blue,
        showscale=False,
        cmin=0,
        cmax=1
    ))
    zmax_val = max(zmax_val, np.max(testing[2]))

    testing = spherical_triangle(O, F, A)
    fig.add_trace(go.Surface(
        z=testing[2],
        x=testing[0],
        y=testing[1],
        colorscale=uniform_blue,
        showscale=False,
        cmin=0,
        cmax=1
    ))
    zmax_val = max(zmax_val, np.max(testing[2]))

    testing = spherical_triangle(F, E, D)
    uniform_green = [[0, 'rgb(0, 255, 0)'], [1, 'rgb(0, 255, 0)']]  # A colorscale that's blue throughout
    fig.add_trace(go.Surface(
        z=testing[2],
        x=testing[0],
        y=testing[1],
        colorscale=uniform_green,
        showscale=False,
        cmin=0,
        cmax=1
    ))
    zmax_val = max(zmax_val, np.max(testing[2]))
 
    testing = spherical_triangle(F, O, D)
    fig.add_trace(go.Surface(
        z=testing[2],
        x=testing[0],
        y=testing[1],
        colorscale=uniform_green,
        showscale=False,
        cmin=0,
        cmax=1
    ))
    zmax_val = max(zmax_val, np.max(testing[2]))

    testing = spherical_triangle(D, C, B)
    uniform_red = [[0, 'rgb(255, 0, 0)'], [1, 'rgb(255, 0, 0)']]  # A colorscale that's blue throughout
    fig.add_trace(go.Surface(
        z=testing[2],
        x=testing[0],
        y=testing[1],
        colorscale=uniform_red,
        showscale=False,
        cmin=0,
        cmax=1
    ))
    zmax_val = max(zmax_val, np.max(testing[2]))

    testing = spherical_triangle(D, O, B)
    fig.add_trace(go.Surface(
        z=testing[2],
        x=testing[0],
        y=testing[1],
        colorscale=uniform_red,
        showscale=False,
        cmin=0,
        cmax=1
    ))
    zmax_val = max(zmax_val, np.max(testing[2]))

    #planes
    uniform_purple = [[0, 'rgb(153, 153, 255)'], [1, 'rgb(153, 153, 255)']]
    planeAB = orthogonal_plane(A, B, zmax_val)
    fig.add_trace(go.Surface(
        z=planeAB[2],
        x=planeAB[0],
        y=planeAB[1],
        colorscale=uniform_purple,
        showscale=False,
        opacity=0.4  # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    planeBC = orthogonal_plane(B, C, zmax_val)
    fig.add_trace(go.Surface(
        z=planeBC[2],
        x=planeBC[0],
        y=planeBC[1],
        colorscale=uniform_purple,
        showscale=False,
        opacity=0.4  # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    uniform_salmon = [[0, 'rgb(255, 153, 153)'], [1, 'rgb(255, 153, 153)']]
    planeCD = orthogonal_plane(C, D, zmax_val)
    surface_color_data = np.full(planeCD[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=planeCD[2],
        x=planeCD[0],
        y=planeCD[1],
        colorscale=uniform_salmon,
        surfacecolor=surface_color_data,
        showscale=False,
        opacity=0.4  # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    planeDE = orthogonal_plane(D, E, zmax_val)
    surface_color_data = np.full(planeDE[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=planeDE[2],
        x=planeDE[0],
        y=planeDE[1],
        colorscale=uniform_salmon,
        surfacecolor=surface_color_data,
        showscale=False,
        opacity=0.4  # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    uniform_lightgreen = [[0, 'rgb(153, 255, 204)'], [1, 'rgb(153, 255, 204)']]
    planeEF = orthogonal_plane(E, F, zmax_val)
    surface_color_data = np.full(planeEF[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=planeEF[2],
        x=planeEF[0],
        y=planeEF[1],
        colorscale=uniform_lightgreen,
        surfacecolor=surface_color_data,
        showscale=False,
        opacity=0.4  # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    planeFA = orthogonal_plane(F, A, zmax_val)
    surface_color_data = np.full(planeFA[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=planeFA[2],
        x=planeFA[0],
        y=planeFA[1],
        colorscale=uniform_lightgreen,
        surfacecolor=surface_color_data,
        showscale=False,
        opacity=0.4  # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    fig.show()

    #pyo.plot(fig, filename='plot2.html', auto_open=False)
    
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

    # Vertices
    A = get_A()
    B = get_B()
    C = get_C()[0]
    D = get_D()[0]
    E = get_E()[0]
    F = get_F()[0]
    O = get_O()[0]

    ball_O = upper_half_to_ball_point([O[0], O[1], 0])
    ball_A = upper_half_to_ball_point([A[0], A[1], 0])
    ball_B = upper_half_to_ball_point([B[0], B[1], 0])
    ball_C = upper_half_to_ball_point([C[0], C[1], 0])
    ball_D = upper_half_to_ball_point([D[0], D[1], 0])
    ball_E = upper_half_to_ball_point([E[0], E[1], 0])
    ball_F = upper_half_to_ball_point([F[0], F[1], 0])

    x_vals = [float(ball_A[0]), float(ball_B[0]), float(ball_F[0]), float(ball_O[0]), float(ball_D[0]), float(ball_E[0]), float(ball_C[0]), 0]
    y_vals = [float(ball_A[1]), float(ball_B[1]), float(ball_F[1]), float(ball_O[1]), float(ball_D[1]), float(ball_E[1]), float(ball_C[1]), 0]
    z_vals = [float(ball_A[2]), float(ball_B[2]), float(ball_F[2]), float(ball_O[2]), float(ball_D[2]), float(ball_E[2]), float(ball_C[2]), 2]  # All points are on the xy-plane
    
    labels = ['A', 'B', 'F', 'O', 'D', 'E', 'C', "inversion center"]

    fig = go.Figure(data=[go.Scatter3d(
        x=x_vals,
        y=y_vals,
        z=z_vals,
        mode='markers+text',
        text=labels,
        textposition="top center",
        marker=dict(size=5, color='black')
    )])

    # Sphere of inversion
    uniform_gray = [[0, 'rgb(128,128,128)'], [1, 'rgb(128,128,128)']]
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    fig.add_trace(go.Surface(z=np.cos(v) + 1,
        x=np.cos(u) * np.sin(v),
        y=np.sin(u) * np.sin(v),
        colorscale=uniform_gray,
        showscale=False,
        opacity=0.2  # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
        ))

    # planes
    planeAB = upper_half_to_ball_surface(orthogonal_plane(A, B, 1000))
    planeBC = upper_half_to_ball_surface(orthogonal_plane(B, C, 1000))
    planeCD = upper_half_to_ball_surface(orthogonal_plane(C, D, 1000))
    planeDE = upper_half_to_ball_surface(orthogonal_plane(D, E, 1000))
    planeEF = upper_half_to_ball_surface(orthogonal_plane(E, F, 1000))
    planeFA = upper_half_to_ball_surface(orthogonal_plane(F, A, 1000))

    uniform_purple = [[0, 'rgb(153, 153, 255)'], [1, 'rgb(153, 153, 255)']]
    # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=planeAB[2],
        x=planeAB[0],
        y=planeAB[1],
        colorscale=uniform_purple,
        showscale=False,
        cmin=0,
        cmax=1
        # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=planeBC[2],
        x=planeBC[0],
        y=planeBC[1],
        colorscale=uniform_purple,
        showscale=False,
        cmin=0,
        cmax=1
        # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    uniform_salmon = [[0, 'rgb(255, 153, 153)'], [1, 'rgb(255, 153, 153)']]
    fig.add_trace(go.Surface(
        z=planeCD[2],
        x=planeCD[0],
        y=planeCD[1],
        colorscale=uniform_salmon,
        showscale=False,
        cmin=0,
        cmax=1
        # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    fig.add_trace(go.Surface(
        z=planeDE[2],
        x=planeDE[0],
        y=planeDE[1],
        colorscale=uniform_salmon,
        showscale=False,
        cmin=0,
        cmax=1
        # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    uniform_lightgreen = [[0, 'rgb(153, 255, 204)'], [1, 'rgb(153, 255, 204)']]
    fig.add_trace(go.Surface(
        z=planeEF[2],
        x=planeEF[0],
        y=planeEF[1],
        colorscale=uniform_lightgreen,
        showscale=False,
        cmin=0,
        cmax=1
        # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))

    fig.add_trace(go.Surface(
        z=planeFA[2],
        x=planeFA[0],
        y=planeFA[1],
        colorscale=uniform_lightgreen,
        showscale=False,
        cmin=0,
        cmax=1
         # Make the plane translucent. Adjust this value (between 0 and 1) to achieve the desired translucency.
    ))
    
    # spheres
    uniform_blue = [[0, 'rgb(0, 0, 255)'], [1, 'rgb(0, 0, 255)']]  # A colorscale that's blue throughout
    face1 = upper_half_to_ball_surface(spherical_triangle(O, B, A))
    face2 = upper_half_to_ball_surface(spherical_triangle(O, F, A))
    surface_color_data = np.full(face1[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=face1[2],
        x=face1[0],
        y=face1[1],
        colorscale=uniform_blue,
        surfacecolor=surface_color_data,
        showscale=False,
        cmin=0,
        cmax=1
    ))

    surface_color_data = np.full(face2[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=face2[2],
        x=face2[0],
        y=face2[1],
        colorscale=uniform_blue,
        surfacecolor=surface_color_data,
        showscale=False,
        cmin=0,
        cmax=1
    ))

    uniform_green = [[0, 'rgb(0, 255, 0)'], [1, 'rgb(0, 255, 0)']]  # A colorscale that's blue throughout
    face3 = upper_half_to_ball_surface(spherical_triangle(F, E, D))
    face4 = upper_half_to_ball_surface(spherical_triangle(F, O, D))
    surface_color_data = np.full(face3[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=face3[2],
        x=face3[0],
        y=face3[1],
        colorscale=uniform_green,
        surfacecolor=surface_color_data,
        showscale=False,
        cmin=0,
        cmax=1
    ))
    
    surface_color_data = np.full(face4[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=face4[2],
        x=face4[0],
        y=face4[1],
        colorscale=uniform_green,
        surfacecolor=surface_color_data,
        showscale=False,
        cmin=0,
        cmax=1
    ))

    face5 = upper_half_to_ball_surface(spherical_triangle(D, C, B))
    face6 = upper_half_to_ball_surface(spherical_triangle(D, O, B))
    uniform_red = [[0, 'rgb(255, 0, 0)'], [1, 'rgb(255, 0, 0)']]  # A colorscale that's blue throughout
    surface_color_data = np.full(face5[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=face5[2],
        x=face5[0],
        y=face5[1],
        colorscale=uniform_red,
        surfacecolor=surface_color_data,
        showscale=False,
        cmin=0,
        cmax=1
    ))

    surface_color_data = np.full(face6[2].shape, 0)  # Create a 2D array with the same shape as your z data, but filled with zeros.
    fig.add_trace(go.Surface(
        z=face6[2],
        x=face6[0],
        y=face6[1],
        colorscale=uniform_red,
        surfacecolor=surface_color_data,
        showscale=False,
        cmin=0,
        cmax=1
    ))
    
    fig.show()
    #pyo.plot(fig, filename='plot3.html', auto_open=False)


vertices_on_complex_plane()
plot_upper_half()
ball_model()

def main():
    pass