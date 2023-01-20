"""
Created in 2019

@author: Anthony Guillaume
"""
#%% Modules
import numpy as np
from quaternion import quaternion 

#%%
def Normalize(Vector):
    return Vector/np.linalg.norm(Vector)

#%%
def VectorPerpendicular(Vector):
    """
    Find a perpendicular 3D vector in some arbitrary direction
    """

    if abs(Vector[0]) < 1e-15:
        return np.array([1, 0, 0])
    if abs(Vector[1]) < 1e-15:
        return np.array([0, 1, 0])
    if abs(Vector[2]) < 1e-15:
        return np.array([0, 0, 1])

    # set arbitrarily a = b =1
    return Normalize(np.array([1, 1, -1.0 * (Vector[0] + Vector[1]) / Vector[2]]))

#%%
def AngleBetweenTwoVectors(U,V):
    """    Return the angle between the vectors U and V ; formula from W.Kahan    """
    u = np.linalg.norm(U) 
    v = np.linalg.norm(V)
    return 2*np.arctan2(np.linalg.norm(U*v - V*u), np.linalg.norm(U*v + V*u))


#%%
def IntersectionLinePlane(A, u, P, n):
    """
    Return the intersection point between a line and A plane. 
    A is a point of the line, u a vector of the line ; P is a point of the plane, n a normal vector
    Line's equation : OM = u*t + OA , t a real
    Plane's equation : n.OP - n.OM = 0
    """
    t = np.dot(n, -A + P) / np.dot(u, n)
    I = u*t + A
    return I


#%%
def SpiralVogel(NbPoint, Radius):
    """
    Return a NbPoint x 2 matrix of 2D points representative of Vogel's spiral with radius Radius
    """
    
    GoldenAngle = np.pi * (3 - np.sqrt(5))
    r = np.sqrt(np.arange(NbPoint) / NbPoint) * Radius
     
    theta = GoldenAngle * np.arange(NbPoint)
     
    Matrix = np.zeros((NbPoint, 2))
    Matrix[:,0] = np.cos(theta)
    Matrix[:,1] = np.sin(theta)
    Matrix = Matrix * r.reshape((NbPoint, 1))
    
    return Matrix

#%%   
def SolverQuadratic(a,b,c):
    """
    Solve the quadratic equation a*x^2 + b*x +c = 0 ; keep only real solutions
    """
    Solution = np.roots([a,b,c])
    RealSolution = []

    for k in range(len(Solution)):
        if abs(Solution[k].imag) < 1e-15:
            RealSolution.append(Solution[k].real)
            
    return RealSolution

#%%
def SolverQuartic(a,b,c,d,e):
    """
    Solve the quartic equation a*x^4 + b*x^3 +c*x^2 + d*x + e = 0 ; keep only real solutions
    """
    Solution = np.roots([a,b,c,d,e])
    SolutionReel = []

    for k in range(len(Solution)):
        if abs(Solution[k].imag) < 1e-15:
            SolutionReel.append(Solution[k].real)
            
    return SolutionReel

#%%   
def KeepPositiveSolution(SolutionList):
    """
    Keep only positive solution (numbers) in the list
    """
    PositiveSolutionList = []
    epsilon = 1e-12
    for k in SolutionList:
        if k > epsilon:
            PositiveSolutionList.append(k)
            
    return PositiveSolutionList

#%%
def ClosestPoint(A, I1, I2):
    """
    Return the closest point I1 or I2 from A
    """
    DistanceSquare1 = np.dot(I1 - A, I1 - A)
    DistanceSquare2 = np.dot(I2 - A, I2 - A)
    if DistanceSquare1 < DistanceSquare2: 
        return I1
    else : 
        return I2

#%%
def FarestPoint(A, I1, I2):
    """
    Return the farthest point I1 or I2 from A
    """
    DistanceSquare1 = np.dot(I1 - A, I1 - A)
    DistanceSquare2 = np.dot(I2 - A, I2 - A)
    if DistanceSquare1 > DistanceSquare2: 
        return I1
    else : 
        return I2

#%%        
def IncludeRectangle(X, Y, Point):
    """
    Check if a 2D Point is a element of the rectangle X x Y
    """
    x = Point[0]
    y = Point[1]
    if abs(x) <= abs(X*0.5) and abs(y) <= abs(Y*0.5):
        return True
    else:
        return False

#%%       
def IncludeDisk(R, Point):
    """
    Check if a 2D Point is a element of the disk of radius R
    """
    x = Point[0]
    y = Point[1]
    if (x**2 + y**2) <= R**2:
        return True
    else:
        return False
    
#%%
def SymmetricalVector(V, SymmetryAxis):
    """
    Return the symmetrical vector to V
    """
    return RotationAroundAxis(SymmetryAxis,np.pi,V)




#%%##############################################################################################
# TRANSLATIONS (very easy but still useful packed into functions) ###############################

def TranslationPoint(Point, T):
    """ Translate Point by Vector T """
    Pointprime = Point + T
    return Pointprime

#%%
def TranslationPointList(PointList, T):
    """ Translate all points in PointList by Vector T """
    PointListPrime = []
    for k in PointList:
        PointListPrime.append(TranslationPoint(k, T))
    return PointListPrime
    
#%%
def TranslationRay(Ray, T):
    """ Translate a Ray by vector T """
    Rayprime = Ray.copy_ray()
    Rayprime.point = Rayprime.point + T
    return Rayprime

#%%
def TranslationRayList(RayList, T):
    """ Translate all rays of RayList by vector T """
    RayListprime = []
    for k in RayList:
        RayListprime.append(TranslationRay(k, T))
    return RayListprime


#%%##############################################################################################
# ROTATIONS using quaternions in one way or another #############################################

def RotationAroundAxis(Axis, Angle, Vector):
    """    Rotate Vector by Angle (in rad) around Axis   """
    rot_axis = np.array([0.] + Axis)
    axis_angle = (Angle*0.5) * rot_axis/np.linalg.norm(rot_axis)
    vec = quaternion(*Vector)
    qlog = quaternion(*axis_angle)
    q = np.exp(qlog)
    VectorPrime = q * vec * np.conjugate(q)
    return VectorPrime.imag 

#%%
def RotationPoint(Point, Axis1, Axis2):
    """  Transform vector Point such that Axis1 in the old coordinate frame becomes vector Axis 2 in the new coordinate system  """
    Angle = AngleBetweenTwoVectors(Axis1,Axis2)
    if Angle < 1e-12:
        Pointprime = Point
    elif Angle > np.pi - 1e-12:
        Pointprime = -Point
    else:
        N = np.cross(Axis1,Axis2)
        Pointprime = RotationAroundAxis(N, Angle, Point)
    return Pointprime

#%%
def RotationPointList(PointList, Axis1, Axis2):
    """  Transform all vectors Point such that Axis1 in the old coordinate frame becomes vector Axis 2 in the new coordinate system  """
    PointListPrime = []
    for k in PointList:
        PointListPrime.append(RotationPoint(k, Axis1, Axis2))
    return PointListPrime

#%%
def RotationRay(Ray, Axis1, Axis2):
    """ Like RotationPoint but with a Ray object """
    Rayprime = Ray.copy_ray()
    OA = Ray.point
    u = Ray.vector
    OB = u + OA
    OAprime = RotationPoint(OA, Axis1, Axis2)
    OBprime = RotationPoint(OB, Axis1, Axis2)
    uprime = OBprime - OAprime
    Rayprime.point = OAprime
    Rayprime.vector = uprime
    return Rayprime

#%%
def RotationRayList(ListeRay, Axis1, Axis2):
    """ Like RotationPointList but with a list of Ray objects """
    ListeRayprime = []
    for k in ListeRay:
        ListeRayprime.append(RotationRay(k, Axis1, Axis2))
    return ListeRayprime

#%%
def RotationAroundAxisRayList(ListeRay, Axis, Angle):
    """ Like RotationAroundAxis but with a list of Ray objects """
    ListeRayPrime = []
    for k in ListeRay:
        Ray = k.copy_ray()
        VectorPrime = RotationAroundAxis(Axis, Angle, Ray.vector)
        Ray.vector = VectorPrime
        ListeRayPrime.append(Ray)
        
    return ListeRayPrime
