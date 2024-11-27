"""
Contains a bunch of useful function for geometric transformations and measurements.
Usually these don't need to be called by users of ART, but they may be useful.

Some general conventions:
- As much as possible, avoid having lists of points or vectors transiting between functions.
    Instead, use the Vector and Point classes defined at the end of this file.
- Functions operating on Rays should preferentially operate on lists of Rays, not individual Rays.
    The reason for that is that it's fairly rare to manipluate a single Ray, and it's easier to
    just put it in a single-element list and call the function that way.

Created in 2019

@author: Anthony Guillaume + Stefan Haessler + André Kalouguine
"""
# %% Modules
import numpy as np
import quaternion
from functools import singledispatch
import logging
import math

logger = logging.getLogger(__name__)

# %% Points and Vectors classes
class Vector(np.ndarray):
    def __new__(cls, input_array):
        input_array = np.asarray(input_array)
        if input_array.ndim > 1:
            return VectorArray(input_array)
        obj = input_array.view(cls)  # Ensure we are viewing as `Vector`
        return obj
    @property
    def norm(self):
        return np.linalg.norm(self)
    def normalized(self):
        return self / self.norm
    def __add__(self, other):
        if isinstance(other, Point):
            return Point(super().__add__(other))
        return Vector(super().__add__(other))
    def __sub__(self, other):
        if isinstance(other, Point):
            return Point(super().__sub__(other))
        return Vector(super().__sub__(other))
    def translate(self, vector):
        return self
    def rotate(self,q):
        return Vector((q*np.quaternion(0,*self)*q.conj()).imag)
    def from_basis(self, r0, r, q):
        return self.rotate(q)
    def to_basis(self, r0, r, q):
        return self.rotate(q.conj())
    def _add_dimension(self, value=0):
        return Vector(np.concatenate((self, [value])))
    def __hash__(self) -> int:
        vector_tuple = tuple(self.reshape(1, -1)[0])
        return hash(vector_tuple)
class VectorArray(np.ndarray):
    def __new__(cls, input_array):
        input_array = np.asarray(input_array)
        if input_array.ndim == 1:
            return Vector(input_array)
        obj = input_array.view(cls)  # Ensure we are viewing as `Vector`
        return obj
    @property
    def norm(self):
        return np.linalg.norm(self, axis=1)
    def __getitem__(self, index):
        return Vector(super().__getitem__(index))
    def normalized(self):
        return self / self.norm[:, np.newaxis]
    def __add__(self, other):
        if isinstance(other, Point):
            return PointArray(super().__add__(other))
        return VectorArray(super().__add__(other))
    def __sub__(self, other):
        if isinstance(other, Point):
            return PointArray(super().__sub__(other))
        return VectorArray(super().__sub__(other))
    def translate(self, vector):
        return self
    def rotate(self,q):
        return VectorArray(quaternion.rotate_vectors(q, self))
    def from_basis(self, r0, r, q):
        return self.rotate(q)
    def to_basis(self, r0, r, q):
        return self.rotate(q.conj())
    def _add_dimension(self, values=0):
        if values == 0:
            values = np.zeros((self.shape[0], 1))
        if not isinstance(values, np.ndarray):
            values = np.array(values)
        return VectorArray(np.concatenate((self, values), axis=1))
    
class Point(np.ndarray):
    def __new__(cls, input_array):
        input_array = np.asarray(input_array)
        if input_array.ndim > 1:
            return PointArray(input_array)
        obj = input_array.view(cls)  # Ensure we are viewing as `Point`
        return obj
    def __add__(self, other):
        return Point(super().__add__(other))
    def __sub__(self, other):
        return Vector(super().__sub__(other))
    def translate(self, vector):
        return Point(super()._add__(vector))
    def rotate(self,q):
        return (self-Origin).rotate(q)
    def from_basis(self, r0, r, q):
        return Point([0,0,0]) + r0 + Vector(self-r0).rotate(q) + r
    def to_basis(self, r0, r, q):
        return Point([0,0,0]) + r0 + Vector(self-r0-r).rotate(q.conj())
    def _add_dimension(self, value=0):
        return Point(np.concatenate((self, [value])))
    def __hash__(self) -> int:
        point_tuple = tuple(self.reshape(1, -1)[0])
        return hash(point_tuple)

class PointArray(np.ndarray):
    def __new__(cls, input_array):
        input_array = np.asarray(input_array)
        if input_array.ndim == 1:
            return Point(input_array)
        obj = input_array.view(cls)  # Ensure we are viewing as `Point`
        return obj
    def __getitem__(self, index):
        return Point(super().__getitem__(index))
    def __add__(self, other):
        return PointArray(super().__add__(other))
    def __sub__(self, other):
        return VectorArray(super().__sub__(other))
    def translate(self, vector):
        return PointArray(super().__add__(vector))
    def rotate(self,q):
        return PointArray(quaternion.rotate_vectors(q, self))
    def from_basis(self, r0, r, q):
        return Point([0,0,0]) + r0 + VectorArray(self-r0).rotate(q) + r
        #return PointArray([Point([0,0,0]) + r0 + (p-r0).rotate(q) + r for p in self])
    def to_basis(self, r0, r, q):
        return Point([0,0,0]) + r0 + VectorArray(self-r0-r).rotate(q.conj())
        #return PointArray([Point([0,0,0]) + r0 + (p-r0-r).rotate(q.conj()) for p in self])
    def _add_dimension(self, values=0):
        if values == 0:
            values = np.zeros((self.shape[0], 1))
        if not isinstance(values, np.ndarray):
            values = np.array(values)
        return PointArray(np.concatenate((self, values), axis=1))
    
Origin = Point([0,0,0])
# %% More traditional vector operations that don't belong in the classes
def Normalize(vector):
    """
    Normalize Vector.
    Obsolete, use use the mgeo.Vector class instead as it has a `normalize` method.
    """
    return vector / np.linalg.norm(vector)

def VectorPerpendicular(vector):
    """
    Find a perpendicular 3D vector in some arbitrary direction
    Undefined behavior, use with caution. There is no unique perpendicular vector to a 3D vector.
    """
    logger.warning("VectorPerpendicular is undefined behavior. There is no unique perpendicular vector to a 3D vector.")
    if abs(vector[0]) < 1e-15:
        return Vector([1, 0, 0])
    if abs(vector[1]) < 1e-15:
        return Vector([0, 1, 0])
    if abs(vector[2]) < 1e-15:
        return Vector([0, 0, 1])

    # set arbitrarily a = b =1
    return Vector([1, 1, -1.0 * (vector[0] + vector[1]) / vector[2]]).normalized()

def AngleBetweenTwoVectors(U, V):
    """
    Return the angle in radians between the vectors U and V ; formula from W.Kahan
    Value in radians between 0 and pi.
    """
    u = np.linalg.norm(U)
    v = np.linalg.norm(V)
    return 2 * np.arctan2(np.linalg.norm(U * v - V * u), np.linalg.norm(U * v + V * u))

def SymmetricalVector(V, SymmetryAxis):
    """
    Return the symmetrical vector to V
    """
    q = QRotationAroundAxis(SymmetryAxis, np.pi)
    return V.rotate(q)

def normal_add(N1, N2):
    """
    Simple function that takes in two normal vectors of a deformation and calculates
    the total normal vector if the two deformations were individually applied.
    Be very careful, this *only* works when the surface is z = f(x,y) and when 
    the deformation is small.
    Might be made obsolete when shifting to a modernised deformation system.
    """
    normal1 = N1.normalized()
    normal2 = N2.normalized()
    grad1 = -normal1[:2] / normal1[2]
    grad2 = -normal2[:2] / normal2[2]
    grad = grad1 + grad2
    total_normal = np.append(-grad, 1)
    return Vector(total_normal).normalized()

# %% Intersection finding
def IntersectionLinePlane(A, u, P, n):
    """
    Return the intersection point between a line and a plane.
    A is a point of the line, u a vector of the line ; P is a point of the plane, n a normal vector
    Line's equation : OM = u*t + OA , t a real
    Plane's equation : n.OP - n.OM = 0
    """
    t = np.dot(n, -A + P) / np.dot(u, n)
    I = u * t + A
    return I

def IntersectionRayListZPlane(RayList, Z=np.array([0])):
    """
    Return the intersection of a list of rays with a different planes with equiations z = Z[i]
    Basically, by default it returns the intersection of the rays with the Z=0 plane but you can 
    give it a few values of Z and it should be faster than calling it multiple times.
    This should let us quickly find the optimal position of the detector as well as trace the caustics.
    If a ray does not intersect the plane... it should replace that point with a NaN. 
    """
    Positions = np.vstack([i.point for i in RayList])
    Vectors = np.vstack([i.vector for i in RayList])
    non_zero = Vectors[:,2] != 0
    Positions = Positions[non_zero]
    Vectors = Vectors[non_zero]
    Z = Z[:, np.newaxis]
    A = Positions[:,2]-Z
    B = -Vectors[:,2]
    #times = (Positions[:,2]-Z)/Vectors[:,2]
    #return A,B
    with np.errstate(divide='ignore', invalid='ignore'):
        times = np.divide(A, B, where=(B != 0), out=np.full_like(A, np.nan))
        #times[times < 0] = np.nan  # Set negative results to NaN
    #return times
    #positive_times = times >= 0
    intersect_positions = Positions[:, :2] + times[:, :, np.newaxis] * Vectors[:, :2]
    result = []
    for i in range(Z.shape[0]):
        # For each plane, we find the intersection points
        #valid_intersections = intersect_positions[i][positive_times[i]]
        valid_intersections = intersect_positions[i]
        result.append(PointArray(valid_intersections))
    return result


# %% Geometrical utilities for plotting
def SpiralVogel(NbPoint, Radius):
    """
    Return a NbPoint x 2 matrix of 2D points representative of Vogel's spiral with radius Radius
    Careful, contrary to most of the code, this is *not* in the 
    ARTcore.Vector or ARTcore.Point format. It is a simple numpy array. 
    The reason is that this is a utility function that can be used both to define directions
    and to generate grids of points.
    """
    GoldenAngle = np.pi * (3 - np.sqrt(5))
    r = np.sqrt(np.arange(NbPoint) / NbPoint) * Radius

    theta = GoldenAngle * np.arange(NbPoint)

    Matrix = np.zeros((NbPoint, 2))
    Matrix[:, 0] = np.cos(theta)
    Matrix[:, 1] = np.sin(theta)
    Matrix = Matrix * r.reshape((NbPoint, 1))

    return Matrix

def find_hull(points):
    """
    Find the convex hull of a set of points using a greedy algorithm.
    This is used to create a polygon that encloses the points.
    """
    # start from leftmost point
    current_point = min(range(len(points)), key=lambda i: points[i][0])
    # initialize hull with current point
    hull = [current_point]
    # initialize list of linked points
    linked = []
    # continue until all points have been linked
    while len(linked) < len(points) - 1:
        # initialize minimum distance and closest point
        min_distance = math.inf
        closest_point = None
        # find closest unlinked point to current point
        for i, point in enumerate(points):
            if i not in linked:
                distance = math.dist(points[current_point], point)
                if distance < min_distance:
                    min_distance = distance
                    closest_point = i
        # add closest point to hull and linked list
        hull.append(closest_point)
        linked.append(closest_point)
        # update current point
        current_point = closest_point
    # add link between last point and first point
    hull.append(hull[0])
    # convert hull to a list of pairs of indices
    indices = [[hull[i], hull[i + 1]] for i in range(len(hull) - 1)]
    return indices


# %% Solvers and utilities for solving equations
def SolverQuadratic(a, b, c):
    """
    Solve the quadratic equation a*x^2 + b*x +c = 0 ; keep only real solutions
    """
    Solution = np.roots([a, b, c])
    RealSolution = []

    for k in range(len(Solution)):
        if abs(Solution[k].imag) < 1e-15:
            RealSolution.append(Solution[k].real)

    return RealSolution


def SolverQuartic(a, b, c, d, e):
    """
    Solve the quartic equation a*x^4 + b*x^3 +c*x^2 + d*x + e = 0 ; keep only real solutions
    """
    Solution = np.roots([a, b, c, d, e])
    RealSolution = []

    for k in range(len(Solution)):
        if abs(Solution[k].imag) < 1e-15:
            RealSolution.append(Solution[k].real)

    return RealSolution


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


def KeepNegativeSolution(SolutionList):
    """
    Keep only positive solution (numbers) in the list
    """
    NegativeSolutionList = []
    epsilon = -1e-12
    for k in SolutionList:
        if k < epsilon:
            NegativeSolutionList.append(k)

    return NegativeSolutionList


# %% Point geometry tools
def ClosestPoint(A: Point, Points: PointArray):
    """
    Given a reference point A and an array of points, return the index of the point closest to A
    """
    distances = (Points-A).norm
    return np.argmin(distances)

def DiameterPointArray(Points: PointArray):
    """
    Return the diameter of the smallest circle (for 2D points) 
    or sphere (3D points) including all the points.
    """
    if len(Points) == 0:
        return None
    return float(np.ptp(Points, axis=0).max())

def CentrePointList(Points):
    """
    Shift all 2D-points in PointList so as to center the point-cloud on the origin [0,0].
    """
    return Points - np.mean(Points, axis=0)

# %% Solid object orientation

def RotateSolid(Object, q):
    """
    Rotate object around basepoint by quaternion q
    """
    Object.q = q*Object.q

def TranslateSolid(Object, T):
    """
    Translate object by vector T
    """
    Object.r = Object.r + T

def RotateSolidAroundInternalPointByQ(Object, q, P):
    """
    Rotate object around P by quaternion q where P is in the object's frame
    """
    pass #TODO

def RotateSolidAroundExternalPointByQ(Object, q, P):
    """Rotate object around P by quaternion q, where P is in the global frame"""
    pass #TODO


# %% Signed distance functions 

def SDF_Rectangle(Point, SizeX, SizeY):
    """Signed distance function for a rectangle centered at the origin"""
    d = np.abs(Point[:2]) - np.array([SizeX, SizeY]) / 2
    return (np.linalg.norm(np.maximum(d, 0)) + np.min(np.max(d, 0)))/2

def SDF_Circle(Point, Radius):
    """Signed distance function for a circle centered at the origin"""
    return np.linalg.norm(Point[:2]) - Radius

def Union_SDF(SDF1, SDF2):
    """Union of two signed distance functions"""
    return np.minimum(SDF1, SDF2)

def Difference_SDF(SDF1, SDF2):
    """Difference of two signed distance functions"""
    return np.maximum(SDF1, -SDF2)

def Intersection_SDF(SDF1, SDF2):
    """Intersection of two signed distance functions"""
    return np.maximum(SDF1, SDF2)



# %% Quaternion calculations
def QRotationAroundAxis(Axis, Angle):
    """
    Return quaternion for rotation by Angle (in rad) around Axis
    """
    rot_axis = Normalize(np.array([0.0] + Axis))
    axis_angle = (Angle * 0.5) * rot_axis
    qlog = np.quaternion(*axis_angle)
    q = np.exp(qlog)
    return q

def QRotationVector2Vector(Vector1, Vector2):
    """
    Return a possible quaternion (among many) that would rotate Vector1 into Vector2.
    Undefined behavior, use with caution. There is no unique quaternion that rotates one vector into another.
    """
    Vector1 = Normalize(Vector1)
    Vector2 = Normalize(Vector2)
    a = np.cross(Vector1, Vector2)
    return np.quaternion(1 + np.dot(Vector1, Vector2), *a).normalized()

def QRotationVectorPair2VectorPair(InitialVector1, Vector1, InitialVector2, Vector2):
    """
    Return the only quaternion that rotates InitialVector1 to Vector1 and InitialVector2 to Vector2.
    Please ensure orthogonality two input and two output vectors.
    """
    Vector1 = Normalize(Vector1)
    Vector2 = Normalize(Vector2)
    Vector3 = Normalize(np.cross(Vector1,Vector2))
    InitialVector1 = Normalize(InitialVector1)
    InitialVector2 = Normalize(InitialVector2)
    InitialVector3 = Normalize(np.cross(InitialVector1,InitialVector2))
    rot_2Initial = np.zeros((3,3))
    rot_2Initial[:,0] = InitialVector1
    rot_2Initial[:,1] = InitialVector2
    rot_2Initial[:,2] = InitialVector3
    rot_2Final = np.zeros((3,3))
    rot_2Final[:,0] = Vector1
    rot_2Final[:,1] = Vector2
    rot_2Final[:,2] = Vector3
    q2Init = quaternion.from_rotation_matrix(rot_2Initial)
    q2Fin = quaternion.from_rotation_matrix(rot_2Final)
    return (q2Fin/q2Init).normalized()


# %% RayList stuff
def RotationRayList(RayList, q):
    """Like RotationPointList but with a list of Ray objects"""
    return [i.rotate(q) for i in RayList]

def TranslationRayList(RayList, T):
    """Translate a RayList by vector T"""
    return [i.translate(T) for i in RayList]
