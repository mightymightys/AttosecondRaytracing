"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%%
import numpy as np
import ART.ModuleGeometry as mgeo


class OpticalElement:
    """ 
    Type is an object of a type-class of element : a mirror or a mask  for the moment
    Position of the center
    Normal is normal axis of the element at center location
    MajorAxis is the direction of the "distinguished axis" of non-rotationally symmetric elements
    like the major axes of toroidal/elliptical mirrors or the off-axis direction of off-axis parabolas.
    That's usually the x-axis in the mirror-coordinate system
    """
    
    def __init__(self, Type, Position, Normal, MajorAxis):
            self._type = Type
            self.position = Position
            self.normal = mgeo.Normalize(Normal) 
            self.majoraxis = mgeo.Normalize(MajorAxis) 
        
    # using property decorator 
    # a getter function 
    @property
    def position(self): 
        return self._position
       
    # a setter function 
    @position.setter 
    def position(self, NewPosition): 
        if type(NewPosition) == np.ndarray and len(NewPosition) == 3: 
              self._position = NewPosition
        else: raise TypeError('Position must be a 3D numpy.ndarray.')
        
        
    @property
    def normal(self): 
        return self._normal 
       
    @normal.setter 
    def normal(self, NewNormal): 
        if type(NewNormal) == np.ndarray and len(NewNormal) == 3 and np.linalg.norm(NewNormal) >0:
           try: 
               if abs(np.dot(mgeo.Normalize(NewNormal),self.majoraxis)) >1e-12: 
              #if the new normal is not perpendicular to the majoraxis, then we rotate the major along with the rotation of the normal vector
                   self._majoraxis = mgeo.RotationAroundAxis(np.cross(self.normal,NewNormal), mgeo.AngleBetweenTwoVectors(self.normal,NewNormal), self.majoraxis)      
           except: pass #this is for the initialization when the majoraxis isn't defined yet and the test above fails
        
           self._normal = mgeo.Normalize(NewNormal)
        else: raise TypeError('Normal must be a 3D numpy.ndarray with finite length.')
        
        
    @property
    def majoraxis(self): 
        return self._majoraxis 
    
    @majoraxis.setter 
    def majoraxis(self, NewMajorAxis): 
        if type(NewMajorAxis) == np.ndarray and len(NewMajorAxis) == 3 and np.linalg.norm(NewMajorAxis) >0:
            if abs(np.dot(self.normal, mgeo.Normalize(NewMajorAxis))) >1e-12: 
                raise ValueError('The normal and major axis of optical elements need to be orthogonal!')
            else:
                self._majoraxis = mgeo.Normalize(NewMajorAxis)
        else: raise TypeError('MajorAxis must be a 3D numpy.ndarray with finite length.')
   
    #make the type property private and providing only a getter method, so it can't be modified after the class instance has been created
    @property
    def type(self): 
        return self._type
 
    
    def __hash__(self):
        position_tuple = tuple(self.position.reshape(1, -1)[0])
        normal_tuple = tuple(self.normal.reshape(1, -1)[0])
        majoraxis_tuple = tuple(self.majoraxis.reshape(1, -1)[0])
        return hash(position_tuple + normal_tuple + majoraxis_tuple) + hash(self.type)
    
    
#%% methods to (mis-)align the OE

    def rotate_pitch_by(self, angle): 
        """
        pitch, i.e. rotation about the axis (NormalAxis x MajorAxis). angle is given in degrees.
        If the plane spanned by NormalAxis and MajorAxis is the incidence plane (normally the case
        in a "clean alignment" situation), then this is simply a modificaiton of the incidence angle by "angle".
        But in general, if the optical element has some odd orientation, there is not a direct correspondence.
        """
        rotation_axis = np.cross(self.normal,self.majoraxis)
        self.normal = mgeo.RotationAroundAxis(rotation_axis, np.deg2rad(angle), self.normal)  
        #the normal.setter function should take care of the majoraxis remaining perpendicular to the normal.
    
    def rotate_roll_by(self, angle): 
        """ roll, i.e. rotation about the MajorAxis. angle is given in degrees. """
        self.normal = mgeo.RotationAroundAxis(self.majoraxis, np.deg2rad(angle), self.normal)  
        
    def rotate_yaw_by(self, angle): 
        """ yaw, i.e. rotation about the NormalAxis. angle is given in degrees. """
        self.majoraxis = mgeo.RotationAroundAxis(self.normal, np.deg2rad(angle), self.majoraxis)  
       
    def rotate_random_by(self, angle): 
        """ rotation about a random axis. angle is given in degrees. """
        self.normal = mgeo.RotationAroundAxis(np.random.random(3), np.deg2rad(angle), self.normal) 
        

    def shift_along_normal(self, distance):
        """ shift along the normal-axis direction. distance is given in mm. """
        self.position += distance*self.normal 
        
    def shift_along_major(self, distance):
        """ shift along the major-axis direction. distance is given in mm. """
        self.position += distance*self.majoraxis
        
    def shift_along_cross(self, distance):
        """ shift along the direction (normal x major). distance is given in mm. """
        self.position += distance*mgeo.Normalize(np.cross(self.normal,self.majoraxis))
        
    def shift_along_random(self, distance):
        """ shift along a direction. distance is given in mm. """
        self.position += distance *mgeo.Normalize(np.random.random(3))
        