"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules
import os
import pickle
#import gzip
import lzma
from time import time
from datetime import datetime
import copy 
import numpy as np
import ART.ModuleGeometry as mgeo
import ART.ModuleMirror as mmirror
import ART.ModuleMask as mmask
import ART.ModuleSource as msource
import ART.ModuleOpticalElement as moe
import ART.ModuleOpticalRay as mray
import ART.ModuleSupport as msupp
import ART.ModuleOpticalChain as moc

#%% 
def hash_list_of_objects(list):
    total_hash = 0
    for x in list:
        total_hash += hash(x)
    return total_hash

def which_indexes(lst):
    """ takes an input-list, and gives you a list of the indeces where the input-list contains another list or a numpy-array """
    indexes = [i for i, x in enumerate(lst) if isinstance(x, (list, np.ndarray))]
    return indexes

#%%
# Save an object to a compressed file
def save_compressed(obj, filename):
    if not type(filename) == str:
        filename = 'kept_data_' + datetime.now().strftime("%Y-%m-%d-%Hh%M")
    
    i = 0
    while os.path.exists(filename + f'_{i}.xz'):
        i += 1
    filename = filename + f'_{i}'
    #with gzip.open(filename + '.gz', 'wb') as f:
    with lzma.open(filename + '.xz', 'wb') as f:
        pickle.dump(obj, f)
    print("Saved results to " + filename +".xz.")
    print("->To reload from disk do: kept_data = mp.load_compressed('" + filename +"')")

# Load an object from a compressed file
def load_compressed(filename):
    #with gzip.open(filename + '.gz', 'rb') as f:
    with lzma.open(filename + '.xz', 'rb') as f:
        obj = pickle.load(f)
    return obj

#%%  tic-toc timer functions
_tstart_stack = []
def tic():
    _tstart_stack.append(time())

def toc(fmt="Elapsed: %s s"):
    print(fmt % (time() - _tstart_stack.pop()))
    
    
#%%
def FindCentralRay(RayList):
    for k in RayList:
        if k.number == 0:
            return k
    return None


def StandardDeviation1D(List):      
    V = 0
    m = np.mean(List)
    for k in List:
        V = V + (k - m)**2
        
    return np.sqrt(V / len(List))


def StandardDeviation2D(List):       
    V = 0
    mx = np.mean([k[0] for k in List])
    my = np.mean([k[1] for k in List])
    for k in List:
        V = V + (k[0] - mx)**2 + (k[1] - my)**2
        
    return np.sqrt(V / len(List))

def StandardDeviation3D(List):       
    V = 0
    mx = np.mean([k[0] for k in List])
    my = np.mean([k[1] for k in List])
    mz = np.mean([k[2] for k in List])
    for k in List:
        V = V + (k[0] - mx)**2 + (k[1] - my)**2 + (k[2] - mz)**2
        
    return np.sqrt(V / len(List))

def StandardDeviation(List):
    if type(List[0]) in [int, float, np.float64]:
        return StandardDeviation1D(List)
    elif len(List[0]) == 2:
        return StandardDeviation2D(List)
    elif len(List[0]) == 3:
        return StandardDeviation3D(List)
    else:
        return None

#%%    
def WeightedStandardDeviation1D(List, Weights):
    V = 0
    sq = 0
    m, s = np.average(List, weights=Weights, returned=True)
    for k in range(len(List)):
        V = V + Weights[k]*(List[k] - m)**2
        sq = sq + Weights[k]**2     
    
    return np.sqrt(V / (s - sq/s) )
    
def WeightedStandardDeviation2D(List, Weights):
 
    V = 0
    sq = 0
    mx, s = np.average([k[0] for k in List], weights=Weights, returned=True)
    my    = np.average([k[1] for k in List], weights=Weights, returned=False)
    for k in range(len(List)):
        V = V + Weights[k]*( (List[k][0] - mx)**2 + (List[k][1] - my)**2  )
        sq = sq + Weights[k]**2     

    return np.sqrt(V / (s - sq/s) )
     
def WeightedStandardDeviation3D(List, Weights):
 
    V = 0
    sq = 0
    mx, s = np.average([k[0] for k in List], weights=Weights, returned=True)
    my    = np.average([k[1] for k in List], weights=Weights, returned=False)
    mz    = np.average([k[2] for k in List], weights=Weights, returned=False)
    
    for k in range(len(List)):
        V = V + Weights[k]*( (List[k][0] - mx)**2 + (List[k][1] - my)**2 + (List[k][2] - mz)**2 )
        sq = sq + Weights[k]**2     

    return np.sqrt(V / (s - sq/s) )

def WeightedStandardDeviation(List, Weights):
    Condition1 = type(List[0]) == int or type(List[0]) == float or type(List[0]) == np.float64
    Condition2 = type(Weights[0]) == int or type(Weights[0]) == float or type(Weights[0]) == np.float64
    Condition3 = len(List) == len(Weights)
    
    if Condition1 and Condition2 and Condition3:
        return WeightedStandardDeviation1D(List, Weights)
    elif len(List[0]) == 2 and Condition2 and Condition3:
        return WeightedStandardDeviation2D(List, Weights)
    elif len(List[0]) == 3 and Condition2 and Condition3:
        return WeightedStandardDeviation3D(List, Weights)
    else:
        return None    

#%%
def DiameterPointList(PointList):
    """
    Return the diameter of the smallest circle (for 2D points) or sphere (3D points) including all the points
    """
    if len(PointList) == 0:
        return None
        
    elif len(PointList[0]) == 2:
        Xmax = PointList[0][0]
        Xmin = PointList[0][0]
        Ymax = PointList[0][1]
        Ymin = PointList[0][1]
        
        for k in PointList:
            x = k[0]
            y = k[1]
            
            Xmax = max(x,Xmax)
            Xmin = min(x,Xmin)
            Ymax = max(y,Ymax)
            Ymin = min(y,Ymin)
            
        X = abs(Xmax - Xmin)
        Y = abs(Ymax - Ymin)
        Diameter = max(X, Y)
        
        return Diameter
    
    elif len(PointList[0]) == 3:
        Xmax = PointList[0][0]
        Xmin = PointList[0][0]
        Ymax = PointList[0][1]
        Ymin = PointList[0][1]
        Zmax = PointList[0][2]
        Zmin = PointList[0][2]
        
        for k in PointList:
            x = k[0]
            y = k[1]
            z = k[2]
            
            Xmax = max(x,Xmax)
            Xmin = min(x,Xmin)
            Ymax = max(y,Ymax)
            Ymin = min(y,Ymin)
            Zmax = max(z,Zmax)
            Zmin = min(z,Zmin)
            
        X = abs(Xmax - Xmin)
        Y = abs(Ymax - Ymin)
        Z = abs(Zmax - Zmin)
        Diameter = max(X, Y)
        Diameter = max(Diameter, Z)
        
        return Diameter
    
    
#%%
def CentrePointList(ListPoint):
    ListCoordX = []
    ListCoordY = []
    appendX = ListCoordX.append
    appendY = ListCoordY.append
    
    for k in ListPoint:
        appendX(k[0])
        appendY(k[1])
    
    CentreX = (np.amax(ListCoordX) + np.amin(ListCoordX)) * 0.5
    CentreY = (np.amax(ListCoordY) + np.amin(ListCoordY)) * 0.5
    
    ListPointCentre = []
    append = ListPointCentre.append
    for k in ListPoint:
        x = k[0] - CentreX
        y = k[1] - CentreY
        append(np.array([x,y]))
        
    return ListPointCentre

#%%
def ReturnNumericalAperture(RayList, RefractiveIndex):
    
    CentralRay = FindCentralRay(RayList)
    if CentralRay is None:
        CentralVector = np.array([0,0,0])
        for k in RayList:
            CentralVector = CentralVector + k.vector
        CentralVector = CentralVector/len(RayList)  
    else:    
        CentralVector = CentralRay.vector
    ListAngleAperture = []
    for k in RayList:
        ListAngleAperture.append(mgeo.AngleBetweenTwoVectors(CentralVector, k.vector))

    return np.sin(np.amax(ListAngleAperture)) * RefractiveIndex

#%%
def ReturnAiryRadius(Wavelength, NumericalAperture):
    if NumericalAperture > 1e-3 and Wavelength is not None: 
        return 1.22*0.5*Wavelength/NumericalAperture
    else:
        return 0  #for very small numerical apertures, diffraction effects becomes negligible and the Airy Radius becomes meaningless
    


#%%
def FindOptimalDistanceBIS(movingDetector, Amplitude, Step, RayList, OptFor, IntensityWeighted):
    ListSizeSpot = []
    ListDuration = []
    ListFitness = []
    if IntensityWeighted:
        Weights = [k.intensity for k in RayList]

    movingDetector.shiftByDistance(-Amplitude)
    n = int(2*Amplitude / Step)
    for i in range(n):
        ListPointDetector2DCentre = movingDetector.get_PointList2DCentre(RayList)
        if OptFor in ["intensity", "spotsize"]:
            if IntensityWeighted:
                SpotSize = WeightedStandardDeviation(ListPointDetector2DCentre,Weights)
            else:
                SpotSize = StandardDeviation(ListPointDetector2DCentre)       
            ListSizeSpot.append(SpotSize)
        
        if OptFor in ["intensity", "duration"]:
            DelayList = movingDetector.get_Delays(RayList)
            if IntensityWeighted:
                Duration = WeightedStandardDeviation(DelayList,Weights)
            else:
                Duration = StandardDeviation(DelayList)            
            ListDuration.append(Duration)
        
        if OptFor == "intensity":
            Fitness = SpotSize**2*Duration 
        elif OptFor == "duration":
            Fitness = Duration 
        elif OptFor == "spotsize":
            Fitness = SpotSize
        ListFitness.append(Fitness)
        
        movingDetector.shiftByDistance(Step)
        
    FitnessMin = min(ListFitness)
    ind = ListFitness.index(FitnessMin)
    if OptFor in ["intensity", "spotsize"]:
        OptSizeSpot = ListSizeSpot[ind]
    else: 
        OptSizeSpot = np.nan
    if OptFor in ["intensity", "duration"]:    
        OptDuration = ListDuration[ind]
    else:
        OptDuration = np.nan
        
    movingDetector.shiftByDistance(-(n-ind)*Step)

    return movingDetector, OptSizeSpot, OptDuration


def FindOptimalDistance(Precision, Detector, RayList, OptFor="intensity", Amplitude=None, IntensityWeighted=False, verbose=False):
    if OptFor not in ["intensity", "size", "duration"]:
        raise NameError("I don`t recognize what you want to optimize the detector distance for. OptFor must be either 'intensity', 'size' or 'duration'.")    
        
    FirstDistance = Detector.get_distance()
    ListPointDetector2DCentre = Detector.get_PointList2DCentre(RayList)
    SizeSpot = 2*StandardDeviation(ListPointDetector2DCentre)
    NumericalAperture = ReturnNumericalAperture(RayList, 1)
    if Amplitude is None:
        Amplitude = min(4*np.ceil(SizeSpot / np.tan(np.arcsin(NumericalAperture))), FirstDistance)
    #Step = Amplitude/10
    Step = Amplitude/5
    
    if verbose: print(f"Searching optimal detector position for *{OptFor}* within [{FirstDistance-Amplitude:.3f}, {FirstDistance+Amplitude:.3f}] mm...", end='', flush=True)
    movingDetector = Detector.copy_detector()
    for k in range(Precision+1):
        movingDetector, OptSizeSpot, OptDuration = FindOptimalDistanceBIS(movingDetector, Amplitude*0.1**k, Step*0.1**k, RayList, OptFor, IntensityWeighted)
    
    if not FirstDistance -Amplitude+10**-Precision < movingDetector.get_distance() < FirstDistance+Amplitude-10**-Precision:
        print("There`s no minimum-size/duration focus in the searched range.") 

    print('\r\033[K', end='', flush=True) #move to beginning of the line with \r and then delete the whole line with \033[K            
    return movingDetector, OptSizeSpot, OptDuration



#%%
def RayTracingCalculation(source_rays, optical_elements):
    
    ez = np.array([0,0,1])
    ex = np.array([1,0,0])
    RayListHistory = []
    
    for k in range(0,len(optical_elements)):
        
        if k==0:  RayList = source_rays
        else:     RayList = RayListHistory[k-1]
    
        #Mirror = optical_elements[k].type
        Position = optical_elements[k].position
        n = optical_elements[k].normal
        m = optical_elements[k].majoraxis
        
        #transform rays into mirror coordinate system
        RayList = mgeo.TranslationRayList(RayList, -Position)
        RayList = mgeo.RotationRayList(RayList, n, ez)
            #FIRST rotate m in the same way as you just did the ray list,
            #THEN rotate the rays again to match the NEW m with the x-axis !
        mPrime = mgeo.RotationPoint(m, n, ez) 
        RayList = mgeo.RotationRayList(RayList, mPrime, ex)
        RayList = mgeo.TranslationRayList(RayList, optical_elements[k].type.get_centre())
        
        #optical element acts on the rays:
        if 'Mirror' in optical_elements[k].type.type:
            RayList = mmirror.ReflectionMirrorRayList(optical_elements[k].type, RayList)
        elif optical_elements[k].type.type == 'Mask':
            RayList = mmask.TransmitMaskRayList(optical_elements[k].type, RayList) 
        else: raise NameError('I don`t recognize the type of optical element ' + optical_elements[k].type.type + '.')
        
        #transform new rays back into "lab frame"
        RayList = mgeo.TranslationRayList(RayList, -optical_elements[k].type.get_centre())
        RayList = mgeo.RotationRayList(RayList, ex, mPrime)
        RayList = mgeo.RotationRayList(RayList, ez, n)
        RayList = mgeo.TranslationRayList(RayList, Position)
        
        RayListHistory.append(RayList)
        
    return RayListHistory



#%%
def singleOEPlacement(SourceProperties: dict, OpticsList: list, DistanceList: list[(int, float)], IncidenceAngleList: list[(int, float)], IncidencePlaneAngleList: list[(int, float)], Description: str):
    Divergence = SourceProperties["Divergence"]
    SourceSize = SourceProperties["SourceSize"]
    RayNumber = SourceProperties["NumberRays"]
    Wavelength = SourceProperties["Wavelength"]

    IncidencePlaneAngleList = [np.deg2rad((i%360)) for i in IncidencePlaneAngleList] # wrap angles into [0,360]deg and convert to radian
    IncidenceAngleList = [np.deg2rad((i%360)) for i in IncidenceAngleList]  # wrap angles into [0,360]deg and convert to radian
    
    # Launch source ray bundle from the origin [x,y,z]=[0,0,0] into the x-direction:
    SourcePosition = np.array([0,0,0])
    SourceDirection = np.array([1,0,0])
    if Divergence == 0:
        if SourceSize ==0:
            try:
                radius = 0.5*min(OpticsList[0].support.dimX,OpticsList[0].support.dimY) # for rect. support
            except:
                radius = OpticsList[0].support.radius # otherwise it must be a round support
        else: 
            radius = SourceSize/2;
        SourceRayList = msource.PlaneWaveDisk(SourcePosition, SourceDirection, radius, RayNumber, Wavelength = Wavelength)
    else:
        if SourceSize ==0:
            SourceRayList = msource.PointSource(SourcePosition, SourceDirection, Divergence, RayNumber, Wavelength = Wavelength)
        else: 
            SourceRayList = msource.ExtendedSource(SourcePosition, SourceDirection, SourceSize, Divergence, RayNumber, Wavelength = Wavelength)
    
    SourceRayList = msource.ApplyGaussianIntensityToRayList(SourceRayList, 1/np.e**2) #Intensity drops to 1/e^2 at edge of ray-bundle
    
    # Now successively create and align optical elements 
    SourceRay = [mray.Ray(SourcePosition, SourceDirection)] # a single source ray in the central direction of the ray-bundle created above
    OpticalElements = []
    OpticalElementCentre = SourcePosition
    CentralVector = SourceDirection 
    RotationAxis = np.array([0,1,0])  #Vector perpendicular to the incidence plane, i.e. initially the incidence plane is the x-z-plane.
    
    for k, Optic in enumerate(OpticsList):
        # for convex mirrors, rotated them by 180° them so we reflect from the "back side"
        if Optic.type == 'SphericalCX Mirror' or Optic.type == 'CylindricalCX Mirror':
            IncidenceAngleList[k] = np.pi - IncidenceAngleList[k] 
        
        # shift OpticalElementCentre-point from that of the preceding OE, by the required distance along the central ray vector
        OpticalElementCentre = CentralVector * DistanceList[k] + OpticalElementCentre
                 
        if abs(IncidencePlaneAngleList[k]-np.pi) < 1e-10:
            RotationAxis = -RotationAxis
        else:
            RotationAxis = mgeo.RotationAroundAxis(CentralVector, -IncidencePlaneAngleList[k], RotationAxis)
                
        OpticalElementNormal = mgeo.RotationAroundAxis(RotationAxis, -np.pi/2 + IncidenceAngleList[k], np.cross(CentralVector,RotationAxis))
        
        OpticalElementMajorAxis = np.cross(RotationAxis,OpticalElementNormal)
                    
        Element = moe.OpticalElement(Optic, OpticalElementCentre, OpticalElementNormal, OpticalElementMajorAxis)
        
        OpticalElements.append(Element)
        auxChain = moc.OpticalChain(SourceRay, OpticalElements)
        
        if 'Mirror' in Optic.type:
            OutRays = auxChain.get_output_rays()
            CentralVector =  OutRays[-1][0].vector
        elif Optic.type == 'Mask' :
            # CentralVector remains unchanged
            # in our auxChain (and *not* in the OpticalElements-list), replace mask by completely tansparent "fake version"
            # to be sure that our CentralVector serving as alignment-guide always passes
            FakeMask = mmask.Mask(msupp.SupportRoundHole(Radius=100, RadiusHole=100, CenterHoleX=0, CenterHoleY=0))
            auxChain.optical_elements[-1] = moe.OpticalElement(FakeMask, OpticalElementCentre, OpticalElementNormal,OpticalElementMajorAxis)     
        else: 
            raise NameError('I don`t recognize the type of optical element ' + OpticalElements[k].type.type + '.')
 
    return moc.OpticalChain(SourceRayList,OpticalElements,Description)


def OEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList=None, Description = '', render: bool =False):
    """ Automatically place optical elements in the "lab frame" according to given distances and incidence angles. Outputs an OpticalChain-object.
    The source is placed at the origin, and points into the x-direction. The optical elements then follow.
    As long as the angles in IncidencePlaneAngleList are 0, the incidence plane remains the x-z plane, i.e. the optical elements are rotated about
    the y-axis to set the desied incidence angle between OE-normal and the direction of incidence of the incoming ray-bundle. 
    Otherwise, the incidence plane gets rotated.
    One of the elements of one of the lists 'DistanceList', 'IncidenceAngleList', or 'IncidencePlaneAngleList' can be a list of numpy-array.
    In that case, a list of OpticalChain-objects is created which can be looped over.
    """    
    if render: from mayavi import mlab
    
    if IncidencePlaneAngleList is None: IncidencePlaneAngleList = np.zeros(len(OpticsList))
    
    nest_indx_distance  = which_indexes(DistanceList)
    nest_indx_incidence = which_indexes(IncidenceAngleList)
    nest_indx_incplane  = which_indexes(IncidencePlaneAngleList)
    
    total_nested = len(nest_indx_incidence + nest_indx_incplane +nest_indx_distance) 
    if total_nested >1:
        raise ValueError('Only one element of one of the lists IncidenceAngleList, IncidencePlaneAngleList, or DistanceList can be a list or array itself. Otherwise things get too tangled...')
    elif total_nested ==1:
        i = (nest_indx_incidence + nest_indx_incplane +nest_indx_distance)[0]
        loop_variable_name = OpticsList[i].type + '_idx_' + str(i)
        if nest_indx_incidence: #if the list is not empty
            loop_variable_name += " incidence angle (deg)"
            loop_variable =   copy.deepcopy(IncidenceAngleList[i])
            loop_list = IncidenceAngleList
        elif nest_indx_distance:
            loop_variable_name += " distance (mm)"
            loop_variable =  copy.deepcopy(DistanceList[i])
            loop_list = DistanceList
        elif nest_indx_incplane:
            loop_variable_name += " incidence-plane angle rotation (deg)"
            loop_variable =   copy.deepcopy(IncidencePlaneAngleList[i])
            loop_list = IncidencePlaneAngleList
   
        OpticalChainList = []
        for x in loop_variable:
            loop_list[i] = x
            ModifiedOpticalChain = singleOEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList, Description)
            ModifiedOpticalChain.loop_variable_name = loop_variable_name
            ModifiedOpticalChain.loop_variable_value = x
            OpticalChainList.append(ModifiedOpticalChain)
            
            if render:
                #to produce and save renderings of the setup for each iteration in the loop:    
                mlab.options.offscreen = True
                fig = ModifiedOpticalChain.render()
                name = loop_variable_name + "=" + "{:.2f}".format(x) +".png"
                name = name.replace(" ", "_")
                print('...saving image...', end='', flush=True) 
                mlab.savefig(name, figure=fig, magnification=3)
                print('\r\033[K', end='', flush=True) #move to beginning of the line with \r and then delete the whole line with \033[K
            
        if render: mlab.options.offscreen = False
        
        return OpticalChainList
    
    elif total_nested ==0:
        OpticalChain = singleOEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList, Description) #all simple
        if render: 
            mlab.options.offscreen = True
            fig = OpticalChain.render()
            name = "OpticalChain_" + datetime.now().strftime("%Y-%m-%d-%Hh%M") + ".png"
            mlab.savefig(name, figure=fig, magnification=3)
            mlab.options.offscreen = False
        
        return OpticalChain

