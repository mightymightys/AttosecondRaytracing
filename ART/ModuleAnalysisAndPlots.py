#%% Data-analysis and plotting modules
"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import ART.ModuleProcessing as mp
import ART.ModuleGeometry as mgeo

#%%
def getDetectorPoints(RayListAnalysed, Detector):   
    DetectorPointList2DCentre = Detector.get_PointList2DCentre(RayListAnalysed)
    DectectorPoint2D_Xcoord = np.array([k[0]*1e3 for k in DetectorPointList2DCentre]) # mm to µm
    DectectorPoint2D_Ycoord = np.array([k[1]*1e3 for k in DetectorPointList2DCentre])

    FocalSpotSize = mp.DiameterPointList(DetectorPointList2DCentre)
    SpotSizeSD = mp.StandardDeviation(DetectorPointList2DCentre)
    return DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, FocalSpotSize, SpotSizeSD


#%%
def getETransmission(RayListIn, RayListOut):
    """ Calculates the energy transmission from RayListIn to RayListOut in percent. """
    ETransmission = 100 * sum(Ray.intensity for Ray in RayListOut) / sum(Ray.intensity for Ray in RayListIn)
    return ETransmission

#%%
def GetResultSummary(Detector, RayListAnalysed, verbose=False):
    """Calculate and return various results for the given detector and RayList."""
    DetectorPointList2DCentre = Detector.get_PointList2DCentre(RayListAnalysed)
    FocalSpotSizeSD = mp.StandardDeviation(DetectorPointList2DCentre)
    DelayList = Detector.get_Delays(RayListAnalysed)
    DurationSD = mp.StandardDeviation(DelayList)

    if verbose:
        FocalSpotSize = mp.DiameterPointList(DetectorPointList2DCentre)
        summarystring = "At the detector distance of " + '{:.3f}'.format(Detector.get_distance()) +" mm we get:\n" \
                        + "Spatial std : " + '{:.3f}'.format(FocalSpotSizeSD*1e3) + " \u03BCm and min-max: " + '{:.3f}'.format(FocalSpotSize*1e3) +" \u03BCm\n" \
                        + "Temporal std : " + '{:.3e}'.format(DurationSD) + " fs and min-max : " + '{:.3e}'.format(max(DelayList)-min(DelayList)) + " fs"

        print(summarystring)
    
    return FocalSpotSizeSD, DurationSD

#%%
def SpotDiagram(RayListAnalysed, Detector, Wavelength, DrawAiryAndFourier, ColorCoded=None):  
    NumericalAperture = mp.ReturnNumericalAperture(RayListAnalysed, 1) #NA determined from final ray bundle
    if DrawAiryAndFourier : 
        AiryRadius = mp.ReturnAiryRadius(Wavelength, NumericalAperture)*1e3 # in µm
    else:
        AiryRadius = 0
   
    
    DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, FocalSpotSize, SpotSizeSD = getDetectorPoints(RayListAnalysed, Detector)
    
    if ColorCoded=='Intensity':
        IntensityList = [k.intensity for k in RayListAnalysed]
        z = np.asarray(IntensityList)
        zlabel = "Intensity (arb.u.)"
        title = "Intensity + Spot Diagram\n press left/right to move detector position"
        addLine=''
    elif ColorCoded=='Incidence':
        IncidenceList = [np.rad2deg(k.incidence) for k in RayListAnalysed] # degree
        z = np.asarray(IncidenceList)
        zlabel = "Incidence angle (deg)"
        title = "Ray Incidence + Spot Diagram\n press left/right to move detector position"
        addLine = ''
    elif ColorCoded=='Delay':
        DelayList = Detector.get_Delays(RayListAnalysed)
        DurationSD = mp.StandardDeviation(DelayList)
        z = np.asarray(DelayList)
        zlabel = "Delay (fs)"
        title = "Delay + Spot Diagram\n press left/right to move detector position"
        addLine = '\n'+'{:.2f}'.format(DurationSD)+' fs SD'
    else:
        z = 'red'
        title = "Spot Diagram\n press left/right to move detector position"
        addLine = ''
    
    
    Dist = Detector.get_distance()
    distStep = min(50, max(0.0005, round(FocalSpotSize/8/np.arcsin(NumericalAperture)*10000)/10000)) #in mm
    
    plt.ion()   
    fig, ax = plt.subplots()
    if DrawAiryAndFourier :
        theta = np.linspace(0, 2*np.pi, 100)
        x = AiryRadius*np.cos(theta)
        y = AiryRadius*np.sin(theta) #
        ax.plot(x,y, c="black")
    
    foo = ax.scatter(DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, c=z, s=15, label='{:.3f}'.format(Dist)+' mm\n'+'{:.1f}'.format(SpotSizeSD*1e3)+' \u03BCm SD'+addLine)

    
    axisLim = 1.1*max(AiryRadius, 0.5*FocalSpotSize*1000)
    ax.set_xlim(-axisLim, axisLim)
    ax.set_ylim(-axisLim, axisLim)
    
    if ColorCoded=='Intensity' or ColorCoded=='Incidence' or ColorCoded=='Delay':
        cbar = fig.colorbar(foo)
        cbar.set_label(zlabel) 
    
    ax.legend(loc='upper right')
    ax.set_xlabel('X (µm)')
    ax.set_ylabel('Y (µm)')
    ax.set_title(title)
    #ax.margins(x=0)

    movingDetector = Detector.copy_detector()
    def press(event):
        nonlocal Dist, distStep, movingDetector, ColorCoded, zlabel, cbar
        if event.key == 'right':   
            movingDetector.shiftByDistance(distStep)
            Dist +=  distStep
        elif event.key == 'left':  
            if Dist > 1.5*distStep:
                movingDetector.shiftByDistance(-distStep)
                Dist -=  distStep
            else:
                movingDetector.shiftToDistance(0.5*distStep)
                Dist = 0.5*distStep
                

        newDectectorPoint2D_Xcoord, newDectectorPoint2D_Ycoord, newFocalSpotSize, newSpotSizeSD = getDetectorPoints(RayListAnalysed, movingDetector)    
        xy = foo.get_offsets()
        xy[:,0] = newDectectorPoint2D_Xcoord
        xy[:,1] = newDectectorPoint2D_Ycoord
        foo.set_offsets(xy)
        
        if ColorCoded=='Delay':
           newDelayList =  np.asarray(movingDetector.get_Delays(RayListAnalysed))
           newDurationSD = mp.StandardDeviation(newDelayList)
           newaddLine = '\n'+'{:.2f}'.format(newDurationSD)+' fs SD'
           foo.set_array(newDelayList)
           foo.set_clim(min(newDelayList),max(newDelayList))
           cbar.update_normal(foo)
            
        foo.set_label('{:.3f}'.format(Dist)+ ' mm\n'+'{:.1f}'.format(newSpotSizeSD*1e3)+' \u03BCm SD'+newaddLine)
        ax.legend(loc='upper right')
        
        axisLim = 1.1*max(AiryRadius, 0.5*newFocalSpotSize*1000)
        ax.set_xlim(-axisLim, axisLim)
        ax.set_ylim(-axisLim, axisLim)
        
        distStep =  min(50,max(0.0005, round(newFocalSpotSize/8/np.arcsin(NumericalAperture)*10000)/10000)) #in mm
        
        fig.canvas.draw_idle() 
    
    fig.canvas.mpl_connect('key_press_event', press)
 
    plt.show()
    
    return fig
    
#%%
def drawDelayGraph(RayListAnalysed, Detector, Distance, Wavelength, DeltaFT, DrawAiryAndFourier, ColorCoded=None, fig=None):
    NumericalAperture = mp.ReturnNumericalAperture(RayListAnalysed, 1) #NA determined from final ray bundle
    AiryRadius = mp.ReturnAiryRadius(Wavelength, NumericalAperture)*1e3 # in µm

    DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, FocalSpotSize, SpotSizeSD = getDetectorPoints(RayListAnalysed, Detector)
    
    DelayList = Detector.get_Delays(RayListAnalysed)
    DurationSD = mp.StandardDeviation(DelayList)
    
    if ColorCoded=='Intensity':
        IntensityList = [k.intensity for k in RayListAnalysed]
    elif ColorCoded=='Incidence':
        IncidenceList = [np.rad2deg(k.incidence) for k in RayListAnalysed] # in degree

    
    plt.ion()   
    if fig is None:
        fig = plt.figure()
    else:
        fig.clear(keep_observers=True)

    ax = Axes3D(fig)
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_zlabel("Delay (fs)")
    
    labelLine = '{:.3f}'.format(Distance)+ ' mm\n'+'{:.1f}'.format(SpotSizeSD*1e3)+' \u03BCm SD\n'+'{:.2f}'.format(DurationSD)+' fs SD'
    if ColorCoded=='Intensity':
        ax.scatter(DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord,DelayList,s=4,c=IntensityList,label=labelLine)
        #plt.title('Delay + Intensity graph\n press left/right to move detector position')
        ax.set_title('Delay + Intensity graph\n press left/right to move detector position')

    elif ColorCoded=='Incidence':
        ax.scatter(DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord,DelayList,s=4,c=IncidenceList,label=labelLine)
        #plt.title('Delay + Incidence graph\n press left/right to move detector position')
        ax.set_title('Delay + Incidence graph\n press left/right to move detector position')
    else: 
        ax.scatter(DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord,DelayList,s=4,c=DelayList,label=labelLine)
        #plt.title('Delay graph\n press left/right to move detector position')
        ax.set_title('Delay graph\n press left/right to move detector position')
        
    ax.legend(loc='upper right')
    
    if DrawAiryAndFourier:
        x = np.linspace(-AiryRadius,AiryRadius,40)
        z = np.linspace(np.mean(DelayList) - DeltaFT*0.5,np.mean(DelayList) + DeltaFT*0.5,40)
        x, z = np.meshgrid(x,z)
        y = np.sqrt(AiryRadius**2 - x**2)
        ax.plot_wireframe(x,y,z,color="grey",alpha=0.1)
        ax.plot_wireframe(x,-y,z,color="grey",alpha=0.1)
    
    axisLim = 1.1*max(AiryRadius, 0.5*FocalSpotSize*1000)
    ax.set_xlim(-axisLim, axisLim)
    ax.set_ylim(-axisLim, axisLim)
    
    plt.show()
    fig.canvas.draw()
    
    return fig, NumericalAperture, AiryRadius, FocalSpotSize

#%%
def DelayGraph(RayListAnalysed, Detector, Wavelength, DeltaFT, DrawAiryAndFourier=False, ColorCoded=None):
    Dist = Detector.get_distance()
    fig, NumericalAperture, AiryRadius, FocalSpotSize = drawDelayGraph(RayListAnalysed, Detector, Dist, Wavelength, DeltaFT, DrawAiryAndFourier, ColorCoded)
        
    distStep =  min(50,max(0.0005, round(FocalSpotSize/8/np.arcsin(NumericalAperture)*10000)/10000)) #in mm

    movingDetector = Detector.copy_detector()
    def press(event):
        nonlocal Dist, distStep, movingDetector, fig
        if event.key == 'right':
            movingDetector.shiftByDistance(distStep)
            Dist +=  distStep
            fig, sameNumericalAperture, sameAiryRadius, newFocalSpotSize = drawDelayGraph(RayListAnalysed, movingDetector, Dist, Wavelength, DeltaFT, DrawAiryAndFourier, ColorCoded, fig)
        elif event.key == 'left':
            if Dist > 1.5*distStep:
                movingDetector.shiftByDistance(-distStep)
                Dist -=  distStep
            else:
                 movingDetector.shiftToDistance(0.5*distStep)
                 Dist = 0.5*distStep
                
            fig, sameNumericalAperture, sameAiryRadius, newFocalSpotSize = drawDelayGraph(RayListAnalysed, movingDetector, Dist, Wavelength, DeltaFT, DrawAiryAndFourier, ColorCoded, fig)
        distStep =  min(50,max(0.0005, round(newFocalSpotSize/8/np.arcsin(NumericalAperture)*10000)/10000)) #in mm
            
    fig.canvas.mpl_connect('key_press_event', press)
    
    return fig

#%% 
def RayRenderGraph(OpticalChain,EndDistance=None,maxRays=150,OEpoints=2000):
    from mayavi import mlab
    RayListHistory = [OpticalChain.source_rays] + OpticalChain.get_output_rays()
    
    if EndDistance == None:
        EndDistance = np.linalg.norm(OpticalChain.source_rays[0].point - OpticalChain.optical_elements[0].position)
    
    print('...rendering image of optical chain...', end='', flush=True) 
    fig = mlab.figure(bgcolor=(1,1,1),size=(1500, 500))
    
    # Ray display
    for k in range(len(RayListHistory)):
        if k != len(RayListHistory)-1:
            knums = list(map(lambda x: x.number, RayListHistory[k])) #make a list of all ray numbers that are still in the game             
            if len(RayListHistory[k+1]) > maxRays:
                rays_to_render = np.random.choice(RayListHistory[k+1], maxRays, replace=False)
            else: rays_to_render = RayListHistory[k+1]

            for j in rays_to_render:
                indx = knums.index(j.number)
                i = RayListHistory[k][indx]
                Point1 = i.point
                Point2 = j.point
                x = np.asarray([Point1[0], Point2[0]])
                y = np.asarray([Point1[1], Point2[1]])
                z = np.asarray([Point1[2], Point2[2]])
                mlab.plot3d(x, y, z, color=(1,0,0),  tube_radius=0.05)

                
        else:
            if len(RayListHistory[k]) > maxRays:
                rays_to_render = np.random.choice(RayListHistory[k], maxRays, replace=False)
            else: rays_to_render = RayListHistory[k]

            for j in rays_to_render:
                Point = j.point
                Vector = j.vector
                x = np.asarray([Point[0], Point[0] + Vector[0]*EndDistance])
                y = np.asarray([Point[1], Point[1] + Vector[1]*EndDistance])
                z = np.asarray([Point[2], Point[2] + Vector[2]*EndDistance])
                mlab.plot3d(x, y, z, color=(1,0,0), tube_radius=0.05) 

    # Optics display
    for OE in OpticalChain.optical_elements:
        OpticPointList = OE.type.get_grid3D(OEpoints) #in the optic's coordinate system
        
        # transform OpticPointList into "lab-frame"
        OpticPointList = mgeo.TranslationPointList(OpticPointList, -OE.type.get_centre())
        MirrorMajorAxisPrime = mgeo.RotationPoint(OE.majoraxis, OE.normal, np.array([0,0,1]))
        OpticPointList = mgeo.RotationPointList(OpticPointList, np.array([1,0,0]), MirrorMajorAxisPrime)
        OpticPointList = mgeo.RotationPointList(OpticPointList, np.array([0,0,1]), OE.normal)
        OpticPointList = mgeo.TranslationPointList(OpticPointList, OE.position)        
        
        # plot 3D-drig of OE as little spheres
        x = np.asarray([i[0] - OE.normal[0]*1.7 for i in OpticPointList])
        y = np.asarray([i[1] - OE.normal[1]*1.7 for i in OpticPointList])
        z = np.asarray([i[2] - OE.normal[2]*1.7 for i in OpticPointList])

        mlab.points3d(x,y,z,mode='sphere',color=(0.66,0.66,0.66),scale_factor=3)
        
        # # or awkwardly transform the OPticPointList into 2D arrays for x,y,z, interpolate,
        # # and then plot a smooth surface ? Actually not faster at all, and not so much nicer either.
        # x_array, y_array, z_array = np.array([]), np.array([]), np.array([])
        # for element in OpticPointList:
        #     x, y, z = element
        #     x_array = np.hstack((x_array, x))
        #     y_array = np.hstack((y_array, y))
        #     z_array = np.hstack((z_array, z))

        # x_positions, y_positions = np.unique(x_array), np.unique(y_array)
        # x_grid, y_grid = np.meshgrid(x_positions, y_positions)

        # x_grid = np.transpose(x_grid)
        # y_grid = np.transpose(y_grid)

        # # Create an empty grid with the same shape as x_grid and y_grid
        # z_grid = np.zeros(x_grid.shape)

        # # Fill the z_grid with the corresponding z values
        # for x, y, z in zip(x_array, y_array, z_array):
        #     x_index = np.where(x_positions == x)[0][0]
        #     y_index = np.where(y_positions == y)[0][0]
        #     z_grid[x_index, y_index] = z

        # z_grid[np.where(z_grid == 0)] = None

        # from scipy.interpolate import griddata

        # # Creating the points and values arrays for the interpolation
        # foopoints = np.column_stack((x_grid[~np.isnan(z_grid)], y_grid[~np.isnan(z_grid)]))
        # foovalues = z_grid[~np.isnan(z_grid)]

        # # Interpolating the data
        # z_grid = griddata(foopoints, foovalues, (x_grid, y_grid), method='linear')
        
        # mlab.mesh(x_grid, y_grid, z_grid, color=(0.66,0.66,0.66))     

    fig.scene._lift()
    mlab.view(azimuth=-90, elevation=90, distance='auto')
    
    print('\r\033[K', end='', flush=True) #move to beginning of the line with \r and then delete the whole line with \033[K
    return fig


#%%
def MirrorProjection(RayListAnalysed, OpticalChain, ReflectionNumber, Detector=None, ColorCoded=None):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    Position = OpticalChain.optical_elements[ReflectionNumber].position
    n = OpticalChain.optical_elements[ReflectionNumber].normal
    m = OpticalChain.optical_elements[ReflectionNumber].majoraxis
    
    #transform rays into the mirror-support reference frame
    #(same as mirror frame but without the shift by mirror-centre)
    RayList = mgeo.TranslationRayList(RayListAnalysed, -Position)
    RayList = mgeo.RotationRayList(RayList, n,  np.array([0,0,1]))
    mPrime = mgeo.RotationPoint(m, n, np.array([0,0,1])) 
    RayList = mgeo.RotationRayList(RayList, mPrime, np.array([1,0,0]))
    
        
    x = np.asarray([k.point[0] for k in RayList])
    y = np.asarray([k.point[1] for k in RayList])
    if ColorCoded=='Intensity':
        IntensityList = [k.intensity for k in RayListAnalysed]
        z = np.asarray(IntensityList)
        zlabel = "Intensity (arb.u.)"
        title = "Ray intensity projected on mirror              "
    elif ColorCoded=='Incidence':
        IncidenceList = [np.rad2deg(k.incidence) for k in RayListAnalysed] # in degree
        z = np.asarray(IncidenceList)
        zlabel = "Incidence angle (deg)"
        title = "Ray incidence projected on mirror              "
    elif ColorCoded=='Delay':
        if Detector is not None:
            z = np.asarray(Detector.get_Delays(RayListAnalysed))
            zlabel = "Delay (fs)"
            title = "Ray delay at detector projected on mirror              "
        else: 
            raise ValueError('If you want to project ray delays, you must specify a detector.')
    else:
        z = 'red'
        title = "Ray impact points projected on mirror"
    
    plt.ion()        
    fig = plt.figure()
    ax = OpticalChain.optical_elements[ReflectionNumber].type.support.ContourSupport(fig)
    p = plt.scatter(x,y,c=z, s=15)
    if ColorCoded=='Delay' or ColorCoded=='Incidence' or ColorCoded=='Intensity':
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(p, cax=cax)
        cbar.set_label(zlabel) 
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    plt.title(title, loc ='right')
    plt.tight_layout()
    
    bbox = ax.get_position()
    bbox.set_points(bbox.get_points()-np.array([[0.01, 0],[0.01, 0]]))
    ax.set_position(bbox)
    plt.show()