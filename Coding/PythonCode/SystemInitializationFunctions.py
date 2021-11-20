import numpy as np
import itertools
from numpy.lib.function_base import append
import scipy.spatial

## This is the function to define the size of somas
def SomaInitialization(NumberofNeurons,RadiusIntervalofSoma):
    Radius = np.random.uniform(RadiusIntervalofSoma[0],RadiusIntervalofSoma[1], size=(NumberofNeurons,1))
    return Radius 

## This is the function to define the space locations of somas
def LocationInitialization(InitializationType,NumberofNeurons,SpaceLimit):
    if InitializationType==1:  ## This is the grid initialization
        LocationMatrix=np.asarray(list(itertools.product(np.arange(0.1*SpaceLimit[0],0.9*SpaceLimit[0],0.8*SpaceLimit[0]/np.around(np.power(NumberofNeurons,1/3)).astype(int)).reshape(np.around(np.power(NumberofNeurons,1/3)).astype(int),1),np.arange(0.1*SpaceLimit[1],0.9*SpaceLimit[1],0.8*SpaceLimit[1]/np.around(np.power(NumberofNeurons,1/3)).astype(int)).reshape(np.around(np.power(NumberofNeurons,1/3)).astype(int),1),np.arange(0.1*SpaceLimit[0],0.9*SpaceLimit[2],0.8*SpaceLimit[2]/np.around(np.power(NumberofNeurons,1/3)).astype(int)).reshape(np.around(np.power(NumberofNeurons,1/3)).astype(int),1))))
    else: ## This is the coordinate initialization
        MiniDistance=0.01*SpaceLimit[0]/np.around(np.power(NumberofNeurons,1/3)).astype(int)
        LocationMatrix=np.hstack((np.random.randint(0.1*SpaceLimit[0], 0.9*SpaceLimit[0], size=(1,1)), np.random.randint(0.1*SpaceLimit[1],0.9*SpaceLimit[1], size=(1,1)),np.random.randint(0.1*SpaceLimit[2],0.9*SpaceLimit[2], size=(1,1))))
        while np.size(LocationMatrix,0)<NumberofNeurons:
            RandomLocation=np.hstack((np.random.randint(0.1*SpaceLimit[0], 0.9*SpaceLimit[0], size=(1,1)), np.random.randint(0.1*SpaceLimit[1],0.9*SpaceLimit[1], size=(1,1)),np.random.randint(0.1*SpaceLimit[2],0.9*SpaceLimit[2], size=(1,1))))
            DistanceMatrix=scipy.spatial.distance.cdist(RandomLocation,LocationMatrix,metric='euclidean')
            if len(np.where(DistanceMatrix<MiniDistance))<1:
                LocationMatrix=np.append(LocationMatrix,RandomLocation,axis=0)
    return LocationMatrix
    
## This is the function to generate the space and somas (we will distinguish between the space inside somas and space outside somas as well)
def SomaandSpaceInitialization(LocationMatrix,Radius,GrainedN):
    InitializedSomaSphereCell = np.empty((np.size(LocationMatrix,0),3),dtype=object)
    RealRadiusCell = np.empty((np.size(LocationMatrix,0),1),dtype=object)
    NeighborCell=np.empty((np.size(LocationMatrix,0),GrainedN*GrainedN),dtype=object)
    for ID in range(np.size(LocationMatrix,0)): ## Here we traverse every neuron
        SomaCenter = LocationMatrix[ID,:] ##This is the center of soma
        SomaRadius = Radius[ID] ## This is the preset radius of the soma. The real radius is generated from this preset value by adding a perturbation 
        u = np.linspace(0, 2 * np.pi, GrainedN)
        v = np.linspace(0, np.pi, GrainedN)
        RealSomaRadiusX = SomaRadius*(0.9+0.2*np.random.rand(GrainedN,GrainedN)) # The perturbation on X direction is in an interval of [0,20%] 
        RealSomaRadiusY = SomaRadius*(0.9+0.2*np.random.rand(GrainedN,GrainedN)) # The perturbation on Y direction is in an interval of [0,20%]
        RealSomaRadiusZ = SomaRadius*(0.9+0.2*np.random.rand(GrainedN,GrainedN)) # The perturbation on Z direction is in an interval of [0,20%]
        InitializedSomaSphereCell[ID,0] = RealSomaRadiusX*np.outer(np.cos(u), np.sin(v)) + SomaCenter[0] # This is the X coordinate
        InitializedSomaSphereCell[ID,1] = RealSomaRadiusY*np.outer(np.sin(u), np.sin(v)) + SomaCenter[1] # This is the Y coordinate
        InitializedSomaSphereCell[ID,2] = RealSomaRadiusZ*np.outer(np.ones(np.size(u)), np.cos(v)) + SomaCenter[2] # This is the Z coordinate
        # The cell named as InitializedSomaSphereCell will contain the coordinates of soma spheres
        XCoordinate=RealSomaRadiusX*np.outer(np.cos(u), np.sin(v)) + SomaCenter[0]
        YCoordinate=RealSomaRadiusY*np.outer(np.sin(u), np.sin(v)) + SomaCenter[1]
        ZCoordinate=RealSomaRadiusZ*np.outer(np.ones(np.size(u)), np.cos(v)) + SomaCenter[2]
        XYZCoordinate=np.append(np.append(np.reshape(XCoordinate,(-1,1)),np.reshape(YCoordinate,(-1,1)),axis=1),np.reshape(ZCoordinate,(-1,1)),axis=1)
        RealRadiusCell[ID,0] =scipy.spatial.distance.cdist([SomaCenter],XYZCoordinate,metric='euclidean')
        # The cell named as RealRadiusCell will contain the real radius values (the distance between each coordinate on the sphere and the soma)
        DistanceBetweenCoordinates=scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(XYZCoordinate,metric='euclidean'))
        DistanceBetweenCoordinates[np.where(DistanceBetweenCoordinates==0)]=np.max(DistanceBetweenCoordinates)
        for IDC in range(np.size(NeighborCell,1)):
            SortedDistanceVector=np.sort(np.unique(DistanceBetweenCoordinates[IDC,:]))
            FoundNeighbors=np.where((DistanceBetweenCoordinates[IDC,:]==SortedDistanceVector[0]) | (DistanceBetweenCoordinates[IDC,:]==SortedDistanceVector[1]))
            Col,Raw=np.array(np.unravel_index(FoundNeighbors, (GrainedN,GrainedN)))
            NeighborCell[ID,IDC]=np.array([Col.T,Raw.T])
    return InitializedSomaSphereCell,RealRadiusCell,NeighborCell

# This is the function to initialize the type of neurons
def NeuralTypeInitialization(NumberofNeurons,NeuralTypePar):
    if NeuralTypePar==1:
        NeuralTypeVector=np.ones((NumberofNeurons,1)) ## All neurons have excitable membranes
    elif NeuralTypePar==2:
        NeuralTypeVector=np.ones((NumberofNeurons,1))*2 ## All neurons have passive membranes
    elif NeuralTypePar==3:
        NeuralTypeVector=np.random.random_integers(1,2, size=(NumberofNeurons,1)) ## The types of neurons are randomized 
    return NeuralTypeVector 

## This is the function to define the space locations of Netrin-1 sources (if there is any)
def NetrinOneFieldInitialization(SpaceLimit,LocationMatrix,NetrinFieldType):
    if NetrinFieldType==1:  ## The Netrin Field is homogeneous everywhere
        SourceLocationMatrix=np.array([]) ## There is no source in the space
    elif NetrinFieldType==2: ## The Netrin Field is non-homogeneous, there is at least one source
        NumofNetrinSource=1 # Number of netrin source
        MiniDistanceBN=100 # This is the minimum distance between any source and neurons
        Min=np.zeros((3,1)) # This is a vector to store the minimum X,Y,Z coordinates
        Max=np.zeros((3,1)) # This is a vector to store the maximum X,Y,Z coordinates
        for IDDim in range(3):
            Min[IDDim]=0.1*SpaceLimit[IDDim]
            Max[IDDim]=0.9*SpaceLimit[IDDim]
            SourceLocationMatrix=np.array([Min[0]+(Max[0]-Min[0])*np.random.rand(1,1),Min[1]+(Max[1]-Min[1])*np.random.rand(1,1),Min[2]+(Max[2]-Min[2])*np.random.rand(1,1)],dtype=np.float16) # Initialize one neuron's coordinates
            while np.size(SourceLocationMatrix,0)<NumofNetrinSource: # while there exist sources without coordinates
                NCoord=np.array([Min[0]+(Max[0]-Min[0])*np.random.rand(1,1),Min[1]+(Max[1]-Min[1])*np.random.rand(1,1),Min[2]+(Max[2]-Min[2])*np.random.rand(1,1)]).reshape(1,3) # Generate one coordinate
                DistanceM=scipy.spatial.distance.cdist(NCoord,LocationMatrix,metric='euclidean') # Calculate the distance from the potential source to somas
                if len(np.where(DistanceM<MiniDistanceBN))<1: 
                    SourceLocationMatrix=np.append(SourceLocationMatrix,np.array(NCoord,dtype=np.float32),axis=0)
    return SourceLocationMatrix

## This function initializes the chemical information on neurons' somas
def ChemicalInitialization(NumberofNeurons,InitializedSomaSphereCell,GrainedN):
    # This function defines the parameters used in chemical systems
    CaOut=1 ## the external calcium concentration fixed at 1 mM
    CaInInitialized=np.array([30,150])/np.power(10,6) ## The internal calcium concentration is initialized at an interval of [50,150] nM, where 1 nM = 10^(-6) mM
    TubulinInitialized=np.array([0,10])/np.power(10,6) ## The internal Tubulin concentration is initialized at an interval of [0,10] nM, where 1 nM = 10^(-6) mM
    UnboundMAP2Initialized=np.array([0,10])/np.power(10,6) ## The internal UMAP-2 concentration is initialized at an interval of [0,10] nM, where 1 nM = 10^(-6) mM
    SubmembraneCell = np.empty((NumberofNeurons,9),dtype=object)
    for ID in range(NumberofNeurons): 
        SubmembraneCell[ID,0] = InitializedSomaSphereCell[ID,0] # This is the X coordinate
        SubmembraneCell[ID,1] = InitializedSomaSphereCell[ID,1] # This is the Y coordinate
        SubmembraneCell[ID,2] = InitializedSomaSphereCell[ID,2] # This is the Z coordinate
        SubmembraneCell[ID,3] = np.array(CaOut*np.ones((GrainedN,GrainedN))) # This is the external calcium concentration outside of a point on sphere
        SubmembraneCell[ID,4] = np.array(CaInInitialized[0]+(CaInInitialized[1]-CaInInitialized[0])*np.random.random_sample((GrainedN,GrainedN))) # This is the internal calcium concentration
        SubmembraneCell[ID,5] = np.array(TubulinInitialized[0]+(TubulinInitialized[1]-TubulinInitialized[0])*np.random.random_sample((GrainedN,GrainedN))) # This is the internal Tubulin concentration
        SubmembraneCell[ID,6] = np.array(UnboundMAP2Initialized[0]+(UnboundMAP2Initialized[1]-UnboundMAP2Initialized[0])*np.random.random_sample((GrainedN,GrainedN))) # This is the internal UMAP-2 concentration
    return SubmembraneCell

# This is the intitialization of development information. Meanwhile, it will determines the data structure
def DevelopInitialization(NumberofNeurons,GrainedN,DurationLength):
    DevelopmentInfoCell = np.empty((NumberofNeurons,5),dtype=object)
    for ID in range(NumberofNeurons): ## Here we traverse every neuron
        DevelopmentInfoCell[ID,0]=np.zeros((GrainedN,GrainedN)) # This is the growth probability information
        DevelopmentInfoCell[ID,1]=np.zeros((GrainedN,GrainedN)) # This is the randomized growth information
        DevelopmentInfoCell[ID,2]=np.empty((GrainedN,GrainedN),dtype=object) # This is the synapes information
        DevelopmentInfoCell[ID,3]=np.empty((GrainedN,GrainedN),dtype=object) # This is the information of growth history
        DevelopmentInfoCell[ID,4]=np.empty((GrainedN,GrainedN),dtype=object) # This is the cache information
        for IDC in range(GrainedN):
            for IDC2 in range(GrainedN):
                DevelopmentInfoCell[ID,2][IDC,IDC2]=np.empty((1,14),dtype=object).reshape(1,-1) # Initialization
                DevelopmentInfoCell[ID,4][IDC,IDC2]=np.empty((0,2),dtype=object) # Initialization
                    # In the above cell, the 1st column saves the coordinates of
       # synaptic components; The 2nd column saves the calcium on synaptic
       # components; The 3rd column saves the Tubulin concentration on
       # synaptic components; The 4th column saves UMAP-2 concentration on
       # synaptic components; The 5th column saves BMAP-2 concentration on
       # synaptic components; The 6th column saves PMAP-2 concentration on
       # synaptic components; The 7th column saves the elongation length on
       # synaptic components; The 8th column saves the branching
       # probability; The 9th column saves the radius of synapse;
       # The 10th column saves the segment type (1 stands for root segment, 
       # 2 stands for intermediate segment, 3 stands for terminal
       # segment); The 11th column saves the previous segment; The
       # 12th column saves the next segment (or segments if there
       # is branching); The 13th column saves the centrifugal
       # order; The 14thcolumn saves temporary information
                for IDC3 in range(1,9):
                    DevelopmentInfoCell[ID,2][IDC,IDC2][0,IDC3]=np.array([0])
                DevelopmentInfoCell[ID,2][IDC,IDC2][0,12]=np.array([0])
                DevelopmentInfoCell[ID,2][IDC,IDC2][0,13]=np.array([1])
                DevelopmentInfoCell[ID,3][IDC,IDC2]=np.zeros((DurationLength,4))
                # the 1th column saves the average elongation rate, the 2nd
                # column saves the variance of elongation rate, the 3rd
                # column saves the number of terminal segments; %% the 4th
                # column saves the normalization base of centrifugal order;
                # Here DevelopmentInfoCell{ID,5} is used to save some
                # information to calculate DevelopmentInfoCell{ID,4}
    return DevelopmentInfoCell