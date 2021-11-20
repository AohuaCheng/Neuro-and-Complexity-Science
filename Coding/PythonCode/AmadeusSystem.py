## This is the main function of Amadeus system. Please ensure that all functions are implemented in the same path.

## Function Import
import scipy.io as io
import numpy as np
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
import sys
import SystemInitializationFunctions
import ChemicalSystemFunctions
import MorphologySystemFunctions
import Visualization

## Parameter Setting
DurationLength=1000 # Here DurationLength is measured in hours. E.g., DurationLength=1000 corresponds to 1000 hours
InitializationType = 2 # Here you can either choose 1 for grid initialization or 2 for coordinate initialization
NeuralTypePar = 1 # There are two types of neurons: neurons with excitable membranes and passive membranes respectively. 
# You can set 1 to let all neurons have excitable membranes or set 2 to let all of them have passive membranes. You can also set 3 to randomly generate membranes.
NumberofNeurons = 1**3 # Here you set the number of neurons, if you choose InitializationType = 1, then NumberofNeurons must be n^3, where n is a natural number
SpaceLimit = [1000,1000,1000] # Here you set the spacial scale (cubic micrometer). Please note that this is a physical parameter, do not set it too outrageous
RadiusIntervalofSoma = [40,60] # Here you set the interval of the radius of each soma. Please note that this is a physical parameter, do not set it too outrageous
SpaceDiscretization=1 # This is the minimum unit for spatial discretization
GrainedN=20 # % Here you need to define the granularity of neuron sphere, "GrainedN", as a natural number. A larger "GrainedN" ensures higher accuracy but requires more computing power 
NetrinFieldType=2 # This is the type of the initialization of netrin-1 field. NetrinFieldType=1 means that the field is homogeneous everywhere. NetrinFieldType=2 means that the netrin-1 field is non-homogeneous

## Initialization for space and soma
Radius=SystemInitializationFunctions.SomaInitialization(NumberofNeurons,RadiusIntervalofSoma) ## Initialize the size value of somas
LocationMatrix=SystemInitializationFunctions.LocationInitialization(InitializationType,NumberofNeurons,SpaceLimit) ## Initialize the space location of somas
InitializedSomaSphereCell,RealRadiusCell,NeighborCell=SystemInitializationFunctions.SomaandSpaceInitialization(LocationMatrix,Radius,GrainedN) ## Initialize the sphere of somas
NeuralTypeVector=SystemInitializationFunctions.NeuralTypeInitialization(NumberofNeurons,NeuralTypePar) ## Initialize the type of neurons
## Initialization for chemical substances
SubmembraneCell=SystemInitializationFunctions.ChemicalInitialization(NumberofNeurons,InitializedSomaSphereCell,GrainedN) ## Initialize the chemical substances on somas
SourceLocationMatrix=SystemInitializationFunctions.NetrinOneFieldInitialization(SpaceLimit,LocationMatrix,NetrinFieldType) ## Initialize the netrin-1 source
## Initialization of Development Information
DevelopmentInfoCell=SystemInitializationFunctions.DevelopInitialization(NumberofNeurons,GrainedN,DurationLength) ## Initialize the development information

## Growth
SubmembraneCellAcrossTime = np.empty((DurationLength,1),dtype=object)
SubmembraneCellAcrossTime[0,0]=SubmembraneCell
for IDT in range(1,DurationLength):
    print("\r", end="")
    print("Neural development: {}%: ".format(int(IDT/DurationLength*100)), "▋" * (int(IDT/DurationLength*100) // 2), end="")
    sys.stdout.flush()
    ## Stage 1 The Development of Lamellipodia
    DevelopmentInfoCell=MorphologySystemFunctions.LamellipodiaGrowthFunction(SubmembraneCell,DevelopmentInfoCell,NumberofNeurons,GrainedN) # Here you will update the Lamellipodia formation and condensation on the soma membrane
    ## Stage 2 The Development of Synapses
    SubmembraneCell,DevelopmentInfoCell=ChemicalSystemFunctions.IterationofChemicalSubstance(LocationMatrix,NeuralTypeVector,SubmembraneCell,DevelopmentInfoCell,RealRadiusCell,NeighborCell,SourceLocationMatrix,NumberofNeurons,GrainedN,IDT) # Here you will update the synaptic formation
    SubmembraneCellAcrossTime[IDT,0]=SubmembraneCell

## Visualization
Visualization.PlotResult(InitializedSomaSphereCell,DevelopmentInfoCell,NumberofNeurons,SpaceLimit,GrainedN,SourceLocationMatrix)