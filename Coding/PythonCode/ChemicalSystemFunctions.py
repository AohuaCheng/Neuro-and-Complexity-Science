import numpy as np
import itertools
from numpy.lib.function_base import append
import scipy
import DiffusionReactionFunctions
import MorphologySystemFunctions

## This is the main function to run the diffusion-reaction processes and morphologyical dynamics
def IterationofChemicalSubstance(LocationMatrix,NeuronType,SubmembraneCell,DevelopmentInfoCell,CellofRealRadius,CellofNeighbors,SourceLocationMatrix,NumberofNeurons,GrainedN,IDT):
    for IDN in range(NumberofNeurons):
        for IDC in range(GrainedN*GrainedN): # Traverse every coordinate on the soma shpere
            Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
            if DevelopmentInfoCell[IDN,1][Coordinate]==1: ## Growth happens on this coordinate
                for IDSeg in range(np.size(DevelopmentInfoCell[IDN,2][Coordinate],0)):
                ## Here is a small tip: although SubmembraneCell and DevelopmentInfoCell will
                ## be updated in this "for-loop", the variable
                ## "np.size(DevelopmentInfoCell[IDN,2][Coordinate]"
                ## counts the size of the original DevelopmentInfoCell that
                ## send into this "for-loop". Therefore, this variable is
                ## constant rather than dynamic
                    SubmembraneCell,DevelopmentInfoCell=DiffusionReactionFunctions.IterationofCa(SubmembraneCell,DevelopmentInfoCell,CellofRealRadius,NeuronType,CellofNeighbors,IDN,IDC,GrainedN,IDSeg)
                    SubmembraneCell,DevelopmentInfoCell=DiffusionReactionFunctions.IterationofTubulin(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,GrainedN,IDSeg,IDT)
                    SubmembraneCell,DevelopmentInfoCell=DiffusionReactionFunctions.IterationofUMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,GrainedN,IDSeg)
                    SubmembraneCell,DevelopmentInfoCell=DiffusionReactionFunctions.IterationofBMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,GrainedN,IDSeg)
                    SubmembraneCell,DevelopmentInfoCell=DiffusionReactionFunctions.IterationofPMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,GrainedN,IDSeg)
                    SubmembraneCell,DevelopmentInfoCell=MorphologySystemFunctions.ElongationAndBranching(LocationMatrix,SubmembraneCell,DevelopmentInfoCell,SourceLocationMatrix,IDN,IDC,GrainedN,IDSeg,IDT)
                    #print([IDSeg,np.size(DevelopmentInfoCell[IDN,2][Coordinate],0),np.size(DevelopmentInfoCell[IDN,2][Coordinate],1)])
                if np.array(DevelopmentInfoCell[IDN,4][Coordinate]).any()!=None:
                    DevelopmentInfoCell[IDN,3][Coordinate][IDT,0]=np.mean(np.array(DevelopmentInfoCell[IDN,4][Coordinate]).reshape(-1,2)[:,0])
                    DevelopmentInfoCell[IDN,3][Coordinate][IDT,1]=np.var(np.array(DevelopmentInfoCell[IDN,4][Coordinate]).reshape(-1,2)[:,0])
                    DevelopmentInfoCell[IDN,3][Coordinate][IDT,2]=np.size(np.array(DevelopmentInfoCell[IDN,4][Coordinate]).reshape(-1,2),0)
                    DevelopmentInfoCell[IDN,3][Coordinate][IDT,3]=np.sum(np.array(DevelopmentInfoCell[IDN,4][Coordinate]).reshape(-1,2)[:,1])
                # the 1th column saves the average elongation rate, the 2nd
                # column saves the variance of elongation rate, the 3rd
                # column saves the number of terminal segments; %% the 4th
                # column saves the normalization base of centrifugal order;
                    DevelopmentInfoCell[IDN,4][Coordinate]=np.array([]) ## Clean the information
            else: ## Growth does not happen on this coordinate
                SubmembraneCell,DevelopmentInfoCell=DiffusionReactionFunctions.IterationofCa(SubmembraneCell,DevelopmentInfoCell,CellofRealRadius,NeuronType,CellofNeighbors,IDN,IDC,GrainedN,-1)
    return SubmembraneCell,DevelopmentInfoCell

