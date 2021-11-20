import numpy as np
import itertools
from numpy.lib.function_base import append # For unknown reason, function "append" can not be imported successfully until we declare it specifically
import scipy
import AxonDirectionFunctions

## This is the function for the out-growth of lamellipodia
def LamellipodiaGrowthFunction(SubmembraneCell,DevelopmentInfoCell,NumberofNeurons,GrainedN):
    MaxGrowingProb=0.00002 # The maximum growing probability per hour
    Alpha=1
    Beta=2
    for IDN in range(NumberofNeurons): ## Here we traverse every neuron
        for IDC in range(GrainedN*GrainedN): ## Here we traverse all coordinates on the soma sphere
            Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
            CoordinateCa=SubmembraneCell[IDN,4][Coordinate]  # This is the submembrane Ca concentration on this coordinate
            MaxCa=np.max(SubmembraneCell[IDN,4]) # This is the maximum submembrane Ca concentration at this moment
            GrowingProb=MaxGrowingProb*(Beta/(Beta-Alpha)*np.power(CoordinateCa/MaxCa,Alpha)-Alpha/(Beta-Alpha)*np.power(CoordinateCa/MaxCa,Beta)) # The growth probability on this coordinate
            DevelopmentInfoCell[IDN,0][Coordinate]=GrowingProb
            RandomizeBase=np.random.binomial(1,GrowingProb,1) # 1 means growth, 0 means non-growth
            DevelopmentInfoCell[IDN,1][Coordinate]=DevelopmentInfoCell[IDN,1][Coordinate]+RandomizeBase>0 # Randomize growth or non-growth
    return DevelopmentInfoCell

## This is the function for the development of synapses
def ElongationAndBranching(LocationMatrix,SubmembraneCell,DevelopmentInfoCell,SourceLocationMatrix,IDN,IDC,GrainedN,IDSeg,IDT):
    TAssembly=0.1 ## Tubulin assembly rate
    TDisassembly=0.0005 ## The disassembly rate of tubulin
    KB=0.01 ## The branching constant
    E=0.5 ## The competition parameter in branching
    S=-0.1 ## Order dependency in branching
    RotationAI=[-np.pi/6,np.pi/6] ## Rotation angle interval in branching
    DisturbanceAI=[-np.pi/8,np.pi/8] ## Rotation angle interval in non-branching
    InitialR=1 ## The radius of root segement
    DecayofR=0.8 ## The decreasing rate of the radius of segements
    InheritanceR=0.9 ## The nheritance rate
    RecordThreshold=10 ## The minimum elongation record threshold
    EnlongationR=1000 ## Enlongation rate
    BranchingR=500 ## Branching rate
    Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
    ## Update the Tubulin concentratio on the synapse (if there is any synapse)
    if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
        EnL=EnlongationR*TAssembly*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,2]*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,4]-TDisassembly ## This is the enlongation length
        if EnL>=RecordThreshold:
            #print(['Growinginginginginginging'])
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,6]=EnL ## Update the  enlongation length
            BaseBranchingProb=KB*(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,5]/(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,5]+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,4])) ## This is the baseline branching probability
            NumofTS=DevelopmentInfoCell[IDN,3][Coordinate][IDT-1,2] ## Number of terminal segments
            Normal=DevelopmentInfoCell[IDN,3][Coordinate][IDT-1,3] ## Normalization base of centrifugal order
            CTerm=1/NumofTS*Normal ## Normalization term of centrifugal order
            BranchingProb=np.array(BranchingR*BaseBranchingProb*np.power(NumofTS,-E)*np.power(2,(-S*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,12]))/CTerm,dtype=np.float32) ## This is the branching probability
            if np.isnan(np.array(BranchingProb,dtype=float)): 
                BranchingProb=np.array([0])
            if BranchingProb>1:
                BranchingProb=np.array([1])
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,7]=list(np.ravel(BranchingProb))[0]
            BranchingOrNot=np.random.binomial(1,DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,7],1) # Generate a random variable for branching or non-branching
            if BranchingOrNot==1: ## Branching
                #print(['Branching'])
                if IDSeg==0: ## This is the root segment
                    MembranceC=[SubmembraneCell[IDN,0][Coordinate],SubmembraneCell[IDN,1][Coordinate],SubmembraneCell[IDN,2][Coordinate]] ## This is the coordinate on membrane
                    BaseVector=MembranceC-LocationMatrix[IDN,:] ## This is the base direction from the soma center
                    BaseVector=BaseVector/np.linalg.norm(BaseVector,ord=2)*EnL ## Here we normalize the direction vector and multiply the enlongation length
                    LVector=AxonDirectionFunctions.Rotation3D(BaseVector,np.array(RotationAI[0]+(RotationAI[1]-RotationAI[0])*np.random.rand(1,3)).reshape(3,1),MembranceC) ## Left branching daughter segment
                    RVector=AxonDirectionFunctions.Rotation3D(BaseVector,np.array(RotationAI[0]+(RotationAI[1]-RotationAI[0])*np.random.rand(1,3)).reshape(3,1),MembranceC) ## Right branching daughter segment
                    LDirection=(LVector-MembranceC)/np.linalg.norm((LVector-MembranceC)) ## The corresponding direction of left segment
                    RDirection=(RVector-MembranceC)/np.linalg.norm((RVector-MembranceC)) ## The corresponding direction of right segment
                    DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,9]=BaseVector/np.linalg.norm(BaseVector) ## Here we save the historical direction information
                    DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,0]=MembranceC ## Here we save the coordinate of the root
                    LSegment=np.empty((1,14),dtype=object) ## This is the cell to save the information of the left branching daughter segment
                    RSegment=np.empty((1,14),dtype=object) ## This is the cell to save the information of the right branching daughter segment
                    for IDL in range(1,4):
                        LSegment[0,IDL]=SubmembraneCell[IDN,IDL+3][Coordinate]*InheritanceR*1/2 ## Initialization of chemical substance
                        RSegment[0,IDL]=SubmembraneCell[IDN,IDL+3][Coordinate]*InheritanceR*1/2 ## Initialization of chemical substance
                    for IDL in range(4,6):
                        LSegment[0,IDL]=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDL]*InheritanceR*1/2 ## Initialization of chemical substance
                        RSegment[0,IDL]=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDL]*InheritanceR*1/2 ## Initialization of chemical substance
                elif IDSeg>0: ## This is not the root segment
                    HistoryD=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,9] ## This is the history direction
                    CurrentC=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,0] ## This is the current coordinate
                    BaseVector=np.ravel(HistoryD)*EnL ## Here we normalize the direction vector and multiply the enlongation length
                    LVector=AxonDirectionFunctions.Rotation3D(BaseVector,np.array(RotationAI[0]+(RotationAI[1]-RotationAI[0])*np.random.rand(1,3)).reshape(3,1),CurrentC) ## Left branching daughter segment
                    RVector=AxonDirectionFunctions.Rotation3D(BaseVector,np.array(RotationAI[0]+(RotationAI[1]-RotationAI[0])*np.random.rand(1,3)).reshape(3,1),CurrentC) ## Right branching daughter segment
                    LDirection=(LVector-CurrentC)/np.linalg.norm((LVector-CurrentC).astype(np.float32)) ## The corresponding direction of left segment
                    RDirection=(RVector-CurrentC)/np.linalg.norm((RVector-CurrentC).astype(np.float32)) ## The corresponding direction of right segment
                    LSegment=np.empty((1,14),dtype=object) ## This is the cell to save the information of the left branching daughter segment
                    RSegment=np.empty((1,14),dtype=object) ## This is the cell to save the information of the right branching daughter segment
                    for IDL in range(1,6):
                        LSegment[0,IDL]=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDL]*InheritanceR*1/2 ## Initialization of chemical substance
                        RSegment[0,IDL]=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDL]*InheritanceR*1/2 ## Initialization of chemical substance
                LSegment[0,0]=LVector ## The coordinate
                RSegment[0,0]=RVector ## The coordinate
                LSegment[0,6]=0  
                RSegment[0,6]=0  
                LSegment[0,7]=0  
                RSegment[0,7]=0  
                LSegment[0,9]=LDirection ## Here we save the historical direction information
                RSegment[0,9]=RDirection ## Here we save the historical direction information
                LSegment[0,10]=IDSeg ## The information of "previous segment"
                RSegment[0,10]=IDSeg ## The information of "previous segment"
                LSegment[0,12]=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,12]+1 ## The centrifugal order
                RSegment[0,12]=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,12]+1 ## The centrifugal order
                LSegment[0,13]=IDT ## The information of time
                RSegment[0,13]=IDT ## The information of time
                LSegment[0,8]=InitialR*np.power(DecayofR,LSegment[0,12]) ## The radius of segment
                RSegment[0,8]=InitialR*np.power(DecayofR,RSegment[0,12]) ## The radius of segment
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]=np.array([np.size(DevelopmentInfoCell[IDN,2][Coordinate],0),np.size(DevelopmentInfoCell[IDN,2][Coordinate],0)+1]) ## Update the information of "next segment"
                #print(np.array([DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11],np.size(DevelopmentInfoCell[IDN,2][Coordinate],0),np.size(DevelopmentInfoCell[IDN,2][Coordinate],0)+1]))
                DevelopmentInfoCell[IDN,2][Coordinate]=np.append(np.append(DevelopmentInfoCell[IDN,2][Coordinate],LSegment,axis=0),RSegment,axis=0) ## Add these two new segments
                ## Here DevelopmentInfoCell{ID,5} is used to save some information to calculate DevelopmentInfoCell{ID,4}
                InformationToRecord=np.append(EnL*np.ones((2,1)),np.append(np.array(np.power(2,(-S*LSegment[0,12]))).reshape(1,1),np.array(np.power(2,(-S*RSegment[0,12]))).reshape(1,1),axis=0),axis=1)
                if np.size(DevelopmentInfoCell[IDN,4][Coordinate])>0:
                    DevelopmentInfoCell[IDN,4][Coordinate]=np.append(np.ravel(DevelopmentInfoCell[IDN,4][Coordinate]),np.ravel(np.array(InformationToRecord)),axis=0)
                else:
                    DevelopmentInfoCell[IDN,4][Coordinate]=np.array(InformationToRecord)
            else: ## Non-branching
                if IDSeg==0: ## This is the root segment
                    DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,6]=EnL ## Update the  enlongation length
                    MembranceC=[SubmembraneCell[IDN,0][Coordinate],SubmembraneCell[IDN,1][Coordinate],SubmembraneCell[IDN,2][Coordinate]] ## This is the coordinate on membrane
                    BaseVector=MembranceC-LocationMatrix[IDN,:] ## This is the base direction from the soma center
                    BaseVector=BaseVector/np.linalg.norm(BaseVector)*EnL ## Here we normalize the direction vector and multiply the enlongation length
                    NVector=AxonDirectionFunctions.Rotation3D(BaseVector,np.array(DisturbanceAI[0]+(DisturbanceAI[1]-DisturbanceAI[0])*np.random.rand(3,1)).reshape(3,1),MembranceC) ## Next segment
                    NDirection=(NVector-MembranceC)/np.linalg.norm(NVector-MembranceC) ## Normalization
                    NVector,NDirection=AxonDirectionFunctions.SynapseGuidancebyNetrinOne(MembranceC,NDirection,EnL,SourceLocationMatrix,IDT) ## Axon guidance
                    DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,9]=BaseVector/np.linalg.norm(BaseVector) ## Here we save the historical direction information
                    DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,0]=MembranceC ## Here we save the coordinate of the root
                    NSegment=np.empty((1,14),dtype=object) ## This is the cell to save the information of the right branching daughter segment
                    for IDL in range(1,4):
                        NSegment[0,IDL]=SubmembraneCell[IDN,IDL+3][Coordinate]*InheritanceR ## Initialization of chemical substance
                    for IDL in range(4,6):
                        NSegment[0,IDL]=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDL]*InheritanceR ## Initialization of chemical substance
                elif IDSeg>0: ## This is not the root segment
                    HistoryD=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,9] ## This is the history direction
                    CurrentC=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,0] ## This is the current coordinate
                    BaseVector=np.ravel(HistoryD)*EnL ## Here we normalize the direction vector and multiply the enlongation length
                    NVector=AxonDirectionFunctions.Rotation3D(BaseVector,np.array(DisturbanceAI[0]+(DisturbanceAI[1]-DisturbanceAI[0])*np.random.rand(3,1)).reshape(3,1),CurrentC) ## Next segment
                    NDirection=(NVector-CurrentC)/np.linalg.norm(np.array(NVector-CurrentC,dtype=np.float32)) ## Normalization
                    NVector,NDirection=AxonDirectionFunctions.SynapseGuidancebyNetrinOne(CurrentC,NDirection,EnL,SourceLocationMatrix,IDT) ## Axon guidance
                    NSegment=np.empty((1,14),dtype=object) ## This is the cell to save the information of the right branching daughter segment
                    for IDL in range(1,6):
                        NSegment[0,IDL]=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDL]*InheritanceR ## Initialization of chemical substance
                NSegment[0,0]=NVector ## The coordinate
                NSegment[0,6]=0 
                NSegment[0,7]=0  
                NSegment[0,9]=(NDirection+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,9])/2 ## Here we save the historical direction information
                NSegment[0,9]=NSegment[0,9]/np.linalg.norm(np.array(NSegment[0,9],dtype=np.float32)) 
                NSegment[0,10]=IDSeg ## The information of "previous segment"
                NSegment[0,12]=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,12] ## The centrifugal order
                NSegment[0,13]=IDT ## The information of time
                NSegment[0,8]=InitialR*np.power(DecayofR,NSegment[0,12]) ## The radius of segment
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]=np.array([np.size(DevelopmentInfoCell[IDN,2][Coordinate],0)]) ## Update the information of "next segment"
                #print(np.array([DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11],np.size(DevelopmentInfoCell[IDN,2][Coordinate],0)+1]))
                DevelopmentInfoCell[IDN,2][Coordinate]=np.append(DevelopmentInfoCell[IDN,2][Coordinate],NSegment,axis=0) ## Add these two new segments
                ## Here DevelopmentInfoCell{ID,5} is used to save some information to calculate DevelopmentInfoCell{ID,4}
                InformationToRecord=[np.array([EnL]),np.power(2,(-S*NSegment[0,12]))]
                if np.size(DevelopmentInfoCell[IDN,4][Coordinate])>0:
                    DevelopmentInfoCell[IDN,4][Coordinate]=np.append(np.ravel(DevelopmentInfoCell[IDN,4][Coordinate]),np.ravel(np.array(InformationToRecord)),axis=0)
                else:
                    DevelopmentInfoCell[IDN,4][Coordinate]=np.array(InformationToRecord)
        elif EnL>=0 and EnL<RecordThreshold: ## We do not need to add new segment, but we need to record the enlongation rate
        ## Here DevelopmentInfoCell{ID,5} is used to save some information to calculate DevelopmentInfoCell{ID,4}
            InformationToRecord=np.array([EnL,np.power(2,(-S*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,12]))]).reshape(1,2)
            DevelopmentInfoCell[IDN,4][Coordinate]=np.append(np.ravel(DevelopmentInfoCell[IDN,4][Coordinate]),np.ravel(np.array(InformationToRecord)),axis=0)
    return SubmembraneCell,DevelopmentInfoCell

