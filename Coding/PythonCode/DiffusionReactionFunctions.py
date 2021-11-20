import numpy as np
import itertools
from numpy.lib.function_base import append
import scipy

## This is the function used for Ca-related diffusion reaction processes
def IterationofCa(SubmembraneCell,DevelopmentInfoCell,CellofRealRadius,NeuronType,CellofNeighbors,IDN,IDC,GrainedN,IDSeg):
    SubmembraneCaDecayInterval=np.array([-10,-50],dtype=np.float32)*np.power(np.array(10,dtype=np.float32),-6) ## The decay of submembrane calcium concentration is [-100,-500] nM, where 1 nM = 10^(-6) mM
    KappaCa=0.1 ## Membrane permeability of Ca
    DCa=36000 ## DCa=10^(-7)cm^2/sec=36000 \mum^2/hour
    OrderDifference=np.power(np.array(10,dtype=np.float32),-4) ## The order difference between the external and submembrane Ca concentration
    CATransport=0.3 ## The active transport rate of Ca
    CDecay=0.01 ## The decay rate of Ca
    IDonS=1 ## The information is save on which column of the synaptic components
    IDonM=4 ## The information is save on which column of the soma membrane
    SingleNeuronType=NeuronType[IDN] ## Load the neural type information
    Neighbors=CellofNeighbors[IDN,IDC] ## Load the neighbor information of each coordinate on the soma sphere
    if SingleNeuronType==1: ## Excitable Membrane
        Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
        if IDSeg==0: ## Growth happens on this coordinate, this is the root segment
            ExternalCa=SubmembraneCell[IDN,3][Coordinate] ## This is the external Ca concentration closed to this coordinate
            NeighborCa=SubmembraneCell[IDN,IDonM][Neighbors[:,0],Neighbors[:,1]] ## This is the neighbor submembrane Ca concentrations
            CoordinateCa=SubmembraneCell[IDN,IDonM][Coordinate]  ## This is the submembrane Ca concentration on this coordinate
            CaDiftoOther=SubmembraneCaDecayInterval[0]+(SubmembraneCaDecayInterval[1]-SubmembraneCaDecayInterval[0])*np.random.rand(1,1) ## This is the decay of Ca concentration to other parts of cytoplasm (far from membrane)
            LaplaceS2=(CoordinateCa+CaDiftoOther+np.sum(NeighborCa)+OrderDifference*ExternalCa+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS])/(1+1+1+np.size(NeighborCa)) ## This is the solution of the Laplace equation by the relaxational method
            ActiveTransport=CATransport*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the active transport term
            Decay=CDecay*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the decay term
            SubmembraneCell[IDN,IDonM][Coordinate]=LaplaceS2-ActiveTransport-Decay ## Update the Ca concentration on shpere
            if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
                LaplaceS3=(SubmembraneCell[IDN,IDonM][Coordinate]+OrderDifference*ExternalCa+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS])/3 ## This is the solution of the Laplace equation by the relaxational method
                ActiveTransport=CATransport*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the active transport term
                Decay=CDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransport-Decay ## Update the tubulin concentration on synaptic segment
            else: ## There is next segment (or segments)
                NextSegments=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).reshape(-1,1)
                TNextSegment=np.zeros((1,np.size(NextSegments))).reshape(-1,1)
                for NS in range(np.size(NextSegments)):
                    TNextSegment[NS]=DevelopmentInfoCell[IDN,2][Coordinate][NextSegments[NS].astype(int),IDonS]
                LaplaceS3=(SubmembraneCell[IDN,IDonM][Coordinate]+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNextSegment)+OrderDifference*ExternalCa)/(3+len(TNextSegment)) ## Update the tubulin concentration on synaptic segment
                ActiveTransportin=CATransport*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the active transport-in term
                ActiveTransportout=CATransport*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the active transport-out term
                Decay=CDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay ## Update the tubulin concentration on synaptic segment
        elif IDSeg>0: ## Growth happens on this coordinate, this is not the root segment
            ExternalCa=SubmembraneCell[IDN,3][Coordinate] ## This is the external Ca concentration closed to this coordinate
            if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
                PreviousSegment=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,10]).reshape(-1,1) ## This is the previous segment
                TNearSegment=np.zeros((1,np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11])))).reshape(-1,1) ## We need to think about the situation where there was a time of branching right before this segment
                for PS in range(np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11]))):
                    TNearSegment[PS]=DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment[PS].astype(int),IDonS]
                LaplaceS3=(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNearSegment)+OrderDifference*ExternalCa)/(len(TNearSegment)+2) ## This is the solution of the Laplace equation by the relaxational method
                ActiveTransport=CATransport*DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,IDonS] ## This is the active transport term
                Decay=CDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransport-Decay ## Update the tubulin concentration on synaptic segment
            else: ## There is next segment (or segments)
                PreviousSegment=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,10]).reshape(-1,1) ## This is the previous segment
                TNearSegment=np.zeros((1,np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11])))).reshape(-1,1) ## We need to think about the situation where there was a time of branching right before this segment
                for PS in range(np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11]))):
                    TNearSegment[PS]=DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment[PS].astype(int),IDonS]
                NextSegments=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).reshape(-1,1) ## We need to think about the situation where there is a time of branching right after this segment
                TNextSegment=np.zeros((1,np.size(NextSegments))).reshape(-1,1)
                for NS in range(np.size(NextSegments)):
                    TNextSegment[NS]=DevelopmentInfoCell[IDN,2][Coordinate][NextSegments[NS].astype(int),IDonS]
                LaplaceS3=(np.sum(TNearSegment)+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNextSegment)+OrderDifference*ExternalCa)/(2+len(TNearSegment)+len(TNextSegment)) ## Update the tubulin concentration on synaptic segment
                ActiveTransportin=1/len(TNearSegment)*CATransport*DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,IDonS]*DevelopmentInfoCell[IDN,3][Coordinate][-1,0] ## This is the active transport-in term
                ActiveTransportout=CATransport*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]*DevelopmentInfoCell[IDN,3][Coordinate][-1,0] ## This is the active transport-out term
                Decay=CDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay ## Update the tubulin concentration on synaptic segment
        elif IDSeg<0: ## Growth does not happen on this coordinate
            ExternalCa=SubmembraneCell[IDN,3][Coordinate] ## This is the external Ca concentration closed to this coordinate
            NeighborCa=np.array(SubmembraneCell[IDN,IDonM])[Neighbors[:,0],Neighbors[:,1]] ## This is the neighbor submembrane Ca concentrations
            CoordinateCa=SubmembraneCell[IDN,IDonM][Coordinate] ## This is the submembrane Ca concentration on this coordinate
            CaDiftoOther=SubmembraneCaDecayInterval[0]+(SubmembraneCaDecayInterval[1]-SubmembraneCaDecayInterval[0])*np.random.rand(1,1) ## This is the decay of Ca concentration to other parts of cytoplasm (far from membrane)
            Decay=CDecay*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the decay term
            SubmembraneCell[IDN,IDonM][Coordinate]=(CoordinateCa+CaDiftoOther+np.sum(NeighborCa)+OrderDifference*ExternalCa)/(1+1+np.size(NeighborCa))-Decay ## Solution of the Laplace equation by the relaxational method
    elif SingleNeuronType==2: ## Passive Membrane
        RealRadiusVector=CellofRealRadius[IDN,0] ## Load the real radius information of each coordinate on the soma sphere
        RealRadius=RealRadiusVector[IDC] ## The real radius of this coor
        Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
        if IDSeg==0: ## Growth happens on this coordinate, this is the root segment
            ExternalCa=SubmembraneCell[IDN,3][Coordinate] ## This is the external Ca concentration closed to this coordinate
            ActiveTransport=CATransport*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the active transport term
            Decay=CDecay*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the decay term
            SubmembraneCell[IDN,IDonM][Coordinate]=(KappaCa/(KappaCa+DCa/RealRadius))*ExternalCa-ActiveTransport-Decay ## Update the tubulin concentration on shpere
            if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
                LaplaceS3=(SubmembraneCell[IDN,IDonM][Coordinate]+OrderDifference*ExternalCa+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS])/3 ## This is the solution of the Laplace equation by the relaxational method
                ActiveTransport=CATransport*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the active transport term
                Decay=CDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransport-Decay ## Update the tubulin concentration on synaptic segment
            else: ## There is next segment (or segments)
                NextSegments=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).reshape(-1,1)
                TNextSegment=np.zeros((1,np.size(NextSegments))).reshape(-1,1)
                for NS in range(np.size(NextSegments)):
                    TNextSegment[NS]=DevelopmentInfoCell[IDN,2][Coordinate][NextSegments[NS].astype(int),IDonS]
                LaplaceS3=(SubmembraneCell[IDN,IDonM][Coordinate]+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNextSegment)+OrderDifference*ExternalCa)/(3+len(TNextSegment)) ## Update the tubulin concentration on synaptic segment
                ActiveTransportin=CATransport*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the active transport-in term
                ActiveTransportout=CATransport*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the active transport-out term
                Decay=CDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay ## Update the tubulin concentration on synaptic segment
        elif IDSeg>0: ## Growth happens on this coordinate, this is not the root segment
            ExternalCa=SubmembraneCell[IDN,3][Coordinate] ## This is the external Ca concentration closed to this coordinate
            if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
                PreviousSegment=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,10]).reshape(-1,1) ## This is the previous segment
                TNearSegment=np.zeros((1,np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11])))).reshape(-1,1) ## We need to think about the situation where there was a time of branching right before this segment
                for PS in range(np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11]))):
                    TNearSegment[PS]=DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment[PS].astype(int),IDonS]
                LaplaceS3=(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNearSegment)+OrderDifference*ExternalCa)/(len(TNearSegment)+2) ## This is the solution of the Laplace equation by the relaxational method
                ActiveTransport=CATransport*DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,IDonS] ## This is the active transport term
                Decay=CDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransport-Decay ## Update the tubulin concentration on synaptic segment
            else: ## There is next segment (or segments)
                PreviousSegment=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,10]).reshape(-1,1) ## This is the previous segment
                TNearSegment=np.zeros((1,np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11])))).reshape(-1,1) ## We need to think about the situation where there was a time of branching right before this segment
                for PS in range(np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11]))):
                    TNearSegment[PS]=DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment[PS].astype(int),IDonS]
                NextSegments=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).reshape(-1,1) ## We need to think about the situation where there is a time of branching right after this segment
                TNextSegment=np.zeros((1,np.size(NextSegments))).reshape(-1,1)
                for NS in range(np.size(NextSegments)):
                    TNextSegment[NS]=DevelopmentInfoCell[IDN,2][Coordinate][NextSegments[NS].astype(int),IDonS]
                LaplaceS3=(np.sum(TNearSegment)+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNextSegment)+OrderDifference*ExternalCa)/(2+len(TNearSegment)+len(TNextSegment)) ## Update the tubulin concentration on synaptic segment
                ActiveTransportin=1/len(TNearSegment)*CATransport*DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,IDonS]*DevelopmentInfoCell[IDN,3][Coordinate][-1,0] ## This is the active transport-in term
                ActiveTransportout=CATransport*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]*DevelopmentInfoCell[IDN,3][Coordinate][-1,0] ## This is the active transport-out term
                Decay=CDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
                DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay ## Update the tubulin concentration on synaptic segment
        elif IDSeg<0: ## Growth does not happen on this coordinate
            ExternalCa=SubmembraneCell[IDN,3][Coordinate] ## This is the external Ca concentration closed to this coordinate
            Decay=CDecay*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the decay term
            SubmembraneCell[IDN,IDonM][Coordinate]=(KappaCa/(KappaCa+DCa/RealRadius))*ExternalCa-Decay ## Update the tubulin concentration on shpere
    return SubmembraneCell,DevelopmentInfoCell

## This is the function used for tubulin-related diffusion reaction processes
def IterationofTubulin(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,GrainedN,IDSeg,IDT):
    TProduction=0.7 ## The production rate of tubulin
    TATransport=0.5 ## The active transport rate of tubulin
    MaxTAT=0.5 ## The maximum active transport rate of tubulin
    TDecay=0.0001 ## The decay rate of tubulin
    TAssembly=0.1 ## The assembly rate of tubulin
    TDisassembly=0.0005 ## The disassembly rate of tubulin
    IDonS=2 ## The information is save on which column of the synaptic components
    IDonM=5 ## The information is save on which column of the soma membrane
    Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
    if IDSeg==0: ## This is the root segment
        LaplaceS2=(SubmembraneCell[IDN,IDonM][Coordinate]+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS])/2 ## This is the solution of the Laplace equation by the relaxational method
        ActiveTransport=np.min(np.array([TATransport*DevelopmentInfoCell[IDN,3][Coordinate][IDT-1,0],MaxTAT]))*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the active transport term
        Decay=TDecay*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the decay term
        SubmembraneCell[IDN,IDonM][Coordinate]=LaplaceS2+TProduction-ActiveTransport-Decay ## Update the tubulin concentration on shpere
        if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
            LaplaceS2=(SubmembraneCell[IDN,IDonM][Coordinate]+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS])/2 ## This is the solution of the Laplace equation by the relaxational method
            Assembly=TAssembly*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,4] ## This is the assembly
            Disassembly=TDisassembly ## This is the disassembly
            Decay=TDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS2+ActiveTransport-Decay-Assembly+Disassembly ## Update the tubulin concentration on synaptic segment
        else: ## There is next segment (or segments)
            NextSegments=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).reshape(-1,1)
            TNextSegment=np.zeros((1,np.size(NextSegments))).reshape(-1,1)
            for NS in range(np.size(NextSegments)):
                TNextSegment[NS]=DevelopmentInfoCell[IDN,2][Coordinate][NextSegments[NS].astype(int),IDonS]
            LaplaceS3=(SubmembraneCell[IDN,IDonM][Coordinate]+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNextSegment))/(2+len(TNextSegment)) ## Update the tubulin concentration on synaptic segment
            ActiveTransportin=np.min(np.array([TATransport*DevelopmentInfoCell[IDN,3][Coordinate][IDT-1,0],MaxTAT]))*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the active transport-in term
            ActiveTransportout=np.min(np.array([TATransport*DevelopmentInfoCell[IDN,3][Coordinate][IDT-1,0],MaxTAT]))*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the active transport-out term
            Decay=TDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay ## Update the tubulin concentration on synaptic segment
    elif IDSeg>0: ## This is not the root segment
        if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
            PreviousSegment=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,10]).reshape(-1,1) ## This is the previous segment
            TNearSegment=np.zeros((1,np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11])))).reshape(-1,1) ## We need to think about the situation where there was a time of branching right before this segment
            for PS in range(np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11]))):
                TNearSegment[PS]=DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment[PS].astype(int),IDonS]
            LaplaceS2=(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNearSegment))/(len(TNearSegment)+1) ## This is the solution of the Laplace equation by the relaxational method
            ActiveTransport=np.min(np.array([TATransport*DevelopmentInfoCell[IDN,3][Coordinate][IDT-1,0],MaxTAT]))*DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,IDonS] ## This is the active transport term
            Assembly=TAssembly*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,4] ## This is the assembly
            Disassembly=TDisassembly ## This is the disassembly
            Decay=TDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS2+ActiveTransport-Decay-Assembly+Disassembly ## Update the tubulin concentration on synaptic segment
        else: ## There is next segment (or segments)
            PreviousSegment=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,10]).reshape(-1,1) ## This is the previous segment
            TNearSegment=np.zeros((1,np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11])))).reshape(-1,1) ## We need to think about the situation where there was a time of branching right before this segment
            for PS in range(np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11]))):
                TNearSegment[PS]=DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment[PS].astype(int),IDonS]
            NextSegments=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).reshape(-1,1) ## We need to think about the situation where there is a time of branching right after this segment
            TNextSegment=np.zeros((1,np.size(NextSegments))).reshape(-1,1)
            for NS in range(np.size(NextSegments)):
                TNextSegment[NS]=DevelopmentInfoCell[IDN,2][Coordinate][NextSegments[NS].astype(int),IDonS]
            LaplaceS3=(np.sum(TNearSegment)+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNextSegment))/(1+len(TNearSegment)+len(TNextSegment)) ## Update the tubulin concentration on synaptic segment
            ActiveTransportin=1/len(TNearSegment)*np.min(np.array([TATransport*DevelopmentInfoCell[IDN,3][Coordinate][IDT-1,0],MaxTAT]))*DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,IDonS] ## This is the active transport-in term
            ActiveTransportout=np.min(np.array([TATransport*DevelopmentInfoCell[IDN,3][Coordinate][IDT-1,0],MaxTAT]))*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the active transport-out term
            Decay=TDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay ## Update the tubulin concentration on synaptic segment
    return SubmembraneCell,DevelopmentInfoCell

## This is the function used for UMAP-related diffusion reaction processes
def IterationofUMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,GrainedN,IDSeg):
    UProduction=0.8 ## The production rate of UMP2
    BindingR=0.4 ## The binding rate
    UnBindingR=0.3 ## The binding rate
    UDecay=0.0003 ## The decay rate of UMP2
    IDonS=3 ## The information is save on which column of the synaptic components
    IDonM=6 ## The information is save on which column of the soma membrane
    Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
    if IDSeg==0: ## This is the root segment
        LaplaceS2=(SubmembraneCell[IDN,IDonM][Coordinate]+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS])/2 ## This is the solution of the Laplace equation by the relaxational method
        Decay=UDecay*SubmembraneCell[IDN,IDonM][Coordinate] ## This is the decay term
        SubmembraneCell[IDN,IDonM][Coordinate]=LaplaceS2+UProduction-Decay ## Update the UMAP concentration on shpere
        if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
            LaplaceS2=(SubmembraneCell[IDN,IDonM][Coordinate]+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS])/2 ## This is the solution of the Laplace equation by the relaxational method
            Decay=UDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
            Binding=BindingR*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the binding term
            Unbinding=UnBindingR*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,4] ## This is the un-binding term
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS2-Decay-Binding+Unbinding ## Update the UMAP concentration on synaptic segment
        else: ## There is next segment (or segments)
            NextSegments=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).reshape(-1,1) ## We need to think about the situation where there is a time of branching right after this segment
            TNextSegment=np.zeros((1,np.size(NextSegments))).reshape(-1,1)
            for NS in range(np.size(NextSegments)):
                TNextSegment[NS]=DevelopmentInfoCell[IDN,2][Coordinate][NextSegments[NS].astype(int),IDonS]
            LaplaceS3=(SubmembraneCell[IDN,IDonM][Coordinate]+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNextSegment))/(2+len(TNextSegment)) ## Update the UMAP concentration on synaptic segment
            Decay=UDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3-Decay ## Update the UMAP concentration on synaptic segment
    elif IDSeg>0: ## This is not the root segment
        if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
            PreviousSegment=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,10]).reshape(-1,1) ## This is the previous segment
            TNearSegment=np.zeros((1,np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11])))).reshape(-1,1) ## We need to think about the situation where there was a time of branching right before this segment
            for PS in range(np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11]))):
                TNearSegment[PS]=DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment[PS].astype(int),IDonS]
            LaplaceS2=(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+sum(TNearSegment))/(len(TNearSegment)+1) ## This is the solution of the Laplace equation by the relaxational method
            Decay=UDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
            Binding=BindingR*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the binding term
            Unbinding=UnBindingR*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,4] ## This is the un-binding term
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS2-Decay-Binding+Unbinding ## Update the UMAP concentration on synaptic segment
        else: ## There is next segment (or segments)
            PreviousSegment=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,10]).reshape(-1,1) ## This is the previous segment
            TNearSegment=np.zeros((1,np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11])))).reshape(-1,1) ## We need to think about the situation where there was a time of branching right before this segment
            for PS in range(np.size(np.size(DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment,11]))):
                TNearSegment[PS]=DevelopmentInfoCell[IDN,2][Coordinate][PreviousSegment[PS].astype(int),IDonS]
            NextSegments=np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).reshape(-1,1) ## We need to think about the situation where there is a time of branching right after this segment
            TNextSegment=np.zeros((1,np.size(NextSegments))).reshape(-1,1)
            for NS in range(np.size(NextSegments)):
                TNextSegment[NS]=DevelopmentInfoCell[IDN,2][Coordinate][NextSegments[NS].astype(int),IDonS]
            LaplaceS3=(np.sum(TNearSegment)+DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]+np.sum(TNextSegment))/(1+len(TNearSegment)+len(TNextSegment)) ## Update the UMAP concentration on synaptic segment
            Decay=UDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
            DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=LaplaceS3-Decay ## Update the UMAP concentration on synaptic segment
    return SubmembraneCell,DevelopmentInfoCell

## This is the function used for BMAP-related diffusion reaction processes
def IterationofBMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,GrainedN,IDSeg):
    BindingR=0.4 ## The binding rate
    UnBindingR=0.3 ## The binding rate
    PR=0.4 ## The phosphorylation rate
    DR=0.3 ## The dephosphorylation rate
    UDecay=0.0003 ## The decay rate of UMP2
    IDonS=4 ## The information is save on which column of the synaptic components
    KF=0.0000005 ## Calcium rate constant for phosphorylation
    KG=0.0000005 ## Calcium rate constant for dephosphorylation
    Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
    if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
        Decay=UDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
        Binding=BindingR*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,3] ## This is the binding term
        Unbinding=UnBindingR*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the un-binding term
        PRL=(np.power(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,1],2))/(KF+(np.power(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,1],2))) ## Phosphorylation rate limit
        DRL=(np.power(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,1],2))/(KG+(np.power(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,1],2))) ## Dephosphorylation rate limit
        Phosphorylation=PR*PRL*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the phosphorylation term
        Dephosphorylation=DR*DRL*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,5] ## This is the dephosphorylation term
        DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=Binding-Unbinding-Phosphorylation+Dephosphorylation-Decay ## Update the BMAP concentration on synaptic segment
    return SubmembraneCell,DevelopmentInfoCell

## This is the function used for PMAP-related diffusion reaction processes
def IterationofPMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,GrainedN,IDSeg):
    PR=0.4 ## The phosphorylation rate
    DR=0.3 ## The dephosphorylation rate
    PDecay=0.0005 ## The decay rate of PMP2
    IDonS=5 ## The information is save on which column of the synaptic components
    KF=0.0000005 ## Calcium rate constant for phosphorylation
    KG=0.0000005 ## Calcium rate constant for dephosphorylation
    Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
    if np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,11]).any()==None: ## There is no next segment
        Decay=PDecay*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the decay term
        PRL=(np.power(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,1],2))/(KF+(np.power(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,1],2))) ## Phosphorylation rate limit
        DRL=(np.power(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,1],2))/(KG+(np.power(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,1],2))) ## Dephosphorylation rate limit
        Phosphorylation=PR*PRL*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,4] ## This is the phosphorylation term
        Dephosphorylation=DR*DRL*DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS] ## This is the dephosphorylation term
        DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,IDonS]=Phosphorylation-Dephosphorylation-Decay ## Update the PMAP concentration on synaptic segment
    return SubmembraneCell,DevelopmentInfoCell







