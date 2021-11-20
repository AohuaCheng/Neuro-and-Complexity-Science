import numpy as np
import math
import itertools
import scipy
from numpy.lib.function_base import append

## This is the 3-d rotation funtion used in synapse direction definition
def Rotation3D(C,T,Origin):
    # T defines rotation angles
    # C is the original coordinate
    # Origin is the origin coordinate of rotation
    ## Define the rotation matrix
    T1=T[0]
    T2=T[1]
    T3=T[2] ## Here T1,T2, and T3 are rotation angles
    Rx=np.matrix([[1,0,0],[0,np.cos(T1),-np.sin(T1)],[0,np.sin(T1),np.cos(T1)]]) ## Rotation in X direction
    Ry=np.matrix([[np.cos(T2),0,np.sin(T2)],[0,1,0],[-np.sin(T2),0,np.cos(T2)]])  ## Rotation in Y direction
    Rz=np.matrix([[np.cos(T3),-np.sin(T3),0],[np.sin(T3),np.cos(T3),0],[0,0,1]])  ## Rotation in Z direction
    R=np.matmul(np.matmul(Rx,Ry),Rz) ## Rotation matrix
    ## Transform the coordinate
    NewC=C*R+Origin 
    return NewC

## This is the function to realize axon guidance by the Netrin-1
def SynapseGuidancebyNetrinOne(CurrentLocation,CurrentDirection,EnL,SourceLocationMatrix,IDT):
    Diffusion=36000 ## DCa=10^(-7)cm^2/sec=36000 \mum^2/hour
    Production=36*np.power(np.array(10,dtype=np.float32),-11) ## Production rate is 10^(-7)nM/sec=36*10^(-5)nM/hour=36*10^(-11)mM/hour 
    DeltaR=EnL ## The width of growth cone is set as 20 mum 
    MinC=0.01*np.power(np.array(10,dtype=np.float32),-15) ## The minimal concentration required for guidance is 0.05*10^(-9)nM or 0.05*10^(-15)mM above the initial concentration
    MinP=0.01 ## The minimal fractional-change of guidance is 0.01
    GuidanceEffect=2 ## The capacity for Netrin-1 to change the axon path
    if np.size(SourceLocationMatrix)>=3: ## The Netrin Field is non-homogeneous
        DisMatrix=scipy.spatial.distance.cdist(np.array(CurrentLocation).reshape(1,3),np.array(SourceLocationMatrix).reshape(-1,3),metric='euclidean') ## This is the distance matrix
        CMatrix=np.multiply(Production*np.power((4*np.pi*Diffusion*DisMatrix),-1),math.erfc(DisMatrix/np.sqrt(4*Diffusion*IDT))) ## Solve the equation of Netrin-1 concentration
        Concen=np.sum(CMatrix) ## This is the concentration of Netrin-1 at current location
        PMatrixTerm1=(SourceLocationMatrix.astype(np.float16).reshape(-1,3)-np.array(CurrentLocation,np.float16))*DeltaR*Production  ## Solve the equation of DeltaC/C (Term 1)
        PMatrixTerm2=np.power(4*np.pi*Diffusion*np.power(DisMatrix,3),-1)*math.erfc(DisMatrix/np.sqrt(4*Diffusion*IDT)) ## Solve the equation of DeltaC/C (Term 2)
        PMatrixTerm3=np.power(4*np.sqrt(np.power(np.pi,3)*np.power(Diffusion,3)*IDT)*np.power(DisMatrix,2),-1)*np.exp(-1*np.power(DisMatrix,2)/(4*Diffusion*IDT)) ## Solve the equation of DeltaC/C (Term 3)
        PMatrix=PMatrixTerm1*(PMatrixTerm2+PMatrixTerm3)/Concen ## Solve the equation of DeltaC/C
        PGradient=np.sum(PMatrix,axis=1) ## This is the DeltaC/C
        ## Decide whether the guidance is possible
        if Concen>=MinC and np.linalg.norm(PGradient)>=MinP: ## The guidance is possible
           Direction=PGradient/np.linalg.norm(PGradient)    
           NDirection=((1-GuidanceEffect*np.linalg.norm(PGradient))*CurrentDirection+GuidanceEffect*np.linalg.norm(PGradient)*Direction)/2 ## Guidance
           NDirection=NDirection/np.linalg.norm(NDirection) ## Normalization of Guidance
           NLocation=EnL*NDirection+CurrentLocation ## Next location
        else: ## The guidance is impossible
           NDirection=CurrentDirection
           NLocation=EnL*CurrentDirection+CurrentLocation ## Next location
    else: ## The Netrin Field is homogeneous
        NDirection=CurrentDirection
        NLocation=EnL*CurrentDirection+CurrentLocation ## Next location
    return NLocation,NDirection
