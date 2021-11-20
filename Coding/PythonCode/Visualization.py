import numpy as np
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
import sys



def PlotResult(InitializedSomaSphereCell,DevelopmentInfoCell,NumberofNeurons,SpaceLimit,GrainedN,SourceLocationMatrix):
    fig=plt.figure()
    ax=fig.gca(projection='3d')
    for IDN in range(NumberofNeurons): ## Here we traverse every neuron
        ax.plot_surface(InitializedSomaSphereCell[IDN,0],InitializedSomaSphereCell[IDN,1],InitializedSomaSphereCell[IDN,2],color='b')
        NumCoor=0
        for IDC in range(GrainedN*GrainedN): # Traverse every coordinate on the soma shpere
            #print("\r", end="")
            #print("Visualization: {}%: ".format(int((IDN*GrainedN*GrainedN+IDC)/(NumberofNeurons*GrainedN*GrainedN)*100)), "▋" * (int((IDN*GrainedN*GrainedN+IDC)/(NumberofNeurons*GrainedN*GrainedN)*100) // 2), end="")
            #sys.stdout.flush()
            Coordinate=np.unravel_index(IDC, (GrainedN,GrainedN)) # Pick a coordinate
            if DevelopmentInfoCell[IDN,1][Coordinate]==1: ## Growth happens on this coordinate
                #print(['Synapses are growing',np.size(DevelopmentInfoCell[IDN,2][Coordinate],0)])
                NumCoor=NumCoor+1
                BaseNum=len(np.where(DevelopmentInfoCell[IDN,1]==1))
                for IDSeg in range(1,np.size(DevelopmentInfoCell[IDN,2][Coordinate],0)):
                    CurrentC=np.ravel(np.array(DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,0]).astype(np.float32))
                    PreviousC=np.ravel(np.array(DevelopmentInfoCell[IDN,2][Coordinate][DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,10],0]).astype(np.float32))
                    #print([np.size(PreviousC,0),np.size(PreviousC,1),np.size(CurrentC,0),np.size(CurrentC,1)])
                    Ratio=3000/SpaceLimit[0]
                    R=DevelopmentInfoCell[IDN,2][Coordinate][IDSeg,8]*Ratio
                    ax.plot([PreviousC[0],CurrentC[0]], [PreviousC[1],CurrentC[1]], [PreviousC[2],CurrentC[2]], color='b', linewidth=R)
    plt.show()
