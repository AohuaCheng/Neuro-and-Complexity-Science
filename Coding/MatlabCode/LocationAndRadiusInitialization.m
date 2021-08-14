function [LocationMatrix,RadiusVector,CellofInitialSP,CellofNeighbors,CellofRealRadius]=LocationAndRadiusInitialization(InitializationType,InitialPar,GrainedN)
RadiusVector=InitialPar.RadiusIntervalofSoma(1) + (InitialPar.RadiusIntervalofSoma(2)-InitialPar.RadiusIntervalofSoma(1)).*rand(InitialPar.NumberofNeurons,1); % This is the vector to store the intial radius
if InitializationType==1 %% This is the grid initialization
   NeuronNinEachDim=(InitialPar.NumberofNeurons).^(1/3); % *dont understand this parameter* 
   Min=zeros(1,3); % This is a vector to store the minimum X,Y,Z coordinates
   Max=zeros(1,3); % This is a vector to store the maximum X,Y,Z coordinates
   Gap=zeros(1,3); % This is a vector to store the gap between two neighbor neurons in each dim
   for IDDim=1:3
       Min(IDDim)=0.1*InitialPar.SpaceLimit(IDDim);
       Max(IDDim)=0.9*InitialPar.SpaceLimit(IDDim);
       Gap(IDDim)=0.8*InitialPar.SpaceLimit(IDDim)/NeuronNinEachDim; % *dont understand this dividing settings*
   end
   [PreX,PreY,PreZ]=meshgrid([Min(1):Gap(1):Max(1)],[Min(2):Gap(2):Max(2)],[Min(3):Gap(3):Max(3)]); % This is coordinate matrix. The raw number equals neuron number and the column number is 3
   LocationMatrix=[PreX(:),PreY(:),PreZ(:)];
elseif InitializationType==2 %% This is the coordinate initialization
   MiniDistanceBN=150; %% This is the minimum distance between any pair of neurons
   Min=zeros(1,3); % This is a vector to store the minimum X,Y,Z coordinates
   Max=zeros(1,3); % This is a vector to store the maximum X,Y,Z coordinates
   for IDDim=1:3
       Min(IDDim)=0.1*InitialPar.SpaceLimit(IDDim);
       Max(IDDim)=0.9*InitialPar.SpaceLimit(IDDim);
   end
   LocationMatrix=[Min(1)+(Max(1)-Min(1))*rand(1,1),Min(2)+(Max(2)-Min(2))*rand(1,1),Min(3)+(Max(3)-Min(3))*rand(1,1)]; %% Initialize one neuron's coordinates
   while size(LocationMatrix,1)<InitialPar.NumberofNeurons %% while there exist neurons without coordinates
        NCoord=[Min(1)+(Max(1)-Min(1))*rand(1,1),Min(2)+(Max(2)-Min(2))*rand(1,1),Min(3)+(Max(3)-Min(3))*rand(1,1)]; %% Generate one coordinate
        DistanceM=pdist([NCoord;LocationMatrix]); %% Calculate the distance between somas
        while any(DistanceM<MiniDistanceBN) 
              NCoord=[Min(1)+(Max(1)-Min(1))*rand(1,1),Min(2)+(Max(2)-Min(2))*rand(1,1),Min(3)+(Max(3)-Min(3))*rand(1,1)];
        end
        disp(['Initialize Neuron-',num2str(size(LocationMatrix,1)/InitialPar.NumberofNeurons)])
        LocationMatrix=[LocationMatrix;NCoord];
   end
end

CellofInitialSP=cell(InitialPar.NumberofNeurons,3); %% Here we save the coordinates on the initialized spheres
CellofNeighbors=cell(InitialPar.NumberofNeurons,1);
CellofRealRadius=cell(InitialPar.NumberofNeurons,1);
[PreX,PreY,PreZ] = sphere(GrainedN); %% Generate a unit sphere
for IDN=1:InitialPar.NumberofNeurons
    CL=LocationMatrix(IDN,:); %% This is the center of soma
    RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
    RadiusX=RadiusVector(IDN)*(0.9+0.2*rand(RealGrainedN,RealGrainedN)); %% Randomly generate the real radius on X direction. Perturbation is at an interval of [0%,20%]
    RadiusY=RadiusVector(IDN)*(0.9+0.2*rand(RealGrainedN,RealGrainedN)); %% Randomly generate the real radius on Y direction. Perturbation is at an interval of [0%,20%]
    RadiusZ=RadiusVector(IDN)*(0.9+0.2*rand(RealGrainedN,RealGrainedN)); %% Randomly generate the real radius on Z direction. Perturbation is at an interval of [0%,20%]
    X = PreX.*RadiusX;
    Y = PreY.*RadiusY;
    Z = PreZ.*RadiusZ;
    CellofInitialSP{IDN,1}=X+CL(1); %% This is the X coordinate
    CellofInitialSP{IDN,2}=Y+CL(2); %% This is the Y coordinate
    CellofInitialSP{IDN,3}=Z+CL(3); %% This is the Z coordinate
    SphereCenterDistance=squareform(pdist([CL;[reshape(CellofInitialSP{IDN,1},[RealGrainedN*RealGrainedN,1]),reshape(CellofInitialSP{IDN,2},[RealGrainedN*RealGrainedN,1]),reshape(CellofInitialSP{IDN,3},[RealGrainedN*RealGrainedN,1])]])); 
    SphereCenterDistance=SphereCenterDistance(1,2:end); %% The real radius (distance from each coordinate to the center)
    CellofRealRadius{IDN,1}=SphereCenterDistance;
    Coordinates=[reshape(CellofInitialSP{IDN,1},[RealGrainedN*RealGrainedN,1]),reshape(CellofInitialSP{IDN,2},[RealGrainedN*RealGrainedN,1]),reshape(CellofInitialSP{IDN,3},[RealGrainedN*RealGrainedN,1])];
    DistanceBetweenCoordinates=squareform(pdist(Coordinates)); %% Here you need to find the neighbors of each coordinate based on the distances
    DistanceBetweenCoordinates(find(DistanceBetweenCoordinates==0))=max(max(DistanceBetweenCoordinates)); 
    Neighbors=cell(size(Coordinates,1),1);
    for IDC=1:size(Coordinates,1)
        SortDistance=sort(unique(DistanceBetweenCoordinates(IDC,:))); %% Sort the distance so as to find the neighbors
        NeighborCoordinates=union(find(DistanceBetweenCoordinates(IDC,:)==SortDistance(1)),find(DistanceBetweenCoordinates(IDC,:)==SortDistance(2))); %% Find the closest neighbors
        [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],NeighborCoordinates);
        TransNC=[Row;Col]';
        Neighbors{IDC,1}=TransNC;
    end
    CellofNeighbors{IDN,1}=Neighbors;
end
    