function [SourceLocationMatrix]=NetrinOneFieldInitialization(InitialPar,LocationMatrix,NetrinFieldType)
if NetrinFieldType==1 %% The Netrin Field is homogeneous everywhere
    SourceLocationMatrix=[]; %% There is no source
elseif NetrinFieldType==2 %% The Netrin Field is non-homogeneous
    NumofNetrinSource=100; %% Number of netrin source
    MiniDistanceBN=150; %% This is the minimum distance between any source and neurons
    Min=zeros(1,3); % This is a vector to store the minimum X,Y,Z coordinates
    Max=zeros(1,3); % This is a vector to store the maximum X,Y,Z coordinates
    for IDDim=1:3
        Min(IDDim)=0.1*InitialPar.SpaceLimit(IDDim);
        Max(IDDim)=0.9*InitialPar.SpaceLimit(IDDim);
    end
    SourceLocationMatrix=[Min(1)+(Max(1)-Min(1))*rand(1,1),Min(2)+(Max(2)-Min(2))*rand(1,1),Min(3)+(Max(3)-Min(3))*rand(1,1)]; %% Initialize one neuron's coordinates
    while size(SourceLocationMatrix,1)<NumofNetrinSource %% while there exist sources without coordinates
        NCoord=[Min(1)+(Max(1)-Min(1))*rand(1,1),Min(2)+(Max(2)-Min(2))*rand(1,1),Min(3)+(Max(3)-Min(3))*rand(1,1)]; %% Generate one coordinate
        DistanceM=pdist([NCoord;LocationMatrix]); %% Calculate the distance between somas
        while any(DistanceM<MiniDistanceBN)
            NCoord=[Min(1)+(Max(1)-Min(1))*rand(1,1),Min(2)+(Max(2)-Min(2))*rand(1,1),Min(3)+(Max(3)-Min(3))*rand(1,1)];
        end
        disp(['Initialize Netrin-1 Source-',num2str(size(SourceLocationMatrix,1)/NumofNetrinSource)])
        SourceLocationMatrix=[SourceLocationMatrix;NCoord];
    end
end


  