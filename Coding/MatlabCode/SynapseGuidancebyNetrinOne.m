function [NLocation,NDirection]=SynapseGuidancebyNetrinOne(CurrentLocation,CurrentDirection,EnL,SourceLocationMatrix,IDT)
%% 
Diffusion=36000; %% DCa=10^(-7)cm^2/sec=36000 \mum^2/hour
Production=36*10^(-11); %% Production rate is 10^(-7)nM/sec=36*10^(-5)nM/hour=36*10^(-11)mM/hour
DeltaR=20; %% The width of growth cone is set as 20 mum 
MinC=0.01*10^(-6); %% The minimal concentration of guidance is 0.01nM or 0.01*10^(-6)mM
MinP=0.01; %% %% The minimal fractional-change of guidance is 0.01

%% Calculate the field of Netrin-1 and find its fractional-change
if ~isempty(SourceLocationMatrix) %% The Netrin Field is non-homogeneous
    DisMatrix=squareform(pdist([CurrentLocation;SourceLocationMatrix])); %% This is the distance matrix
    DisMatrix=DisMatrix(1,2:end); %% We only need the distance between each coordinate and the location of source
    CMatrix=Production*(4*pi*Diffusion*DisMatrix).^(-1).*erfc(DisMatrix./sqrt(4*Diffusion*IDT)); %% This is the concentration of Netrin-1
    PMatrix=-DeltaR*DisMatrix.^(-1).*(1+DisMatrix./sqrt(pi*Diffusion*IDT).*exp(-DisMatrix.^2./(4*Diffusion*IDT))./erfc(DisMatrix./sqrt(4*Diffusion*IDT))); %% This is the matrix of DeltaC/C
    %% Decide whether the guidance is possible
    PossibleGuidance=find((CMatrix>=MinC)&(PMatrix>=MinP)); %% Find the location where the gudiance is possible 
    if ~isempty(PossibleGuidance) %% The guidance is possible
        Direction=(SourceLocationMatrix(PossibleGuidance,:)-CurrentLocation);   
        for IDD=1:size(Direction,1)
            Direction(IDD,:)=Direction(IDD,:)/norm(Direction(IDD,:)); %% Normalization
        end
        NDirection=(CurrentDirection+sum(Direction))/(1+size(Direction,1));  %% Guidance
        NDirection=NDirection/norm(NDirection); %% Normalization
        NLocation=EnL*NDirection+CurrentLocation;  %% Next location
    else %% The guidance is possible
        NDirection=CurrentDirection;
        NLocation=EnL*CurrentDirection+CurrentLocation; %% Next location
    end
else %% The Netrin Field is homogeneous
    NDirection=CurrentDirection;
    NLocation=EnL*CurrentDirection+CurrentLocation; %% Next location
end
