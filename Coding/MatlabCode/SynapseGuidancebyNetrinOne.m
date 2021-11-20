function [NLocation,NDirection]=SynapseGuidancebyNetrinOne(CurrentLocation,CurrentDirection,EnL,SourceLocationMatrix,IDT)
%% 
Diffusion=36000; %% DCa=10^(-7)cm^2/sec=36000 \mum^2/hour
Production=36*10^(-11); %% Production rate is 10^(-7)nM/sec=36*10^(-5)nM/hour=36*10^(-11)mM/hour
DeltaR=EnL; %% The width of growth cone is set as 20 mum 
MinC=0.01*10^(-15); %% The minimal concentration required for guidance is 0.05*10^(-9)nM or 0.05*10^(-15)mM above the initial concentration
MinP=0.01; %% The minimal fractional-change of guidance is 0.01
GuidanceEffect=2; %% The capacity for Netrin-1 to change the axon path
%% Calculate the field of Netrin-1 and find its fractional-change
if ~isempty(SourceLocationMatrix) %% The Netrin Field is non-homogeneous
    DisMatrix=squareform(pdist([CurrentLocation;SourceLocationMatrix])); %% This is the distance matrix
    DisMatrix=DisMatrix(1,2:end); %% We only need the distance between each coordinate and the location of source
    CMatrix=Production*(4*pi*Diffusion*DisMatrix).^(-1).*erfc(DisMatrix./sqrt(4*Diffusion*IDT)); %% Solve the equation of Netrin-1 concentration
    Concen=sum(CMatrix); %% This is the concentration of Netrin-1 at current location
    PMatrix=(SourceLocationMatrix-CurrentLocation)'.*DeltaR*Production.*((4*pi*Diffusion.*DisMatrix.^3).^(-1).*erfc(DisMatrix./sqrt(4*Diffusion*IDT))+(4*sqrt(pi^3*Diffusion^3*IDT)*DisMatrix.^2).^(-1).*exp(-1*DisMatrix.^2/(4*Diffusion*IDT)))/Concen; %% Solve the equation of DeltaC/C
    PGradient=sum(PMatrix,2)'; %% This is the DeltaC/C
    %% Decide whether the guidance is possible
    if (Concen>=MinC)&&(norm(PGradient)>=MinP) %% The guidance is possible
%         disp('Guidance!')
        Direction=PGradient/norm(PGradient);   
        NDirection=((1-GuidanceEffect*norm(PGradient))*CurrentDirection+GuidanceEffect*norm(PGradient)*Direction)/2;  %% Guidance
        NDirection=NDirection/norm(NDirection); %% Normalization of Guidance
        NLocation=EnL*NDirection+CurrentLocation;  %% Next location
    else %% The guidance is impossible
        NDirection=CurrentDirection;
        NLocation=EnL*CurrentDirection+CurrentLocation; %% Next location
    end
else %% The Netrin Field is homogeneous
    NDirection=CurrentDirection;
    NLocation=EnL*CurrentDirection+CurrentLocation; %% Next location
end
