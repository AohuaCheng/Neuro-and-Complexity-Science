%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter Setting Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose the duration length of neural development
DurationLength=2000; % Here DurationLength is measured in hours. E.g., DurationLength=1000 corresponds to 1000 hours
%% Choose the initialization type
InitializationType=2; % Here you can either choose 1 for grid initialization, or 2 for coordinate initialization
%% Choose the neuron type setting
NeuronTypepar=1; % There are two types of neurons: neurons with excitable membranes and passive membranes respectively. 
% You can set 1 to let all neurons have excitable membranes or set 2 to let all of them have passive membranes. You can also set 3 to randomly generate membranes.
%% Define the granularity for neuron sphere
GrainedN=20;  % Here you need to define "GrainedN" as a natural number. GrainedN=20 is recommended as a MATLAB default. A larger "GrainedN" ensures higher accuracy but requires more computing power
%% Define the initialization parameter
InitialPar=InitializationParDefinition(InitializationType,NeuronTypepar,DurationLength); % Here you will get a struct of initialization parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of space and somas 
[LocationMatrix,RadiusVector,CellofInitialSP,CellofNeighbors,CellofRealRadius]=LocationAndRadiusInitialization(InitializationType,InitialPar,GrainedN); % Here you will get an initialized space
%% Initialization of Chemical Substances
[SubmembraneCell]=ChemicalInitialization(GrainedN,InitialPar,CellofInitialSP); % Here you will get the chemical substances distributions
NetrinFieldType=1;
[SourceLocationMatrix]=NetrinOneFieldInitialization(InitialPar,LocationMatrix,NetrinFieldType); 
%% Initialization of Development Information
[DevelopmentInfoCell]=DevelopInitialization(GrainedN,InitialPar); % Here you initialize the development information

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Growth Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubmembraneCellAcrossTime=cell(DurationLength,1); % This is the cell to save membrane information
SubmembraneCellAcrossTime{1,1}=SubmembraneCell;
for IDT=2:DurationLength
    SubmembraneCell=SubmembraneCellAcrossTime{IDT-1,1};
    %% Stage 1 The Development of Lamellipodia
    disp(['Lamellipodia formation and condensation-',num2str(IDT/DurationLength*100),'%'])
    [DevelopmentInfoCell]=LamellipodiaGrowthFunction(SubmembraneCell,DevelopmentInfoCell,InitialPar,GrainedN); % Here you will update the Lamellipodia formation and condensation on the soma membrane
    %% Stage 2 Growth of Synapses
    disp(['Interation of chemical and development information-',num2str(IDT/DurationLength*100),'%'])
    [SubmembraneCell,DevelopmentInfoCell]=IterationofChemicalSubstance(LocationMatrix,SubmembraneCell,DevelopmentInfoCell,CellofRealRadius,CellofNeighbors,SourceLocationMatrix,InitialPar,GrainedN,IDT); % Here you will update the chemical and development information
    SubmembraneCellAcrossTime{IDT,1}=SubmembraneCell;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotResult(CellofInitialSP,DevelopmentInfoCell,InitialPar,GrainedN)

PlotTimeDependentResult(CellofInitialSP,DevelopmentInfoCell,InitialPar,GrainedN,DurationLength)