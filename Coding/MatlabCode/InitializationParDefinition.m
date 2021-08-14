function InitialPar=InitializationParDefinition(InitializationType,NeuronTypepar,DurationLength)
if  InitializationType==1
    InitialPar.NumberofNeurons = 1^3; % Here you set the number of neurons, if you choose InitializationType = 1, then NumberofNeurons must be n^3, where n is a natural number
elseif InitializationType==2
    InitialPar.NumberofNeurons = 1^3; % Here you set the number of neurons, if you choose InitializationType = 2, then NumberofNeurons can be any natural number
elseif (InitializationType~=1) && (InitializationType~=2)
    disp('Please Input a Correct Type (Either 1 or 2)')
end
InitialPar.SpaceLimit = [1500,1500,1500]; % Here you set the spacial scale (cubic micrometer). Please note that the metric is micron (mum)
InitialPar.RadiusIntervalofSoma = [40,60]; % Here you set the interval of the radius of each soma. Please note that the metric is micron (mum)
InitialPar.SpaceDiscretization=1; % This is the minimum unit for spatial discretization
InitialPar.DurationLength=DurationLength; % This is the development duration

if NeuronTypepar==1 %% All neurons have excitable membranes
   InitialPar.NeuronType=ones(InitialPar.NumberofNeurons,1);
elseif NeuronTypepar==2 %% All neurons have passive membranes
   InitialPar.NeuronType=ones(InitialPar.NumberofNeurons,1)*2;
elseif NeuronTypepar==3 %% The types of neurons are randomized 
   InitialPar.NeuronType=randi([1,2],InitialPar.NumberofNeurons,1);
end

