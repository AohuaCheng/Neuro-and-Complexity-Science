function [SubmembraneCell]=ChemicalInitialization(GrainedN,InitialPar,CellofInitialSP)
%% Since there are planty neurons distributed in a large space, it is impossible to generate chemical gradient on the whole space.
%% Otherwise, the fluid mechanics will cost large computing power, preventing the subsequent processes to run.
%% Here we assume that the neurons are distributed in well-mixed media with small perturbation only on microscales. This keeps consistency
%% with the real situations in neural systems. Therefore, we mainly deal with the chemical gradient near soma membrane.

RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
CaOut=1; % the external calcium concentration fixed at 1 mM
CaInInitialized=[30,150]*10^(-6); % The internal calcium concentration is initialized at an interval of [30,150] nM, where 1 nM = 10^(-6) mM
TubulinInitialized=[0,10]*10^(-6); % The internal Tubulin concentration is initialized at an interval of [0,10] nM, where 1 nM = 10^(-6) mM
UnboundMAP2Initialized=[0,10]*10^(-6); % The internal UMAP-2 concentration is initialized at an interval of [0,10] nM, where 1 nM = 10^(-6) mM
SubmembraneCell = cell(InitialPar.NumberofNeurons,7);
for ID=1:InitialPar.NumberofNeurons
    SubmembraneCell{ID,1} = CellofInitialSP{ID,1};
    SubmembraneCell{ID,2} = CellofInitialSP{ID,2};
    SubmembraneCell{ID,3} = CellofInitialSP{ID,3};
    SubmembraneCell{ID,4} = CaOut*ones(RealGrainedN,RealGrainedN); % This is the external calcium concentration outside of a point on sphere
    SubmembraneCell{ID,5} = CaInInitialized(1)+(CaInInitialized(2)-CaInInitialized(1))*rand(RealGrainedN,RealGrainedN); % This is the internal calcium concentration
    SubmembraneCell{ID,6} = TubulinInitialized(1)+(TubulinInitialized(2)-TubulinInitialized(1))*rand(RealGrainedN,RealGrainedN); % This is the internal Tubulin concentration
    SubmembraneCell{ID,7} = UnboundMAP2Initialized(1)+(UnboundMAP2Initialized(2)-UnboundMAP2Initialized(1))*rand(RealGrainedN,RealGrainedN); % This is the internal UMAP-2 concentration
end
