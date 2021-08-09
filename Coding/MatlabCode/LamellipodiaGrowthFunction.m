function [DevelopmentInfoCell]=LamellipodiaGrowthFunction(SubmembraneCell,DevelopmentInfoCell,InitialPar,GrainedN)
MaxGrowingProb=0.00002; %% The maximum growing probability per hour
Alpha=1;
Beta=2;
RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
for IDN=1:InitialPar.NumberofNeurons
    for IDC=1:RealGrainedN*RealGrainedN %% Traverse every coordinate on the soma shpere
        [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate 
        CoordinateCa=SubmembraneCell{IDN,5}(Row,Col);  %% This is the submembrane Ca concentration on this coordinate
        MaxCa=max(max(SubmembraneCell{IDN,5})); %% This is the maximum submembrane Ca concentration at this moment
        GrowingProb=MaxGrowingProb*((Beta/(Beta-Alpha))*(CoordinateCa/MaxCa)^Alpha-(Alpha/(Beta-Alpha))*(CoordinateCa/MaxCa)^Beta); %% The growth probability on this coordinate
        DevelopmentInfoCell{IDN,1}(Row,Col)=GrowingProb;
        NonGrowingProb=(1-GrowingProb)*(1-GrowingProb>0); %% The probability of remaining stationary on this coordinate
        RandomizeBase=binornd(1,GrowingProb); %% 1 means growth, 0 means non-growth
        DevelopmentInfoCell{IDN,2}(Row,Col)=DevelopmentInfoCell{IDN,2}(Row,Col)+RandomizeBase>0; %% Randomize growth or non-growth
    end
end