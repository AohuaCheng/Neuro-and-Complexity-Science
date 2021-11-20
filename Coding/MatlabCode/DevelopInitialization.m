function [DevelopmentInfoCell]=DevelopInitialization(GrainedN,InitialPar)
RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
DevelopmentInfoCell = cell(InitialPar.NumberofNeurons,5);
for ID=1:InitialPar.NumberofNeurons
    DevelopmentInfoCell{ID,1} = zeros(RealGrainedN,RealGrainedN); %% This is the growth probability information
    DevelopmentInfoCell{ID,2} = zeros(RealGrainedN,RealGrainedN); %% This is the randomized growth information
    DevelopmentInfoCell{ID,3} = cell(RealGrainedN,RealGrainedN); %% This is the synapes information
    DevelopmentInfoCell{ID,4} = cell(RealGrainedN,RealGrainedN); %% This is the information of growth history
    DevelopmentInfoCell{ID,5} = cell(RealGrainedN,RealGrainedN); %% This is the cache information
    for IDC=1:RealGrainedN
        for IDC2=1:RealGrainedN
                DevelopmentInfoCell{ID,3}{IDC,IDC2}=cell(1,14); % Initialization
                % In the above cell, the 1st column saves the coordinates of
                % synaptic components; The 2nd column saves the calcium on synaptic
                % components; The 3rd column saves the Tubulin concentration on
                % synaptic components; The 4th column saves UMAP-2 concentration on
                % synaptic components; The 5th column saves BMAP-2 concentration on
                % synaptic components; The 6th column saves PMAP-2 concentration on
                % synaptic components; The 7th column saves the elongation length on
                % synaptic components; The 8th column saves the branching
                % probability; The 9th column saves the radius of synapse;
                % The 10th column saves the segment type (1 stands for root segment, 
                % 2 stands for intermediate segment, 3 stands for terminal
                % segment); The 11th column saves the previous segment; The
                % 12th column saves the next segment (or segments if there
                % is branching); The 13th column saves the centrifugal
                % order; The 14thcolumn saves temporary information
                for IDC3=2:9
                    DevelopmentInfoCell{ID,3}{IDC,IDC2}{1,IDC3}=0;
                end
                DevelopmentInfoCell{ID,3}{IDC,IDC2}{1,13}=0;
                DevelopmentInfoCell{ID,3}{IDC,IDC2}{1,14}=1;
                DevelopmentInfoCell{ID,4}{IDC,IDC2}=zeros(InitialPar.DurationLength,4); %% the 1th column saves the average elongation rate, the 2nd
                % column saves the variance of elongation rate, the 3rd
                % column saves the number of terminal segments; %% the 4th
                % column saves the normalization base of centrifugal order;
                % Here DevelopmentInfoCell{ID,5} is used to save some
                % information to calculate DevelopmentInfoCell{ID,4}
        end
    end
end