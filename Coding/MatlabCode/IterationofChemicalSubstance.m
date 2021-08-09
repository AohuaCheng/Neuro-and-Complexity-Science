function [SubmembraneCell,DevelopmentInfoCell]=IterationofChemicalSubstance(LocationMatrix,SubmembraneCell,DevelopmentInfoCell,CellofRealRadius,CellofNeighbors,SourceLocationMatrix,InitialPar,GrainedN,IDT)
RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
for IDN=1:InitialPar.NumberofNeurons
    for IDC=1:RealGrainedN*RealGrainedN %% Traverse every coordinate on the soma shpere
        [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
        if DevelopmentInfoCell{IDN,2}(Row,Col)==1 %% Growth happens on this coordinate
            for IDSeg=1:size(DevelopmentInfoCell{IDN,3}{Row,Col},1) %% Here is a small tip: although SubmembraneCell and DevelopmentInfoCell will
                % be updated in this "for-loop", the variable
                % "size(DevelopmentInfoCell{ID,3}{Col(IDS),Raw(IDS)},1)"
                % counts the size of the original DevelopmentInfoCell that
                % send into this "for-loop". Therefore, this variable is
                % constant rather than dynamic
                [SubmembraneCell,DevelopmentInfoCell]=IterationofCa(SubmembraneCell,DevelopmentInfoCell,CellofRealRadius,InitialPar,CellofNeighbors,IDN,IDC,RealGrainedN,IDSeg);
                [SubmembraneCell,DevelopmentInfoCell]=IterationofTubulin(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,RealGrainedN,IDSeg,IDT);
                [SubmembraneCell,DevelopmentInfoCell]=IterationofUMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,RealGrainedN,IDSeg);
                [SubmembraneCell,DevelopmentInfoCell]=IterationofBMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,RealGrainedN,IDSeg);
                [SubmembraneCell,DevelopmentInfoCell]=IterationofPMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,RealGrainedN,IDSeg);
                [SubmembraneCell,DevelopmentInfoCell]=ElongationAndBranching(LocationMatrix,SubmembraneCell,DevelopmentInfoCell,SourceLocationMatrix,IDN,IDC,RealGrainedN,IDSeg,IDT);
            end
            if  ~isempty(DevelopmentInfoCell{IDN,5}{Row,Col}) 
                DevelopmentInfoCell{IDN,4}{Row,Col}(IDT,1)=mean(DevelopmentInfoCell{IDN,5}{Row,Col}(:,1));
                DevelopmentInfoCell{IDN,4}{Row,Col}(IDT,2)=var(DevelopmentInfoCell{IDN,5}{Row,Col}(:,1));
                DevelopmentInfoCell{IDN,4}{Row,Col}(IDT,3)=size(DevelopmentInfoCell{IDN,5}{Row,Col},1);
                DevelopmentInfoCell{IDN,4}{Row,Col}(IDT,4)=sum(DevelopmentInfoCell{IDN,5}{Row,Col}(:,2));
                % the 1th column saves the average elongation rate, the 2nd
                % column saves the variance of elongation rate, the 3rd
                % column saves the number of terminal segments; %% the 4th
                % column saves the normalization base of centrifugal order;
                DevelopmentInfoCell{IDN,5}{Row,Col}=[]; % Clean the information
            end
        else %% Growth does not happen on this coordinate
            [SubmembraneCell,DevelopmentInfoCell]=IterationofCa(SubmembraneCell,DevelopmentInfoCell,CellofRealRadius,InitialPar,CellofNeighbors,IDN,IDC,RealGrainedN,-1);
        end
    end
end