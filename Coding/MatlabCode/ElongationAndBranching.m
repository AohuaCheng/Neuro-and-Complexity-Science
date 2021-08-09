function [SubmembraneCell,DevelopmentInfoCell]=ElongationAndBranching(LocationMatrix,SubmembraneCell,DevelopmentInfoCell,SourceLocationMatrix,IDN,IDC,RealGrainedN,IDSeg,IDT)
TAssembly=0.1; %% Tubulin assembly rate
TDisassembly=0.0005; %% The disassembly rate of tubulin
KB=0.01; %% The branching constant
E=0.5; %% The competition parameter in branching
S=-0.1; %% Order dependency in branching
RotationAI=[-pi/6,pi/6]; %% Rotation angle interval in branching
DisturbanceAI=[-pi/8,pi/8]; %% Rotation angle interval in non-branching
InitialR=2; %% The radius of root segement
DecayofR=0.8; %% The decreasing rate of the radius of segements
InheritanceR=0.9; %% The nheritance rate
RecordThreshold=10; %% The minimum elongation record threshold

EnlongationR=1000;  %% Enlongation rate
BranchingR=500;  %% Branching rate

[Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
%% Update the Tubulin concentratio on the soma sphere
% In the cell "DevelopmentInfoCell{IDN,3}{Col(IDS),Raw(IDS)}{xxxxx,yyyyy}",  the 1st column saves the coordinates of
% synaptic components; The 2nd column saves the calcium on synaptic
    % components; The 3rd column saves the Tubulin concentration on
    % synaptic components; The 4th column saves UMAP-2 concentration on
    % synaptic components; The 5th column saves BMAP-2 concentration on
    % synaptic components; The 6th column saves PMAP-2 concentration on
    % synaptic components; The 7th column saves the elongation length on
    % synaptic components; The 8th column saves the branching
    % probability; The 9th column saves the radius of synapse;
    % The 10th column saves the historical growth direction; 
    % The 11th column saves the previous segment; The
    % 12th column saves the next segment (or segments if there
    % is branching); The 13th column saves the centrifugal order 
            % The 14th columns save temporary information
%% Update the Tubulin concentratio on the synapse (if there is any synapse)
if isempty(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}) %% There is no next segment
    EnL=EnlongationR*TAssembly*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,3}*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,5}-TDisassembly; %% This is the enlongation length
    if EnL>=RecordThreshold
        DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,7}=EnL; %% Update the  enlongation length
        BaseBranchingProb=KB*(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,6}/(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,6}+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,5})); %% This is the baseline branching probability
        NumofTS=DevelopmentInfoCell{IDN,4}{Row,Col}(IDT-1,3); %% Number of terminal segments
        Normal=DevelopmentInfoCell{IDN,4}{Row,Col}(IDT-1,4); %% Normalization base of centrifugal order
        CTerm=1/NumofTS*Normal; %% Normalization term of centrifugal order
        BranchingProb=BranchingR*BaseBranchingProb*NumofTS^(-E)*2^(-S*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,13})/CTerm; %% This is the branching probability
        if ~isnan(BranchingProb)
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,8}=BranchingProb;
        end
        BranchingOrNot=binornd(1,DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,8}); %% Generate a random variable for branching or non-branching
        if BranchingOrNot==1 %% Branching
            if IDSeg==1 %% This is the root segment
                MembranceC=[SubmembraneCell{IDN,1}(Row,Col),SubmembraneCell{IDN,2}(Row,Col),SubmembraneCell{IDN,3}(Row,Col)]; %% This is the coordinate on membrane
                BaseVector=MembranceC-LocationMatrix(IDN,:); %% This is the base direction from the soma center
                BaseVector=BaseVector/norm(BaseVector)*EnL; %% Here we normalize the direction vector and multiply the enlongation length
                LVector=Rotation3D(BaseVector,RotationAI(1)+(RotationAI(2)-RotationAI(1))*rand(1,3),MembranceC); %% Left branching daughter segment
                RVector=Rotation3D(BaseVector,RotationAI(1)+(RotationAI(2)-RotationAI(1))*rand(1,3),MembranceC); %% Right branching daughter segment
                LDirection=(LVector-MembranceC)/norm((LVector-MembranceC)); %% The corresponding direction of left segment
                RDirection=(RVector-MembranceC)/norm((RVector-MembranceC)); %% The corresponding direction of right segment
                DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,10}=BaseVector/norm(BaseVector); %% Here we save the historical direction information
                DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,1}=MembranceC; %% Here we save the coordinate of the root
                LSegment=cell(1,14); %% This is the cell to save the information of the left branching daughter segment
                RSegment=cell(1,14); %% This is the cell to save the information of the right branching daughter segment
                for IDL=2:4
                    LSegment{1,IDL}=SubmembraneCell{IDN,IDL+3}(Row,Col)*InheritanceR*1/2; %% Initialization of chemical substance
                    RSegment{1,IDL}=SubmembraneCell{IDN,IDL+3}(Row,Col)*InheritanceR*1/2; %% Initialization of chemical substance
                end
                for IDL=5:6
                    LSegment{1,IDL}=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDL}*InheritanceR*1/2; %% Initialization of chemical substance
                    RSegment{1,IDL}=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDL}*InheritanceR*1/2; %% Initialization of chemical substance
                end
            elseif IDSeg>1 %% This is not the root segment
                HistoryD=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,10}; %% This is the history direction
                CurrentC=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,1}; %% This is the current coordinate
                BaseVector=HistoryD*EnL; %% Here we normalize the direction vector and multiply the enlongation length
                LVector=Rotation3D(BaseVector,RotationAI(1)+(RotationAI(2)-RotationAI(1))*rand(1,3),CurrentC); %% Left branching daughter segment
                RVector=Rotation3D(BaseVector,RotationAI(1)+(RotationAI(2)-RotationAI(1))*rand(1,3),CurrentC); %% Right branching daughter segment
                LDirection=(LVector-CurrentC)/norm((LVector-CurrentC)); %% The corresponding direction of left segment
                RDirection=(RVector-CurrentC)/norm((RVector-CurrentC)); %% The corresponding direction of right segment
                LSegment=cell(1,14); %% This is the cell to save the information of the left branching daughter segment
                RSegment=cell(1,14); %% This is the cell to save the information of the right branching daughter segment
                for IDL=2:6
                    LSegment{1,IDL}=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDL}*InheritanceR*1/2; %% Initialization of chemical substance
                    RSegment{1,IDL}=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDL}*InheritanceR*1/2; %% Initialization of chemical substance
                end
            end
            LSegment{1,1}=LVector; %% The coordinate
            RSegment{1,1}=RVector; %% The coordinate
            LSegment{1,7}=0; 
            RSegment{1,7}=0; 
            LSegment{1,8}=0; 
            RSegment{1,8}=0; 
            LSegment{1,10}=LDirection; %% Here we save the historical direction information
            RSegment{1,10}=RDirection; %% Here we save the historical direction information
            LSegment{1,11}=IDSeg; %% The information of "previous segment"
            RSegment{1,11}=IDSeg; %% The information of "previous segment"
            LSegment{1,13}=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,13}+1; %% The centrifugal order
            RSegment{1,13}=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,13}+1; %% The centrifugal order
            LSegment{1,14}=IDT; %% The information of time
            RSegment{1,14}=IDT; %% The information of time
            LSegment{1,9}=InitialR*DecayofR^LSegment{1,13}; %% The radius of segment
            RSegment{1,9}=InitialR*DecayofR^RSegment{1,13}; %% The radius of segment
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}=[size(DevelopmentInfoCell{IDN,3}{Row,Col},1)+1,size(DevelopmentInfoCell{IDN,3}{Row,Col},1)+2]; %% Update the information of "next segment"
            DevelopmentInfoCell{IDN,3}{Row,Col}=[DevelopmentInfoCell{IDN,3}{Row,Col};LSegment;RSegment]; %% Add these two new segments
            % Here DevelopmentInfoCell{ID,5} is used to save some information to calculate DevelopmentInfoCell{ID,4}
            InformationToRecord=[EnL*ones(2,1),[2^(-S*LSegment{1,13});2^(-S*RSegment{1,13})]];
            DevelopmentInfoCell{IDN,5}{Row,Col}=[DevelopmentInfoCell{IDN,5}{Row,Col};InformationToRecord];
        else %% Non-branching
            if IDSeg==1 %% This is the root segment
                DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,7}=EnL; %% Update the  enlongation length
                MembranceC=[SubmembraneCell{IDN,1}(Row,Col),SubmembraneCell{IDN,2}(Row,Col),SubmembraneCell{IDN,3}(Row,Col)]; %% This is the coordinate on membrane
                BaseVector=MembranceC-LocationMatrix(IDN,:); %% This is the base direction from the soma center
                BaseVector=BaseVector/norm(BaseVector)*EnL; %% Here we normalize the direction vector and multiply the enlongation length
                NVector=Rotation3D(BaseVector,DisturbanceAI(1)+(DisturbanceAI(2)-DisturbanceAI(1))*rand(1,3),MembranceC); %% Next segment
                NDirection=(NVector-MembranceC)/norm(NVector-MembranceC); %% Normalization
                [NVector,NDirection]=SynapseGuidancebyNetrinOne(MembranceC,NDirection,EnL,SourceLocationMatrix,IDT); %% Axon guidance
                DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,10}=BaseVector/norm(BaseVector); %% Here we save the historical direction information
                DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,1}=MembranceC; %% Here we save the coordinate of the root
                NSegment=cell(1,14); %% This is the cell to save the information of the left branching daughter segment
                for IDL=2:4
                    NSegment{1,IDL}=SubmembraneCell{IDN,IDL+3}(Row,Col)*InheritanceR; %% Initialization of chemical substance
                end
                for IDL=5:6
                    NSegment{1,IDL}=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDL}*InheritanceR; %% Initialization of chemical substance
                end
            elseif IDSeg>1 %% This is not the root segment
                HistoryD=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,10}; %% This is the history direction
                CurrentC=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,1}; %% This is the current coordinate
                BaseVector=HistoryD*EnL; %% Here we normalize the direction vector and multiply the enlongation length
                NVector=Rotation3D(BaseVector,DisturbanceAI(1)+(DisturbanceAI(2)-DisturbanceAI(1))*rand(1,3),CurrentC); %% Next segment
                NDirection=(NVector-CurrentC)/norm(NVector-CurrentC); %% Normalization
                [NVector,NDirection]=SynapseGuidancebyNetrinOne(CurrentC,NDirection,EnL,SourceLocationMatrix,IDT); %% Axon guidance
                NSegment=cell(1,13); %% This is the cell to save the information of the left branching daughter segment
                for IDL=2:6
                    NSegment{1,IDL}=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDL}*InheritanceR; %% Initialization of chemical substance
                end
            end
            NSegment{1,1}=NVector; %% The coordinate
            NSegment{1,7}=0;
            NSegment{1,8}=0; 
            NSegment{1,10}=(NDirection+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,10})/2; %% Here we save the historical direction information
            NSegment{1,10}=NSegment{1,10}/norm(NSegment{1,10});
            NSegment{1,11}=IDSeg; %% The information of "previous segment"
            NSegment{1,13}=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,13}; %% The centrifugal order
            NSegment{1,14}=IDT; %% The information of time
            NSegment{1,9}=InitialR*DecayofR^NSegment{1,13}; %% The radius of segment
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}=size(DevelopmentInfoCell{IDN,3}{Row,Col},1)+1; %% Update the information of "next segment"
            DevelopmentInfoCell{IDN,3}{Row,Col}=[DevelopmentInfoCell{IDN,3}{Row,Col};NSegment]; %% Add these two new segments
            % Here DevelopmentInfoCell{ID,5} is used to save some information to calculate DevelopmentInfoCell{ID,4}
            InformationToRecord=[EnL,2^(-S*NSegment{1,13})];
            DevelopmentInfoCell{IDN,5}{Row,Col}=[DevelopmentInfoCell{IDN,5}{Row,Col};InformationToRecord];
        end
    elseif (EnL>=0)&(EnL<RecordThreshold)%% We do not need to add new segment, but we need to record the enlongation rate
        % Here DevelopmentInfoCell{ID,5} is used to save some information to calculate DevelopmentInfoCell{ID,4}
            InformationToRecord=[EnL,2^(-S*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,13})];
            DevelopmentInfoCell{IDN,5}{Row,Col}=[DevelopmentInfoCell{IDN,5}{Row,Col};InformationToRecord];
    end
end