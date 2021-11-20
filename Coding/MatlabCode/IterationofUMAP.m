function [SubmembraneCell,DevelopmentInfoCell]=IterationofUMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,RealGrainedN,IDSeg)
UProduction=0.8; %% The production rate of UMP2
BindingR=0.4; %% The binding rate
UnBindingR=0.3; %% The binding rate
UDecay=0.0003; %% The decay rate of UMP2
IDonS=4; %% The information is save on which column of the synaptic components
IDonM=7; %% The information is save on which column of the soma membrane

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
% The 10th column saves the historical growth direction; The 11th column 
% saves the previous segment; The
% 12th column saves the next segment (or segments if there
% is branching); The 13th column saves the centrifugal order 
        % The 14th columns save temporary information
%% Update the Tubulin concentratio on the synapse (if there is any synapse)
if IDSeg==1 %% This is the root segment
    LaplaceS2=(SubmembraneCell{IDN,IDonM}(Row,Col)+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS})/2; %% %% This is the solution of the Laplace equation by the relaxational method
    Decay=UDecay*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the decay term
    SubmembraneCell{IDN,IDonM}(Row,Col)=LaplaceS2+UProduction-Decay; %% Update the tubulin concentration on shpere
    if isempty(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}) %% There is no next segment
        LaplaceS2=(SubmembraneCell{IDN,IDonM}(Row,Col)+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS})/2; %% %% This is the solution of the Laplace equation by the relaxational method
        Decay=UDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
        Binding=BindingR*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the binding term
        Unbinding=UnBindingR*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,5}; %% This is the un-binding term
        DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS2-Decay-Binding+Unbinding; %% Update the UMAP concentration on synaptic segment
    else %% There is next segment (or segments)
        NextSegments=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12};
        TNextSegment=zeros(1,length(NextSegments));
        for NS=1:length(NextSegments)
            TNextSegment(NS)=DevelopmentInfoCell{IDN,3}{Row,Col}{NextSegments(NS),IDonS};
        end
        LaplaceS3=(SubmembraneCell{IDN,IDonM}(Row,Col)+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}+sum(TNextSegment))/(2+length(TNextSegment)); %% Update the UMAP concentration on synaptic segment
        Decay=UDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
        DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3-Decay; %% Update the UMAP concentration on synaptic segment
    end
elseif IDSeg>1 %% This is not the root segment
    if isempty(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}) %% There is no next segment
        PreviousSegment=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,11}; %% This is the previous segment
        TNearSegment=zeros(1,length(length(DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,12}))); %% We need to think about the situation where there was a time of branching right before this segment
        for PS=1:length(length(DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,12}))
            TNearSegment(PS)=DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment(PS),IDonS};
        end
        LaplaceS2=(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}+sum(TNearSegment))/(length(TNearSegment)+1); %% %% This is the solution of the Laplace equation by the relaxational method
        Decay=UDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
        Binding=BindingR*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the binding term
        Unbinding=UnBindingR*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,5}; %% This is the un-binding term
        DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS2-Decay-Binding+Unbinding; %% Update the tubulin concentration on synaptic segment
    else %% There is next segment (or segments)
        PreviousSegment=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,11}; %% This is the previous segment
        TNearSegment=zeros(1,length(length(DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,12}))); %% We need to think about the situation where there was a time of branching right before this segment
        for PS=1:length(length(DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,12}))
            TNearSegment(PS)=DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment(PS),IDonS};
        end
        NextSegments=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}; %% We need to think about the situation where there is a time of branching right after this segment
        TNextSegment=zeros(1,length(NextSegments));
        for NS=1:length(NextSegments)
            TNextSegment(NS)=DevelopmentInfoCell{IDN,3}{Row,Col}{NextSegments(NS),IDonS};
        end
        LaplaceS3=(sum(TNearSegment)+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}+sum(TNextSegment))/(1+length(TNearSegment)+length(TNextSegment)); %% Update the UMAP concentration on synaptic segment
        Decay=UDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
        DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3-Decay; %% Update the UMAP concentration on synaptic segment
    end
end
        

