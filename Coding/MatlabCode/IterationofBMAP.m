function [SubmembraneCell,DevelopmentInfoCell]=IterationofBMAP(SubmembraneCell,DevelopmentInfoCell,IDN,IDC,RealGrainedN,IDSeg)
BindingR=0.4; %% The binding rate
UnBindingR=0.3; %% The binding rate
PR=0.4; %% The phosphorylation rate
DR=0.3; %% The dephosphorylation rate
UDecay=0.0003; %% The decay rate of UMP2
IDonS=5; %% The information is save on which column of the synaptic components
KF=0.0000005; %% Calcium rate constant for phosphorylation
KG=0.0000005; %% Calcium rate constant for dephosphorylation

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
    Decay=UDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
    Binding=BindingR*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,4}; %% This is the binding term
    Unbinding=UnBindingR*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the un-binding term
    PRL=(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,2}^2)/(KF+(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,2}^2));  %% Phosphorylation rate limit
    DRL=(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,2}^2)/(KG+(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,2}^2));  %% Dephosphorylation rate limit
    Phosphorylation=PR*PRL*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the phosphorylation term
    Dephosphorylation=DR*DRL*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,6}; %% This is the dephosphorylation term
    DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=Binding-Unbinding-Phosphorylation+Dephosphorylation-Decay; %% Update the BMAP concentration on synaptic segment
end
