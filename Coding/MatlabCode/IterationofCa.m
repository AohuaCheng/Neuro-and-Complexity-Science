function [SubmembraneCell,DevelopmentInfoCell]=IterationofCa(SubmembraneCell,DevelopmentInfoCell,CellofRealRadius,InitialPar,CellofNeighbors,IDN,IDC,RealGrainedN,IDSeg)
SubmembraneCaDecayInterval=[-10,-50]*10^(-6); % The decay of submembrane calcium concentration is [-100,-500] nM, where 1 nM = 10^(-6) mM
KappaCa=0.1; %% Membrane permeability of Ca
DCa=36000; %% DCa=10^(-7)cm^2/sec=36000 \mum^2/hour
OrderDifference=10^(-4); %% The order difference between the external and submembrane Ca concentration
CATransport=0.3; %% The active transport rate of Ca
CDecay=0.01; %% The decay rate of Ca
IDonS=2; %% The information is save on which column of the synaptic components
IDonM=5; %% The information is save on which column of the soma membrane

SingleNeuronType=InitialPar.NeuronType(IDN);
Neighbors=CellofNeighbors{IDN,1}; % Load the neighbor information of each coordinate on the soma sphere
if SingleNeuronType==1 %% Excitable Membrane
    [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
    if IDSeg==1 %% Growth happens on this coordinate, this is the root segment
        ExternalCa=SubmembraneCell{IDN,4}(Row,Col); %% This is the external Ca concentration closed to this coordinate
        NeighborCa=SubmembraneCell{IDN,IDonM}(Neighbors{IDC,1}(:,1),Neighbors{IDC,1}(:,2)); %% This is the neighbor submembrane Ca concentrations
        CoordinateCa=SubmembraneCell{IDN,IDonM}(Row,Col);  %% This is the submembrane Ca concentration on this coordinate
        CaDiftoOther=SubmembraneCaDecayInterval(1)+(SubmembraneCaDecayInterval(2)-SubmembraneCaDecayInterval(1))*rand(1,1); %% This is the decay of Ca concentration to other parts of cytoplasm (far from membrane)
        LaplaceS2=(CoordinateCa+CaDiftoOther+sum(sum(NeighborCa))+OrderDifference*ExternalCa+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS})/(1+1+1+numel(NeighborCa)); %% %% This is the solution of the Laplace equation by the relaxational method
        ActiveTransport=CATransport*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the active transport term
        Decay=CDecay*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the decay term
        SubmembraneCell{IDN,IDonM}(Row,Col)=LaplaceS2-ActiveTransport-Decay; %% Update the Ca concentration on shpere
        if isempty(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}) %% There is no next segment
            LaplaceS3=(SubmembraneCell{IDN,IDonM}(Row,Col)+OrderDifference*ExternalCa+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS})/3; %% %% This is the solution of the Laplace equation by the relaxational method
            ActiveTransport=CATransport*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the active transport term
            Decay=CDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3+ActiveTransport-Decay; %% Update the tubulin concentration on synaptic segment
        else %% There is next segment (or segments)
            NextSegments=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12};
            TNextSegment=zeros(1,length(NextSegments));
            for NS=1:length(NextSegments)
                TNextSegment(NS)=DevelopmentInfoCell{IDN,3}{Row,Col}{NextSegments(NS),IDonS};
            end
            LaplaceS3=(SubmembraneCell{IDN,IDonM}(Row,Col)+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}+sum(TNextSegment)+OrderDifference*ExternalCa)/(3+length(TNextSegment)); %% Update the tubulin concentration on synaptic segment
            ActiveTransportin=CATransport*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the active transport-in term
            ActiveTransportout=CATransport*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the active transport-out term
            Decay=CDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay; %% Update the tubulin concentration on synaptic segment
        end
    elseif IDSeg>1 %% Growth happens on this coordinate, this is not the root segment
        ExternalCa=SubmembraneCell{IDN,4}(Row,Col); %% This is the external Ca concentration closed to this coordinate
        if isempty(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}) %% There is no next segment
            PreviousSegment=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,11}; %% This is the previous segment
            TNearSegment=zeros(1,length(length(DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,12}))); %% We need to think about the situation where there was a time of branching right before this segment
            for PS=1:length(length(DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,12}))
                TNearSegment(PS)=DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment(PS),IDonS};
            end
            LaplaceS3=(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}+sum(TNearSegment)+OrderDifference*ExternalCa)/(length(TNearSegment)+2); %% %% This is the solution of the Laplace equation by the relaxational method
            ActiveTransport=CATransport*DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,IDonS}; %% This is the active transport term
            Decay=CDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3+ActiveTransport-Decay; %% Update the tubulin concentration on synaptic segment
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
            LaplaceS3=(sum(TNearSegment)+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}+sum(TNextSegment)+OrderDifference*ExternalCa)/(2+length(TNearSegment)+length(TNextSegment)); %% Update the tubulin concentration on synaptic segment
            ActiveTransportin=1/length(TNearSegment)*CATransport*DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,IDonS}*DevelopmentInfoCell{IDN,4}{Row,Col}(end,1); %% This is the active transport-in term
            ActiveTransportout=CATransport*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}*DevelopmentInfoCell{IDN,4}{Row,Col}(end,1); %% This is the active transport-out term
            Decay=CDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay; %% Update the tubulin concentration on synaptic segment
        end
    elseif IDSeg<0 %% Growth does not happen on this coordinate
        ExternalCa=SubmembraneCell{IDN,4}(Row,Col); %% This is the external Ca concentration closed to this coordinate
        NeighborCa=SubmembraneCell{IDN,IDonM}(Neighbors{IDC,1}(:,1),Neighbors{IDC,1}(:,2)); %% This is the neighbor submembrane Ca concentrations
        CoordinateCa=SubmembraneCell{IDN,IDonM}(Row,Col);  %% This is the submembrane Ca concentration on this coordinate
        CaDiftoOther=SubmembraneCaDecayInterval(1)+(SubmembraneCaDecayInterval(2)-SubmembraneCaDecayInterval(1))*rand(1,1); %% This is the decay of Ca concentration to other parts of cytoplasm (far from membrane)
        Decay=CDecay*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the decay term
        SubmembraneCell{IDN,IDonM}(Row,Col)=(CoordinateCa+CaDiftoOther+sum(sum(NeighborCa))+OrderDifference*ExternalCa)/(1+1+numel(NeighborCa))-Decay; %% Solution of the Laplace equation by the relaxational method
    end
elseif SingleNeuronType==2 %% Passive
    RealRadiusVector=CellofRealRadius{IDN,1}; % Load the real radius information of each coordinate on the soma sphere
    RealRadius=RealRadiusVector(IDC); %% The real radius of this coor
    [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
    if IDSeg==1 %% %% Growth happens on this coordinate, this is the root segment
        ExternalCa=SubmembraneCell{IDN,4}(Row,Col); %% This is the external Ca concentration closed to this coordinate
        ActiveTransport=CATransport*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the active transport term
        Decay=CDecay*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the decay term
        SubmembraneCell{IDN,IDonM}(Row,Col)=(KappaCa/(KappaCa+DCa/RealRadius))*ExternalCa-ActiveTransport-Decay; %% Update the tubulin concentration on shpere
        if isempty(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}) %% There is no next segment
            LaplaceS3=(SubmembraneCell{IDN,IDonM}(Row,Col)+OrderDifference*ExternalCa+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS})/3; %% %% This is the solution of the Laplace equation by the relaxational method
            ActiveTransport=CATransport*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the active transport term
            Decay=CDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3+ActiveTransport-Decay; %% Update the tubulin concentration on synaptic segment
        else %% There is next segment (or segments)
            NextSegments=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12};
            TNextSegment=zeros(1,length(NextSegments));
            for NS=1:length(NextSegments)
                TNextSegment(NS)=DevelopmentInfoCell{IDN,3}{Row,Col}{NextSegments(NS),IDonS};
            end
            LaplaceS3=(SubmembraneCell{IDN,IDonM}(Row,Col)+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}+sum(TNextSegment)+OrderDifference*ExternalCa)/(3+length(TNextSegment)); %% Update the tubulin concentration on synaptic segment
            ActiveTransportin=CATransport*SubmembraneCell{IDN,IDonM}(Row,Col); %% This is the active transport-in term
            ActiveTransportout=CATransport*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the active transport-out term
            Decay=CDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay; %% Update the tubulin concentration on synaptic segment
        end
    elseif IDSeg>1 %% %% Growth happens on this coordinate, this is not the root segment
        ExternalCa=SubmembraneCell{IDN,4}(Row,Col); %% This is the external Ca concentration closed to this coordinate
        if isempty(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,12}) %% There is no next segment
            PreviousSegment=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,11}; %% This is the previous segment
            TNearSegment=zeros(1,length(length(DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,12}))); %% We need to think about the situation where there was a time of branching right before this segment
            for PS=1:length(length(DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,12}))
                TNearSegment(PS)=DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment(PS),IDonS};
            end
            LaplaceS3=(DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}+sum(TNearSegment)+OrderDifference*ExternalCa)/(length(TNearSegment)+2); %% %% This is the solution of the Laplace equation by the relaxational method
            ActiveTransport=CATransport*DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,IDonS}; %% This is the active transport term
            Decay=CDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3+ActiveTransport-Decay; %% Update the tubulin concentration on synaptic segment
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
            LaplaceS3=(sum(TNearSegment)+DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}+sum(TNextSegment)+OrderDifference*ExternalCa)/(2+length(TNearSegment)+length(TNextSegment)); %% Update the tubulin concentration on synaptic segment
            ActiveTransportin=1/length(TNearSegment)*CATransport*DevelopmentInfoCell{IDN,3}{Row,Col}{PreviousSegment,IDonS}*DevelopmentInfoCell{IDN,4}{Row,Col}(end,1); %% This is the active transport-in term
            ActiveTransportout=CATransport*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}*DevelopmentInfoCell{IDN,4}{Row,Col}(end,1); %% This is the active transport-out term
            Decay=CDecay*DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}; %% This is the decay term
            DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,IDonS}=LaplaceS3+ActiveTransportin-ActiveTransportout-Decay; %% Update the tubulin concentration on synaptic segment
        end
    elseif IDSeg<0 %% Growth does not happen on this coordinate
        ExternalCa=SubmembraneCell{IDN,4}(Row,Col); %% This is the external Ca concentration closed to this coordinate
        Decay=CDecay*SubmembraneCell{IDN,IDonM}(Col,Col); %% This is the decay term
        SubmembraneCell{IDN,IDonM}(Row,Col)=(KappaCa/(KappaCa+DCa/RealRadius))*ExternalCa-Decay; %% Update the tubulin concentration on shpere
    end
end