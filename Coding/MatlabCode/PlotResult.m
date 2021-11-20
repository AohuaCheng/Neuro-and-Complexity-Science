function PlotResult(CellofInitialSP,DevelopmentInfoCell,InitialPar,GrainedN,SourceLocationMatrix)

if ~isempty(SourceLocationMatrix) %% The Netrin Field is non-homogeneous
    Fig=figure('Color',[1,1,1]);
    hold on;
    for IDN=1:InitialPar.NumberofNeurons
        s=surf(CellofInitialSP{IDN,1},CellofInitialSP{IDN,2},CellofInitialSP{IDN,3},DevelopmentInfoCell{IDN,1});
        s.EdgeColor = 'none';
        RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
        NumCoor=0;
        for IDC=1:RealGrainedN*RealGrainedN %% Traverse every coordinate on the soma shpere
            [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
            if DevelopmentInfoCell{IDN,2}(Row,Col)==1 %% Growth happens on this coordinate
                NumCoor=NumCoor+1;
                BaseNum=length(find(DevelopmentInfoCell{IDN,2}==1));
                for IDSeg=2:size(DevelopmentInfoCell{IDN,3}{Row,Col},1)
                    disp(['Ploting-','Neuron-',num2str(IDN/InitialPar.NumberofNeurons*100),'%','-Coord-',num2str(NumCoor/length(find(DevelopmentInfoCell{IDN,2}==1))*100),'%','-Seg-',num2str(IDSeg/size(DevelopmentInfoCell{IDN,3}{Row,Col},1)*100),'%'])
                    CurrentC=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,1};
                    PreviousC=DevelopmentInfoCell{IDN,3}{Row,Col}{DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,11},1};
                    Ratio=3000/InitialPar.SpaceLimit(1);
                    R=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,9}*Ratio;
                    plot3([PreviousC(1),CurrentC(1)],[PreviousC(2),CurrentC(2)],[PreviousC(3),CurrentC(3)],'Color',[0, 0.8-0.6/BaseNum*NumCoor, 0.2+0.6/BaseNum*NumCoor],'LineWidth',R)
                end
            end
        end
    end
    scatter3(SourceLocationMatrix(:,1),SourceLocationMatrix(:,2),SourceLocationMatrix(:,3),600,'LineWidth',3,'MarkerEdgeColor',[0.3010 0.7450 0.9330])
    view(30,45);
    colormap winter
    TitleNeeded='Multiple Neurons';
    title(TitleNeeded);
    set(gca,'FontSize',36,'Fontname', 'Calibri');
    set(gca,'linewidth',3)
    set(Fig, 'position', get(0,'ScreenSize'));
else %% The Netrin Field is homogeneous
    Fig=figure('Color',[1,1,1]);
    hold on;
    for IDN=1:InitialPar.NumberofNeurons
        s=surf(CellofInitialSP{IDN,1},CellofInitialSP{IDN,2},CellofInitialSP{IDN,3},DevelopmentInfoCell{IDN,1});
        s.EdgeColor = 'none';
        RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
        NumCoor=0;
        for IDC=1:RealGrainedN*RealGrainedN %% Traverse every coordinate on the soma shpere
            [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
            if DevelopmentInfoCell{IDN,2}(Row,Col)==1 %% Growth happens on this coordinate
                NumCoor=NumCoor+1;
                BaseNum=length(find(DevelopmentInfoCell{IDN,2}==1));
                for IDSeg=2:size(DevelopmentInfoCell{IDN,3}{Row,Col},1)
                    disp(['Ploting-','Neuron-',num2str(IDN/InitialPar.NumberofNeurons*100),'%','-Coord-',num2str(NumCoor/length(find(DevelopmentInfoCell{IDN,2}==1))*100),'%','-Seg-',num2str(IDSeg/size(DevelopmentInfoCell{IDN,3}{Row,Col},1)*100),'%'])
                    CurrentC=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,1};
                    PreviousC=DevelopmentInfoCell{IDN,3}{Row,Col}{DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,11},1};
                    Ratio=3000/InitialPar.SpaceLimit(1);
                    R=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,9}*Ratio;
                    plot3([PreviousC(1),CurrentC(1)],[PreviousC(2),CurrentC(2)],[PreviousC(3),CurrentC(3)],'Color',[0, 0.8-0.6/BaseNum*NumCoor, 0.2+0.6/BaseNum*NumCoor],'LineWidth',R)
                end
            end
        end
    end
    view(30,45);
    colormap winter
    TitleNeeded='Multiple Neurons';
    title(TitleNeeded);
    set(gca,'FontSize',36,'Fontname', 'Calibri');
    set(gca,'linewidth',3)
    set(Fig, 'position', get(0,'ScreenSize'));
end