figure('Color',[1,1,1]);
for IDN=1
    disp(IDN)
RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
NumCoor=0;
for IDC=1:RealGrainedN*RealGrainedN %% Traverse every coordinate on the soma shpere
        [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
        if DevelopmentInfoCell{IDN,2}(Row,Col)==1 %% Growth happens on this coordinate
           NumCoor=NumCoor+1;
           GrowthStartM=DevelopmentInfoCell{IDN,3}{Row,Col}{2,14}; %% This is the moment when the growth start
           MaximumCO=zeros(1,DurationLength);
           for IDT=1:DurationLength
               GrowthLength=cell2mat(DevelopmentInfoCell{IDN,3}{Row,Col}(:,7));
               GrowthTime=cell2mat(DevelopmentInfoCell{IDN,3}{Row,Col}(:,14));
               HistoryInfo=GrowthLength(find(GrowthTime<=IDT));
               if ~isempty(HistoryInfo)
                   MaximumCO(IDT)=sum(HistoryInfo);
               end
           end
           subplot(2,3,IDN)
           hold on;
           plot(smooth(smooth(MaximumCO)),'LineWidth',3,'Color',[0, 0.2+0.6/length(find(DevelopmentInfoCell{IDN,2}==1))*NumCoor, 0.8-0.6/length(find(DevelopmentInfoCell{IDN,2}==1))*NumCoor])
           TitleNeeded=['SL of Neuron-',num2str(IDN)];
           title(TitleNeeded);
           xlabel('Hours')
           ylabel('\mum')
           set(gca,'FontSize',36,'Fontname', 'Calibri');
        end
end
end

figure('Color',[1,1,1]);
for IDN=1:6
    disp(IDN)
RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
NumCoor=0;
for IDC=1:RealGrainedN*RealGrainedN %% Traverse every coordinate on the soma shpere
        [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
        if DevelopmentInfoCell{IDN,2}(Row,Col)==1 %% Growth happens on this coordinate
            NumCoor=NumCoor+1;
           subplot(2,3,IDN)
           hold on;
           h=histfit(cell2mat(DevelopmentInfoCell{IDN,3}{Row,Col}(:,13)),5);
           h(1).FaceColor = [0, 0.2+0.6/length(find(DevelopmentInfoCell{IDN,2}==1))*NumCoor, 0.8-0.6/length(find(DevelopmentInfoCell{IDN,2}==1))*NumCoor];
           h(1).FaceAlpha=0.4;
           h(1).LineStyle='none';
           h(2).LineWidth = 3;
           h(2).Color = [.2 .2 .2];
           TitleNeeded=['CO of Neuron-',num2str(IDN)];
           title(TitleNeeded);
           set(gca,'FontSize',36,'Fontname', 'Calibri');
        end
end
end

figure('Color',[1,1,1]);
for IDN=1:6
    disp(IDN)
RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
NumCoor=0;
for IDC=1:RealGrainedN*RealGrainedN %% Traverse every coordinate on the soma shpere
        [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
        if DevelopmentInfoCell{IDN,2}(Row,Col)==1 %% Growth happens on this coordinate
           NumCoor=NumCoor+1;
           subplot(2,3,IDN)
           hold on;
           plot(DevelopmentInfoCell{IDN,4}{Row,Col}(:,1),'LineWidth',3,'Color',[0, 0.2+0.6/length(find(DevelopmentInfoCell{IDN,2}==1))*NumCoor, 0.8-0.6/length(find(DevelopmentInfoCell{IDN,2}==1))*NumCoor])
           TitleNeeded=['GR of Neuron-',num2str(IDN)];
           title(TitleNeeded);
           xlabel('Hours')
           ylabel('\mum/Hours')
           set(gca,'FontSize',36,'Fontname', 'Calibri');
        end
end
end
