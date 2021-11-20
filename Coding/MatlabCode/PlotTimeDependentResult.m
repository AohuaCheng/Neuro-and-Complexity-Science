function PlotTimeDependentResult(CellofInitialSP,DevelopmentInfoCell,InitialPar,GrainedN,DurationLength)

Fig=figure('Color',[1,1,1]);
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(Fig,'JavaFrame');
pause(0.1);					
set(jFrame,'Maximized',1);	
pause(0.1);					
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	

TimeDis=50;
for IDT=0:TimeDis:DurationLength
    disp(IDT)
    hold on;
    for IDN=1:InitialPar.NumberofNeurons
        s=surf(CellofInitialSP{IDN,1},CellofInitialSP{IDN,2},CellofInitialSP{IDN,3},DevelopmentInfoCell{IDN,1});
        s.EdgeColor = 'none';
        RealGrainedN=GrainedN+1; %% When you set GrainedN=x, there are x+1 parts on the sphere of somas
        NumCoor=0;
        for IDC=1:RealGrainedN*RealGrainedN %% Traverse every coordinate on the soma shpere
            [Row,Col]=ind2sub([RealGrainedN,RealGrainedN],IDC); %% Pick a coordinate
            if DevelopmentInfoCell{IDN,2}(Row,Col)==1 %% Growth happens on this coordinate
                GrowthTime=cell2mat(DevelopmentInfoCell{IDN,3}{Row,Col}(:,14)); %% This is the time step of the growth
                NumCoor=NumCoor+1;
                BaseNum=length(find(DevelopmentInfoCell{IDN,2}==1));
                HistoryInfo=find((GrowthTime<=IDT)&(GrowthTime>IDT-1*TimeDis));
                HistoryInfo = setdiff(HistoryInfo,1); 
                if length(HistoryInfo)>1
                    for PreIDSeg=1:length(HistoryInfo)
                        IDSeg=HistoryInfo(PreIDSeg);
%                         disp(['Ploting-','Neuron-',num2str(IDN/InitialPar.NumberofNeurons*100),'%','-Coord-',num2str(NumCoor/length(find(DevelopmentInfoCell{IDN,2}==1))*100),'%','-Seg-',num2str(IDSeg/size(DevelopmentInfoCell{IDN,3}{Row,Col},1)*100),'%'])
                        CurrentC=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,1};
                        PreviousC=DevelopmentInfoCell{IDN,3}{Row,Col}{DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,11},1};
                        Ratio=3000/InitialPar.SpaceLimit(1);
                        R=DevelopmentInfoCell{IDN,3}{Row,Col}{IDSeg,9}*Ratio;
                        plot3([PreviousC(1),CurrentC(1)],[PreviousC(2),CurrentC(2)],[PreviousC(3),CurrentC(3)],'Color',[0, 0.8-0.6/BaseNum*NumCoor, 0.2+0.6/BaseNum*NumCoor],'LineWidth',R)
                    end
                end
            end
        end
    end
    hold off;
    view(30,45);
    xlim([0,6000])
    ylim([0,6000])
    zlim([0,6000])
    xlabel('\mum')
    ylabel('\mum')
    zlabel('\mum')
    colormap winter
    TitleNeeded=['Multiple Neurons at ',num2str(IDT),' Hours'];
    title(TitleNeeded);
    set(gca,'FontSize',36,'Fontname', 'Calibri');
    set(gca,'linewidth',3)
    drawnow
    frame = getframe(Fig); 
    Im = frame2im(frame); 
    [Imind,Cm] = rgb2ind(Im,256);
    filename = 'D:\Files\InformationDynamics\Amadeus\ResultImage\Gif\Multi-1.gif';
    if IDT == 0       
        imwrite(Imind,Cm,filename,'gif','WriteMode','overwrite', 'Loopcount',inf);
    else 
        imwrite(Imind,Cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end
end