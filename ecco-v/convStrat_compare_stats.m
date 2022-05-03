% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

showPlot='on';

startTime=datetime(2015,7,4,0,5,30);
endTime=datetime(2015,7,4,0,8,0);

plotCases=0;

dataDir='/scr/snow2/rsfdata/projects/cset/cloudSat/GEOPROF/hdf/';
dataDir2='/scr/snow2/rsfdata/projects/cset/cloudSat/PRECIP-COLUMN/hdf/';

figdir=['/scr/sci/romatsch/other/convStratPaperHCR/'];

allFiles=dir([dataDir,'*.hdf']);

EccoCS1Dall=[];

colmapSC=[0,0.1,0.6;
            0.38,0.42,0.96;
            0.65,0.74,0.86;
            0.32,0.78,0.59;
            1,0,0;
            1,0,1;
            1,1,0;
            0.99,0.77,0.22;
            0.7,0,0];

%% Loop through files
for ii=1:length(allFiles)

    infile=allFiles(ii).name;
    infile2=[infile(1:19),'_CS_2C-PRECIP-COLUMN_GRANULE_P1_R05_E06_F00.hdf'];

    disp(['File ',num2str(ii),' of ',num2str(length(allFiles)),': ',infile]);

    %% Read data

    data=[];

    yearIn=infile(1:4);
    dayIn=infile(5:7);
    hourIn=infile(8:9);
    minIn=infile(10:11);
    secIn=infile(12:13);

    timeIn=hdfread([dataDir,infile],'Profile_time');

    yearStart=datetime(str2num(yearIn),1,1);
    timeStart=yearStart+days(str2num(dayIn)-1)+hours(str2num(hourIn))+minutes(str2num(minIn))+seconds(str2num(secIn));

    data.time=timeStart+seconds(timeIn{:});

    longitude=hdfread([dataDir,infile],'Longitude');
    data.longitude=double(longitude{:});
    latitude=hdfread([dataDir,infile],'Latitude');
    data.latitude=double(latitude{:});
    binSize=hdfread([dataDir,infile],'Vertical_binsize');
    binSize=binSize{:};

    DBZ=hdfread([dataDir,infile],'Radar_Reflectivity');
    DBZ(DBZ==-8888)=nan;
    data.DBZ=DBZ./100;
    data.DBZ=double(data.DBZ)';
    data.FLAG=hdfread([dataDir,infile],'CPR_Cloud_mask');
    data.FLAG(data.FLAG==-9)=nan;
    data.FLAG=double(data.FLAG');
    TOPO=hdfread([dataDir,infile],'DEM_elevation');
    data.TOPO=double(TOPO{:});
    data.TOPO(data.TOPO==-9999)=0;

    %% Read conv_strat flag

    csFlag=hdfread([dataDir2,infile2],'Conv_strat_flag');
    csFlag=double(csFlag{:});
    csFlag(csFlag<1)=nan;
    meltAlt=hdfread([dataDir2,infile2],'Freezing_level');
    meltAlt=double(meltAlt{:});
    meltAlt(meltAlt<0)=nan;

    divAlt=meltAlt+4;
    %% Prepare data

    % Flag non-cloud echo
    data.DBZ(data.FLAG<30)=nan;
    data.DBZ(106:end,:)=[];

    % Remove data over land
    data.DBZ(:,data.TOPO>0)=nan;

    nanCols=any(~isnan(data.DBZ),1);
    dataCols=sum(nanCols);

    if dataCols==0
        continue
    end

    % Create asl
    data.asl=0:binSize:104*binSize;
    data.asl=repmat(data.asl,length(data.time),1);
    data.asl=double(flipud(data.asl'));

    % Create melting layer
    data.MELTING_LAYER=nan(size(data.DBZ));
    data.MELTING_LAYER(data.asl>=meltAlt.*1000)=20;
    data.MELTING_LAYER(data.asl<meltAlt.*1000)=10;

    % Create fake temperature profile
    data.TEMP=nan(size(data.DBZ));
    data.TEMP(data.asl>=divAlt.*1000)=-30;
    data.TEMP(data.asl<divAlt.*1000)=10;

    ylimUpper=(max(data.asl(~isnan(data.DBZ)))./1000)+0.5;

    %% Texture from reflectivity and velocity

    disp('Calculating reflectivity texture ...');

    pixRadDBZ=3; % Radius over which texture is calculated in pixels. Default is 3.
    dbzBase=-10; % Reflectivity base value which is subtracted from DBZ. Default is -10.

    dbzText=f_reflTexture(data.DBZ,pixRadDBZ,dbzBase);

    %% Convectivity

    % Convectivity
    upperLimDBZ=20;
    convDBZ=1/upperLimDBZ.*dbzText;

    %% Basic classification

    disp('Basic classification ...');

    stratMixed=0.4; % Convectivity boundary between strat and mixed.
    mixedConv=0.5; % Convectivity boundary between mixed and conv.

    classBasic=f_classBasic(convDBZ,stratMixed,mixedConv,data.MELTING_LAYER);

    %% Sub classification

    disp('Sub classification ...');

    classSub=f_classSub(classBasic,data.asl,data.TOPO,data.MELTING_LAYER,data.TEMP);

    classSubPlot=classSub;
    classSubPlot(classSub==14)=1;
    classSubPlot(classSub==16)=2;
    classSubPlot(classSub==18)=3;
    classSubPlot(classSub==25)=4;
    classSubPlot(classSub==30)=5;
    classSubPlot(classSub==32)=6;
    classSubPlot(classSub==34)=7;
    classSubPlot(classSub==36)=8;
    classSubPlot(classSub==38)=9;

    stratConv1D=max(classSubPlot,[],1);

    EccoCS1D=cat(2,stratConv1D',csFlag');
    EccoCS1D(any(isnan(EccoCS1D),2),:)=[];

    EccoCS1Dall=cat(1,EccoCS1Dall,EccoCS1D);

    if plotCases
        diffCols=diff(nanCols);
        startInds=find(diffCols==1);
        endInds=find(diffCols==-1);

        if startInds(1)>endInds(1)
            startInds=cat(2,1,startInds);
        end
        if length(startInds)>length(endInds)
            endInds=cat(2,endInds,size(data.DBZ,2));
        end

        % 1D ecco
        stratConv1D=max(classSubPlot,[],1);
        time1D=data.time(~isnan(stratConv1D));
        stratConv1D=stratConv1D(~isnan(stratConv1D));

        col1D=colmapSC(stratConv1D,:);

        % 1D cloudsat
        csFlag(csFlag<1)=nan;

        time1DcloudSat=data.time(~isnan(csFlag));
        cloudSatFlag1D=csFlag(~isnan(csFlag));

        colmapCloudSat=[1,0,0;
            0,0,1;
            1,1,0];

        col1DcloudSat=colmapCloudSat(cloudSatFlag1D,:);

        for jj=1:length(startInds)
            if endInds(jj)-startInds(jj)>300
                close all

                wi=15;
                hi=10;

                fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
                fig1.PaperPositionMode = 'manual';
                fig1.PaperUnits = 'inches';
                fig1.Units = 'inches';
                fig1.PaperPosition = [0, 0, wi, hi];
                fig1.PaperSize = [wi, hi];
                fig1.Resize = 'off';
                fig1.InvertHardcopy = 'off';

                set(fig1,'color','w');

                s1=subplot(4,1,1);

                colormap jet

                hold on
                surf(data.time(startInds(jj):endInds(jj)),data.asl(:,startInds(jj):endInds(jj))./1000,data.DBZ(:,startInds(jj):endInds(jj)),'edgecolor','none');
                view(2);
                ylabel('Altitude (km)');
                caxis([-35 25]);
                ylim([0 ylimUpper]);
                xlim([data.time(startInds(jj)),data.time(endInds(jj))]);
                set(gca,'XTickLabel',[]);
                cb1=colorbar;
                grid on
                box on

                s2=subplot(4,1,2);

                hold on
                surf(data.time(startInds(jj):endInds(jj)),data.asl(:,startInds(jj):endInds(jj))./1000,convDBZ(:,startInds(jj):endInds(jj)),'edgecolor','none');
                view(2);
                ylabel('Altitude (km)');
                caxis([0 1]);
                ylim([0 ylimUpper]);
                xlim([data.time(startInds(jj)),data.time(endInds(jj))]);
                cb2=colorbar;
                set(gca,'XTickLabel',[]);
                grid on
                box on

                s5=subplot(30,1,28);

                hold on
                scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
                set(gca,'clim',[0,1]);
                set(gca,'YTickLabel',[]);
                set(gca,'XTickLabel',[]);
                s5.Colormap=colmapSC;
                xlim([data.time(startInds(jj)),data.time(endInds(jj))]);
                grid on
                box on

                s6=subplot(30,1,30);

                hold on
                scat1=scatter(time1DcloudSat,ones(size(time1DcloudSat)),10,col1DcloudSat,'filled');
                set(gca,'clim',[0,4]);
                set(gca,'YTickLabel',[]);
                s6.Colormap=colmapCloudSat;
                xlim([data.time(startInds(jj)),data.time(endInds(jj))]);
                grid on
                box on
                text(datetime(2015,7,4,0,8,31),3,'Convective','edgecolor','r','Margin',0.1)
                text(datetime(2015,7,4,0,8,31),0,'Stratiform','edgecolor','b','Margin',0.1)
                text(datetime(2015,7,4,0,8,31),-3,'Shallow','edgecolor','y','Margin',0.1)
               
                s4=subplot(4,1,3);

                hold on
                surf(data.time(startInds(jj):endInds(jj)),data.asl(:,startInds(jj):endInds(jj))./1000,classSubPlot(:,startInds(jj):endInds(jj)),'edgecolor','none');
                view(2);
                ylabel('Altitude (km)');
                caxis([0 10]);
                ylim([0 ylimUpper]);
                xlim([data.time(startInds(jj)),data.time(endInds(jj))]);
                s4.Colormap=colmapSC;
                caxis([0.5 9.5]);
                cb4=colorbar;
                cb4.Ticks=1:9;
                cb4.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
                    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
                set(gca,'XTickLabel',[]);
                grid on
                box on

                s1pos=s1.Position;
                s5pos=s5.Position;
                s5.Position=[s5pos(1),s5pos(2),s1pos(3),s5pos(4)];
                s6pos=s6.Position;
                s6.Position=[s6pos(1),s6pos(2),s1pos(3),s6pos(4)];
             
            end
        end
    end
end

%% Prepare for plotting

convCS=EccoCS1Dall(EccoCS1Dall(:,2)==1,:);
stratCS=EccoCS1Dall(EccoCS1Dall(:,2)==2,:);
shallowCS=EccoCS1Dall(EccoCS1Dall(:,2)==3,:);

edges=0.5:9.5;
convDist=histcounts(convCS(:,1),edges);
stratDist=histcounts(stratCS(:,1),edges);
shallowDist=histcounts(shallowCS(:,1),edges);

colmapSC=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    1,0,0;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    0.7,0,0];

%% Plot
close all

wi=9;
hi=7;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(3,1,1);

hold on
convPerc=convDist./size(convCS,1).*100;
b=bar(convPerc,1);
for jj=1:length(convDist)
    b.FaceColor='flat';
    b.CData(jj,:)=colmapSC(jj,:);
end
ylabel('Percent (%)');
xlim([0.5,9.5]);
set(gca,'XTick',1:9);
set(gca,'XTickLabel',{'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'});
set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on
title(['(a) CloudSat convective, strat: ',num2str(sum(convPerc(1:3))),...
    ', mixed: ',num2str(convPerc(4)),...
    ', conv: ',num2str(sum(convPerc(5:9)))],'FontSize',11,'FontWeight','bold');

s2=subplot(3,1,2);

hold on
stratPerc=stratDist./size(stratCS,1).*100;
b=bar(stratPerc,1);
for jj=1:length(stratDist)
    b.FaceColor='flat';
    b.CData(jj,:)=colmapSC(jj,:);
end
ylabel('Percent (%)');
xlim([0.5,9.5]);
set(gca,'XTick',1:9);
set(gca,'XTickLabel',{'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'});
set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on
title(['(b) CloudSat stratiform, strat: ',num2str(sum(stratPerc(1:3))),...
    ', mixed: ',num2str(stratPerc(4)),...
    ', conv: ',num2str(sum(stratPerc(5:9)))],'FontSize',11,'FontWeight','bold');

s3=subplot(3,1,3);

hold on
shallowPerc=shallowDist./size(shallowCS,1).*100;
b=bar(shallowPerc,1);
for jj=1:length(shallowDist)
    b.FaceColor='flat';
    b.CData(jj,:)=colmapSC(jj,:);
end
ylabel('Percent (%)');
xlim([0.5,9.5]);
set(gca,'XTick',1:9);
set(gca,'XTickLabel',{'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'});
set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on
title(['(c) CloudSat shallow, strat: ',num2str(sum(shallowPerc(1:3))),...
    ', mixed: ',num2str(shallowPerc(4)),...
    ', conv: ',num2str(sum(shallowPerc(5:9)))],'FontSize',11,'FontWeight','bold');


set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'cloudSatStats.png'],'-dpng','-r0')

save([figdir,'cloudSatComparison.mat'],'convPerc','stratPerc','shallowPerc','colmapSC');