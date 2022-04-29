% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

showPlot='on';

startTime=datetime(2015,7,4,0,5,30);
endTime=datetime(2015,7,4,0,8,0);

meltAlt=3.5; % Estimated altitude of melting layer in km
divAlt=6; % Estimated altitude of divergence level in km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir='/scr/snow2/rsfdata/projects/cset/cloudSat/GEOPROF/hdf/';
dataDir2='/scr/snow2/rsfdata/projects/cset/cloudSat/PRECIP-COLUMN/hdf/';

infile2='2015184230321_48842_CS_2C-PRECIP-COLUMN_GRANULE_P1_R05_E06_F00.hdf';

figdir=['/scr/sci/romatsch/other/convStratPaperHCR/'];

allFiles=dir([dataDir,'*.hdf']);

EccoCS1Dall=[];

%% Loop through files
for ii=1:length(allFiles)

    infile=allFiles(ii).name;
    infile2=[infile(1:19),'_CS_2C-PRECIP-COLUMN_GRANULE_P1_R05_E06_F00.hdf'];

    disp(infile);

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
    
    %% Prepare data

    % Flag non-cloud echo
    data.DBZ(data.FLAG<30)=nan;
    data.DBZ(106:end,:)=[];

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

    pixRadDBZ=10; % Radius over which texture is calculated in pixels. Default is 50.
    dbzBase=-10; % Reflectivity base value which is subtracted from DBZ.

    dbzText=f_reflTexture(data.DBZ,pixRadDBZ,dbzBase);

    %% Convectivity

    % Convectivity
    upperLimDBZ=25;
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
end

%% Prepare for plotting

convCS=EccoCS1Dall(EccoCS1Dall(:,2)==1,:);
stratCS=EccoCS1Dall(EccoCS1Dall(:,2)==2,:);
shallowCS=EccoCS1Dall(EccoCS1Dall(:,2)==3,:);

edges=0:10;
convDist=histcounts(convCS(:,1));
stratDist=histcounts(stratCS(:,1));
shallowDist=histcounts(shallowCS(:,1));

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
b=bar(convDist./size(convCS,1).*100,1);
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
text(1,32,'(a) CloudSat convective','FontSize',11,'FontWeight','bold');

s2=subplot(3,1,2);

hold on
b=bar(stratDist./size(stratCS,1).*100,1);
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
text(1,32,'(b) CloudSat stratiform','FontSize',11,'FontWeight','bold');

s3=subplot(3,1,3);

hold on
b=bar(shallowDist./size(shallowCS,1).*100,1);
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
text(1,32,'(c) CloudSat shallow','FontSize',11,'FontWeight','bold');


% s1.Position=[0.049 0.715 0.82 0.27];
% s2.Position=[0.049 0.432 0.82 0.27];
% s4.Position=[0.049 0.15 0.82 0.27];

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'cloudSatStats.png'],'-dpng','-r0')
