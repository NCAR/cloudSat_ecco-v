% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

showPlot='on';

startTime=datetime(2015,7,4,0,5,30);
endTime=datetime(2015,7,4,0,8,30);

meltAlt=3.5; % Estimated altitude of melting layer in km
divAlt=6; % Estimated altitude of divergence level in km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir='/scr/snow2/rsfdata/projects/cset/cloudSat/hdf/';
infile='2015184230321_48842_CS_2B-GEOPROF_GRANULE_P1_R05_E06_F00.hdf';

figdir=['/scr/sci/romatsch/other/convStratCloudSat/'];

%% Read data

yearIn=infile(1:4);
dayIn=infile(5:7);
hourIn=infile(8:9);
minIn=infile(10:11);
secIn=infile(12:13);

timeIn=hdfread([dataDir,infile],'Profile_time');

yearStart=datetime(str2num(yearIn),1,1);
timeStart=yearStart+days(str2num(dayIn)-1)+hours(str2num(hourIn))+minutes(str2num(minIn))+seconds(str2num(secIn));

timeAll=timeStart+seconds(timeIn{:});
[min1,firstInd]=min(abs(timeAll-startTime));
[min2,lastInd]=min(abs(timeAll-endTime));

longitude=hdfread([dataDir,infile],'Longitude');
longitude=longitude{:};
latitude=hdfread([dataDir,infile],'Latitude');
latitude=latitude{:};
binSize=hdfread([dataDir,infile],'Vertical_binsize');
binSize=binSize{:};

DBZ=hdfread([dataDir,infile],'Radar_Reflectivity');
DBZ(DBZ==-8888)=nan;
DBZ=DBZ./100;
FLAG=hdfread([dataDir,infile],'CPR_Cloud_mask');
FLAG(FLAG==-9)=nan;
TOPO=hdfread([dataDir,infile],'DEM_elevation');
TOPO=TOPO{:};
TOPO(TOPO==-9999)=0;

%% Get right times
data.time=timeAll(firstInd:lastInd);
data.longitude=longitude(firstInd:lastInd);
data.latitude=latitude(firstInd:lastInd);
data.DBZ=double(DBZ(firstInd:lastInd,:))';
data.FLAG=FLAG(firstInd:lastInd,:)';
data.TOPO=double(TOPO(firstInd:lastInd));

%% Prepare data

% Flag non-cloud echo
data.DBZ(data.FLAG<30)=nan;
data.DBZ(106:end,:)=[];

% Create asl
data.asl=0:binSize:104*binSize;
data.asl=repmat(data.asl,length(data.time),1);
data.asl=flipud(data.asl');

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

%% Plot strat conv

disp('Plotting conv/strat ...');

close all

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

% 1D
stratConv1D=max(classSubPlot,[],1);
time1D=data.time(~isnan(stratConv1D));
stratConv1D=stratConv1D(~isnan(stratConv1D));

colmapSC=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    1,0,0;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    0.7,0,0];

col1D=colmapSC(stratConv1D,:);

close all

f1 = figure('Position',[200 500 1500 1100],'DefaultAxesFontSize',12,'visible',showPlot);

s1=subplot(4,1,1);

colormap jet

hold on
surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-35 25]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
colorbar
grid on
box on
title('Reflectivity (dBZ)')
s1pos=s1.Position;
s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];

s2=subplot(4,1,2);

hold on
surf(data.time,data.asl./1000,convDBZ,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 1]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
colorbar
grid on
box on
title('Convectivity')
s2pos=s2.Position;
s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];

s3=subplot(4,1,3);

hold on
surf(data.time,data.asl./1000,classBasic,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0.5 3.5]);
colormap(s3,[0,0,1;0.32,0.78,0.59;1 0 0]);
cb=colorbar;
cb.Ticks=[1,2,3];
cb.TickLabels=cat(2,{'Stratiform','Mixed','Convective'});
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
title('Basic stratiform/convective partitioning')
grid on
box on
s3pos=s3.Position;
s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];

s5=subplot(30,1,30);

hold on
scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
set(gca,'clim',[0,1]);
set(gca,'YTickLabel',[]);
s5.Colormap=colmapSC;
xlim([data.time(1),data.time(end)]);
grid on
box on
s5pos=s5.Position;
s5.Position=[s5pos(1),s5pos(2)-0.023,s1pos(3),s5pos(4)];

s4=subplot(4,1,4);

hold on
surf(data.time,data.asl./1000,classSubPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 10]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
s4.Colormap=colmapSC;
caxis([0.5 9.5]);
cb=colorbar;
cb.Ticks=1:9;
cb.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
set(gca,'XTickLabel',[]);
grid on
box on
title('Stratiform/convective partitioning')
s4pos=s4.Position;
s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];

linkaxes([s1 s2 s3 s4],'xy');

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'convStrat_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
