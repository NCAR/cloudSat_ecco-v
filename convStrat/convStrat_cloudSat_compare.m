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
infile='2015184230321_48842_CS_2B-GEOPROF_GRANULE_P1_R05_E06_F00.hdf';

dataDir2='/scr/snow2/rsfdata/projects/cset/cloudSat/PRECIP-COLUMN/hdf/';
infile2='2015184230321_48842_CS_2C-PRECIP-COLUMN_GRANULE_P1_R05_E06_F00.hdf';

figdir=['/scr/sci/romatsch/other/convStratPaperHCR/'];

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

%% Read conv_strat flag

csFlag=hdfread([dataDir2,infile2],'Conv_strat_flag');
csFlag=double(csFlag{:});
csFlag=csFlag(firstInd:lastInd);

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

%% Prepare for plotting

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

% 1D ecco
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

% 1D cloudsat
csFlag(csFlag<1)=nan;

time1DcloudSat=data.time(~isnan(csFlag));
cloudSatFlag1D=csFlag(~isnan(csFlag));

colmapCloudSat=[1,0,0;
    0,0,1;
    1,1,0];

col1DcloudSat=colmapCloudSat(cloudSatFlag1D,:);

% Text
textAlt=11.3;
textDate=datetime(2015,7,4,0,5,32);

%% Plot
disp('Plotting conv/strat ...');

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

s1=subplot(4,1,1);

colormap jet

hold on
surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-35 25]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
set(gca,'XTickLabel',[]);
cb1=colorbar;
grid on
box on
text(textDate,textAlt,'(a) Reflectivity (dBZ)','FontSize',11,'FontWeight','bold');

s2=subplot(4,1,2);

hold on
surf(data.time,data.asl./1000,convDBZ,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 1]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
cb2=colorbar;
set(gca,'XTickLabel',[]);
grid on
box on
text(textDate,textAlt,'(b) Convectivity','FontSize',11,'FontWeight','bold');

s5=subplot(30,1,28);

hold on
scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
set(gca,'clim',[0,1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
s5.Colormap=colmapSC;
xlim([data.time(1),data.time(end)]);
grid on
box on

s6=subplot(30,1,30);

hold on
scat1=scatter(time1DcloudSat,ones(size(time1DcloudSat)),10,col1DcloudSat,'filled');
set(gca,'clim',[0,4]);
set(gca,'YTickLabel',[]);
s6.Colormap=colmapCloudSat;
xlim([data.time(1),data.time(end)]);
grid on
box on
text(datetime(2015,7,4,0,8,31),3,'Convective','edgecolor','r','Margin',0.1)
text(datetime(2015,7,4,0,8,31),0,'Stratiform','edgecolor','b','Margin',0.1)
text(datetime(2015,7,4,0,8,31),-3,'Shallow','edgecolor','y','Margin',0.1)
text(textDate,3.5,'(d) Echo type CloudSat','FontSize',11,'FontWeight','bold');

s4=subplot(4,1,3);

hold on
surf(data.time,data.asl./1000,classSubPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 10]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
s4.Colormap=colmapSC;
caxis([0.5 9.5]);
cb4=colorbar;
cb4.Ticks=1:9;
cb4.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
set(gca,'XTickLabel',[]);
grid on
box on
text(textDate,textAlt,'(c) Echo type','FontSize',11,'FontWeight','bold');

s1.Position=[0.049 0.715 0.82 0.27];
s2.Position=[0.049 0.432 0.82 0.27];
s4.Position=[0.049 0.15 0.82 0.27];
s5.Position=[0.049 0.12 0.82 0.02];
s6.Position=[0.049 0.057 0.82 0.02];

cb1.Position=[0.875,0.715,0.02,0.27];
cb2.Position=[0.875,0.432,0.02,0.27];
cb4.Position=[0.875,0.15,0.02,0.27];

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'cloudSatCompare.png'],'-dpng','-r0')
