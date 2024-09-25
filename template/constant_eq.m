
%stripdir='/fs/project/howat.4/EarthDEM/region*/strips_unf/2m/'; %Directory of strip DEM files.
%imagedir='~/blue/data/ArcticDEMdata/earthquakes/turkeyeq/planet/'; %Directory of image files.
imagedir='/blue/chunlidai/data/imagery/WorldView/'; %Directory of image files.
codedir='/blue/chunlidai/apps/horizontalv/code/';  %Directory of codes.
cosidir='/Users/chunli/Desktop/Osuworkmac/asar/elliot/landsat/'; %%directory of files as input of COSI-Corr

flagtype=8; %image data types: 1, Landast 8; 2, Aster; 3, Planet; 4, PGC images; 5, Landsat 7 striped images; 6, subset; 8 Worldview images
imageores=1; %image resolution; 15 m landsat; 3 OR 4 m Planet; 1 OR 2 m Wolrdview;
flagtool=4; % tools to use: 1, COSI-COrr; 2, SETSM SDM; 3, Seongsu's MIMIC2; 4, COSI-COrr Plus

%
flagscarp=0; % 0 no scarp shapefile (run all grids);
	   % Otherwise, only run on scarps; 1 scarp in maxtrix in image coordinates; 2 scarp in lon lat as shapefile polygon;
scarpfile='scarp.shp';
flagplot=0; %Recommend 0. %1 only when using scarp shapefile
%
%'5N'	'7V'	'6V'	'7V'	'6V'	'6V'	'6V'	'7V'
% Allen GrandPlateau Portage Chisana Nassau Sanctuary Serpentine Nunatak

projgdal='epsg:32637'; % epsg:3413  epsg:32606 
projstrin='UTM zone 37 north'; %'polar stereo north'; % 'polar stereo south'; 'UTM zone 45 north';

mons=7;mone=8; %mon>=5&mon<=10; %mons, start month of snow-free seasons. %mone, end month of snow free months.
		%suggest to use all seasons; For latitude 80N, use July and August for summer;
		%For latitude 70N, use June to September for summer.
mons=1;mone=12; %all season.

algorithmin=2;%Recommend 2. Fit model for time series. 
	%e.g.  1 linear (ice melting); 2 constant (landslides) ; 3 constant + linear (Okmok volcano)

resr=60; %imageores*4 *5 ; %60; %30*5; % Recommend 60 m ( Old 150 m). Output spatial resolution in meters. 2; %4;%40.;
resr_cosi=imageores*4; %15*4; %Recommend 60 m. (step 4; output pixel from COSI-Corr); for plotcosicorscarp_auto.m
			% make sure it's consistent with the step size in cosiplus.json.

%plot scale
scale1=5e-3*5;
scalebar=20; %plot a scale bar of 20 m
scalexy=[384.79401       5649.894 ]*1e3; %location to plot the scale

year_start=0;year_end=9999; % min max year of measurement; year >= years and year <= yeare %Recommend no limit [0, 9999]
flagoutput=0; %Recommend 0. if 1 output lots of figures and data for plotting, 0 save space;

mindate=datenum('2020/02/06'); %datenum('2015/01/01');%0; %datenum('2013/01/01');%0;% %delete data befor this date
maxdate=9e9;%datenum('2025/10/17'); %9e9;delete data after this date
deldate=[];%deldate=[datenum('1999/11/03'),datenum('2001/07/28'),datenum('2015/05/24')]; %hi ; deldate=[datenum('2010/05/31'),datenum('2012/10/25'),datenum('2015/06/17'),datenum('2013/06/06'),datenum('2017/02/03')];

currentdir=pwd;
sdmlist=[currentdir,'/sdmlist.txt']; %input file list for SDM.
corrfilelist=[currentdir,'/corrfilelist.txt']; %a list of absolute file names of displacement map results from image correlation.

poolsize=1;

flag_vel_iw=1; %1, calculate velocity on all areas including land, ice, water;
	%0, calculate velocity only on land, NOT on ice or water;
	%2, calculate velocity on land and ice, Not on water. (Not functioning)
flag_pts=0; % 1, only calculate time series of given points; 0, calculate velocity for all grid points.
runstep=2; % 1, only run steps until cosi-corr; 2, run all steps.

flageq=1; %1, earthquake, no cropping, reserved image coverage; 0, landslide, can crop
flagcropin=0; %1 apply subsetting (default); 0 do not apply the subsetting.
coverthres=2; % the selected images should cover the tile > 60% for landslide; > 20% for earthquake.
flagwarp=1;
eqepoch(1)=datenum('20230206011734','yyyymmddHHMMSS');timefix=1;
eqepoch(2)=datenum('20230206102448','yyyymmddHHMMSS');

flagvolc=1; %duration of event time;
cosipluslist=[currentdir,'/cosipluslist.txt'];
subtiledx=25; %10 km, or 25 km, divisible by 50 km. Size of subtiles for dividing 50 km tile.

currentdir=pwd;
imagesubdir=[currentdir,'/imagesubdir/'];  %directory of subset images

corrresultdir=[currentdir,'/results/'];  %directory for correlation results.
corrcosiresultdir=[currentdir,'/cosiresults/'];  %directory for correlation results.

if ~exist(imagesubdir,'dir')
  mkdir(imagesubdir)
end
if ~exist(corrresultdir,'dir')
  mkdir(corrresultdir)
end
if ~exist(corrcosiresultdir,'dir')
  mkdir(corrcosiresultdir)
end

picdir='./pics/';
if ~exist(picdir,'dir')
  mkdir(picdir)
end

if ~exist('para','dir')
  mkdir para
end
