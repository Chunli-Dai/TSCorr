% main program for getting coastline for each ArcticDEM tile
% Requirements: gdal software
% %%%% inputs needed

constant 

macdir='/Users/chunlidai/surge/';
macdir=[];

currentdir=pwd;

addpath(genpath([codedir]));

shpname='./GSHHS/GSHHS_f_L1.shp';% a priori coastline shapefile

% %%%% control parameters
width=2e3; %buffer width of the a priori coastline, e.g., 2km.

if ~exist('mat0.mat','file') %readding boundary; time 1 hour for 90751 strip files and 10130 mono xml files

% %%% Preparation: get the list of strip files and boundries
filename='imagelistall.dat';
if ~exist(filename,'file')
   %flagtype=1; %image data types: 1, Landast 8; 2, Aster; 3, Planet; 4, PGC images; 5, Landsat 7 striped images;
   if flagtype==4
   %str=sprintf('find  %s -name ''*mdf.txt'' > %s',deblank(stripdir),filename);
   %str=sprintf('find  %s -name ''*meta.txt'' > %s',deblank(imagedir),filename);
   elseif flagtype==1 || flagtype==6 %
   str=sprintf('find  %s -name ''L*_B8.TIF'' > %s',deblank(imagedir),filename);
   elseif flagtype==2
   %AST_L1T_00302092005192338_20150508061919_118525_V.tif
   str=sprintf('find  %s -name ''AST_L1T_*_V.tif'' > %s',deblank(imagedir),filename);
   elseif flagtype==3 %Planet images
% 20180928_merge.tif
   %str=sprintf('find  %s -name ''*_merge.tif'' > %s',deblank(imagedir),filename);
   str=sprintf('find  %s -name ''*AnalyticMS*.tif'' > %s',deblank(imagedir),filename);
   elseif flagtype==8 % Worldview images
	  str=sprintf('find  %s -name ''W*.tif'' > %s',deblank(imagedir),filename);
   end
	str
  [status, cmdout]=system(str);
   % filename need to be sorted by date
   % sort -n -t'_' -k19 corrfilelist -o corrfilelist_sorted
end

fprintf ('\n Step 0: geting the boundary for all files in the region.\n')
%READ INPUT PARAMETERS; getting the boundaries for all files
% filename='boundaries_reg31.dat';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
range=zeros(n,4);XYbg=cell(n,1);projgdalg=cell(n,1);
for i=1:n
   ifile=[macdir,fgetl(fid)];
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
   %satname=f{i}(1:4);

   % get the boundary from xml file
   try
   [XYbi,rangei,projgdali]=imagebd(ifile);
   catch e
        fprintf(['\n Tilemain_nov.m There was an error! The message was:',e.message,'; image file:',ifile]);
	XYbi=[0 0;];    rangei=[0 0 0 0]; projgdali='epsg:32637';
   end

   range(i,1:4)=rangei;XYbg{i}=XYbi;projgdalg{i}=projgdali;
end

save mat0.mat -v7.3

else 
load mat0.mat
end

%% nov, Collect number of overlapping images.

%lat lon
x0=-180; y0=30;
xe=-100; ye=75;
resrc=0.02; % 0.01 degree =1km;
ofile='nov.tif'; %hi
if exist(ofile,'file')
nov=readGeotiff(ofile);
nx=length(nov.x);ny=length(nov.y);
else %40m resolution
% integer type, 0 for nan.
%resrc=400;
nov.x=x0:resrc:xe;
nov.y=ye:(-resrc):y0;
nx=length(nov.x);ny=length(nov.y);
nov.z=uint16(zeros(ny,nx));
end

resrc=mean(diff(nov.x));
%M=logical(size(nov.z));

novt=nov; % for this region
novt.z=uint16(zeros(ny,nx));

for i=1:n
        XYbi=XYbg{i};
        Xb=XYbi(:,1);Yb=XYbi(:,2);

%       if strcmp(projgdal(1:7),'epsg:32')
                %xy to lat lon
                projgdalj=projgdalg{i};
                xj=Xb;yj=Yb;
                [latj,lonj]=xy2latlon(xj,yj,projgdalj);
                Xb=lonj;Yb=latj;
%       end

        idx=round((Xb-nov.x(1))/resrc)+1;
        idy=round((Yb-nov.y(1))/(-resrc))+1;
        Mb = poly2mask(idx,idy, ny,nx); % build polygon mask
        novt.z=novt.z + uint16(Mb);
end

%update nov with the new novt when nov is zero and novt is larger.
M=novt.z>nov.z;
nov.z(M)=novt.z(M);

projstr='polar stereo north'; %actually lat lon
%projstr=projstrin;
writeGeotiff(ofile,nov.x,nov.y,uint16(nov.z),12,255,projstr)

%save testnov.mat -v7.3

exit %stop here to get the mat0.mat

