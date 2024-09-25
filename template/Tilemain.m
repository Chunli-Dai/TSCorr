% main program for getting coastline for each ArcticDEM tile
% Requirements: gdal software
% %%%% inputs needed

constant 

macdir='/Users/chunlidai/surge/';
macdir=[];

%addpath(genpath([codedir]));

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
   str=sprintf('find  %s -name ''*meta.txt'' > %s',deblank(imagedir),filename);
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
	  %/blue/chunlidai/data/imagery/WorldView//WV01_20140101082940_1020010029E44B00_14JAN01082940-P1BS-500096395080_01_P002_u16ns32637.tif
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
   [XYbi,rangei,projgdali]=imagebd(ifile);
   range(i,1:4)=rangei;XYbg{i}=XYbi;projgdalg{i}=projgdali;
end

save mat0.mat -v7.3

else 
load mat0.mat
end

dx=100e3;x0=-4000e3;y0=-4000e3;
dxs=dx/2/10; %5km
dxs=dx/50; %2km
%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;

% ArcticDEM mosaic tile grids
%Eureka lon lat: -85.9418555555556	79.988633
%32_33_2_2_5m_v2.0 ; 
%32_34_1_2_10_10 ; %work on tile, size of 5 km by 5 km.
%name convention; yid_xid_xids_yids_xidss_yidss
% inputtype=3;
%get the input from input.txt
currentdir=pwd;
filename=[currentdir,'/input.txt'] 
fid = fopen(filename) 
inputtype=fscanf(fid, '%d', [1, 1])';
if inputtype ==1
    str=fgetl(fid);tilefile=fgetl(fid);
elseif inputtype==2 || inputtype==3
   latlon=fscanf(fid, '%f', [2, 1])';
   lateq=latlon(1);loneq=latlon(2);
elseif inputtype ==4
   rang0=fscanf(fid, '%f', [4, 1])';
end

            
            switch inputtype
              case 1 %based on input xid etc.
%                   xid=33;yid=32;xids=2;yids=2;xidss=1;yidss=1;
%                   tilefile=sprintf('%02d_%02d_%01d_%01d_%02d_%02d.tif',yid,xid,xids,yids,xidss,yidss);
                  tilefile=deblank(tilefile);
              case 2 %find the block based on input coordinates; %not compatible with earthdem tiles.
%                   loneq=-85.9418555555556; lateq=79.988633; %Eureka
%                   lateq=61.143; loneq= -148.16;
                  [xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
                  %x=x0+(xid-1)*dx+(xids-1)*dx/2+(xidss-1)*dxs;y=y0+(yid-1)*dx+(yids-1)*dx/2+(yidss-1)*dxs;
                  xid=floor((xeq-x0)/dx)+1; yid=floor((yeq-y0)/dx)+1;
                  xids=floor((xeq-(x0+(xid-1)*dx))/(dx/2))+1;
                  yids=floor((yeq-(y0+(yid-1)*dx))/(dx/2))+1;
                  xidss=floor((xeq-(x0+(xid-1)*dx+(xids-1)*dx/2))/dxs)+1;
                  yidss=floor((yeq-(y0+(yid-1)*dx+(yids-1)*dx/2))/dxs)+1;
                  tilefile=sprintf('%02d_%02d_%01d_%01d_%02d_%02d.tif',yid,xid,xids,yids,xidss,yidss) ;
              case 3 %Use a rectangle box around the input coordinates;
                % get the data boundary, rang0, of this DEM tile 
%                   loneq=-151.1055; lateq=59.399; % 59.399, -151.1055; Kinnikinnick Landslide, which happened between Aug 23 and Sept 4 2017
		  [xeq,yeq]=latlon2xy(lateq,loneq,projgdal);
		% if strcmp(projgdal,'epsg:3413')
                %   [xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
	        % else
		%   [xeq,yeq]=mfwdtran(utmstruct,lateq,loneq);
		% end
                  exb=10e3; %10e3/2; %2.5e3;
%		  rang0=[xeq-exb xeq+exb yeq-exb yeq+exb ];
		  %resr=30; ranget=round(rang0/resr)*resr;rang0=ranget;
		  resr=15; %consistent with Landsat 15 m resolution images;
		  rang0=round([xeq-exb yeq-exb xeq+exb yeq+exb]/resr)*resr;
		  rang0=rang0+resr/2; %%consistent with Landsat grids.

		  %rang0=[490-3.5 490+3.5 6671-4  6674]*1e3; %7 by 7 km

                  tilefile=rang0  
                % [Co]=Landslide(rang0,range,XYbg,f,fdir);
	      case 4 %use boundary
                 % xeq=(rang0(1)+rang0(2))/2;yeq=(rang0(3)+rang0(4))/2;
		 % [lateq,loneq]=polarstereo_inv(xeq,yeq,[], [],70,-45);
                  tilefile=rang0;

              otherwise
                  disp('Unknown input zone! Try again!') 
            end    
            

            tic
%           [Co]=Slowmoving(tilefile,range,XYbg,f,fdir);
            [Co]=Slowmoving(tilefile,XYbg,projgdalg,f,fdir);
            fprintf(['\n Tile ',tilefile])
            toc

	    if flagoutput==0 %remove all files except output
 	      % [status, cmdout]=system('rm -fr [a-Z]*');
%          [status, cmdout]=system('rm -f [h-Z]*'); %save dX4Sg.mat
%          [status, cmdout]=system('rm -f [c-C]*');
%          [status, cmdout]=system('rm -f e*');
%          [status, cmdout]=system('rm -rf outcoregdem');
	    end % flagoutput
            
            exit

