function [Co]=Slowmoving(tilefile,XYbg,projgdalg,f,fdir)
%usage:             [Co]=slowmoving(tilefile,range,XYbg,f,fdir);
constant

Co=[];

%Grids of ArcticDEM Tiles
% 54_06_2_2_5m_v2 yid_xid_xids_yids
% 1,2 ; 2,2
% 1,1 ; 2,1
%x=x0+(xid-1)*dx; %left edge of a box;
dx=100e3;x0=-4000e3;y0=-4000e3;%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;
%dxs=dx/2/10; %5km
dxs=dx/50; %2km
% %%% Preparation: get rang0 from tile name; Need tilefile (instead of rang0) as the input for creating output filenames.
if isa(tilefile,'char')
    rang0=getbox(tilefile); 
    odir=[tilefile];
else %given rang0
    rang0=tilefile;
    odir=['givensite']; %notice odir should be all characters;
end

%create file 'imagelist.dat' which include lists of images only within the selected tile.
[XYbg,projgdalg,f,fdir]=creatsublist(rang0,XYbg,projgdalg,f,fdir); % need rang0 as the input

%Plot for time series
%Read files including all points (lon lat) for time series;
%For a 5 km block, check four points on regular grids (at e.g. +(1.25, 1.25), (3.75, 1.25), (1.25, 3.75), (3.75, 3.75) km pixels ).
dxq1=abs(rang0(2)-rang0(1))/2;dyq1=abs(rang0(4)-rang0(3))/2;
xq1grid=[rang0(1)+dxq1/2; rang0(1)+dxq1/2*3; rang0(1)+dxq1/2; rang0(1)+dxq1/2*3;];
yq1grid=[rang0(3)+dyq1/2; rang0(3)+dyq1/2; rang0(3)+dyq1/2*3; rang0(3)+dyq1/2*3;];
%[latq1grid,lonq1grid]=polarstereo_inv(xq1grid,yq1grid,[], [],70,-45); %[60.1768700940536 -141.184899979958]

%Creat a file xyq1.dat, which selecting all points from Higlist2023.gmt that are within the tile.
%getxyq(rang0); %to be created. Output: x, y coordinates in the given projection, point id index (line number).
[Higx,Higy,Higid] = getxyq(rang0);

filenameq1='xyq1.dat'; % n by 3; x, y id
if exist(filenameq1,'file')
fidq1 = fopen(filenameq1);
nq1 = linecount(fidq1);
fidq1 = fopen(filenameq1);
xyq1=fscanf(fidq1, '%f', [3, nq1])';

[nj,ni]=size(xyq1);
if nj~=nq1
             fprintf(['\n The format of ',filenameq1,' is wrong, should be n by 3.\n']);
	     fprintf(['\n Use only points from Higlist.\n']);
	     xyq1=[Higx(:), Higy(:), Higid(:)];
else
             fprintf(['\n Adding ',num2str(length(nj)),' more points.\n']);
	     Higxy=[Higx(:), Higy(:), Higid(:)];
	     xyq1=[Higxy;xyq1;];
end

xq1=[xyq1(:,1);xq1grid(:);];
yq1=[xyq1(:,2);yq1grid(:);];
else
xq1=[xq1grid(:);];
yq1=[yq1grid(:);];
end
[latq1,lonq1]=xy2latlon(xq1,yq1,projgdal);
nq1=length(xq1);

fprintf(['\n Total of ',num2str(nq1),' points to check time series. \n'])
fprintf(['\n Input points: '])
for i=1:nq1-length(xq1grid)
    fprintf(['\n ',num2str(i), ' (lon lat x y): ',num2str([lonq1(i), latq1(i), xq1(i), yq1(i)]), '; Higlist index:',num2str(xyq1(i,3))])
end
fprintf(['\n\n Regular grid points: '])
for i=nq1-length(xq1grid)+1:nq1
    fprintf(['\n ',num2str(i), ' (lon lat x y): ',num2str([lonq1(i), latq1(i), xq1(i), yq1(i)])])
end

    
fprintf('\n Step 1: Data Preparation \n')

fprintf('\n Step 1.1: get rock mask. \n')
if ~exist('icewater.mat','file')
Co=prep_rangerock(rang0);
end

fprintf('\n Step 1.2: prepare subset of images. \n')

Co=preparebatchcosicorr(rang0,XYbg,projgdalg,f,fdir);

%runstep=1; % 1, only run steps until cosi-corr; 2, run all steps.
if runstep==2

fprintf('\n Step 2: get displacements for all pairs of images \n')

Co=getvel();

fprintf('\n Step 3: retrieve time series \n')

[Co]=plotcosicorscarp_auto(rang0,xq1,yq1,odir);
end

return
end
