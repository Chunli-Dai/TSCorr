function [Co]=plotcosicorscarp_auto(rang0,xq1,yq1,odir)
%Plot the output of all pairs of displacement from cosicor 
% and get the average displacement over the scarp area.
% Modified from /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/plotcosicorscarp_main.m
% Changes: automated search on regular grids.

%addpath(genpath(['/Users/chunlidai/surge/home/dai.56/arcticdemapp/landslide/code1/']));
%addpath('/Users/chunlidai/ArchivedBoxSync/OpticalImages/MIMC2/MIMC2_demo_package_release_170412/envi/');
constant

Co=[];

nq1=length(xq1);

if flageq==1; %earthquake
	namestr='disp'; %displacement in meter
elseif flageq==0; %0, landslide
	namestr='rate'; %rate in meter/year
end


OutName=[odir,'_',namestr,'.tif'];
if exist(OutName,'file')
	fprintf(['\n Final result file exist, skip this tile:',OutName,' \n'])
	return
end


%flagtaan=1; %0 not taan landslie, 1 Taan Landslide

% change: now use absolute directory in corrfilelist
%filename='corrfilelist';%landsat-7 and 8; /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/
filename=corrfilelist;
%corrresultdir='/Users/chunlidai/share/sdm/runsite22/landsatcosicorr/results/';
landsatdir=imagesubdir; %imagedir;%directory of input images for plotting background
icewaterfile=[currentdir,'/icewater.mat'];
load(icewaterfile)
clear tag water2 %free space

%Taan see /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/plotcosicorscarp_main_taan.m

fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
for i=1:n
   ifile=[fgetl(fid)];
   [diri,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[diri,'/'];%working on two regions 
end
fclose(fid);

%Only work on subset of all images around the scarp;
%resr=30*5; %output grid
%resr_cosi=15*4; %(step 4; output pixel from COSI-Corr);
%rang0=[490-3.5 490+3.5 6671-4  6674]*1e3; %7 by 7 km

xout=rang0(1):resr:rang0(2);yout=rang0(4):-resr:rang0(3);
rate=zeros(length(yout),length(xout));ratestd=rate;
ratex=rate;ratey=rate;ratexstd=rate;rateystd=rate;
%winsize=resr_cosi*7;%window size in meters for getting the motion rate:7by7 pixels (output pixel from COSI-Corr);
% rang0=rang0+[-4e3 +4e3 -4e3 +4e3];
rang0=rang0+[-3e3 +3e3 -3e3 +3e3];
scarp.x=rang0(1):resr:rang0(2);scarp.y=rang0(4):-resr:rang0(3);scarp.z=zeros(length(scarp.y),length(scarp.x));

%Get the control points 
% load /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/icewater_sv1barry.mat
% load /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/plotdata/icewaterbarry_sv2_1994a2009.mat %High resoluton;from RGI 1994 to 2009
clear control
control=icewater; control.z=~icewater.z;

% For output pixel selection, flag_vel_iw
controlo.x=xout;controlo.y=yout;
controlo.z=interp2(control.x,control.y,double(control.z), xout,yout','*nearest',0);

% load /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/outlinexy.mat

figure;imagesc(icewater.x*1e-3,icewater.y*1e-3,icewater.z);colorbar
xlabel('x (km)');ylabel('y (km)');
axis equal;
% axis([431 445 6772 6784]);
view(0,-90);
% hold on;plot(x*1e-3,y*1e-3,'b.-') %scarp outline in outlinexy.mat
colormap gray
title('Ice and water surfaces')
saveas(gcf,'icewater','fig')
clear icewater %free space
close all
 
% demdir='/Users/chunlidai/share/sdm/runsite22/landsatcosicorr/results/';
% landsatdir='/Applications/bda/LandsatBarry/';

%load all data for faster computation
fprintf('\n Loading all data:')
tic

Dall=cell(n,1);

for i=1:n %get the displacement of all pairs of displacements;
    
%ifile=[corrresultdir,f{i}];
ifile=[fdir{i},f{i}];
if ~exist(ifile,'file')
    continue
end

%flagtool=2; % tools to use: 1, COSI-COrr; 2, SETSM SDM; 3, Seongsu's MIMIC2.
if flagtool==1
hdrfile=[strtrim(ifile),'.hdr'];
if ~exist(hdrfile,'file')
    fprintf(['\n File does not exist:', hdrfile,'; i=',num2str(i)]);
end
[D,info]=enviread(ifile); 

%double check if the map projection is consistent with the one in constant.m
str1=info.coordinate_system_string;
[C,matches] = strsplit(str1,{'UTM_Zone_','"'});
if ~strcmpi(z1,C{2})
   fprintf(['\n Error: the map projection of file is not consistent with the given projection :', ifile,'; i=',num2str(i)]);
end

elseif flagtool==2   %2, SETSM SDM;
%construct a D and info. D matrix: n by m by 3; info.x info.y
% output: ./results/correlation_YYYYMMDDHHMMSS_LC08_B8_sub.TIF_vs_YYYYMMDDHHMMSS_LC08_B8_sub.TIF
filestr=f{i};

[C,matches] = strsplit(filestr,{'correlation_','_vs_'});
% date1=C{2}(1:14);date2=C{3}(1:14);

filename=C{2};filename2=C{3};
date1=filename2date(filename);
date2=filename2date(filename2);

%dt=datenum(date2,'yyyymmddHHMMSS')-datenum(date1,'yyyymmddHHMMSS'); %date difference in days
str1=date1;
    if length(strtrim(str1)) ==8
       datefmt='yyyymmdd';
    elseif length(strtrim(str1)) ==14
       datefmt='yyyymmddHHMMSS';
    else
       fprintf(['\n plotcosicorscarp_auto.m error: date format is not yyyymmddHHMMSS or yyyymmdd: ',str1,' \n'])
    end
dt=datenum(date2,datefmt)-datenum(date1,datefmt); %date difference in days

ifiley=strrep(ifile,'dx.tif','dy.tif');
datax=readGeotiff(ifile); % runL8_dx.tif
datay=readGeotiff(ifiley); % runL8_dy.tif

%reduce resolution from 15 m to 60 m,  imageores to resr_cosi
dsr=1/(resr_cosi/imageores);  %15/60
datax.z = imresize(datax.z,dsr);
datax.x = imresize(datax.x,dsr);
datax.y = imresize(datax.y,dsr);
datay.z = imresize(datay.z,dsr);
datay.x = imresize(datay.x,dsr);
datay.y = imresize(datay.y,dsr);
%

vx=datax.z*dt;vy=datay.z*dt; % vx=vx.z*dt;%m/day -> m
M=vx==0|vy==0;
vx(M)=nan;vy(M)=nan; %set 0 to void
[ny,nx]=size(vx);
D=zeros(ny,nx,2);

info.x=datax.x;info.y=datax.y;
D(:,:,1)=vx;D(:,:,2)=vy; %meter

elseif flagtool==3 %see script_demo.m

ifile=[param.path.vmap,'/',filename_vmap];
load(ifile); 
ppfile=strrep(ifile,'vmap_','pp_result_');ppfile=strrep(ppfile,'.mat','_raw.mat');
load(ppfile);
info.x=vmap.x;info.y=vmap.y;
vx=vmap.vx*vmap.datediff/365.25;%m/yr -> m
vy=vmap.vy*vmap.datediff/365.25;
[ny,nx]=size(vx);
D=zeros(ny,nx,2);
D(:,:,1)=vx;D(:,:,2)=vy; %meter

elseif flagtool==4 % 4, COSI-COrr Plus

try
    data=readGeotiff(ifile); 
catch e
    fprintf(['\n Error reading ifile:',ifile,' The message was:\n%s',e.message]);
    continue
end

    info.x=data.x;info.y=data.y;
    [ny,nx,nz]=size(data.z);
    D=nan(ny,nx,2);
    D(:,:,1)=data.z(:,:,1);
    D(:,:,2)=data.z(:,:,2); %meter

end
Dall{i}.data=D;Dall{i}.info=info;
end % for i=1:n 
toc

save t1.mat -v7.3

% Coregistration; 17 hours for 26129 files.
coregfile=[currentdir,'/coreg.mat'];
if exist(coregfile,'file')
	fprintf(['\n Load Coregistration: ',coregfile,'.\n'])
	load(coregfile) 
else

fprintf('\n Start Coregistration.\n')
% For earthquakes, there is the risk of ZERO displacements on only one side
% of the fault. Maybe delete these pairs.

%initialize empty file
% file='scarpdxdy.txt';
% fid10 = fopen('scarpdxdy.txt','w');
scarpdxdy(n)=struct('i',[],'datem',[],'dates',[],'txcstd',[],'tycstd',[],'vxc',[],'vyc',[],'Mf',[],'dx',[],'dy',[],'ratio',[],'x',[],'y',[]); %ratio  

for i=1:n %get the displacement of all pairs of displacements;
    if isempty(Dall{i})
        continue
    end
    
%ifile=[demdir,f{i}];
ifile=[fdir{i},f{i}];
%Fig. 2 displacement map
% ifile='/Users/chunlidai/share/sdm/runsite22/landsatcosicorr/2011vs2015ws32step4.tif';%ASTER; i=110
% ifile='/Users/chunlidai/share/sdm/runsite22/landsatcosicorr/20150617vs20170615ws32step4.tif';%Landsat;i=361
% infile='/Applications/bda/LandsatBarry/LC08_L1TP_067017_20170615_20170629_01_T1_B8.TIF';
% ifile='/Users/chunlidai/share/sdm/runsite22/monocosicorr/run1/correlation_IK01_20020402213700_2002040221371390000010000699_po_945612_pan_0000000_u16ns3413_sub_utm2m.tif_vs_WV01_20100427213038_102001000DA0E100_10APR27213038-P1BS-052121841010_04_P007_u16ns3413_sub_utm2m.tif_18'; %mono images 
% ifile='/Users/chunlidai/share/sdm/runsite22/landsatcosicorr/result4aster/correlation_AST_L1T_00309272011211825_20150607194411_118314_V.tif_vs_AST_L1T_00307092017211916_20170710102105_13989_V.tif_113';
% ifile='/Users/chunlidai/share/sdm/runsite22/landsatcosicorr/20100526vs20110927ws32step4.tif';
% ifile='/Users/chunlidai/share/sdm/runsite22/landsatcosicorr/20110927vs20170709ws32step4.tif';
% ifile='/Users/chunlidai/share/sdm/runsite22/landsatcosicorr/20100526vs20140926ws32step4.tif';
% ifile='/Users/chunlidai/share/sdm/runsite22/landsatcosicorr/20140926vs20170709ws32step4.tif';
%planet
% ifile='/Users/chunlidai/share/sdm/runsite22/planet/results/correlation_20160929_202103_0e0e_3B_Visual.tif_vs_20170722_202631_1044_3B_Visual.tif_1';

if 0 %flagplot==1  %&& flagplotsv==1
    
%plot the master image
% infile='/Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/site22landsat/LC08_L1TP_068017_20150617_20170226_01_T1_B8.TIF';
% infile='/Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/plotdata/BarryArmPlanet/copy1/20160929_202103_0e0e_3B_Visual_geo.tif';
[~,name,~] =fileparts([strtrim(ifile)]);
% text1=[name(11:18),' - ',name(1:8)];

% standard name convention; PGC monoimages and planet images
% imagesubdir, 20131028000000_LC08_B8_sub.TIF, from correlation_20131028000000_LC08_B8_sub.TIF_vs_20201024000000_LC08_B8_sub.TIF_dx.tif
    %id=findstr(name,'correlation_');id2=findstr(name,'_vs_');
    %infile=[landsatdir,name(id(1)+12:(id2(1)-1))];
    [C,matches] = strsplit(name,{'correlation_','_vs_'});
     infile=[landsatdir,'/',C{2}];

data=readGeotiff(infile,'map_subset',rang0);
if contains(name,'AST')
%data.z=double(data.z(:,:,2));%ASTER band 2
end
figure;
set(gcf,'Color','white')
set(gca,'FontSize', 14);
set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
imagesc(data.x*1e-3,data.y*1e-3,data.z);colorbar
xlabel('x (km)');ylabel('y (km)');
% caxis([6000 13000])
view(0,-90)
% axis equal;axis([310 320 1250 1260]);
axis equal;
% axis([433 443 6776 6785]);
% axis([431 445 6772 6784]);%default zone
% axis([435 442 6775.5 6782.5])%%zoom in 
view(0,90);
%hold on;plot(x*1e-3,y*1e-3,'g.-') %scarp outline in outlinexy.mat
hold on;plot(scarpx*1e-3,scarpy*1e-3,'m.-') %scarp outline in outlinexy.mat
colormap gray
title([name(13:59),'; i=',num2str(i)])
box on

% caxis([0 50])
freezeColors

end

%scale1=5e-3*5;
%scalebar=20; %plot a scale bar of 20 m
%scalexy=[434.29501      6782.9702]*1e3; %location to plot the scale

% figure
% [Co]=plotcosicorr(ifile,scale1,scalebar,scalexy,control,scarp,i,Dall{i});
[Co]=plotcosicorr(ifile,scale1,scalebar,scalexy,control,[],i,Dall{i});
% pause
scarpdxdy(i)=Co; %Must have the exact same elements.
end % for i=1:n 
% fclose(fid10);
fprintf('\n Done with Coregistration.\n')
toc
%save coreg.mat scarpdxdy -v7.3
save(coregfile,'scarpdxdy', '-v7.3');
end % if exist

clear Dall ; %free space

%  numlabs; %get the number of workers
nx=length(xout);
ratecell=cell(nx,1);
flagscarp=flagscarp;
sz = getenv('SLURM_NTASKS');
sz=str2num(sz);
fprintf(['\n ',num2str(sz),' worker(s) allocated in job.slurm.\n'])
[status, cmdout]=system('free -h')
whos
poolobj=parpool(sz);
parfor io=1:nx 

    [status, cmdout]=system('free -h')

    [rate_io, ratestd_io, ratex_io, ratexstd_io, ratey_io, rateystd_io]= getratesub(io,xout,yout,flagscarp,controlo,xq1,yq1,scarp,info,control,scarpdxdy);
    ratecell_io=cell(6,1);
    ratecell_io{1}=rate_io;
    ratecell_io{2}=ratestd_io;
    ratecell_io{3}=ratex_io;
    ratecell_io{4}=ratexstd_io;
    ratecell_io{5}=ratey_io;
    ratecell_io{6}=rateystd_io;
    ratecell{io}=ratecell_io;
end %io
delete(poolobj)

%save t2.mat -v7.3

for io=1:nx
    %rate(:,io)=rate_io(:,io);
    %ratestd(:,io)=ratestd_io(:,io);
    %ratex(:,io)=ratex_io(:,io);
    %ratexstd(:,io)=ratexstd_io(:,io);
    %ratey(:,io)=ratey_io(:,io);
    %rateystd(:,io)=rateystd_io(:,io);
    rate(:,io)=ratecell{io}{1};
    ratestd(:,io)=ratecell{io}{2};
    ratex(:,io)=ratecell{io}{3};
    ratexstd(:,io)=ratecell{io}{4};
    ratey(:,io)=ratecell{io}{5};
    rateystd(:,io)=ratecell{io}{6};
end

figure;
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
imagesc(xout*1e-3,yout*1e-3,rate);title('rate');view(0,-90);colorbar;colormap jet;
dn=1;scale1=5e-3*5*10;
hold on;hold all;quiver(xout(1:dn:end)*1e-3,yout(1:dn:end)*1e-3,(ratex(1:dn:end,1:dn:end))*scale1,(ratey(1:dn:end,1:dn:end))*scale1,0,'k','linewidth',2);
view(0,90)

saveas(gcf,'rate','fig')

figure;
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]);
imagesc(xout*1e-3,yout*1e-3,ratestd);title('rate std');view(0,-90);colorbar;colormap jet
saveas(gcf,'ratestd','fig')

save rate.mat xout yout rate ratestd ratex ratexstd ratey rateystd

%projstr=projgdal;
projstr=projstrin;
OutName=[odir,'_',namestr,'.tif'];
writeGeotiff(OutName,xout,yout,double(rate),4,nan,projstr); 
%OutName=[odir,'_ratestd.tif'];
OutName=[odir,'_',namestr,'std.tif'];
writeGeotiff(OutName,xout,yout,double(ratestd),4,nan,projstr); 
OutName=[odir,'_',namestr,'x.tif'];
writeGeotiff(OutName,xout,yout,double(ratex),4,nan,projstr); 
OutName=[odir,'_',namestr,'xstd.tif'];
writeGeotiff(OutName,xout,yout,double(ratexstd),4,nan,projstr); 
OutName=[odir,'_',namestr,'y.tif'];
writeGeotiff(OutName,xout,yout,double(ratey),4,nan,projstr); 
OutName=[odir,'_',namestr,'ystd.tif'];
writeGeotiff(OutName,xout,yout,double(rateystd),4,nan,projstr); 

if 0
[X,Y]=meshgrid(xout,yout);
[LAT,LON]=polarstereo_inv(X,Y,[],[],70,-45);
out=[LON(:),LAT(:),rate(:),ratestd(:)];
save -ascii rate.txt out
end

%z1 ='7N';
[X,Y]=meshgrid(xout,yout);

[LAT,LON]=xy2latlon(X(:),Y(:),projgdal);
%[LAT,LON]=minvtran(utmstruct,X,Y);
out=[LON(:),LAT(:),rate(:),ratestd(:),ratex(:),ratexstd(:),ratey(:),rateystd(:)];
save -ascii rate.txt out

    close all

return
end
