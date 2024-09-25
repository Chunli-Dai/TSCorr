function [Co]=plotcosicorr(ifile,scale1,scalebar,scalexy,control,scarp,iorder,Dall)
%load cosicorr results and plot the vector map
%Modified from /Users/chunlidai/ArchivedBoxSync/OpticalImages/MIMC2/MIMC2_demo_package_release_170412/script_demo.m
% addpath('/Users/chunlidai/ArchivedBoxSync/OpticalImages/MIMC2/MIMC2_demo_package_release_170412/envi/');
constant
%flagplot=0; %1;
flagscreen=0; %1 output results at screen; 0 no output
Co=[];

[~,name,ext] =fileparts([strtrim(ifile)]);
name=[name,ext]; %sep 2023.

%standarized filenames e.g.: */correlation_20131028000000_LC08_B8_sub.TIF_vs_20201024000000_LC08_B8_sub.TIF_dx.tif
[C,matches] = strsplit(name,{'correlation_','_vs_'});
%datem=C{2}(1:14);dates=C{3}(1:14);
filename=C{2};filename2=C{3};
datem=filename2date(filename);
dates=filename2date(filename2);

text1=[dates,' - ',datem];
texto=[dates,'m',datem];
    
if exist('Dall','var') %load data
    D=Dall.data;info=Dall.info;
else
    [D,info]=enviread(ifile);
end
txc=nanmean(nanmean(D(:,:,1))); %east/west
tyc=nanmean(nanmean(D(:,:,2))); %north/south
Dx=D(:,:,1);Dy=D(:,:,2);Ds=sqrt(Dx.^2+Dy.^2);
Ms=Ds<10;
txc=nanmean(nanmean(Dx(Ms))); %east/west
tyc=nanmean(nanmean(Dy(Ms))); %north/south
% vxc=D(:,:,1)-txc;
% vyc=D(:,:,2)-tyc;
vxc=D(:,:,1);
vyc=D(:,:,2);
vs=sqrt(vxc.^2+vyc.^2); %magnitude of velocity (meter)
vxcorg=vxc;vycorg=vyc;

%filter is the angle of displacement over a neighber 9 pixels is larger
%than a threshold , remove pixels;
%refers to ~/arcticdemapp/landslide/code1/Landslidefilter.m 
% using convolution to compute std, faster,
angle=atan2(vyc,vxc)*180/pi;
kernel_size = 3; % dx=11*resr;
k=ones(kernel_size)/(kernel_size.^2);
uxmn = nanconv(angle, k, 'noedge');
uxstd=sqrt( abs(nanconv(angle.^2,k,'noedge') - uxmn.^2 )); % obtain standard deviation
Mf=uxstd>=60; %small magnitude of displacement should still have same direction;
% Mf=uxstd>=60&vs>5;% big, but with random directions; displacement angle  homogeneous;
% figure;imagesc(info.x*1e-3,info.y*1e-3,uxstd);colorbar

if flagscreen==1
fprintf(['[tx, ty (m)]=',num2str([txc,tyc]),'\n'])
end

if exist('control','var')
    %if control points map exist; apply the mask to get tx ty, and std;
    % control points good points percentage;
    Mrock=interp2(control.x,control.y,double(control.z), info.x,info.y','*nearest',0);
    
    %1 within the scarp
    if ~isempty(scarp) % exist('scarp','var')
    Mscarp=interp2(scarp.x,scarp.y,double(scarp.z), info.x,info.y','*nearest',0);
%     Mscarp=logical(Mscarp); %old
    Mscarp=logical(Mscarp)&logical(Mrock); %new August 2020; only calculate the motion of scarp area within rock area, excluding ice/water.
    else
       Mscarp=false(size(Mrock));     
    end
    Mrock=logical(Mrock);
    
%     Mcontrol=Mrock&~Mscarp;
    Mcontrol=Mrock&~Mscarp&~isnan(vs);
%     save Mcontrol.mat info Mcontrol
    
%     Mcsel=~Mf&Ms&Mcontrol;%used this for landsat
    Mcsel=~Mf&Mcontrol;
    txc=nanmedian(nanmedian(Dx(Mcsel))); %east/west
    tyc=nanmedian(nanmedian(Dy(Mcsel))); %north/south
    txcstd=nanstd((Dx(Mcsel)));
    tycstd=nanstd((Dy(Mcsel)));
    
    %selected points for control area (second iteration): exclude outliers
    Mcsel=Mcsel&abs(Dx-txc)<3*txcstd&abs(Dy-tyc)<3*tycstd;

    Mcsel2=~Mf&abs(Dx-txc)<3*txcstd&abs(Dy-tyc)<3*tycstd; %without Mcontrol, good points in all areas.

    txc=nanmedian(nanmedian(Dx(Mcsel))); %east/west
    tyc=nanmedian(nanmedian(Dy(Mcsel))); %north/south
    txcstd=nanstd((Dx(Mcsel)));
    tycstd=nanstd((Dy(Mcsel)));

    %apply threshold for good and bad cases
    nptc=sum(Mcontrol(:)); %all points within Mcontrol
    ncgood=sum(Mcsel(:)); %good points selected for txc calculation
    ratioc=ncgood/nptc*100;
    %good case: 54% 55% 33% 44% 25% 12% 45%
    %bad case: 38% 28%
    
    %ratioc has to > 3%
    if flagscreen==1
    fprintf(['\n [tx, ty, txstd, tystd (m), ratioc]=',num2str([txc,tyc,txcstd,tycstd, round(ratioc)]),'%% \n'])
    end
    
    vxc=vxc-txc;vyc=vyc-tyc;% apply shift;
    vxcorg=vxcorg-txc;vycorg=vycorg-tyc;% apply shift;
%     vxc(~Mcontrol&~Mscarp)=nan;vyc(~Mcontrol&~Mscarp)=nan;%map out ice and water but keep scarp
    vxc_out=vxc; vyc_out=vyc; %output displacement maps without masking out control areas.
    
    %filter out random big displacements; but keep scarp
    % vxc(Mf&~Mscarp)=nan;vyc(Mf&~Mscarp)=nan;
    vxc_out(~Mcsel2&~Mscarp)=nan;vyc_out(~Mcsel2&~Mscarp)=nan;%mask out Mf (random direction) and outliers, but including ice/water.
    vxc(~Mcsel&~Mscarp)=nan;vyc(~Mcsel&~Mscarp)=nan;%map out non selected points;
    vxc2=vxc;vyc2=vyc;

    %within the scarp, map out displacement that has direction 20 away from
    %mean direction;
    if 1
        angle=atan2(vyc,vxc)*180/pi;%update the angle after txc tyc correction!!
    meanangle=nanmedian(angle(Mscarp&~Mf));
    stdangle=nanstd(angle(Mscarp));
    Mscarpbad=~(abs(angle-meanangle)<=20)&Mscarp;  %;label nans as bad ones too; abs(angle-meanangle)>20&Mscarp; 
    vxc(Mscarpbad)=nan;vyc(Mscarpbad)=nan;%
    end

end

if ~isfield(info, 'x') 
    fprintf(['\n Warning: no geographic coordinates, construct the map coordinates.\n'])
    [m,n]=size(vxc);
    npixel=16;
    info.x=[1:n]*1e3*npixel;info.y=[1:m]*1e3*npixel;
end

if ~exist('scalexy','var')
    scalexy(1)=(info.x(1)+(info.x(end)-info.x(1))*0.05);
    scalexy(2)=(info.y(1)+info.y(end))/2;
end

if flagplot==1
dn=2;%figure

if exist('vxc2','var')
    %plot all vectors including ones within the scarp;
    hold all;
    %all vectors
    quiver(info.x(1:dn:end)*1e-3,info.y(1:dn:end)*1e-3,(vxcorg(1:dn:end,1:dn:end))*scale1,(vycorg(1:dn:end,1:dn:end))*scale1,0,'g','linewidth',2);
    %selected vectors with good direction, and ALL vectors within scarp
    quiver(info.x(1:dn:end)*1e-3,info.y(1:dn:end)*1e-3,(vxc2(1:dn:end,1:dn:end))*scale1,(vyc2(1:dn:end,1:dn:end))*scale1,0,'b','linewidth',2);
    axis equal
end

    if 0  %clean figure;
         M1=(sqrt(vxc.^2+vyc.^2))<20&(Mrock&~Mscarp);vxc(~M1)=nan;vyc(~M1)=nan;%control area
%          M1=~(abs(angle-meanangle)<=6)&Mscarp;vxc(M1)=nan;vyc(M1)=nan;%scarp area;
         M1=~(abs(angle-meanangle)<=20)&Mscarp;vxc(M1)=nan;vyc(M1)=nan;%scarp area; %same as Mscarpbad
    end

hold all;
% scale1=5e-3;
% scale1=2e-3;%glacierWVwsize512zone50cm1.tif
% scale1=15e-2;%glacierasar.tif'
% scalexy=[(info.x(1)+(info.x(end)-info.x(1))*0.05)*1e-3,(info.y(1)+info.y(end))/2*1e-3,(0)*scale1,(10)*scale1];

%selected vectors with good directions, and GOOD vectors within scarp
quiver(info.x(1:dn:end)*1e-3,info.y(1:dn:end)*1e-3,(vxc(1:dn:end,1:dn:end))*scale1,(vyc(1:dn:end,1:dn:end))*scale1,0,'r','linewidth',2);
% quiver(-2357.2,621,(0)*scale1,(2)*scale1,0,'r','linewidth',3) %2m scale
% hold all;quiver((info.x(1)+(info.x(end)-info.x(1))*0.05)*1e-3,(info.y(1)+info.y(end))/2*1e-3,(0)*scale1,(10)*scale1,0,'r','linewidth',3) %2m scale
hold all;quiver(scalexy(1)*1e-3,scalexy(2)*1e-3,(scalebar)*scale1,(0)*scale1,0,'k','linewidth',2,'MaxHeadSize',5) %2m scale
axis equal
title(['displacement vector: ',text1])

if exist('iorder','var')
    title(['displacement vector: ',text1,';i=',num2str(iorder)])
end
% title(['displacement vector: ',ifile])
saveas(gcf,[texto,'vector'],'fig')


end

%use all velocity for 2d map; plot the bad ones;
% vxc=D(:,:,1);
% vyc=D(:,:,2);

if 0% && flagplot==1
    %Filtering criteria: control surface: std of angle< 60degree &(displacement - its mean )<3std; 
    %scarp area keep displacement that has (angle - mean angle)<20;
figure;imagesc(info.x*1e-3,info.y*1e-3,sqrt(vxc.^2+vyc.^2),'alphadata',~isnan(vxc));%,'alphadata',Ds<30);
colorbar;colormap jet
axis equal
view(0,-90);title(['displacement magnitude: ',ifile])
xlabel('x (km)');ylabel('y (km)');
end

if 0 %1 %flagplot==1
dn=2;
%for paper; 
%Filtering criteria: constrol surface: ds < 20m; scarp area: same as above, (angle - mean angle)<20; 
figure
M1=(sqrt(vxcorg.^2+vycorg.^2))<20&(Mrock&~Mscarp);%control area has to < 20m
hold all;imagesc(info.x*1e-3,info.y*1e-3,sqrt(vxcorg.^2+vycorg.^2),'alphadata',((Mscarp&~Mscarpbad)|M1));%,'alphadata',Ds<30);
colorbar;colormap jet
axis equal
view(0,-90);title(['displacement magnitude: ',ifile])
xlabel('x (km)');ylabel('y (km)');
caxis([0 5])
ofile=['./pics/',texto,'vs'];
print('-dpng','-r400',ofile)
%filter vector the same as M2
M2=((Mscarp&~Mscarpbad)|M1);
vycorg2=vycorg;vxcorg2=vxcorg;vxcorg2(~M2)=nan;vycorg2(~M2)=nan;
quiver(info.x(1:dn:end)*1e-3,info.y(1:dn:end)*1e-3,(vxcorg2(1:dn:end,1:dn:end))*scale1,(vycorg2(1:dn:end,1:dn:end))*scale1,0,'r','linewidth',2);
saveas(gcf,['./pics/',texto,'vs'],'fig')

end

%all data over scarp and control surface, no filter;
if 0%flagplot==1
figure;imagesc(info.x*1e-3,info.y*1e-3,sqrt(vxcorg.^2+vycorg.^2),'alphadata',Mcontrol|Mscarp);%,'alphadata',Ds<30);
colorbar;colormap jet
axis equal
view(0,-90);title(['displacement magnitude: ',ifile])
xlabel('x (km)');ylabel('y (km)');
caxis([0 60])
% axis([-3065 -3062.5 1037.2 1041]);view(0,-90)%fig.3
% caxis([0 20])

% figure;histogram(sqrt(vxc.^2+vyc.^2))
% xlabel('displacement (m)');title(ifile)
end

% out.x=info.x; out.y=info.y; out.v=sqrt(vxc.^2+vyc.^2);
% save displacement.mat out

%save it to the structure instead
Co.i=iorder;
Co.datem=(datem); Co.dates=(dates);
Co.txcstd=txcstd; Co.tycstd=tycstd; 
Co.vxc=vxc_out; Co.vyc=vyc_out; Co.Mf=Mf;
Co.x=info.x;Co.y=info.y;

if ~isempty(scarp)% exist('scarp','var')

    %1 within the scarp
%     Mscarp=interp2(scarp.x,scarp.y,double(scarp.z), info.x,info.y','*nearest',0);

    %ratio of good points within the scarp;
    npts=sum(Mscarp(:));
    ngood=sum(~isnan(vxc(Mscarp)));%all good points plotted; %sum(sum(~Mf&Mscarp&~Mscarpbad));
    ratio=ngood/npts*100;
    
    %apply threshold for good and bad cases
    %good case: 69% 65% 30% 61% 10% 24% 36%
    %30% 18%
  
    %bad case: 
    if ratio<30
        multi=3; %larger uncertainty
    else
        multi=1;
    end
    
    %estimate the mean dx dy over the scarp area; 
    dx=nanmedian(vxc(Mscarp)); %nanmean(vxc(Mscarp))?
    dy=nanmedian(vyc(Mscarp)); %nanmean(vyc(Mscarp))
    
    if flagscreen==1
    fprintf(['\n [dx, dy, stdangle, ratioscarp]=',num2str([dx,dy,round(stdangle),round(ratio)]),'%% \n'])
    end
    
    if 0 % save it to file
    file='scarpdxdy.txt';
    fid10 = fopen('scarpdxdy.txt','a');
    % master image date, slave image date, dx, dy, ratio of good points
    % within the scarp;
    fprintf(fid10,'%d %d %12.6f %12.6f %12.6f %12.6f %12.6f \n',str2double(datem), str2double(dates), dx, dy, txcstd, tycstd, ratio);
%     fprintf(fid3,'%12.6f %12.6f  %23.15e %23.15e   %23.15e %23.15e\n',LAT(jy,jx),LON(jy,jx),trend,trest,rate,ratestd);

    fclose(fid10);
    end
    
     %save it to the structure instead
    Co.dx=dx; Co.dy=dy; Co.ratio=ratio;
    
%     pause
else
    Co.dx=[]; Co.dy=[]; Co.ratio=[];

end


close all

return
end
