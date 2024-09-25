%simple linear fit function, modified from /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/plotcosicorscarp_main.m
function [rate,ratestd]=linearfit(time,ds,dse,flagplot)
%modified from 

% flagplot=0;
if nargin==3
    flagplot=0;
elseif nargin==4
    %use the given flagplot
end

% it can be 1e-9; 
%if dse(1)==0;dse(1)=median(dse); end%assign the median weight to correct the zero weight for the first measurements
%if dse(1)<1;dse(1)=median(dse); end%assign the median weight to correct the zero weight for the first measurements
M=dse<1;dse(M)=median(dse); %some of the point in the middle can have very small std, which skew the linear fit.

% load dstime_aster.mat; dataast=datao;
% load dstimelandsat_sv3woflagwo20000412.mat; dataland=datao;
% 
% ds2=interp1(dataland.time,dataland.ds,dataast.time,'linear',NaN);
% dx2=interp1(dataland.time,dataland.dx,dataast.time,'linear',NaN);
% dy2=interp1(dataland.time,dataland.dy,dataast.time,'linear',NaN);
% Mt=dataast.time>datenum('2002/03/10')&dataast.time<datenum('2014/09/26');
% ds2(Mt)=nan;dx2(Mt)=nan;dy2(Mt)=nan;
% % hold all;plot(dataast.time,ds2,'b+')
% offsetds=nanmean(dataast.ds-ds2); %-7.55m;
% offsetdx=nanmean(dataast.dx-dx2); %-9.04m;
% offsetdy=nanmean(dataast.dy-dy2); %0.03m;
% % dataast.ds=dataast.ds-offsetds;
% dataast.dx=dataast.dx-offsetdx;
% dataast.dy=dataast.dy-offsetdy;
% dataast.ds=sqrt(dataast.dx.^2+dataast.dy.^2); %see /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/getdstime.m

% time=[dataast.time(:);dataland.time(:)];ds=[dataast.ds(:);dataland.ds(:)];
% dse=[dataast.dse(:);dataland.dse(:)];

%linear fit
time1=datenum('2008/06/21');%seg1, t<=time1
time2=datenum('2013/09/15');%seg2, time1<t<=time2;%seg3 t> time2;
time3=datenum('2016/10/25');
for seg=1:1 %1:4
    if seg==1 %'31-Jul-1999' and '21-Jun-2008'
%         M=time<=time1;
        M=true(size(time));
    elseif seg==2  %'26-May-2010' '15-Sep-2013'
        M=time>time1&time<=time2;
    elseif seg==3  % '19-Mar-2014'  '25-Oct-2016'
        M=time> time2&time<=time3;
    elseif seg==4  %'11-Mar-2017'  '10-Mar-2020'
        M=time>time3;
    end
    epochorg=time(M);T6=ds(M);
    T6std=dse(M); %ones(size(T6));%dse(M); %mitigate outliers;
%     T6std=ones(size(T6));%dse(M);
    
    %fit y=a+bt;
    yr=365.25;
    formatSpec = '%6.1f';

    mp=2;
    lenf=length(epochorg);
    yobs=T6;
    epoch=epochorg/yr;
    tm=mean(epoch);
    P=diag(T6std(:).^-2); 


    AM=zeros(lenf,mp);  %
    AM(:,1)=1.;
    AM(:,2)=epoch-tm; %unit year

    var=inv(AM'*P*AM);
    est=var*AM'*P*yobs;
    etilde=yobs-AM*est;
    sigma02hat=etilde'*P*etilde/(lenf-mp); %mp=2
    fprintf(['\n seg, sqrt(sigma02hat):',num2str([seg, sqrt(sigma02hat)]),'\n'])
    %4m, 8m,5m, 3m
    var=var*sigma02hat;
      
    rate=est(mp);ratestd=sqrt(var(mp,mp));
    fit=AM*est;

    if 0 %flagplot==1
    hold all;
    epochp=epoch*yr; % change unit to day
    plot(epoch*yr,fit,'r-','linewidth',2,'markersize',14)
    text(max(epochp)-(max(epochp)-min(epochp))*0.3,max(T6)-(max(T6)-min(T6))*0.4 ,['Trend=',num2str(rate,formatSpec),'\pm',num2str(ratestd,formatSpec),'m/yr'],'FontSize',12)
    end
end

ratei=rate;

    if flagplot==1
epochp=epoch*yr; % change unit to day
figure;
set(gcf,'Color','white')
set(gca,'FontSize', 14);
set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
errorbar(time,ds,dse,'r.')
% legend('Aster','Landsat-7 and 8','Location','Northwest')
datetick('x','mm/yyyy')
% eptick=datenum({'2000/07/01','2002/07/01','2013/07/01','2015/07/01','2017/07/01','2019/07/01'});
% eptick=datenum({'2013/07/01','2015/07/01','2017/07/01','2019/07/01','2021/07/01'});
% set(gca,'XTick',eptick)
% set(gca,'XTickLabel',eptick,'FontSize',12,'XTickLabelMode','manual')
% datetick('x','mm/yyyy','keepticks')
box on
ylabel('Displacement (m)')
plot(epoch*yr,fit,'r-','linewidth',2,'markersize',14)
text(max(epochp)-(max(epochp)-min(epochp))*0.3,max(T6)-(max(T6)-min(T6))*0.4 ,['Trend=',num2str(rate,formatSpec),'\pm',num2str(ratestd,formatSpec),'m/yr'],'FontSize',12)

if 0
figure;
hold all
plot(datao.time,datao.dx,'g.-',datao.time,datao.dy,'b.-',datao.time,datao.ds,'r.-')
shadedErrorBar(datao.time,datao.dx,datao.dxe,'g',1)
shadedErrorBar(datao.time,datao.dy,datao.dye,'b',1)
shadedErrorBar(datao.time,datao.ds,datao.dse,'r',1)
legend('dx','dy','ds')
datetick('x','mm/yyyy')
eptick=datenum({'2000/07/01','2002/07/01','2013/07/01','2015/07/01','2017/07/01','2019/07/01'});
set(gca,'XTick',eptick)
set(gca,'XTickLabel',eptick,'FontSize',12,'XTickLabelMode','manual')
datetick('x','mm/yyyy','keepticks')
end

    end
    
return
end
