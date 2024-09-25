function [rate, ratestd, ratex, ratexstd, ratey, rateystd]= getratesub(io,xout,yout,flagscarp,controlo,xq1,yq1,scarp,info,control,scarpdxdy);
%control is used.

constant

flageq2=0; % default, not earthquake
if exist('flageq', 'var') == 1
    flageq2=flageq;
end

nq1=length(xq1);
resr=median(diff(xout));
% winsize=resr_cosi*7;%window size in meters for getting the motion rate:7by7 pixels (output pixel from COSI-Corr);
% winsize=resr*7;%window size in meters for getting the motion rate:7by7 pixels (output pixel from COSI-Corr);
winsize=7; %window size in pixels (pixel size of scarp.x).
n=length(scarpdxdy);

%rate=zeros(length(yout),length(xout));
rate=zeros(length(yout),1);
ratestd=rate;
ratex=rate;ratey=rate;ratexstd=rate;rateystd=rate;

    for jo=1:length(yout)
%         rate(j,i)=;
%target pixel
            xp=xout(io);yp=yout(jo);

            if flagscarp==3
            jq=scarpid;
            xp=xq1(jq); yp=yq1(jq);
            end

            if flag_vel_iw==1 
                % calculate all points
            elseif flag_vel_iw==0 % calculate velocity only on land, NOT on ice or water;
                if controlo.z(jo,io)~=1 && (flagscarp==0)  %only work on rock surfaces; skip ice/water
                    continue;
                end
            end

            %an array of requested locations for plotting time series
            flagxp=0;
            for jq=1:nq1
                jxq1=find(abs(xout-xq1(jq))<=resr/2);jyq1=find(abs(yout-yq1(jq))<=resr/2);
                if ismember(io,jxq1)&&ismember(jo,jyq1) %jx=jxq1(1), jy=jyq1(2); not good
                   fprintf(['\n xq1(jq) with jq:',num2str(jq),'; [lenx leny]=',num2str([length(jxq1) length(jyq1)]),';'])
                   flagxp=1;break;
                end
            end

            if flagxp==1 %plot
                flagplotsv=1;
            else %not plot
                flagplotsv=0;
                jq=0;
            end %0 %plotting

            % flag_pts=1; % 1, only calculate time series of given points; 0, calculate velocity for all grid points.
            if flag_pts==1 && flagxp ==0
                    continue;
            end

%Get a dense grid with the same resolution as scarpdxdy(i).x or resr_cosi, but smaller extent than xout yout.
xouthigh=max([xout(1),xp-50*resr_cosi]):resr_cosi:min([xout(end),xp+50*resr_cosi]);
youthigh=min([yout(1),yp+50*resr_cosi]):-resr_cosi:max([yout(end),yp-50*resr_cosi]);

%Get the scarp area; 
% add option B: read input shapefile
if flagscarp==1 || flagscarp==2 %~=0 %use scarp shapefile
%[scarp,scarpx,scarpy]=getscarp(scarp.x,scarp.y);
resrhigh=3; % 3 m
rang0=[min(scarp.x) max(scarp.x) min(scarp.y) max(scarp.y)];
scarphigh.x=rang0(1):resrhigh:rang0(2); scarphigh.y=rang0(4):-resrhigh:rang0(3); %use a high resolution grid for given shapefiles
[scarp,scarpx,scarpy]=getscarp(scarphigh.x,scarphigh.y);
else
%recreate scarp
scarp.x=xouthigh; scarp.y=youthigh; 
scarp.z=zeros(length(scarp.y),length(scarp.x));
resr_scarp=median(diff(scarp.x));

iop=find(abs(scarp.x-xp)<=resr_scarp/2);jop=find(abs(scarp.y-yp)<=resr_scarp/2);
scarp.z(jop,iop)=1;
%buffer by window size
% M=(imdilate(scarp.z,ones(round(winsize/resr_cosi))));
M=(imdilate(scarp.z,ones(round(winsize))));
scarp.z=M;

%scarp outline 
B = bwboundaries(M); 
k=1;scarpx=scarp.x(B{k}(:,2)); scarpy=scarp.y(B{k}(:,1));
hold on;plot(scarpx*1e-3,scarpy*1e-3,'m.-') %scarp outline in outlinexy.mat
end

% The scarp area &rock area need to have: >= 10 pixels
%info=Dall{1}.info; %Give errors if Dall{1} is empty.
% Mscarp=interp2(scarp.x,scarp.y,double(scarp.z), info.x,info.y','*nearest',0);
Mscarp=interp2(scarp.x,scarp.y,double(scarp.z), xouthigh,youthigh','*nearest',0);
% Mrock=interp2(controlo.x,controlo.y,double(controlo.z),xout,yout','*nearest',0); %coarse grids
Mrock=interp2(control.x,control.y,double(control.z), xouthigh,youthigh','*nearest',0);

if flag_vel_iw==1 
    % calculate all points
    Mscarp=logical(Mscarp);
    %To avoid edge effect from the smoothing window in this step.
    %if the target pixel is on land, select only pixels on land within the window for next steps.
    %if the target pixel is on ice/water, select only pixels on ice/water.
    if controlo.z(jo,io) ==1
        Mscarp=logical(Mscarp)&logical(Mrock); 
    elseif controlo.z(jo,io) ==0
        Mscarp=logical(Mscarp)&logical(~Mrock); 
    end
elseif flag_vel_iw==0 % calculate velocity only on land, NOT on ice or water;
    Mscarp=logical(Mscarp)&logical(Mrock); %new August 2020; only calculate the motion of scarp area within rock area, excluding ice/water.
end

if 0
%remove the Barry Arm landslide area;
Mscarpt1=interp2(scarp.x,scarp.y,double(scarp.z), control.x,control.y','*nearest',0);
control.z(Mscarpt1==1)=0;
end


if sum(Mscarp(:)) <10 %only work on rock surfaces, skip window boxes with little valid pixels. 
continue
end

fprintf(['\n Processing io jo:',num2str([io jo]),' /(',num2str([length(xout) length(yout)]),'). \n'])
tic

%load data from scarpdxdy
% scarpdxdy_io: same as scarpdxdy, except the dx dy ratio is for grid io.
scarpdxdy_io(n)=struct('i',[],'datem',[],'dates',[],'txcstd',[],'tycstd',[],'vxc',[],'vyc',[],'Mf',[],'dx',[],'dy',[],'ratio',[],'x',[],'y',[]); %ratio  

for i=1:n
    if isempty(scarpdxdy(i).i)||isnan(scarpdxdy(i).txcstd)
        continue;  %scarpdxdy_io(i).i=[];
    end
    % If the correlation map does not cover the target pixel (xp, yp), skip it.
    infox=scarpdxdy(i).x;infoy=scarpdxdy(i).y;
    if ~(min(infox)<xp && xp<max(infox) && min(infoy)<yp && yp < max(infoy))
        continue
    end

    % Get the matrix subset to speed up interp2. Reduce computation time for each point from 10 sec to 0.5 sec.
    Mx=infox>=(min(xouthigh)-resr)&infox<=(max(xouthigh)+resr);
    My=infoy>=(min(youthigh)-resr)&infoy<=(max(youthigh)+resr);

%     vxc=scarpdxdy(i).vxc; vyc=scarpdxdy(i).vyc; Mf=scarpdxdy(i).Mf;
    vxc=interp2(infox(Mx),infoy(My),scarpdxdy(i).vxc(My,Mx), xouthigh,youthigh','*linear',nan);
    vyc=interp2(infox(Mx),infoy(My),scarpdxdy(i).vyc(My,Mx), xouthigh,youthigh','*linear',nan);
    Mf=interp2(infox(Mx),infoy(My),scarpdxdy(i).Mf(My,Mx), xouthigh,youthigh','*nearest',1); % Mf 1 bad, 0 good.

    scarpdxdy_io(i)=scarpdxdy(i); 
    scarpdxdy_io(i).vxc=[];scarpdxdy_io(i).vyc=[]; scarpdxdy_io(i).Mf=[]; %release space
    %within the scarp, map out displacement that has direction 20 away from
    %mean direction;
    if 1
    angle=atan2(vyc,vxc)*180/pi;%update the angle after txc tyc correction!!
    meanangle=nanmedian(angle(Mscarp&~Mf));
    stdangle=nanstd(angle(Mscarp));
    Mscarpbad=~(abs(angle-meanangle)<=20)&Mscarp;  %;label nans as bad ones too; abs(angle-meanangle)>20&Mscarp; 
    vxc(Mscarpbad)=nan;vyc(Mscarpbad)=nan;%
    end
    
        %ratio of good points within the scarp;
    npts=sum(Mscarp(:));
    ngood=sum(~isnan(vxc(Mscarp)));%all good points plotted; %sum(sum(~Mf&Mscarp&~Mscarpbad));
    ratio=ngood/npts*100;
      
    %estimate the mean dx dy over the scarp area; 
    dx=nanmedian(vxc(Mscarp)); %nanmean(vxc(Mscarp))?
    dy=nanmedian(vyc(Mscarp)); %nanmean(vyc(Mscarp))
    
    %save it to the structure 
    scarpdxdy_io(i).dx=dx; scarpdxdy_io(i).dy=dy; scarpdxdy_io(i).ratio=ratio;

end


     if 0 % flagplotsv ==1
        ofile=['scarpdxdy_jq',num2str(jq),'.txt'];
        [status, cmdout]=system(['cp scarpdxdy.txt ',ofile]);
%        save(ofile,'demp','dempmt', 'epochorg','timefix','eqepoch','flagplotsv');%let jq=0
      end

try
file='scarpdxdy.txt';
[datao]=getdstime(file,scarpdxdy_io);

%get the rate;
time=[datao.time(:);];ds=[datao.ds(:);];
dse=[datao.dse(:);];
if flageq2==0
%[ratei,rateistd]=linearfit(time,ds,dse);
[rateix,rateixstd]=linearfit(datao.time,datao.dx,datao.dxe);
[rateiy,rateiystd]=linearfit(datao.time,datao.dy,datao.dye);
elseif flageq2==1
    timefix=1;
    epochorg=datao.time;dempmt=ones(size(datao.time));
    demp=datao.dx;
    [oflag,rateix,rateixstd,ratet1,ratestdt1,eqm,eqs,eqe]=timeseries2(demp,dempmt,epochorg,timefix,eqepoch,flagplotsv,jq);
    demp=datao.dy;
    [oflag,rateiy,rateiystd,ratet1,ratestdt1,eqm,eqs,eqe]=timeseries2(demp,dempmt,epochorg,timefix,eqepoch,flagplotsv,jq);
end
ratei=sqrt(rateix^2+rateiy^2);
dfdx=rateix./ratei;
dfdy=rateiy./ratei;
rateistd=sqrt(dfdx.^2.*rateixstd.^2+dfdy.^2.*rateiystd.^2); %assume no correlation between dx dy

%rate(jo,io)=ratei;ratestd(jo,io)=rateistd;
%ratex(jo,io)=rateix;ratexstd(jo,io)=rateixstd;
%ratey(jo,io)=rateiy;rateystd(jo,io)=rateiystd;
rate(jo)=ratei;ratestd(jo)=rateistd;
ratex(jo)=rateix;ratexstd(jo)=rateixstd;
ratey(jo)=rateiy;rateystd(jo)=rateiystd;
catch e
    datao=[];time=[];
    fprintf('getdstime.m There was an error! The message was:\n%s',e.message);
end

fclose all;
close all;
toc


if (flagplotsv ==1 || flagscarp~=0) &&~isempty(datao)
    ofile=['dstime_jq',num2str(jq),'.mat'];
    %save dstime_mono.mat datao
    save(ofile,'datao');
    % out=[datestr(datao.time,26), datao.ds,datao.dse];
    file=['dstime_jq',num2str(jq),'.txt'];
    fid10 = fopen(file,'w');
    for i=1:length(time)
    fprintf(fid10,'\n %s %12.6f %12.6f',datestr(datao.time(i),26), datao.ds(i),datao.dse(i));
    end
    fclose(fid10);

    figure;
    hold all
    plot(datao.time,datao.dx,'g.-',datao.time,datao.dy,'b.-',datao.time,datao.ds,'r.-')
    shadedErrorBar(datao.time,datao.dx,datao.dxe,'g',1)
    shadedErrorBar(datao.time,datao.dy,datao.dye,'b',1)
    shadedErrorBar(datao.time,datao.ds,datao.dse,'r',1)
    legend('dx','dy','ds')
    datetick('x','mm/yyyy')
    %eptick=datenum({'2000/07/01','2002/07/01','2013/07/01','2015/07/01','2017/07/01','2019/07/01'});
    %set(gca,'XTick',eptick)
    %set(gca,'XTickLabel',eptick,'FontSize',12,'XTickLabelMode','manual')
    %datetick('x','mm/yyyy','keepticks')
    box on
    ylabel('Cumulative Displacement (m)')
    ofile=['dstime_jq',num2str(jq)];
    saveas(gcf,ofile,'fig')
    print('-dpng','-r400',ofile)

    %error bar plots
    [rateix,rateixstd]=linearfit(datao.time,datao.dx,datao.dxe,1);
    ylabel('Cumulative Displacement x (m)')
    ofile=['dxtime_jq',num2str(jq)];
    saveas(gcf,ofile,'fig')

    [rateiy,rateiystd]=linearfit(datao.time,datao.dy,datao.dye,1);
    ylabel('Cumulative Displacement y (m)')
    ofile=['dytime_jq',num2str(jq)];
    saveas(gcf,ofile,'fig')

    [ratei,rateistd]=linearfit(time,ds,dse,1);
    ylabel('Cumulative Displacement (m)')
    ofile=['dstime_jq_errbar',num2str(jq)];
    saveas(gcf,ofile,'fig')
    
    close all

end  % if flagplotsv ==1
%%

    if flagscarp~=0 % Only run on scarp
        return
    end

    end%jo
    
return
end
