function [Co]=prep_rangerock(rang0)
constant

Co=[];
%addpath(genpath(['/home/dai.56/arcticdemapp/landslide/code1/']));

%get rang0
%Only work on subset of all images around the scarp;
%rang0=[431 445 6772 6784]*1e3; %UTM '6N' coordinates;
rang0=rang0+[-4e3 +4e3 -4e3 +4e3]; 
%z1 ='6N';

%utm to polar stereographic coordinates;

x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];

%Not for antarctica
if strcmp(projgdal,'epsg:3413')
	rang0_polar=rang0;
else
	fprintf(['\n xy2latlon using gdal calculating ',num2str(length(x0)),' points takes'])
	tic
%[lat,lon]=minvtran(utmstruct,x0,y0);
[lat,lon]=xy2latlon(x0,y0,projgdal);
toc
[Xb,Yb]=polarstereo_fwd(lat,lon,[],[],70,-45);
rang0_polar=round([min(Xb) max(Xb) min(Yb) max(Yb)]/40)*40;
end

%get rock mask 
%refer to /fs/project/howat.4/dai.56/chunliScripts/scripts/plotsrtmvs.m

[tag]=getrgi(rang0_polar); %rock 1 ; ice 0;

[water]=getocean(rang0_polar); %1 water; 0 land


if strcmp(projgdal,'epsg:3413')
	data=tag;
	water2.z=interp2(water.x,water.y,double(water.z),tag.x,tag.y','*linear',1);  %water;
else
fprintf('\n prep_rangerock.m Polar stereographic coordinates to UTM. \n')
tic
%data.x=rang0(1):40:rang0(2);data.y=rang0(4):-40:rang0(3);
data.x=rang0(1):10:rang0(2);data.y=rang0(4):-10:rang0(3);
[X,Y]=meshgrid(data.x,data.y);
size(X)
%[lat,lon]=minvtran(utmstruct,X,Y);
[lat,lon]=xy2latlon(X,Y,projgdal);
toc
[X,Y]=polarstereo_fwd(lat,lon,[],[],70,-45);
data.z=interp2(tag.x,tag.y,double(tag.z), X,Y,'*linear',0); %0 ice
water2.z=interp2(water.x,water.y,double(water.z), X,Y,'*linear',1); %1 water;
tag=data; %1 rock, 0 ice;

end

%get water mask;

icewater=data;
icewater.z=~tag.z|water2.z; % 1 ice or water; 0 rock
icewaterfile=[currentdir,'/icewater.mat'];
%save  icewater.mat   icewater tag water2   -v7.3  
save(icewaterfile, 'icewater','tag','water2','-v7.3');

%% write output
%projstr='polar stereo north';
projstr=projstrin; 
%projstr='UTM';%'Transverse_Mercator'; %'UTM_Zone_6N';
%zone='6N';
OutName=['icewater.tif']; % %1 is good (landslide); 0 is bad
writeGeotiff(OutName,icewater.x,icewater.y,int32(icewater.z),3,0,projstr)

return
end
