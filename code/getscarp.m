function [scarp,scarpx,scarpy]=getscarp(xout,yout)
% Usage [scarp,scarpx,scarpy]=getscarp(scarp.x,scarp.y);
% read input scarp shapefile or matrix
% output scarp matrix, scarp.x scarp.y scarp.z (1 inside scarp, 0 outside).
% scarpx, scarpy: outline

constant

scarp.x=xout;scarp.y=yout;
scarp.z=zeros(length(scarp.y),length(scarp.x));


%flagscarp=2; % 0 no scarp shapefile; 1 scarp in maxtrix in image coordinates; 2 scarp in lon lat as shapefile polygon;
%scarpfile='scarp.shp';

if flagscarp==1 % read scarp matrix
% /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/scarp.mat
scarpin=load(scarpfile);
scarpin=scarpin.scarp;
scarp.z=interp2(scarpin.x,scarpin.y,double(scarpin.z), scarp.x,scarp.y','*nearest',0);

elseif flagscarp==2 % read scarp shapefile
S2 = shaperead(scarpfile);

% figure;mapshow(S2);
ns=length(S2);

wm.x=xout; wm.y=yout;
nwy=length(wm.y);nwx=length(wm.x);
Msv=false(nwy,nwx);%1 land, 0, ocean
for j=1:ns
        XYbi=[S2(j).X(:),S2(j).Y(:)];
        lon1=XYbi(:,1);lat1=XYbi(:,2);

	[Xb,Yb]=latlon2xy(lat1,lon1,projgdal);
%       if strcmp(projgdal,'epsg:3413')
%       	[Xb,Yb]=polarstereo_fwd(lat1,lon1,[],[],70,-45);
%	else
%                [Xb,Yb]=mfwdtran(utmstruct,lat1,lon1);
%	end
        
        % check whether this polygon intersect with the block
	resr=mean((diff(wm.x)));
        idx=round((Xb-wm.x(1))/resr)+1;
        idy=round((Yb-wm.y(1))/(-resr))+1;
        
        %Oct 10, 2018:fix bug 14; separate polygons using the NaN;
        M=[idx(:),idy(:)];
        idx = any(isnan(M),2);
        idy = 1+cumsum(idx);
        idz = 1:size(M,1);
        C = accumarray(idy(~idx),idz(~idx),[],@(r){M(r,:)});

        for k=1:length(C)
            idx=C{k}(:,1);idy=C{k}(:,2);   
            Mb=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
            Msv=Msv|Mb;%lost the holes 
        end        
        
end
    
scarp.z=double(Msv); % 1 inside scarp, 0 outside

end % if flagscarp


%scarp outline 
M=scarp.z;
B = bwboundaries(M);
k=1;scarpx=scarp.x(B{k}(:,2)); scarpy=scarp.y(B{k}(:,1));
hold on;plot(scarpx*1e-3,scarpy*1e-3,'m.-') %scarp outline in outlinexy.mat

return
end
