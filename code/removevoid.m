function [M1]=removevoid(rang0,XYbg,projgdalg,f,fdir);
% find the index ind of images that have all zeros.
% This may happen when using gdal to find the boundary of images, which may include zeros parts of polygons.

constant

M1=[];

ratio1=zeros(length(fdir),1);

x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
[lat0,lon0]=xy2latlon(x0,y0,projgdal);

for j=1:length(fdir)
   	ifile=[fdir{j},'/',f{j}];
	projgdalj=projgdalg{j}; % e.g.,'epsg:32637'

        [x0j,y0j]=latlon2xy(lat0,lon0,projgdalj);
	rang0j=[min(x0j) max(x0j) min(y0j) max(y0j)];
	data=readGeotiff(ifile,'map_subset',rang0j);
	M=data.z==0|isnan(data.z);
	ratio1(j)=sum(M(:))./length(data.z(:))*100; % percentage of void pixels.

end

thres=20;

if exist('coverthres','var')
    thres=coverthres;
end


M1=((100-ratio1)<thres); % image must have > 20% good pixels; find bad images;

fprintf(['\n removevoid.m ',num2str(sum(M1)),' out of ',num2str(length(M1)),' images have too much void pixels within the input tile.\n'])

return
end
