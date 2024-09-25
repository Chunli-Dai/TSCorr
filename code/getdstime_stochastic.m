function [datao]=getdstime(filename)
%get the time series of displacement from all pairs of relative
%displacement.
%file: file that have the relative displacement of all pairs of images;
% refers to adjustOffsets.m

mindate=0; %datenum('2013/01/01');%0;% %delete data befor this date
deldate=datenum('2000/04/12');%datenum('2008/06/21'); %date to be deleted
thres_del=1; %1 ; if ratio < thres_del, delete this data;
thres_less=30; % 30; if ratio < thres_less, give less weight to this data;
% thres_del=30; thres_less=50; %v4

fid10 = fopen(filename,'r'); %file='scarpdxdy.txt';
n = linecount(fid10);

datem=zeros(n,1);dates=zeros(n,1);
dx=zeros(n,1);dy=zeros(n,1);
dxe=zeros(n,1);dye=zeros(n,1);
flag=ones(n,1);

%read mannually marked flag file
%recommend: delete pairs with ratio < 30; bad weight for 30< ratio < 50;
%good for ratio >50;
filename2='corrflag.txt';
if exist(filename2,'file')
fid2=fopen(filename2);
ifile=fgetl(fid2);
for k=1:n
    data1=fscanf(fid2, '%f', [4, 1])';ifile=fgetl(fid2);
    flag(k)=data1(2); %stdangle(k)=data1(3);ratio2(k)=data1(4);
    ind(k)=data1(1);
end
fclose(fid2)
else %do nothing
end

fid10 = fopen(filename);
idd=[];
for k=1:n
    data1=fscanf(fid10, '%f', [7, 1])';ifile=fgetl(fid10);
    str1=num2str(data1(1));
    year=str1(1:4);mon=str1(5:6);day=str1(7:8);
    datem(k)=datenum([year,'/',mon,'/',day]);%master
    str1=num2str(data1(2));
    year=str1(1:4);mon=str1(5:6);day=str1(7:8);
    dates(k)=datenum([year,'/',mon,'/',day]); %slave
    ratio=data1(7);
    
    %bad case: 
%     if ratio<30||flag(k)==0.5
    if ratio<thres_less||flag(k)==0.5
        multi=100; %3; %larger uncertainty
    else
        multi=1;
    end
    
    dx(k)=data1(3);dy(k)=data1(4);dxe(k)=data1(5)*multi;dye(k)=data1(6)*multi;
    
    if ratio<thres_del || isnan(data1(3)) || isnan(data1(4))||flag(k)==0 ||(datem(k)<=mindate||dates(k)<=mindate)||(datem(k)==deldate||dates(k)==deldate)
%     if ratio<1 || isnan(data1(3)) || isnan(data1(4))||flag(k)==0 ||(datem(k)<=mindate||dates(k)<=mindate)
        %ratio<1
        idd=[idd(:);k];
    end

end
fclose(fid10);

datem(idd)=[];dates(idd)=[];dx(idd)=[];dy(idd)=[];dxe(idd)=[];dye(idd)=[];

time=unique([datem(:);dates(:)]);

%construct offsets
%                     i: [n<D7>1 double] image 1 index
%                     j: [n<D7>1 double] image 2 index
%                    dx: [n<D7>1 double] x displacement (dem 2 - dem 1) 
%                    dy: [n<D7>1 double] y displacement (dem 2 - dem 1) 
%                   dxe: [n<D7>1 double] x formal uncertainty
%                   dye: [n<D7>1 double] y formal uncertainty

count=0; %record the count of  pairs;
clear offsets
for k=1:length(dx)
    i=find(time==datem(k));
    j=find(time==dates(k));
    
    %save results to offsets
    count=count+1;
    offsets.i(count)=i; %dem 1 index
    offsets.j(count)=j; %dem 2 index
    offsets.dx(count)=dx(k); % [n<D7>1 double] x offset (image 2 - image 1) 
    offsets.dy(count)=dy(k); % [n<D7>1 double] y offset (image 2 - image 1)
    offsets.dxe(count)=dxe(k); % [n<D7>1 double] x uncertianty
    offsets.dye(count)=dye(k); % [n<D7>1 double] y uncertianty

end

%% The following is from adjustOffsets.m

% make sure all fields are column vectors
offsets=structfun( @(x) x(:),offsets,'UniformOutput',false);

% get the number of DEMs in the stack from the max j of i-j pairs
% Ndems= max(offsets.j);%
Ndems= max([offsets.j(:);offsets.i(:)]); %chunli use this.-> if the input
% has i>j combination (the above line causes error); This works better.

Npairs = length(offsets.dx);

% find DEMs with no acceptable pairs
i_missing = setdiff(1:Ndems,unique(offsets.i));
j_missing = setdiff(1:Ndems,unique(offsets.j));
n_missing = intersect(i_missing,j_missing);

% Build design and weight matrices
A = zeros(Npairs,Ndems); % initialize design matrix

linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsets.i);
A(linearInd) = -1;
linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsets.j);
A(linearInd) = +1;

% remove filtered dems
A(:,n_missing) = [];

% add dx(t1)=0 lines
[na,ma]=size(A);
At1=zeros(1,ma);At1(1)=1;
% At1=zeros(1,Ndems);At1(1)=1;
A = [A;At1];

dx = [offsets.dx;0];
dy = [offsets.dy;0];

dxe = [offsets.dxe; 1];%suggest to change uncertainty of first measurements from 1 to 10? -> doesn't matter? 
dye = [offsets.dye; 1];

% wx = 1./dxe.^2;
% wy = 1./dye.^2;
 %chunli: January 2020, use Pc=sqrt(P), then the results would be equivalent to inv(A'PA)*A'Py;
wx = 1./dxe;
wy = 1./dye;

n = 1:Ndems;
n(n_missing) = [];

dX = nan(Ndems,1);
dY = nan(Ndems,1);

dX(n) = (wx.*A)\(wx.*dx);
dY(n) = (wy.*A)\(wy.*dy);

% The uncertainties
lenf=length(dx);mp=length(n);
P=diag(1./dxe.^2);
etilde=dx-A*dX;
ml=1;%rank of constraint 
sigma02hat_x=etilde'*P*etilde/(lenf-mp+ml); 
var_x=inv(A'*P*A)*sigma02hat_x;

P=diag(1./dye.^2);
etilde=dy-A*dY;
sigma02hat_y=etilde'*P*etilde/(lenf-mp+ml); 
var_y=inv(A'*P*A)*sigma02hat_y;

%uncertainties of estimated parameters;
dXe = nan(Ndems,1);
dYe = nan(Ndems,1);

dXe(n)=sqrt(diag(var_x));
dYe(n)=sqrt(diag(var_y));


datao.time=time;
datao.dx=dX;
datao.dy=dY;
datao.ds=sqrt(dX.^2+dY.^2);
datao.dxe=dXe;
datao.dye=dYe;

% dse1=sqrt(dXe.^2+dYe.^2); %to check;wrong
dfdx=dX./datao.ds;
dfdy=dY./datao.ds;
datao.dse=sqrt(dfdx.^2.*dXe.^2+dfdy.^2.*dYe.^2); %assume no correlation between dx dy

return
end