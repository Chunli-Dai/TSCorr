function [datao]=getdstime_add(filename)
%get the time series of displacement from adding pairs of relative
%displacement.
%file: file that have the relative displacement of all pairs of images;
% refers to adjustOffsets.m

mindate=0;%datenum('2015/01/01');%0; %datenum('2013/01/01');%0;% %delete data befor this date
maxdate=9e9;%datenum('2015/01/01');;%datenum('2015/10/17'); %9e9;delete data after this date
deldate=[datenum('2010/05/31'),datenum('2012/10/25'),datenum('2015/06/17'),datenum('2013/06/06'),datenum('2017/02/03')]; 
%0;datenum('2000/04/12');%datenum('2008/06/21'); %date to be deleted
deldate=[datenum('13-Feb-2003')];
% deldate=[datenum('2018/10/31')]; %for compare_dstime.m
deldate=[];
thres_del=1;%30;%6; %1 ; if ratio < thres_del, delete this data;
thres_less=30;%50;%30; % 30; if ratio < thres_less, give less weight to this data;
% thres_del=30; thres_less=50; %v4

fid10 = fopen(filename,'r'); %file='scarpdxdy.txt';
n = linecount(fid10);
fclose(fid10);

datem=zeros(n,1);dates=zeros(n,1);
dx=zeros(n,1);dy=zeros(n,1);
dxe=zeros(n,1);dye=zeros(n,1);
flag=ones(n,1);

%read mannually marked flag file
%recommend: delete pairs with ratio < 30; bad weight for 30< ratio < 50;
%good for ratio >50;
filename2='corrflag.txt';
if 0    % exist(filename2,'file')
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
    multi=1;
    
    dx(k)=data1(3);dy(k)=data1(4);dxe(k)=data1(5)*multi;dye(k)=data1(6)*multi;
    
    if ratio<thres_del || isnan(data1(3)) || isnan(data1(4))||flag(k)==0 ||(datem(k)<=mindate||dates(k)<=mindate)||(datem(k)>=maxdate||dates(k)>=maxdate)
%     if ratio<1 || isnan(data1(3)) || isnan(data1(4))||flag(k)==0 ||(datem(k)<=mindate||dates(k)<=mindate)
        %ratio<1
        idd=[idd(:);k];
    end
    
    for mi=1:length(deldate)
        if (datem(k)==deldate(mi)||dates(k)==deldate(mi))
            idd=[idd(:);k];
        end
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

if 1
%delete these pairs that only one image j for each image i;
threspairj=1; %if the total number of pairs for images j  <= threspairj, delete these pairs;
Ndems=max([offsets.j(:);offsets.i(:)]);
iall=1:Ndems;
idd=[];
for k=1:length(iall)
    M=(offsets.i==iall(k))|(offsets.j==iall(k));%this dem appears as master or slave;
    if sum(M)<=threspairj %delete these pairs
        idd=[idd;find(M)];
    end
end
offsets.i(idd)=[];offsets.j(idd)=[];offsets.dx(idd)=[];
offsets.dy(idd)=[];offsets.dxe(idd)=[];offsets.dye(idd)=[];

%delete rebundant pairs; only use one pair i, i+1, delete all other pairs
idd=[]; %ones to keep
for k=1:length(iall)
    M=(offsets.i==iall(k))&(offsets.j==iall(k)+1);%this dem appears as master or slave;
    
    if 0 %manually delete  '26-Mar-2013' '31-Oct-2018' for compare_dstime.m
         %not a good way, better to use deldate
        if k==7; M=(offsets.i==7 & offsets.j==9);end
        if k==8; M=[];end
        if k==28; M=(offsets.i==28 & offsets.j==30);end
        if k==29; M=[];end
    end 
    
    idd=[idd;find(M)];
end
iall2=1:length(offsets.i);
idd=iall2(~ismember(iall2,idd));
% idd(end)=[]; %manually save one more pair for least-sqaures

offsets.i(idd)=[];offsets.j(idd)=[];offsets.dx(idd)=[];
offsets.dy(idd)=[];offsets.dxe(idd)=[];offsets.dye(idd)=[];

end

%% The following is from adjustOffsets.m

% make sure all fields are column vectors
offsets=structfun( @(x) x(:),offsets,'UniformOutput',false);

% get the number of DEMs in the stack from the max j of i-j pairs
% Ndems= max(offsets.j);%
Ndems= max([offsets.j(:);offsets.i(:)]); %chunli use this.-> if the input
% has i>j combination (the above line causes error); This works better.

Npairs = length(offsets.dx);

% find DEMs with no acceptable pairs; parameters
i_missing = setdiff(1:Ndems,unique(offsets.i));
j_missing = setdiff(1:Ndems,unique(offsets.j));
n_missing = intersect(i_missing,j_missing);

% Build design and weight matrices
A = zeros(Npairs,Ndems); % initialize design matrix

linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsets.i);
A(linearInd) = -1;
linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsets.j);
A(linearInd) = +1;

% remove filtered dems;; parameters
A(:,n_missing) = [];

% add dx(t1)=0 lines
[na,ma]=size(A);
At1=zeros(1,ma);At1(1)=1;
% A = [A;At1];

dx = [offsets.dx];
dy = [offsets.dy];

dxe = [offsets.dxe;]; 
dye = [offsets.dye;];

% wx = 1./dxe.^2;
% wy = 1./dye.^2;
 %chunli: January 2020, use Pc=sqrt(P), then the results would be equivalent to inv(A'PA)*A'Py;
% wx = 1./dxe;
% wy = 1./dye;

n = 1:Ndems;
n(n_missing) = [];

dX = nan(Ndems,1);
dY = nan(Ndems,1);

% dX(n) = (wx.*A)\(wx.*dx);
% dY(n) = (wy.*A)\(wy.*dy);

% The uncertainties
lenf=length(dx);mp=length(n);
for ij=1:2
    if ij==1
        P=diag(1./dxe.^2);
        y=dx;
    elseif ij==2
        P=diag(1./dye.^2);
        y=dy;
    end

    %fixed constraint;
    k0=0; K=At1;
    ml=1;%rank of constraint 

    Nmat=A'*P*A;KK=K'*K;c=A'*P*y;
    inv1=inv(Nmat+KK);
    ksi=inv1*c+inv1*K'*inv(K*inv1*K')*(k0-K*inv1*c);
    etilde=y-A*ksi;
    sigma02hat=1;%etilde'*P*etilde/(lenf-mp+ml); 
    Dksi=sigma02hat*(inv1-inv1*K'*inv(K*inv1*K')*K*inv1);

    if ij==1
    dX(n) = ksi;
    var_x=Dksi;
    elseif ij==2
    dY(n) = ksi;
    var_y=Dksi;
    end

end

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