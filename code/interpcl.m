function [clx2,cly2,S2,id12]=interpcl(clx,cly,dc)
%Interpolation of centerline to any given interval along centerline distance.
% clx, cly input centerline nodes, in meter.
% dc input, requested centerline interval in meter
% or dc is a coordinate of a given point, which is assumed to be on the
% line (clx, cly).
% clx2, cly2, S2: output centerline nodes' coordinates and distance in meter
% id12, the index of first and second points that are closet to the target in input line (clx, cly);

clx=clx(:);cly=cly(:);% be a column rather than an array.

S = [0; cumsum(sqrt(diff(clx(:)).^2+diff(cly(:)).^2))];
lencl=max(S);% total centerline length

clx2=[];cly2=[];S2=[];id12=[];

[nrow,ncol]=size(dc);
if ncol==1
    
p.dx=dc*1e-3; %km
p.x=0:p.dx:lencl*1e-3; %km
p.x=p.x*1e3;%m

%For each new node, find the two original nodes that includes this new
%node.
kpre=1;clx3=[];cly3=[];idk3=[];
clx2=zeros(length(p.x),1);cly2=clx2;
for i=1:length(p.x)
    
    %find k that S(k)<= p.x(i) <S(k+1)
    [k, ~] = find(S <= p.x(i), 1, 'last'); % find index of middle centerline node

    if ~isempty(k)
        if k< length(S) && k>=1 
            pt1=[clx(k) cly(k)];pt2=[clx(k+1) cly(k+1)];

            ds12=sqrt((pt1(1)-pt2(1))^2+(pt1(2)-pt2(2))^2);
            ds12p=p.x(i)-S(k); %has to be >=0;
            if ds12p <0;printf(['Warning: interpcl.m: ds12p < 0, find the wrong nodes! \n']);end

            ratio=ds12p/ds12;
            dx=(clx(k+1)-clx(k))*ratio;
            dy=(cly(k+1)-cly(k))*ratio;
            clx2(i)=clx(k)+dx;cly2(i)=cly(k)+dy;
        else %k=length(S)
            if S(k)==p.x(i)
                clx2(i)=clx(k);cly2(i)=cly(k);
            else
                printf(['Warning: interpcl.m: k is not within a reasonable range! \n'])
            end
        end
    
    else
        printf(['Warning: interpcl.m: k is empty! \n'])
    end
    
    %for validation
    clx3=[clx3(:);clx(kpre+1:k);clx2(i)];
    cly3=[cly3(:);cly(kpre+1:k);cly2(i)];  
    kpre=k;
    idk3=[idk3;length(clx3)]; %id of the interpolated nodes;
end

S2b = [0; cumsum(sqrt(diff(clx2(:)).^2+diff(cly2(:)).^2))];
S3 = [0; cumsum(sqrt(diff(clx3(:)).^2+diff(cly3(:)).^2))];
S2=p.x(:); %cumulative centerline distance of new nodes along original centerline.
fprintf(['\n Check the new centerline node, maximum centerline distance difference S3-p.x:',num2str(max(abs(S3(idk3)-p.x(:)))),'\n']);

elseif ncol==2
    %given any point, find it's distance along the given profile
    x0=dc(1,1);y0=dc(1,2);
    clx2=x0;cly2=y0;
    
    % make clx to be ascending
    if mean(diff(clx)) >=0 %ascending order
        clxa=clx;clya=cly;Sa=S;
        flaga=1;
    else %descending order
       clxa=flip(clx(:));clya=flip(cly(:));
       Sa=S(end)-flip(S(:));
       flaga=0;
    end
        
    
    %find k that clx(k)<= x0 <clx(k+1) 
    [k, ~] = find(clxa <= x0, 1, 'last'); % find index of middle centerline node
    if isempty(k) %
        pt1=[clxa(1) clya(1)];
        ds12p=sqrt((pt1(:,1)-x0).^2+(pt1(:,2)-y0).^2); 
        s0=0-ds12p;
        id12=[1, 2]; %index of first and second point that are closet to target
    else
    if k<= length(clxa) && k>=1 
        pt1=[clxa(k) clya(k)];%pt2=[clxa(k+1) clya(k+1)];

        ds12p=sqrt((pt1(:,1)-x0).^2+(pt1(:,2)-y0).^2); 
        if (x0-pt1(1)) <0;printf(['Warning: interpcl.m: ds12p < 0, find the wrong nodes! \n']);end
        s0=Sa(k)+ds12p; 
        
        id12=[k, k+1];
    end
    end
    
    if flaga==1
        S2=s0;
    elseif flaga==0 %change back to original descending order to 
        S2=S(end)-s0;
        id12=length(S)-id12+1;
    end
%     figure;plot(clx,S,'.-',x0,S2,'>');title('clx S original')
%     figure;plot(clxa,Sa,'.-',x0,s0,'>');title('clx S ascending order')

end

return
end
