
%xx: column 1:3: x, y, dx
load('frontposition_barry.mat')
%plot front position on the right axis
% left_color = [0 0 0]; %black
% right_color = [0 0 1]; %b
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
figure
hold all;
% yyaxis right
plot(idates,xx(p,3)-xx(p(1),3),'-ob','linewidth',2)
% set(gca,'fontsize',16)
% datetick('x','mm/yy','keeplimits');
% set(gcf,'color','w')
ylabel('Glacier Calving Front Retreat (m)')
ylabel('Glacier Calving Front Position (m)')
% hold on;plot([min(idates) max(idates)],[dists,dists],'g-','linewidth',1.5)
% hold on;plot([min(idates) max(idates)],[distn,distn],'m-','linewidth',1.5)
hold on;h=area([[min(idates)-1000 max(idates)+1000]],[distn distn],dists);
h.FaceColor=[230 230 230]./255;
plot(idates,xx(p,3)-xx(p(1),3),'-ob','linewidth',2)
legend('Glacier calving front retreat','Location','Northwest')
eptick=datenum({'2000/01/01','2005/01/01','2010/01/01','2015/01/01','2020/01/01'});
set(gca,'XTick',eptick)
set(gca,'XTickLabel',eptick,'FontSize',12,'XTickLabelMode','manual')
datetick('x','mm/yyyy','keepticks')
% axis([[min(idates) max(idates) -500 4000]])
axis([730180.488574868 738470.1005229958 -500 4000])
box on

% % Get the correlation with time ds in /Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/plotcosicorscarp_main.m
gtime=idates;gdis=xx(p,3)-xx(p(1),3);
[timer,I]=sort(time); dsr=ds(I);
hold on;plot(timer,dsr,'+')
gdisi=interp1(gtime,gdis,timer,'linear',NaN);
hold on;plot(timer,gdisi,'+')
R = corrcoef(dsr,gdisi); %0.9906


ofile1='glacierfront.shp';
shp1 = struct('Geometry', 'Polyline', 'X', xx(:,1), 'Y', xx(:,2));
shapewrite(shp1, ofile1);

%get the position of scarp outline along centerline
load outlinexy.mat
psouth=[437895,6778500];
pnorth=[438870,6780510];
[lats,lons]=minvtran(utmstruct,psouth(1),psouth(2)); %-148.1531   61.1347
[lats,lons]=minvtran(utmstruct,pnorth(1),pnorth(2)); %-148.1357   61.1529

addpath('/Users/chunlidai/Box/ISSMCoastline/ContourAlgorithm/poly_poly_dist/');
[d_min, x_d_min, y_d_min,is_vertex, idx_c]=p_poly_dist(psouth(1), psouth(2), xx(:,1), xx(:,2), false); 
dist=d_min; 
inds=idx_c; %not tested; %index of the closest point on xref;
[d_min, x_d_min, y_d_min,is_vertex, indn]=p_poly_dist(pnorth(1), pnorth(2), xx(:,1), xx(:,2), false); 


dists=xx(inds,3)-xx(p(1),3);distn=xx(indn,3)-xx(p(1),3);
fprintf(['\nCenterline distance of landslide south boundary: ',num2str(dists),' \n']) %1351.5 m
fprintf(['\nCenterline distance of landslide north boundary: ',num2str(distn),' \n']) %3499.2 m

hold on;plot([min(idates) max(idates)],[dists,dists],'g-','linewidth',1.5) 
hold on;plot([min(idates) max(idates)],[distn,distn],'m-','linewidth',1.5)


file='glacierfrontposition.txt';
fid10 = fopen(file,'w');
for i=1:length(idates)
fprintf(fid10,'\n %s %12.2f %12.2f %12.2f ',datestr(idates(i),26),xx(p(i),3)-xx(p(1),3),xx(p(i),1),xx(p(i),2));
end
fclose(fid10)


%plot the location of scarp outline, and glacier center line.
figure;
set(gcf,'Color','white')
set(gca,'FontSize', 18);
set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
hold all;
plot(xx(:,1)*1e-3,xx(:,2)*1e-3,'.-')
plot(x*1e-3,y*1e-3,'b.')
hold on;plot(psouth(1)*1e-3,psouth(2)*1e-3,'g>');
hold on;plot(pnorth(1)*1e-3,pnorth(2)*1e-3,'m>');
xlabel('x coordinates (km, UTM zone 6N)');ylabel('y coordinates (km)')
axis equal
box on
hold on;plot(xx(inds,1)*1e-3,xx(inds,2)*1e-3,'g>')
hold on;plot(xx(indn,1)*1e-3,xx(indn,2)*1e-3,'m>')


% The attached matfile contains the ice front position along the glacier centerline, also plotted.
% 
% idates are the image dates (datenum)
% xx is an array [dx x y] where dx is the distance along a flow line transect from an arbitrary point in front of the glacier and x and y are the transect coordinates.
% p is the index location of glacier front at the idate.
% so. To get the attached plot you do:
% 
% >> load('frontposition_barry.mat')
% >> plot(idates,xx(p)-xx(p(1)),'-ok','linewidth',2)
% >> set(gca,'fontsize',16)
% >> datetick('x','mm/yy','keeplimits');
% >> set(gcf,'color','w')
% >> ylabel('Glacier Calving Front Retreat (m)')
% 
% The rapid retreat is textbook unstable tidewater retreat, where the glacier front was up on a bathymtric high (actually and island) before 2008 and then rapidly retreated through deeper bathymetry, probably an overdeepening, unltil it hit shallower water and stabilitzed in 2015.
% 
% Ian