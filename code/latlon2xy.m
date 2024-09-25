function [xj,yj]=latlon2xy(latj,lonj,projgdalj)

	   if strcmp(projgdalj,'epsg:3413'); % Arctic
	      [xj,yj]=polarstereo_fwd(latj,lonj,[],[],70,-45);
	   elseif strcmp(projgdalj,'epsg:3031') %Antarctica
	      [xj,yj]=polarstereo_fwd(latj,lonj,[],[],-71,0);
	   else %utm
	   
	   % e.g., UTM Zone 2S: EPSG 32702; UTM Zone 3N: EPSG 32603;
	       %Using gdal, e.g., echo $lon $lat | gdaltransform -s_srs EPSG:4326 -t_srs EPSG:32719
	       crs_epsg=str2double(projgdalj(6:end));%e.g.32636;
	       p1 = projcrs(crs_epsg); %projcrs, Since R2020b
	       [xj,yj] = projfwd(p1, latj,lonj);

   	   end

return
end
