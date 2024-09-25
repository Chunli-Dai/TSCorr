function [XYbi,rangei,projgdali]=imagebd(ifile)
%get the Polygon boundary, and the rectangle range of a mono image from .xml file
% e.g. /data1/pgc_projects/dai_aleutians_multi_mono/imagery//WV02/WV02_20130923221601_1030010026BD6B00_13SEP23221601-M1BS-500127380110_01_P001.xml 
% Or from tif file using gdalinfo
	constant
	projgdali=''; % initialize; %need to update the retrieval of projgdal from xml or meta.txt mdf.txt.

    vstr={'<ULLON>','<ULLAT>','<URLON>','<URLAT>','<LRLON>','<LRLAT>','<LLLON>','<LLLAT>'};
%   flagfmt1=1; %flag for format 1; default format
    
    vstr2={'<upperLeft>','<upperRight>','<lowerRight>','<lowerLeft>'};%mono image xml
%   flagfmt2=1; %

	vstr3={'X:','Y:'}; %3 strip meta file contain vstr4 also.
	vstr4={'Upper left coordinates'};%4 scene meta file, tif_results
	
	if ~exist(ifile,'file')
		warning(['File does not exist:',ifile])
        	XYbi=[0 0]; 	rangei=[0 0 0 0];
		return;
	end

	[demdir,name,ext] =fileparts([strtrim(ifile)]);
        name1=[name(end-3:end),ext];

       if contains(name1,'.tif')||contains(name1,'.TIF')
          %skip 
       else

	metafile=ifile;
	c=textread(metafile,'%s','delimiter','\n');
	
    %check the format 
	[demdir,name,ext] =fileparts([strtrim(ifile)]);
	name1=[name(end-3:end),ext];
    i=1;str=vstr{i};
    r1=find(~cellfun(@isempty,strfind(c,str)));
    i=1;str=vstr2{i};
    r2=find(~cellfun(@isempty,strfind(c,str)));

    	end %if

    if strcmp(name1,'meta.txt') %3 (strip) or 4 (tif_results)
        i=1;str=vstr3{i};
        r3=find(~cellfun(@isempty,strfind(c,str)));
        i=1;str=vstr4{i};
        r4=find(~cellfun(@isempty,strfind(c,str)));
        if ~isempty(r3)
           flagfmt=3;
        else
           flagfmt=4;
        end
    elseif strcmp(name1,'_mdf.txt')
        flagfmt=5;
    elseif contains(name1,'.tif')||contains(name1,'.TIF')
        flagfmt=6;
    elseif ~isempty(r1)
        flagfmt=1;
    elseif ~isempty(r2)
        flagfmt=2;
    end
    

% Get the Footprint Vertices X, Y, close the loop
    if flagfmt==1 
% /fs/byo/howat-data5/pgc_deliv/chunli/DaiTerHyd/deliv3_AlaskaRiverGaugesOrthoImagery_20190315/ortho_imagery_non-max_ona/QB02_20020222212303_1010010000215600_02FEB22212303-M1BS-052800669020_01_P001_u16ns3413.xml     
%                         <ULLON>-1.501389837000000e+02</ULLON>
%                         <ULLAT>6.126182419000000e+01</ULLAT>
    
	n=length(vstr)/2;
	lon=zeros(n+1,1);lat=zeros(n+1,1);
	for i=1:n*2
	str=vstr(i);
	r=find(~cellfun(@isempty,strfind(c,str)));
	if isempty(r)        
		warning(['xml file is different as anticipated.',ifile])
        XYbi=[0 0]; 	rangei=[0 0 0 0];
		return;
	end %
	c2=c{r(1)};
	r1=strfind(c2,'>');r2=strfind(c2,'</');
    c2([1:r1(1),r2(1):end])='';
	z = sscanf(c2, '%g', 1);
	j=ceil(i/2);
	if mod(i,2)  %1 odd number, 0 even
           lon(j)=z;
	else
	   lat(j)=z;
	end
	end % if i

	%add projgdali
    str='<EPSG_CODE>';
    r=find(~cellfun(@isempty,strfind(c,str)));
    c2=c{r(1)};
        pattern = '<EPSG_CODE>(\d+)</EPSG_CODE>';
        matches = regexp(c2, pattern, 'tokens');
        % Check if there is a match
        if ~isempty(matches)
            % Extracted number as a string
            projgdali= matches{1}{1};
            projgdali=['epsg:',projgdali];

    else
            fprintf(['\n EPSG code not found in:',ifile,'\n']);
    end


	lon(n+1)=lon(1);lat(n+1)=lat(1); %close the loop
%[Xb,Yb]=polarstereo_fwd(lat,lon,[], [],70,-45);
	[Xb,Yb]=latlon2xy(lat,lon,projgdali);

    elseif flagfmt ==2
        %/data4/EarthDEM/alaska_2018oct26/qb_wv_alaska_metadata/QB02_20030408203848_1010010001C95401_03APR08203848-P1BS-000000073404_01_P001.xml 
%           <upperLeft>
%                 <latitude>61.21</latitude>
%                 <longitude>-145.416111111</longitude>
%           </upperLeft>
        
    n=length(vstr2);
	lon=zeros(n+1,1);lat=zeros(n+1,1);
	for i=1:n
	str=vstr2(i);
	r=find(~cellfun(@isempty,strfind(c,str)));
	if isempty(r)        
		warning(['xml file is different as anticipated.',ifile])
        XYbi=[0 0]; 	rangei=[0 0 0 0];
		return;
	end %
	c2=c{r(1)+1};
	r1=strfind(c2,'>');r2=strfind(c2,'</');
    c2([1:r1(1),r2(1):end])='';
	zlat = sscanf(c2, '%g', 1);
    
    c2=c{r(1)+2};
	r1=strfind(c2,'>');r2=strfind(c2,'</');
    c2([1:r1(1),r2(1):end])='';
	zlon = sscanf(c2, '%g', 1);

       lon(i)=zlon;
	   lat(i)=zlat;
	end % for i        
	lon(n+1)=lon(1);lat(n+1)=lat(1); %close the loop
	[Xb,Yb]=polarstereo_fwd(lat,lon,[], [],70,-45);
	%[Xb,Yb]=latlon2xy(lat,lon,projgdali);

    elseif flagfmt ==3
	%/data3/ArcticDEM/region_34_alaska_north/strips/2m/WV02_20160806_103001005A286A00_103001005A1A8300_seg1_2m_meta.txt
    %X: -2060408 -2058776 -2058248 -2057592 -2057376 -2054626 -2054282 -2054186 -2052908 -1962418 -1961946 -1957764 -1957764 -1959332 -1960252 -1962620 -1962630 -1963390 -2060336 -2060464 -2060480 -2060408 
    % Y: 740374 745446 747054 748996 749476 755420 756140 756332 756412 723102 722926 720648 720622 715246 712112 704224 704224 704488 739328 739376 739384 740374 
    
        r=find(~cellfun(@isempty,strfind(c,vstr3(1))));
        Xbs=deblank(strrep(c{r(1)},vstr3{1},''));
	if ~isempty(Xbs)
        Xb=strread(Xbs,'%d');
        r=find(~cellfun(@isempty,strfind(c,vstr3(2))));
        Ybs=deblank(strrep(c{r(1)},vstr3{2},''));
        Yb=strread(Ybs,'%d');   
	else
          warning(['flagfmt=3; xml file is different as anticipated.',ifile])
          XYbi=[0 0;]; 	rangei=[0 0 0 0];return
	end
    elseif flagfmt ==4
	%/home/dai.56/data2/ArcticDEM/region_31_alaska_south/tif_results/2m/GE01_20161023_1050010006C6CC00_1050010006C6CA00_501020521050_01_P003_501020521060_01_P002_2_meta.txt
	% chunli/scripts/run_overlap.sh  gives the boundary of image including edges; to use matchtag to exclude edges.
%     Output dimensions=3621  8841
%     Upper left coordinates=-3065560.000000  827080.000000

    if 1 %method with matchtag, too slow, 5 sec; 1 sec if (imresize with 0.01 instead of 0.1)
	tic  %too slow, 5 sec
        infile= strrep(metafile,'meta.txt','matchtag.tif');
        fprintf(['flagfmt=4; use matchtag for image boundary:',ifile])
        if exist(infile,'file') 
        data=readGeotiff(infile); %0 is good; 1 is edge; > 1 bad data.
        else
        fprintf([infile,' not exist.'])
           XYbi=[0 0]; 	rangei=[0 0 0 0];
		return;
        end
        datar.x= imresize(data.x,0.01);datar.y= imresize(data.y,0.01);
        datar.z= imresize(data.z,0.01);
        [X,Y]=meshgrid(datar.x,datar.y);
        idx=X(datar.z==1);idy=Y(datar.z==1);
        cf = 0.1;%0.5; %boundary curvature factor (0= point boundary, 1 =conv hull)
        k = boundary(idx,idy,cf);
        Xb=idx(k);Yb=idy(k);
	toc
	elseif 0 %method using xml file; cons: xml file may not be provided. 2 sec or 1 sec
	tic
	  [~,filename,~]=fileparts(metafile);
	  satname=filename(1:4);
	  imagestr={'Image 1=','Image 2='};
	  for i=1:2
            %c1=c;r1=strfind(c1,'Image 1=');if(isempty(r1)); warning('Image 1 not found') ; end
            r1=find(~cellfun(@isempty,strfind(c,imagestr{i})));
	    if(isempty(r1)); warning([imagestr{i},' not found']) ;XYbi=[0 0];  rangei=[0 0 0 0];return; end
            clear c1;c1=c{r1(1)};
            r1=strfind(c1,'/');c1(1:(r1(end)))='';
            if strcmp(satname,'WV01')
                mfile{i}=deblank(c1);
            else
                mfile{i}=deblank(strrep(c1,'P1BS','M1BS'));%e.g.WV02_20160304214247_1030010052B75A00_16MAR04214247-M1BS-500641617080_01_P009.tif
            end
	    xmlfile=strrep(mfile{i},'.tif','.xml');
	    [status, cmdout ]=system(['find -L ',multidir,' -name ',xmlfile]); %returns empty for cmdout when not found. status=0 always
	    if  ~isempty(cmdout) %if .ntf file is found, tif file is produced and stored in orthworkdir.
    		mfile{i}=deblank(cmdout);
	    end
	  end % i
	  [XYbi1,rangei1]=imagebd(mfile{1});
	  [XYbi2,rangei2]=imagebd(mfile{2});
	  %get the common area boundary
	  poly1 = polyshape(XYbi1(:,1),XYbi1(:,2));
	  poly2 = polyshape(XYbi2(:,1),XYbi2(:,2));
	  polyout = intersect(poly1,poly2);
	  if ~isempty(polyout.Vertices)
          %Xb=polyout.Vertices(:,1);Yb=polyout.Vertices(:,2);
	  Xb=[polyout.Vertices(:,1);polyout.Vertices(1,1)];Yb=[polyout.Vertices(:,2);polyout.Vertices(1,2);];
	  else
	    XYbi=[0 0;];    rangei=[0 0 0 0]; 
	    return
	  end
          toc
	else % use meta.txt, which covers the entire map, not excluding edge. 0.2 sec, fast but not accurate.
	  % chunli/scripts/run_overlap.sh  gives the boundary of image including edges; to use matchtag to exclude edges.
            r1=find(~cellfun(@isempty,strfind(c,'Output Resolution=')));
	    c1=c{r1(1)};
	    r1=strfind(c1,'=');c1(1:r1)='';
            res = sscanf(c1, '%g', 1);
	
            r1=find(~cellfun(@isempty,strfind(c,'Output dimensions=')));
	    c1=c{r1(1)};
	    r1=strfind(c1,'=');c1(1:r1)='';
            [nxy]= sscanf(c1, '%g', 2);
	    nx=nxy(1);ny=nxy(2);

            r1=find(~cellfun(@isempty,strfind(c,'Upper left coordinates=')));
	    c1=c{r1(1)};
	    r1=strfind(c1,'=');c1(1:r1)='';
            xy= sscanf(c1, '%g', 2);
	    xl=xy(1);yu=xy(2);
	    xr= xl + ( nx - 1 ) * res ;
	    yl=yu - ( ny - 1 ) * res ;
	    Xb=[xl xr xr xl xl];Xb=Xb(:);
	    Yb=[yl yl yu yu yl];Yb=Yb(:);
	    if isempty(Xb)||isempty(Yb)
	        XYbi=[0 0;];    rangei=[0 0 0 0];
                return;
	    end
	end %method
    
    elseif flagfmt == 5
        %  '/Users/chunlidai/share/landslide/site1Eureka/stripdata/SETSM_GE01_20090804_10504100049D6800_1050410002F81C00_seg1_2m_v3.0_mdf.txt'
%         X1 = -699068.0;
%         Y1 = -842660.0;
%         X2 = -699052.0;
%         Y2 = -842636.0;
	str1='X1 =';
        r1=find(~cellfun(@isempty,strfind(c,str1)));
	str='horizontalCoordSysOGCWKT';
        r2=find(~cellfun(@isempty,strfind(c,str)));

	n=round((r2-r1)/2);
	lon=zeros(n+1,1);lat=zeros(n+1,1);
	for i=1:n*2
        j=ceil(i/2);
        if mod(i,2)  %1 odd number, 0 even
           stri=deblank(strrep(str1,'X1',['X',num2str(j)]));
        else
           stri=deblank(strrep(str1,'X1',['Y',num2str(j)]));
        end

	str=stri;
	r=find(~cellfun(@isempty,strfind(c,str)));
	if isempty(r)        
		warning(['xml file is different as anticipated.',ifile])
        XYbi=[0 0]; 	rangei=[0 0 0 0];
		return;
	end %
	c2=c{r(1)};

	Ybs=deblank(strrep(c{r(1)},stri,''));
    Ybs=deblank(strrep(Ybs,';',''));
    z = sscanf(Ybs, '%g', 1);

%	r1=strfind(c2,'=');r2=strfind(c2,';');
%    c2([1:r1(1),r2(1):end])='';
%	z = sscanf(c2, '%g', 1);
	if mod(i,2)  %1 odd number, 0 even
           lon(j)=z;
	else
	   lat(j)=z;
	end
	end % if i
	lon(n+1)=lon(1);lat(n+1)=lat(1); %close the loop
	Xb=lon;Yb=lat; 
		

    elseif flagfmt == 6 %Geotiff files, using gdalinfo # to be written 
	str1=['gdalinfo ',ifile];
	[status, cmdout]=system([str1]);
	
	% get Xb Yb
	% Specify the target string to match
        targetString ='Upper Left  (';
	values=getvaluematch(targetString,cmdout);
	ulx=values(1);uly=values(2);

	targetString='Lower Left  (';
	values=getvaluematch(targetString,cmdout);
	llx=values(1);lly=values(2);

	targetString='Upper Right (';
	values=getvaluematch(targetString,cmdout);
	urx=values(1);ury=values(2);

	targetString='Lower Right (';	
        values=getvaluematch(targetString,cmdout);
        lrx=values(1);lry=values(2);
	Xb=[llx lrx urx ulx llx]'; Yb=[lly lry ury uly lly]';

	%Double check if the projection is consistent with the given one in constant.m
	%find the projection number
        pattern='ID\["EPSG",(\d+)';
        matches=regexp(cmdout,pattern,'match');
        if isempty(matches);
		fprintf(['\n Error: gdalinfo does not give the projection info: ',ifile,'\n']);
	end

        value=regexp(matches{end},'\d+','match');  %e.g.,value{1}='32637';
	projgdali=['epsg:',value{1}];
	% force projection to a fixed value
% 	projgdali=projgdal;

%       if ~strcmp(projgdali,projgdal)
	%  fprintf(['\n Image has different projection than that in the constant.m: ',ifile,', ' , projgdali,' \n']);
	%  fprintf(['\n Converting the projection for the boundary. \n']);
%	end

    % Unkown format
    else 
        		warning(['xml file is different as anticipated.',ifile])
        XYbi=[0 0;]; 	rangei=[0 0 0 0];
		return;
    end
    
        XYbi=[Xb(:),Yb(:)]; %n by 2
	rangei=[min(Xb) max(Xb) min(Yb) max(Yb)];
return
end

