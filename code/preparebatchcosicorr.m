function [Co]=preparebatchcosicorr(rang0,XYbg,projgdalg,f,fdir)
%prepare the files for Batch correlation in Cosi-corr;
%get the boundary of Master image;
%get the targets images;


%addpath(genpath(['/Users/chunlidai/surge/home/dai.56/arcticdemapp/landslide/code1/']));
Co=[];
constant

flageq2=0; % default, not earthquake
if exist('flageq', 'var') == 1
    flageq2=flageq;
end

flagcrop=1; %1 apply subsetting (default); 0 do not apply the subsetting.
if exist('flagcropin', 'var') == 1
    flagcrop=flagcropin;
end

if 0
%Taan 
filename='filelistTaanASTER'; %need to be sorted by date
dir='/Users/chunli/Desktop/Osuworkmac/asar/tyndall/TaanASTER/'; %directory of files as input of COSI-Corr
filename='/Users/chunlidai/Box/NSFpermafrost2019/manuscript/runsdm_cosicorr/cosicorrinputTaanLandsat/filelistTaanLandsat'; %need to be sorted by date
dir='/Users/chunli/Desktop/Osuworkmac/asar/tyndall/landsat/'; %directory of files as input of COSI-Corr
end

dir=cosidir;

% old
if ~exist('fdir', 'var') 
filename=[currentdir,'/imagelist.dat'];
%imagesubdir;
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
fdir=cell(n,1);f=cell(n,1);
range=zeros(n,4);XYbg=cell(n,1);projgdalg=cell(n,1);
for i=1:n
   ifile=[fgetl(fid)];
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
      % get the boundary from xml file
   [XYbi,rangei,projgdali]=imagebd(ifile);
   range(i,1:4)=rangei;XYbg{i}=XYbi;projgdalg{i}=projgdali;
end
end %exist fdir


%Only work on subset of all images in the region of interest

% Name the subset (e.g. 20150617000000_LC08_L1TP_068017_20150617_20170226_01_T1_B8_sub.TIF) to be compatible with MIMC2.
%  Shorten the name: YYYYMMDDHHMMSS_LC08_B8_sub.TIF
% creating the result file name structure as the output of Cosi-Corr: e.g. correlation_LE07_L1TP_063018_19990728_20161003_01_T1_B8.TIF_vs_LE07_L1TP_062018_20020713_20160929_01_T1_B8.TIF_22

bandi=1;

% if flagtype==2 %ASTER
% bandi=2; %RED band
% else
% bandi=1;
% end

n=length(fdir);
for i=1:n
ifile=[fdir{i},'/',f{i}];
filei=f{i};
% flagtype=1; %image data types: 1, Landast 8; 2, Aster; 3, Planet; 4, PGC images; 5, Landsat 7 striped images;
if flagtype==1
%LC08_L1TP_048025_20201024_20201106_02_T1_B8.TIF
ymdstr=[filei(18:25),'000000']; %eg 20150617000000
mission=['_',filei(1:4)]; %e.g. _LC08, _AST0, _PLAN, _WV02,
bandstr='_B8_sub.TIF'; % _B8_sub.TIF; Band 8 for landsat, band 2 for aster, band 1 for planet and worldview.
bandi=1; %8;
elseif flagtype==2
%AST_L1T_00305072020193649_20200508100034_14305_V.tif
ymdstr=[filei(16:19),filei(12:15),filei(20:25)];
mission='_AST0';
bandstr='_B2_sub.TIF';
bandi=2; %RED band
elseif flagtype==3
%20180924_merge.tif
%ymdstr=[filei(1:8),'000000'];
%mission='_PLAN';
%bandstr='_B1_sub.TIF';
%20230120_075156_40_247b_3B_AnalyticMS_SR_clip.tif
ymdstr=[filei(1:8),filei(10:15)];
mission='_PLAN';
bandstr='_AB_sub.TIF'; %all bands, Visual products which are 3 band “RGB” NIR.

elseif flagtype==4 
%WV02_20111008_103001000EA22900_103001000D580C00_2m_lsf_seg1_meta.txt
%SETSM_s2s041_W1W1_20081112_1020010005BC6D00_1020010005C6A700_2m_lsf_seg3_ortho.tif
filei=strrep(filei,'SETSM_s2s041_','');
ymdstr=[filei(6:13),'000000'];
mission=['_',filei(1:4)];
bandstr='_B1_sub.TIF';
elseif flagtype==5
%LE07_L1TP_063018_19990728_20161003_01_T1_B8.TIF
ymdstr=[filei(18:25),'000000']; 
mission=['_',filei(1:4)]; 
bandstr='_B8_sub.TIF'; 
bandi=1; %8;
elseif flagtype==6
%20130505000000_LC08_B8_sub.TIF
ymdstr=filei(1:14);
mission=filei(15:19);
%bandstr=filei(20:30);
bandstr=filei(20:end);
bandi=1;
elseif flagtype==8
	%WV01_20121101085030_102001001D815E00_12NOV01085030-P1BS-505057038040_01_P002_u16ns32637.tif
	ymdstr=[filei(6:19)];
	mission=['_',filei(1:4)];
	bandstr='_B1_sub.TIF';
end
fsub{i}=[ymdstr,mission,bandstr];
fdirsub{i}=imagesubdir;
ofile=[fdirsub{i},'/',fsub{i}];

% gdal_translate -projwin 431000 6784000 445000 6772000  LC08_L1TP_068017_20150617_20170226_01_T1_B8.TIF 20150617000000_LC08_L1TP_068017_20150617_20170226_01_T1_B8_sub.TIF
% gdalwarp $infile $ofile3 -t_srs epsg:3413 -tr $ores $ores -te -3429122 2374856 -3406392 2396770
testr=num2str([rang0(1) rang0(3) rang0(2)  rang0(4)]); %;[-te xmin ymin xmax ymax]

%Get image subsets within the given boundary.
if flagcrop==1
% Input image should be unsigned integer 16 bits for SDM: UInt16
% gdalwarp $infile $ofile3 -t_srs  $proj -te $region -r bilinear -tr 15 15 -ot UInt16  #interpolate 
str=['gdalwarp ',ifile,' ',ofile,' -t_srs ',projgdal,' -te ', testr,' -r bilinear  -tr ',num2str([imageores imageores]),' -ot UInt16'] 
[status, cmdout]=system(str) 

if flagtype==2 % % ||flagtype==5 || flagtype==1 %ASTER %get only band 2
fprintf(['\n preparebatchcosicorr.m select band ', num2str(bandi),' from the input image:', ifile])
str=['mv ',ofile,' t1.tif'];
[status, cmdout]=system(str);
str=['gdal_translate',' -b ', num2str(bandi),' t1.tif ',' ',ofile,' -t_srs ',projgdal,' -te ', testr,' -r bilinear  -tr ',num2str([imageores imageores]),' -ot UInt16'];
[status, cmdout]=system(str);
end %ASTER

else % do not crop, but make sure images are in same projection and same output pixel size.

	projgdali=projgdalg{i};
	if ~exist('flagwarp','var'); flagwarp=0;end %default 0

	%using gdal to convert projection
	if flagwarp==1 || ~strcmp(projgdali,projgdal)
           fprintf(['\n Image has different projection than that in the constant.m: ',ifile,', ' , projgdali,' \n']);
           fprintf(['\n Converting the projection for the boundary, but no cropping. \n']);
	   % the output grid should match the base images. The coordinates should be divisible by imageores.
   	   XYbj=XYbg{i};
           xj=XYbj(:,1);yj=XYbj(:,2);
        %  xy to lat lon
	   [latj,lonj]=xy2latlon(xj,yj,projgdali);
        %  lat lon to xy
	   [xj,yj]=latlon2xy(latj,lonj,projgdal);
	   rangib=[min(xj) max(xj) min(yj) max(yj)];

	   ranget=round(rangib/imageores)*imageores;rangib=ranget;
	   testrb=num2str([rangib(1) rangib(3) rangib(2)  rangib(4)]); %;[-te xmin ymin xmax ymax]

	   str=['gdalwarp ',ifile,' ',ofile,' -t_srs ',projgdal,' -te ', testrb,' -r bilinear  -tr ',num2str([imageores imageores]),' -ot UInt16'] 
	   [status, cmdout]=system(str) 

	else  %link the images to imagesubdir
           fprintf(['\n Link image to imagesubdir: ',ifile,', ' , projgdali,' \n']);
	   str=['ln -fs ',ifile,'  ',ofile];
	   [status, cmdout]=system(str);
        end

end %flagcrop

end %for i


%input for COSI-Corr
fid10 = fopen('bases.txt','w');
fid11 = fopen('targets.txt','w');   
%sdmlist='sdmlist.txt';
fid12 = fopen(sdmlist,'w');   
fid14 = fopen(cosipluslist,'w');   

ores=resr_cosi ; %; output resolution in meters.

% input option for COSI-Corr
option='2';%'6a';  %'6b';%'1';
        %1, 1-element: filename (full path);
        %2, 2-elements: filename (full path) + band number (first band = 1).
        %6, 6-elements: filename (full path) + band number + spatial extent (Xstart, Xend, Ystart, Yend).

ndem=n;
bandc=1; %bandi; % input band number for COSI Corr; 

for idem1=1:ndem %index of DEM 1

    %read DEM 1
    i=idem1;
    metafile=[fdirsub{i},'/',fsub{i}]; %metafileref=metafile;

%   if flagcrop==1 %use the cropped images image1 image2, or linked images.
        basefile=fsub{i};
        image1=[fdirsub{idem1},'/',basefile];  %full path name
%   else  %use the file in the original folder instead of the new folder imagesubdir.
%       basefile=f{i};
%       image1=[fdir{idem1},'/',basefile];
%   end
    baseofile=fsub{i}; %to construct short output name
    
    if ~strcmp(option,'1')&&~strcmp(option,'2') %need the spatial extent.
    data=readGeotiff(metafile);

    idx=find(data.x>=rang0(1)&data.x<=rang0(2));
    idy=find(data.y>=rang0(3)&data.y<=rang0(4));
    % parameters for base file list in CosiCorr-Guide.pdf:spatial extent (Xstart, Xend, Ystart, Yend).
    para=[min(idx) max(idx) min(idy) max(idy)];
    % fprintf(num2str(para))
    end
    
    for idem2=(idem1+1):ndem
        j=idem2;

%       if flagcrop==1 %use the cropped images image1 image2, or linked images. 
            targetfile=[fsub{j}];
            image2=[fdirsub{idem2},'/',targetfile];
%       else
%           targetfile=[f{j}];
%           image2=[fdir{idem2},'/',targetfile];
%       end
        targetofile=fsub{j}; %to construct short output name

        % For earthqaukes, big study area. Check overlapping rate.
        % Skip pairs with no overlap.
        if flageq2==1
            [~,~,overlap_ratiom]=getoverlap(XYbg{i},XYbg{j});
	    thres=20; %default
	    if exist('coverthres','var')
 		   thres=coverthres;
	    end
            if overlap_ratiom<=thres %20 % skip this pair
               continue 
            end
        end % flageq2
        % End earthqaukes
        
	%e.g. /Users/chunli/Desktop/Osuworkmac/asar/tyndall/landsat/LE07_L1TP_063018_19990728_20161003_01_T1_B8.TIF 1  10315  10781   7028   7494 
        if strcmp(option,'1')
            str1=[dir,basefile];
        elseif strcmp(option,'2')
            %str1=[dir,basefile,'  1'];
            str1=[dir,basefile,'  ',num2str(bandc)];
        elseif strcmp(option,'6a')
            %str1=[dir,basefile,' ',num2str([ 1 para])]; %pan band for   Landsat
            str1=[dir,basefile,' ',num2str([ bandc para])]; 
        elseif strcmp(option,'6b')
%            str1=[dir,basefile,' ',num2str([ 2 para])]; %green for aster -> Band 2 is Red; Band 2 is named green in COSICOrr, but its wavelength is Red.
            str1=[dir,basefile,' ',num2str([ bandc para])]; 
        end
	%e.g. /Users/chunli/Desktop/Osuworkmac/asar/tyndall/landsat/LE07_L1TP_062018_20020408_20160928_01_T1_B8.TIF 
        str2=[dir,targetfile]; % band 2 green band for salves images too.

        
        %% construct files for SETSM

%str=['time /fs/project/howat.4/SETSM/setsm -SDM 2 -image /home/dai.56/chunliwork/sdm/runsite22/landsatcosicorr/run4/20150617000000_LC08_L1TP_068017_20150617_20170226_01_T1_B8_sub.TIF -image /home/dai.56/chunliwork/sdm/runsite22/landsatcosicorr/run4/20170615000000_LC08_L1TP_067017_20170615_20170629_01_T1_B8_sub.TIF -outpath /home/dai.56/chunliwork/sdm/runsite22/landsatcosicorr/run4/runL8ores100 -sdm_as 0.08 -sdm_days 729 -outres 100'];

% 	date1=basefile(1:14);date2=targetfile(1:14);
    date1=filename2date(baseofile);date2=filename2date(targetofile);
% % creating the result file name structure as the output of Cosi-Corr: e.g. correlation_LE07_L1TP_063018_19990728_20161003_01_T1_B8.TIF_vs_LE07_L1TP_062018_20020713_20160929_01_T1_B8.TIF_22
	outpath=[corrresultdir,'/correlation_',baseofile,'_vs_',targetofile];

    if length(strtrim(date1)) == 8
       datefmt='yyyymmdd';
    elseif length(strtrim(date1)) ==14
       datefmt='yyyymmddHHMMSS';
    else
       fprintf(['\n preparebatchcosicorr.m error: date format is not yyyymmddHHMMSS or yyyymmdd: ',str1,' \n'])
    end
    dt=datenum(date2,datefmt)-datenum(date1,datefmt); %date difference in days
    %dt=datenum(date2,'yyyymmddHHMMSS')-datenum(date1,'yyyymmddHHMMSS'); %date difference in days

	str3=['time setsm -SDM 2 -image ',image1,' -image ',image2,' -outpath ',outpath,' -sdm_as 0.08 -sdm_days ',num2str(dt),' -outres ',num2str(imageores)];
	%SETSM SDM works best when using the same resolution as input image.

    %% construct files for COSI-COrr Plus
        
    % python correlate_cli.py correlate para/test_planeti.json
    jsonfile=['para/cosiplus',num2str(i),'vs',num2str(j),'.json'];
    str4a=['cp cosiplus.json ',jsonfile,'; sed -i ''s|base.tif|',image1,'|g'' ',jsonfile,'; sed -i ''s|target.tif|',image2,'|g'' ',jsonfile,'; sed -i ''s|correlation_base_vs_target|',outpath,'|g'' ',jsonfile];
    [status, cmdout]=system(str4a);

    str4b=['python correlate_cli.py correlate ',jsonfile];

        fprintf(fid10,'%s \n',str1);
        fprintf(fid11,'%s \n',str2);
        fprintf(fid12,'%s \n',str3);
        fprintf(fid14,'%s \n',str4b);

    end
end
fclose(fid10);
fclose(fid11);
fclose(fid12);
fclose(fid14);

%input for SETSM SDM; ignore this
if flagtool==2
  ores=resr_cosi; %; output resolution in meters.
end


return
end


