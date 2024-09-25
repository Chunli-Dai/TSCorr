function [Co]=getvel()
% given pairs of images; run sdm; write list of output files;
constant
Co=[];

%ores=resr_cosi ; %; output resolution in meters.

%sdmlist='sdmlist.txt'; %input file list for SDM.


%create corrfilelist % 'corrfilelist'
%corrfilelist='corrfilelist'; %a list of absolute file names of displacement map results from image correlation.
%flagtool=2; % tools to use: 1, COSI-COrr; 2, SETSM SDM; 3, Seongsu's MIMIC2.
if flagtool==2

    fid12 = fopen(sdmlist,'r');
    n = linecount(fid12);

    strg=cell(n,1);
    strog=cell(n,1);
    fid12 = fopen(sdmlist);
    for i=1:n
      str=[fgetl(fid12)];
      strg{i}=str;
    end
    fclose(fid12);
    
    %poolobj=parpool(poolsize);
    for i=1:n
        str=strg{i};
        
        %str=['time /fs/project/howat.4/SETSM/setsm -SDM 2 -image /home/dai.56/chunliwork/sdm/runsite22/landsatcosicorr/run4/20150617000000_LC08_L1TP_068017_20150617_20170226_01_T1_B8_sub.TIF -image /home/dai.56/chunliwork/sdm/runsite22/landsatcosicorr/run4/20170615000000_LC08_L1TP_067017_20170615_20170629_01_T1_B8_sub.TIF -outpath /home/dai.56/chunliwork/sdm/runsite22/landsatcosicorr/run4/runL8ores100 -sdm_as 0.08 -sdm_days 729 -outres 100'];
        
        fprintf([str,'\n'])
        [status, cmdout]=system(str) 
        
        %refer to ~/arcticdemapp/landslide/code1/strip2date.m 
        %newStr = extractBetween(str,"quick "," fox")
        outpath=extractBetween(str,"-outpath","-sdm_as");
        outpath=strtrim(outpath{1});
        [~,name,ext] =fileparts([strtrim(outpath)]);
        name=[name,ext];
        % output: /home/dai.56/chunliwork/sdm/runsite22/landsatcosicorr/run4/runL8ores100/runL8ores100_dx.tif
        strog{i}=[outpath,'/',name,'_dx.tif'];
    
    end %for i
    %delete(poolobj)
    
    fid13 = fopen(corrfilelist,'w');
    for i=1:n
    fprintf(fid13,'%s \n',strog{i});
    end
    fclose(fid13);

elseif flagtool==1 %Use CosiCorr
%corrcosiresultdir 
    str=sprintf('ls  %s/correlation_*.TIF_*[0-9]  > %s',deblank(corrcosiresultdir),corrfilelist);
    [status, cmdout]=system(str);

elseif flagtool==4 %Use CosiCorr Plus

    coregfile=[currentdir,'/coreg.mat'];
    if exist(coregfile,'file') % it tooks 14 hours to run 963 pairs with 10 CPUs.
        fprintf('\n getvel.m: coreg.mat and results exist; no need to rerun cosicorr plus.\n')
    else

    fid12 = fopen(cosipluslist,'r');
    n = linecount(fid12);

    strg=cell(n,1);
    fid12 = fopen(cosipluslist);
    for i=1:n
      str=[fgetl(fid12)];
      strg{i}=str;
    end
    fclose(fid12);

    numNodesStr = getenv('SLURM_NNODES');
    sn=str2num(numNodesStr)

    sz = getenv('SLURM_NTASKS');
    sz=str2num(sz);

%   sz=sn*sz;
    fprintf(['\n ',num2str(sz),' worker(s) allocated in job.slurm.\n'])

    [status, cmdout]=system('free -h')

    for i=1:n
        str=strg{i};
	%skip if the output file exists.

%	str1=[str];
	str1=['time ',str,' & '];
        fprintf([str1,'\n'])
        [status, cmdout]=system(str1);  %hang jobs with &  
                
        % python correlate_cli.py correlate para/test_planeti.json &
	if mod(i,sz)==0 || i==n % wait for the command to finish when i is mulitple times of ncpus.
    	   [status, cmdout]=system('free -h')
    	   [status, cmdout]=system('ps ux |grep correlate_cli.py |wc -l');
    	   ncmd=str2num(cmdout)-2;  % number of command lines that contain correlate_cli.py
	   fprintf(['\n Number of cosi-corr+ pairs running:',num2str(ncmd),'; i/n=',num2str([i,n]),' \n']);
	   while ncmd >0
        	pause(30); %30 second;
        	[status, cmdout]=system('ps ux |grep correlate_cli.py |wc -l');
        	ncmd=str2num(cmdout)-2; % number of command lines that contain correlate_cli.py
	   	fprintf(['\n Number of cosi-corr+ pairs running:',num2str(ncmd),'; i/n=',num2str([i,n]),' \n']);
   	   end
     	 %  pause(300); %30 second;
	end
    end %for i

    [status, cmdout]=system('ps ux ')
    [status, cmdout]=system('free -h')

    end % exist

    % corrfilelist
    str=sprintf('ls  %s/correlation_*_vs_*.TIF  > %s',deblank(corrresultdir),corrfilelist);
    [status, cmdout]=system(str);

end


return
end
