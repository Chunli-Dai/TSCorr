function datestr=filename2date(filename)
% Input: the filename in string 
% Output: date in string format yyyymmddHHMMSS

%refers to plotcosicorr.m; May 2023
% [C,matches] = strsplit(name,{'correlation_','_vs_'});
% filename=C{2};filename2=C{3};
% datem=filename2date(filename);

ifile=filename;

[~,name,ext] =fileparts([strtrim(ifile)]);
name=[name,ext];
if contains(name,'AST')
    % correlation_AST_L1T_00306242000213418_20150409231130_50053_V.tif_vs_AST_L1T_00309272011211825_20150607194411_118314_V.tif_1
    % AST_L1T_00304292006212432_20150514043317_89008_V.tif
    id=findstr(name,'AST');
    datem=[name(id(1)+(15:18)),name(id(1)+(11:14))];
%         dates=[name(id(2)+(15:18)),name(id(2)+(11:14))];
elseif contains(name,'sub_utm2m') %PGC monoimages
    %correlation_IK01_20020402213700_2002040221371390000010000699_po_945612_pan_0000000_u16ns3413_sub_utm2m.tif_vs_WV01_20100427213038_102001000DA0E100_10APR27213038
    %IK01_20020402213700_2002040221371390000010000699_po_945612_pan_0000000_u16ns3413_sub_utm2m.tif
%     id=findstr(name,'correlation_');id2=findstr(name,'_vs_');
%     infile1=[name(id(1)+12:(id2(1)-1))];
%     infile2=[name((id2(1)+4):end)];
      infile1=name;
      datem=infile1(6:13);
%         dates=infile2(6:13);
elseif contains(name,'LE07_L1') ||contains(name,'LC08_L1') || contains(name,'LC09_L1') || contains(name,'L1TP') || contains(name,'L1GT') %Landsat
    %correlation_LC08_L1TP_068017_20150617_20170226_01_T1_B8.TIF_vs_LC08_L1TP_067017_20160831_20170222_01_T1_B8.TIF_1
    %correlation_LE07_L1TP_063018_19990728_20161003_01_T1_B8.TIF_vs_LE07_L1TP_062018_20020713_20160929_01_T1_B8.TIF_22
    %correlation_20131028000000_LC08_B8_sub.TIF_vs_20201024000000_LC08_B8_sub.TIF_dx.tif
    %LC09_L1TP_060019_20220604_20220604_02_T1_B8.TIF
    datem=name(18:25); %name(30:37);
%         dates=name(81:88);
elseif contains(name,'sub.TIF') || contains(name,'sub.tif') % subset ; %
    %correlation_20131028000000_LC08_B8_sub.TIF_vs_20201024000000_LC08_B8_sub.TIF_dx.tif
    %20131028000000_LC08_B8_sub.TIF
    datem=name(1:14); 

elseif contains(name,'Visual')||contains(name,'AnalyticMS') %planet monoimages
%     correlation_20160929_202103_0e0e_3B_Visual.tif_vs_20170722_202631_1044_3B_Visual.tif_1
%   20160929_202103_0e0e_3B_Visual.tif
%         id=findstr(name,'correlation_');id2=findstr(name,'_vs_');
%         infile1=[name(id(1)+12:(id2(1)-1))];
%         infile2=[name((id2(1)+4):end)];
%         datem=infile1((6:13)-5);
%         dates=infile2((6:13)-5);
    datem=name(1:8);
elseif contains(name,'vmap_') %MIMC outputs; 
        %'pp_result_20130315000000_20170310000000_raw.mat'
        %   'vmap_20130315000000_20170310000000.mat'
        datem=name(6:13);
%         dates=name(21:28);
elseif contains(name,'WV')||contains(name,'GE')
        datem=name(6:19);
       % WV01_20130310081719_102001001F1B5D00_13MAR10081719-P1BS-500058845020_01_P004.tif
elseif contains(name,'W1')||contains(name,'W2')||contains(name,'W3')
%SETSM_s2s041_W1W1_20081112_1020010005BC6D00_1020010005C6A700_2m_lsf_seg3_ortho.tif
	name=strrep(name,'SETSM_s2s041_','');
	datem=name(6:13);
else
    fprintf(['\n filename2date.m error: filename is not compatible with our code, ',name]);
    datestr='00000000000000';
    return;
end

%output format yyyymmddHHMMSS
str1=datem;
if length(strtrim(str1)) ==8
%    datefmt='yyyymmdd';
   datestr=[datem,'000000'];
elseif length(strtrim(str1)) ==14
%    datefmt='yyyymmddHHMMSS';
   datestr=datem;
else
   fprintf(['\n filename2date.m error: date format is not yyyymmddHHMMSS or yyyymmdd: ',str1,' \n'])
end


return
end
