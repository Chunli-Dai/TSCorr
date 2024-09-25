#/bin/sh -f

if [ $# -lt 1 ]
then
  echo "usage: $0 filelist "
  exit
fi

filelist=$1
ores=$2  #8m

odir='./subset/';
#odir='./orthosubset2mfixgrid/';
#for line in `cat filelistpgc`; do
#for line in `cat nadirviewlist.txt`; do
mkdir $odir
for line in ` cat $filelist `; do
infile=$line
dir1=$(dirname $line)
#odir=$dir1

filename=$(basename $line)
#"${firstString/Suzi/$secondString}
ofile1=$odir/${filename/.tif/_sub.tif}
ofile2=$odir/${filename/.tif/_sub_utm2m.tif}

ofile3=$odir/${filename/.tif/_utm.tif}
#LC08_L1TP_048025_20201024_20201106_02_T1_B8.TIF 20131028000000_LC08_B8_sub.TIF
satstr=${filename:0:4}
#ofile3=$odir/${filename:17:8}"000000_LC08_B8_sub.TIF"
ofile3=$odir/${filename:17:8}"000000_"$satstr"_B8_sub.TIF"

if [ -f "$ofile3" ]; then
ofile3=${ofile3/sub.TIF/sub2.TIF}
fi

echo $infile $ofile3

#gdal_translate -projwin -3113600 730300 -3102200 722000 /home/dai.56/chunliwork/sdm/datasite22/qb4202_Barry_Glacier_Pan_Mono_Ortho_Imagery_2020may08/ortho/WV03_20170203215639_10400100271A1300_17FEB03215639-P1BS-501172500040_01_P008_u16ns3413.tif WV03_20170203215639_10400100271A1300_17FEB03215639-P1BS-501172500040_01_P008_u16ns3413_sub.tif
#gdalwarp WV02_20150516_1030010041649400_103001004243B900_2m_lsf_seg4_ortho_prep.tif WV02_20150516utm_2m.tif -t_srs epsg:32606 -tr 2 2
#-te xmin ymin xmax ymax -te 432936 6772572 443098 6785034

#gdal_translate -projwin -3113600 730300 -3102200 722000 $infile $ofile1
#gdalwarp $ofile1 $ofile2 -t_srs epsg:32606 -tr 2 2 -te 432936 6772572 443098 6785034

#boundary of gdalinfo W2W2_20180611_103001007E7F1900_103001007F843600_2m_lsf_seg1_ortho.tif 
#-3429122 -3406392 2374856 2396770

#hi gdalwarp $infile $ofile3 -t_srs epsg:3413 -tr $ores $ores -te -3429122 2374856 -3406392 2396770
#gdalwarp $infile $ofile3 -t_srs epsg:3413 -tr $ores $ores 

#str=['gdalwarp ',ifile,' ',ofile,' -t_srs ',projgdal,' -te ', testr,' -tr ',num2str([imageores imageores]),' -ot UInt16'] 
#gdalwarp /fs/project/howat.4/images4sdm/Elliot/Orthophotos//20180804_merge.tif /fs/project/howat.4/dai.56/chunliwork/sdm/runelliot/runplanet/site1/imagesubdir//20180804000000_PLAN_B1_sub.TIF -t_srs epsg:32610 -te 382170  5643360   392160  5653350 -tr 4  4 -ot UInt16
#gdalwarp $infile $ofile3 -t_srs epsg:32610 -te 382170  5643360   392160  5653350 -tr 15 15 -ot UInt16
gdalinfo $infile
proj='epsg:32607'
region='426000	6870000	446000	6890000' #Chisana
region='323422.5         6533467.5          343417.5         6553462.5' #Grand Plateau
region='668077.5	6533977.5	688072.5	6553987.5' #Grand Plateau epsg:32607
#region='537037.5         7504957.5          557047.5         7524952.5'  #Allen River
region='391177.5         6726112.5          411187.5         6746107.5'  #Portage
#region='412117.5         6669997.5          432127.5         6689992.5' #Nassau
#region='369277.5         7044772.5          389287.5         7064767.5'  #Sanctuary
#region='420202.5         6765427.5          440212.5         6785422.5' #Serpentine
region='600907.5         6629287.5          620902.5         6649282.5'  #Nunatak
#gdalwarp $infile $ofile3 -te $region -r bilinear -tr 15 15 -ot UInt16 #interpolation makes it blurry
#gdalwarp $infile $ofile3 -t_srs  $proj -te $region -tr 15 15 -ot UInt16  #interpolate with nearest point, yielding shifts which can be coregistered. But itmay also show shifts in different directions.
gdalwarp $infile $ofile3 -t_srs  $proj -te $region -r bilinear -tr 15 15 -ot UInt16  #interpolate 

#rm $ofile1
#exit
done
