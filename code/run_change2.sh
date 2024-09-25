#/bin/sh -f
# 
# Modified from /home/dai.56/chunli/scripts/run_overlap_strip.sh 

inputtype=3; 
#inputtype=1; 
# 1 %based on input xid etc.
# 2 %find the block based on input coordinates;
# 3 %Use a rectangle box around the input coordinates;

if [ $inputtype -eq 2 -o $inputtype -eq 3  ] ; then
#loneq=;lateq=;workdir=site1;
file='aoi.txt'
nline=`wc -l < $file`
echo $nline
#5 11 12 16 17 18 22 24 25
for (( i=1; i<=$nline; i++ )) #1:80
#for i in 5 11 12 16 17 18 22 24 25
do 
latloneq=`sed -n ''$i'p' aoi.txt |awk  -F' ' '{ print $1, $2 }'`

workdir=site$i
if [ ! -d $workdir ] 
then
    mkdir $workdir
fi
cd $workdir
ln -fs ../mat0.mat .
ln -fs ../*.m .
ln -fs ../job*pbs .

echo $inputtype > input.txt
echo $latloneq >> input.txt
qsub jobpar.pbs
cd ../
#exit
done #i

else #inputtype=1
#50_08_2_1_2m_v3.0 
for (( xid=8; xid<=8; xid++ )) #1:80
do 
for (( yid=50; yid<=50; yid++ )) #1:80
do 
for (( xids=2; xids<=2; xids++ )) #1:2
do 
for (( yids=1; yids<=1; yids++ ))  #1:2
do 
for (( xidss=1; xidss<=1; xidss++ )) #1:25
do 
for (( yidss=1; yidss<=25; yidss++ )) #1:25
do 

xidc=`printf "%2.2d" "$xid"`
yidc=`printf "%2.2d" "$yid"`
xidsc=`printf "%1.1d" "$xids"`
yidsc=`printf "%1.1d" "$yids"`
xidssc=`printf "%2.2d" "$xidss"`
yidssc=`printf "%2.2d" "$yidss"`
workdir="$yidc"_"$xidc"_"$xidsc"_"$yidsc"_"$xidssc"_"$yidssc"
echo $workdir
if [ ! -d $workdir ] 
then
    mkdir $workdir
fi
cd $workdir
ln -fs ../mat0.mat .
ln -fs ../*.m .
ln -fs ../job*pbs .

echo $inputtype > input.txt
echo $workdir >> input.txt
qsub jobpar.pbs
cd ../

exit

done #yidss
done #xidss
done #yids
done #xids
done #yid
done  #xid

fi

#xr=`awk '$1!="X:"{next};!i++{min=$2;max=$2;}{for(j=2;j<=NF;++j){min=(min<$j)?min:$j;max=(max>$j)?max:$j}}END{printf " %d\n", max}' $infile`
#echo x range: $xl $xr '; y range:' $yl $yu
#  if [ $yu -ge $yl0 -a $yu -le $yu0 -o $yl -ge $yl0 -a $yl -le $yu0 ] ; then
#  fi

