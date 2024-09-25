#/bin/sh -f
# 
# From /home/chunlidai/blue/chunlidai/landslide/arcticdem_34_alaska_north/run_change_group_pgc_par.sh

inputtype=5;
# 1 %based on input xid etc.
# 2 %find the block based on input coordinates;
# 3 %Use a rectangle box around the input coordinates;
# 5 %use a list of tile names

subtiledx=$(awk -F'[^0-9]+' '/subtiledx=/{print $2}' constant.m)
nx=$((50 / subtiledx))
ny=$((50 / subtiledx))

if [ $inputtype -eq 5  ] ; then

#use a list of tile names
#tilelist=('39_17_2_2' '40_17_1_2' '40_17_1_1' '39_17_1_2');
#tilelist=('39_17_1_2');

# get tilelist from a filelist
# ls -d ??_??_?_? > tilelist
i=0
for line in `cat tilelist`
do
i=$i+1;
tilelist[$i-1]=$line;
done

nlist=${#tilelist[@]}
echo Total number of tiles: $nlist.

echo ${tilelist[*]}

shrundir=`pwd`; #e.g. /u/sciteam/chunli/scratch/chunli/arcticdem_08_canada_baffin/
count=0
countj=1;

#PGC tilename convention "$yidc"_"$xidc"_"$yidsc"_"$xidsc"

for (( i=1; i<=$nlist; i++ ))
do
tilename=${tilelist[$i-1]} #39_17_2_2 "$yidc"_"$xidc"_"$yidsc"_"$xidsc"
			   # or utm10n_47_04_1_2

workdirtile1=$tilename
echo workdirtile1 $workdirtile1
if [ ! -d $workdirtile1 ] 
then
    mkdir $workdirtile1
fi

# 39_17_2_2
cd $workdirtile1
pwd
#pwdsv=`pwd`;
dir_tile50km=`pwd`;

jobidg=(); # reset $jobid matrix

# subtiles

for (( yidss=1; yidss<=ny; yidss++ )) #1:5
do
for (( xidss=1; xidss<=nx; xidss++ )) #1:5
do
yidssc=`printf "%2.2d" "$yidss"`
xidssc=`printf "%2.2d" "$xidss"`
subtile="$tilename"_"$yidssc"_"$xidssc" #utm37n_42_02_1_1_01_01 #"$yidc"_"$xidc"_"$yidsc"_"$xidsc"_"$yidssc"_"$xidssc"
workdir=$dir_tile50km/$subtile # full path
echo yidss $yidssc workdir $workdir
if [ ! -d $workdir ]
then
    mkdir $workdir
fi
cd $workdir
pwd

ln -fs $shrundir/mat0.mat .
cp  $shrundir/job.slurm .
cp  $shrundir/compile.sh .
cp  $shrundir/*.m .
cp  $shrundir/correlate_cli.py .
cp  $shrundir/cosiplus.json .

#ln -fs $shrundir/run_Tilemain.sh .
ln -fs $shrundir/Tilemain .
cp /home/chunlidai/blue/apps/horizontalv/template/run_Tilemain.sh .

flageq=$(awk -F'[=;]' '/flageq=/{print $2}' constant.m)
if [ $flageq -eq 1 ]; then
    namestr='disp'  # earthquake, displacement in meter
elif [ $flageq -eq 0 ]; then
    namestr='rate'  # landslide, rate in meter/year
fi

jumpfile=$subtile"_"$namestr".tif" #e.g., utm37n_42_02_1_1_01_01_dispx.tif
if [ -f "$jumpfile" ]; then
    echo "$jumpfile exists. Continue to the next tile. "
else
    echo "$jumpfile does not exist. Carry out the calculation. "

#Update projection information in constant.m based on tile names (e.g., utm10n_47_04_2_2).
projgdal=$(awk -F"'" '/projgdal/{print $2}' constant.m)
str1=${subtile:0:3} # e.g.,'utm'
if [ "$str1" = "utm" ]; then
    echo "UTM projection"
    utmid=${subtile:3:3}; #e.g.,10n
    #sed -i "s|37N|$utmid|g" constant.m
    nsstr=${subtile:5:1}; # 'n' or 's'
    #Northern hemisphere, e.g., UTM Zone 3N: EPSG 32603;
    if [[ "${nsstr,,}" == 'n' ]]; then
    sed -i "s|epsg:32637|epsg:326${utmid:0:2}|g" constant.m
    sed -i "s|UTM zone 37 north|UTM zone ${utmid:0:2} north|g" constant.m
     #Southern hemisphere, e.g., UTM Zone 2S: EPSG 32702;
    elif [[ "${nsstr,,}" == 's' ]]; then
    sed -i "s|epsg:32637|epsg:327${utmid:0:2}|g" constant.m
    sed -i "s|UTM zone 37 north|UTM zone ${utmid:0:2} south|g" constant.m
    fi
elif [ $projgdal = "epsg:3413" ]; then  
    echo "Polar Stereographic North projection"
#    sed -i "s|projgdal='epsg:32637';|projgdal='epsg:3413;'|g" constant.m
    # projstrin update
    sed -i "s|UTM zone 37 north|polar stereo north|g" constant.m

elif [ $projgdal = "epsg:3031" ]; then
    echo "Antarctic Polar Stereographic projection, epsg:3031"
    sed -i "s|UTM zone 37 north|polar stereo south|g" constant.m
else  
	echo Projection is not as expected. Double check projgdal projstrin in constant.m
fi

echo 1 > input.txt
echo $subtile >> input.txt

sbatch job.slurm

#wait if the current number of jobs in queue is >=450
while true
do
sleep 1s #wait 5 seconds
njobs=`squeue --me | wc -l`
maxjobs=2; #300 #6; #16; # max: 9*40 for 5 years; 10* 180/40 *40 until august 2023
let "njobs-=1"
if [[ $njobs -lt $maxjobs ]]
then
   echo The number of current jobs is $njobs ', < '$maxjobs . 
   break
else
   echo The number of current jobs is $njobs ', wait to submit new jobs until it is <' $maxjobs .
   sleep 10m #wait 1 minute
fi
done #while
# wait

fi #if jumpfile exist

cd $dir_tile50km/ #../  #41_16_1_1/

done #xidss
done #yidss

cd $shrundir/ #../  #$shrundir

echo list all jobs ${jobidg[@]} #list all jobs


done # tile list for (( i=1; i<=$nlist; i++ ))

fi


