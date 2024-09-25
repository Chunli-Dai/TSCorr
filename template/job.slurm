#!/bin/bash
#SBATCH -t 95:59:00
#SBATCH --mem=70GB 
# #900GB
#SBATCH -N 1
#SBATCH -n 10 
#SBATCH -e %j.err
#SBATCH -o %j.out
#SBATCH --priority=TOP

#burst
#SBATCH --account=chunlidai
#SBATCH --qos=chunlidai-b

module load conda
conda info --envs
conda activate CosiCorr

# subtile=#utm37n_42_02_1_1_01_01
file=input.txt
subtile=$(awk 'NR==2' $file)
projgdal=$(awk -F"'" '/projgdal/{print $2}' constant.m)
str1=${subtile:0:3} # e.g.,'utm'
if [ "$str1" = "utm" ]; then
    echo "UTM projection"
echo Start compiling files
rm -f run_Tilemain.sh mccExcludedFiles.log Tilemain
./compile.sh 
cp /home/chunlidai/blue/apps/horizontalv/template/run_Tilemain.sh .
echo Compiling files end.

else
    echo "Not UTM Projection; Use the compiled Tilemain files; No need to re-compile files."
fi

# User variables
#matlab < text1.m
#./run_Tilemain.sh  /apps/matlab/r2020a
#./run_Tilemain.sh  /apps/matlab/r2020b
./run_Tilemain.sh  /apps/matlab/r2023b
 # projcrs, only works for matlab Since R2020b

#remove files
pwdstr=`pwd`
pwdstrlast=`echo $(basename $pwdstr)`
#if [ "$pwdstrlast" = "$workdir1" ]; then
if true
then
echo Remove files in this directory `pwd`.
rm -f mat0.mat t*.mat 
rm -f core.*
rm -f core
rm -f coreg.mat  Tilemain Tilemain.ctf*
rm -rf imagesubdir/*
rm -rf results/*
# rm -f *.fig
fi
