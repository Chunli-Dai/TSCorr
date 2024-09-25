
echo Start compiling files
rm -f run_Tilemain.sh mccExcludedFiles.log Tilemain

module unload matlab
module load matlab/2023b # # projcrs, only works for matlab Since R2020b
mcc -m Tilemain.m -a /blue/chunlidai/apps/horizontalv/code/ -a constant.m

cp /home/chunlidai/blue/apps/horizontalv/template/run_Tilemain.sh .
echo Compiling files end.

