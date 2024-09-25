
Steps:

step 1: download Landsat7/8 data from https://earthexplorer.usgs.gov/
selection criteria: collection 2, level 1 (L1TP), Band 8, cloud cover  0 to 10%; sun elevation: 20 to 90 deg 

step 2 (This step replaced by preparebatchcosicorr.m): 
subset to target area; and select pan band.
Prepare for ortho2utm2m.sh: matlab< prepare_te.m
./ortho2utm2m.sh

step 3: see ~/blue/apps/horizontalv/template/steps.txt
Start from Slowmoving.m or Tilemain.m
