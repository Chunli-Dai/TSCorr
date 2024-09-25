# TSCorr
Time-Series image Correlation software (TSCorr)

Steps:
1\ cp /home/chunlidai/blue/apps/horizontalv/template/* .
2\ modify tilelist
3\ Update imagedir in constant.m
   Update projgdal, projstrin if you are running ArcticDEM tiles.
4\ #Get mat0.mat # This step only needs to be done once for each folder of images.
sbatch job_nov.slurm #run Tilemain_nov.m to get mat0.mat.
5\ #bundle all jobs
./compile.sh
nohup ./run_change_group_pgc_par.sh > outrun1 &



Reference:
Dai, C., Higman, B., Lynett, P.J., Jacquemart, M., Howat, I.M., Liljedahl, A.K., Dufresne, A., Freymueller, J.T., Geertsema, M., Ward Jones, M. and Haeussler, P.J., 2020. Detection and assessment of a large and potentially tsunamigenic periglacial landslide in Barry Arm, Alaska. Geophysical Research Letters, 47(22), p.e2020GL089800.
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020GL089800 
