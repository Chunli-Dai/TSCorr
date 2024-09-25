#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MATLAB Runtime environment for the current $ARCH and executes 
# the specified command.
#
exe_name=$0
exe_dir=`dirname "$0"`

# Change the path to deployable archive
#cachedir=/home/chunlidai/blue/cache/$RANDOM
cachedir=/blue/chunlidai/cache/$(whoami)
mkdir -p $cachedir
export MCR_CACHE_ROOT=$cachedir

echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT="$1"
  echo ---
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
  export LD_LIBRARY_PATH;

    # User added environmen variables; For SETSM
  export PATH=$PATH:/home/chunlidai/blue/apps/landslide/code1/SETSM
  export LD_LIBRARY_PATH=/home/chunlidai/blue/apps/landslide/code1/SETSM/tools/tiff-4.0.3/lib:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/home/chunlidai/blue/apps/landslide/code1/SETSM/tools/libgeotiff-1.4.2/lib:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/home/chunlidai/blue/apps/landslide/code1/SETSM/tools/proj-5.1.0/lib:$LD_LIBRARY_PATH
    # For CosiCorr Plus
  export LD_LIBRARY_PATH=~/blue/apps/CCPlus/Geospatial-COSICorr3D/geoCosiCorr3D/lib:$LD_LIBRARY_PATH
  module load conda 
  conda info --envs
  conda activate CosiCorr

  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  shift 1
  args=
  while [ $# -gt 0 ]; do
      token=$1
      args="${args} \"${token}\"" 
      shift
  done
  eval "\"${exe_dir}/Tilemain\"" $args
fi

#rm $cachedir -rf

exit

