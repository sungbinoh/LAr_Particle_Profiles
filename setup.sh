export LArProf_WD=`pwd`
export LArProf_LIB_PATH=$LArProf_WD/lib/
mkdir -p $LArProf_LIB_PATH
mkdir -p $LArProf_WD/tar

export LArProfV="v1"
mkdir -p $LArProf_WD/data/$LArProfV
export DATA_DIR=$LArProf_WD/data/$LArProfV

#### USER INFO ####
export LArProfLogEmail='sungbin.oh@cern.ch'
export LArProfLogWeb=''
export LArProfLogWebDir=''

#### Setup Directories ####
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup root v6_22_08d -q e20:p392:prof

echo "@@@@ Setup output directories"
export LArProfRunlogDir=$LArProf_WD/output/log
export LArProfOutputDir=$LArProf_WD/output/result

alias pdout="cd $LArProfOutputDir"

export MYBIN=$LArProf_WD/bin/
export PYTHONDIR=$LArProf_WD/python/
export PATH=${MYBIN}:${PYTHONDIR}:${PATH}

export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$LArProf_WD/ProfileMakers/include/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LArProf_LIB_PATH

source $LArProf_WD/bin/BashColorSets.sh

## submodules ##
#source bin/CheckSubmodules.sh

## Todo list ##
python python/PrintToDoLists.py
source $LArProf_WD/tmp/ToDoLists.sh
rm $LArProf_WD/tmp/ToDoLists.sh

## Log Dir ##
echo "* Your Log Directory Usage"
#du -sh $LArProfRunlogDir
echo "-----------------------------------------------------------------"
CurrentGitBranch=`git branch | grep \* | cut -d ' ' -f2`
printf "> Current LAr_Particle_Profiles branch : "${BRed}$CurrentGitBranch${Color_Off}"\n"
