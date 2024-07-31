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
source /cvmfs/fermilab.opensciencegrid.org/packages/common/setup-env.sh ## -- For Alma 9
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load root@6.28.12

echo "@@@@ Setup output directories"
export LArProfRunlogDir=$LArProf_WD/output/log
export LArProfOutputDir=$LArProf_WD/output/result

alias pdout="cd $LArProfOutputDir"

export MYBIN=$LArProf_WD/bin/
export PYTHONDIR=$LArProf_WD/python/
export PATH=${MYBIN}:${PYTHONDIR}:${PATH}

export ROOT_INCLUDE_PATH=/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/include/:$ROOT_INCLUDE_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PDSPAna_LIB_PATH:/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/lib/

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
