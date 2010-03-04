#!/bin/bash

usage='Usage: -n <name>'

args=`getopt r: -- "$@"`
if test $? != 0
     then
         echo $usage
         exit 1
fi

eval set -- "$args"
for i 
  do
  case "$i" in
      -n) shift; name=$2;shift;;
  esac      
done

echo 'Checking CRAB status for ' ${name} 

if [ "X"${name} == "X" ]
    then
    echo "INVALID DIR NAME! Please give a valid one!"
    echo $usage
    exit 
fi

# setup crab environment
source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh;
eval `scramv1 runtime -sh`;
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh;

cd ${name};

crab -status;
nrun=`crab -status 2>&1 | grep -c RUN`;
npend=`crab -status 2>&1 | grep -c PEND`;
ndone=`crab -status 2>&1 | grep -c DONE`;

if [ "${nrun}" == "0" ] && [ "${npend}" == "0" ]
then
    echo "Run "${name} "is done..." "run:" $nrun "pend:" $npend "done:" $ndone
    crab -get
else
    echo "Run "${name} "NOT yet done..." "run:" $nrun "pend:" $npend "done:" $ndone
fi

#cd -;
