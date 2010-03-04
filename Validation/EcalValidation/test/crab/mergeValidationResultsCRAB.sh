#!/bin/bash

usage='Usage: -n <name> -d <crab_subdirectory> -o <output_dir> -w <work_dir> -m <merge_dir>'

args=`getopt rdowm: -- "$@"`
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
      -d) shift; crab_dir=$2;shift;;
      -o) shift; output_dir=$2;shift;;
      -w) shift; work_dir=$2;shift;;
      -m) shift; merge_dir=$2;shift;;
  esac      
done

#this_dir=`pwd`;
#echo "this is your working dir"
#echo `pwd`

if [ "X"${name} == "X" ]
    then
    echo "INVALID DIR NAME! Please give a valid one!"
    echo $usage
    exit
fi

if [ "X"${work_dir} == "X" ]
    then
    work_dir=`pwd`;
    echo " using default work dir" ${work_dir}
else
    echo " using work dir "${work_dir}
fi

if [ "X"${merge_dir} == "X" ]
    then
    merge_dir=${name}
    echo " using default merge dir" ${merge_dir}
else
    echo " using merge dir "${merge_dir}
fi

if [ "X"${crab_dir} == "X" ]
    then
    crab_dir=`\ls -rt1 ${work_dir}/${merge_dir} | grep "crab_" | tail -1 | awk '{print $NF}'`;
    echo " using default output dir" ${crab_dir}
else
    echo " using output dir "${crab_dir}
fi


echo 'Merging CRAB output ' ${name} 'crab_dir' ${crab_dir}

cd ${work_dir}/${merge_dir}/${crab_dir}/res;
#pwd;

# check root files
nroot=`\ls EcalValidation_*root | grep -vc ${name}`;
nmergeroot=`\ls EcalValidation_*root | grep -c ${name}`;
#echo $nmergeroot

if [ "${nroot}" == "0" ] && [ "${nmergeroot}" == "0" ]
then
    echo " NO root files" $nroot ".. exiting"
    exit
else
    echo " $nroot root files, $nmergeroot merged files"
fi

#rm -f EcalValidation_${name}.root
# now to hadd
hadd -f EcalValidation_${name}.root EcalValidation_*root

# check log files
nstdout=`\ls -l CMSSW*stdout | wc | awk '{print $1}'`;
nmergeout=`\ls EcalValidation_*log | grep -c ${name}`;
if [ "${nstdout}" == "0" ] #&& [ "${nmergeout}" == "0" ]
then
    echo " NO stdout files" $nstdout 
    #exit
else
    echo " $nstdout stdout files, $nmergeout merged files"
    cat CMSSW*stdout > EcalValidation_${name}.log
fi

gzip -f EcalValidation_${name}.log

if [ "X"${output_dir} == "X" ]
    then
    output_dir=/castor/cern.ch/user/c/ccecal/VALIDATION
    echo " use default output dir" $output_dir
    rfmkdir ${output_dir}
    rfchmod 775 ${output_dir}
    rfcp EcalValidation_${name}.log.gz ${output_dir}
    rfcp EcalValidation_${name}.root ${output_dir}
    rfchmod 775 ${output_dir}/EcalValidation_${name}.log.gz
    rfchmod 775 ${output_dir}/EcalValidation_${name}.root
else
    echo " using output dir "${output_dir}
    rfmkdir ${output_dir}
    rfchmod 775 ${output_dir}
    rfcp EcalValidation_${name}.log.gz ${output_dir}
    rfcp EcalValidation_${name}.root ${output_dir}
    rfchmod 775 ${output_dir}/EcalValidation_${name}.log.gz
    rfchmod 775 ${output_dir}/EcalValidation_${name}.root
#    rm -f EcalValidation_${name}.log
fi

size=`\ls -l EcalValidation_${name}.root | grep -c ${name}`;
echo $size

if [ "${size}" == "0" ]
then
    echo " Warning: your merged file has size ZERO.. will not delete root files."
else
    echo
#    echo " Warning: deleting root files"
#    \ls EcalValidation_*root | grep -v ${name} | awk '{print "rm -f "$NF}' |sh
fi

#rm -rf CMSSW*stdout
#rm -rf CMSSW*stderr
#rm -rf crab*xml

#\ls -l;

mv EcalValidation_${name}.root ${work_dir}/${merge_dir}
