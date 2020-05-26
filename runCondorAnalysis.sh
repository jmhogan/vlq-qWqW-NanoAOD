#!/bin/bash

hostname
date

infilename=${1}
outfilename=${2}
outputDir=${3}
skipevents=${4}
maxevents=${5}

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
scramv1 project CMSSW CMSSW_10_2_16
echo "Made cmssw"
cd CMSSW_10_2_16

echo "unpacking tar"
tar -xf ../qWqW_nano.tar
rm ../qWqW_nano.tar

cd src/vlq-qWqW-NanoAOD/

echo "cmsenv"
eval `scramv1 runtime -sh`

XRDpath=root://cmsxrootd.fnal.gov/

echo "creating ${outfilename} by reading ${infilename}"
root -l -b -q runAnalysis.C\(\"${infilename}\",\"${outfilename}\",${skipevents},${maxevents}\)

echo "ROOT Files:"
ls -l *.root
rm puppiCorr.root
# copy output to eos

echo "xrdcp output for condor"
for FILE in *.root
do
  echo "xrdcp -f ${FILE} root://cmseos.fnal.gov/${outputDir}/${FILE}"
  xrdcp -f ${FILE} root://cmseos.fnal.gov/${outputDir}/${FILE} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    #rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done

echo "done"
