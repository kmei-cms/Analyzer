#!/bin/tcsh

set filelists=root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/filelists/
set filelists_Kevin=root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/filelists_Kevin/

set SCRIPTSDIR=${CMSSW_BASE}/src/Framework/Framework/scripts
if ! ( $PATH:q =~ *${SCRIPTSDIR}:q* ) then
    setenv  PATH ${SCRIPTSDIR}:$PATH
    echo "adding ${SCRIPTSDIR} to the path"
endif

# set up the top tagger libraries etc
echo "|-----------------------------|"
echo "|      Set up top tagger      |"
echo "|-----------------------------|"
echo ""
cmsenv
source ${CMSSW_BASE}/src/TopTagger/TopTagger/test/taggerSetup.csh
echo "sourced taggerSetup"

# get most up to date samples config
if (! -f sampleSets.cfg) then
    echo ""
    echo "|--------------------------------------|"
    echo "|     Soft linking the sample configs  |"
    echo "|--------------------------------------|"
    getSamplesCfg.sh
endif

## Check repos for updates
#echo ""
#echo "|--------------------------------------|"
#echo "|      Checking repos for updates      |"
#echo "|--------------------------------------|"
#status.sh

## Copy over filelists if they are not present or have changed
#echo""
#echo "|-----------------------------|"
#echo "|      Copying filelists      |"
#echo "|-----------------------------|"
#echo""
#if (! -d condor/filelists ) then
#    echo "No filelists found, copying from eos"
#    mkdir condor/filelists/
#    xrdcp -r ${filelists} condor/filelists/
#    ln -s condor/filelists filelists
#else
#    echo "You already have the filelists. To get the current version, delete condor/filelists and run this again"
#endif
#
#if (! -d condor/filelists_Kevin ) then
#    echo "No filelists_Kevin found, copying from eos"
#    mkdir condor/filelists_Kevin/
#    xrdcp -r ${filelists_Kevin} condor/filelists_Kevin/
#    ln -s condor/filelists_Kevin filelists_Kevin
#else
#    echo "You already have the filelists_Kevin. To get the current version, delete condor/filelists_Kevin and run this again"
#endif
