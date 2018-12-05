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

# get scale factor root files (Should fix this)
if (! -f allInOne_BTagEff.root) then
    echo ""
    echo "|--------------------------------------|"
    echo "|  Copying scale factor files (fix me) |"
    echo "|--------------------------------------|"
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/allInOne_BTagEff.root .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/PileupHistograms_0121_69p2mb_pm4p6.root .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/allInOne_leptonSF_Moriond17.root .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/CSVv2_Moriond17_B_H.csv .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/allInONe_HtSFDist_2016.root .
endif

# Check repos for updates
if ("$1" == "-s") then
    echo ""
    echo "|--------------------------------------|"
    echo "|      Checking repos for updates      |"
    echo "|--------------------------------------|"
    echo "If it asks for your password too many times you can do something like the following:"
    echo "         ssh-add ~/.ssh/id_rsa"
    status.sh
endif
