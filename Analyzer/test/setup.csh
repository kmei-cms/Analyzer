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
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/allInOne_BTagEff.root .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/allInOne_leptonSF_2016.root .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/allInOne_leptonSF_2017.root .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/CSVv2_Moriond17_B_H.csv .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/allInOne_SFMean.root .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/PileupHistograms_0121_69p2mb_pm4p6.root .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/pu_ratio.root .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/L1prefiring_jetpt_2017BtoF.root .

    #Need to add all root files to this "FullRun2" directory when they are produced: for now keep the old files from above (SUS-19-004_Final)
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/DeepCSV_102XSF_V1.csv . 
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/DeepCSV_2016LegacySF_V1.csv .  
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/DeepCSV_94XSF_V4_B_F.csv .
    xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/PileupHistograms_2018_69mb_pm5.root .
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
