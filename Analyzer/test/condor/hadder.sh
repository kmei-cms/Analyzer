#!/bin/bash

dir=${CMSSW_BASE}/src/Analyzer/Analyzer/test/condor/output-files-slim-FULL

function strip {
    local STRING=${1#$"$2"}
    echo ${STRING%$"$2"}
}

# Make sure directory ends with "/*"
if [[ $dir != */ ]]
then
    d="$dir/*"
else
    d="$dir*"
fi

echo "${dir}"

# Loop over output-files
for f in $d
do
    # Only interested in directories
    [ -d "${f}" ] || continue

    # Get the name of the sample collection
    sc=$(strip "$f" "$dir/")
    
    if [ "${sc}" = "AllSignal" ] 
    then 
        continue
    fi

    #echo "hadd ${sc}.root *.root"
    hadd -fk haddfiles/${sc}.root ${f}/*.root

    #if [ "${sc}" = "Rare" ]
    #then
    #    #echo "hadd ${sc}.root *.root"
    #    hadd ${f}/${sc}.root ${f}/*.root
    #else
    #    #echo "hadd ${sc}.root *${sc}*.root"
    #    hadd ${f}/${sc}.root ${f}/*${sc}*.root
    #fi

done
#cp output-files-slim/AllSignal/MySlimAnalysis_rpv_stop_350_0.root haddfiles/rpv_stop_350.root
#cp output-files-slim/AllSignal/MySlimAnalysis_rpv_stop_450_0.root haddfiles/rpv_stop_450.root
#cp output-files-slim/AllSignal/MySlimAnalysis_rpv_stop_550_0.root haddfiles/rpv_stop_550.root
#cp output-files-slim/AllSignal/MySlimAnalysis_rpv_stop_650_0.root haddfiles/rpv_stop_650.root
#cp output-files-slim/AllSignal/MySlimAnalysis_rpv_stop_750_0.root haddfiles/rpv_stop_750.root
#cp output-files-slim/AllSignal/MySlimAnalysis_rpv_stop_850_0.root haddfiles/rpv_stop_850.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_850_SYY_0.root haddfiles/stealth_stop_850_SYY.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_750_SYY_0.root haddfiles/stealth_stop_750_SYY.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_650_SYY_0.root haddfiles/stealth_stop_650_SYY.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_550_SYY_0.root haddfiles/stealth_stop_550_SYY.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_450_SYY_0.root haddfiles/stealth_stop_450_SYY.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_350_SYY_0.root haddfiles/stealth_stop_350_SYY.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_350_SHuHd_0.root haddfiles/stealth_stop_350_SHuHd.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_450_SHuHd_0.root haddfiles/stealth_stop_450_SHuHd.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_550_SHuHd_0.root haddfiles/stealth_stop_550_SHuHd.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_650_SHuHd_0.root haddfiles/stealth_stop_650_SHuHd.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_750_SHuHd_0.root haddfiles/stealth_stop_750_SHuHd.root
#cp output-files-slim/AllSignal/MySlimAnalysis_stealth_stop_850_SHuHd_0.root haddfiles/stealth_stop_850_SHuHd.root
