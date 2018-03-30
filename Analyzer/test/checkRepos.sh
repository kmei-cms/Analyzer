#!/bin/bash

dir=${CMSSW_BASE}/src

# Make sure directory ends with "/"
if [[ $dir != */ ]]
then
    dir="$dir/*"
else
    dir="$dir*"
fi

# Loop all sub-directories
for f in $dir
do
    # Only interested in directories
    [ -d "${f}" ] || continue

    # Check if directory is a git repository
    if [ -d "$f/.git" ]
    then
        echo -en "\033[0;35m"
        echo "${f}"
        echo -en "\033[0m"

        mod=0
        cd $f
        
        # Fetch and check for changes
        if [ $(git remote | grep test -c) -ne 0 ]
        then
            git fetch test 
            git status -uno
        else
            echo "Need to:       git remote add test repo"
            echo "Repo options:" 
            echo "               https://github.com/StealthStop/Framework.git"
            echo "               https://github.com/StealthStop/Analyzer.git"
            echo "               https://github.com/susy2015/TopTagger.git"
            echo "               https://github.com/susy2015/TopTaggerTools.git"
        fi

        # Check for modified files
        if [ $(git status | grep modified -c) -ne 0 ]
        then
            mod=1
            echo -en "\033[0;31m"
            echo "Modified files"
            echo -en "\033[0m"
        fi

        # Check for untracked files
        if [ $(git status | grep Untracked -c) -ne 0 ]
        then
            mod=1
            echo -en "\033[0;31m"
            echo "Untracked files"
            echo -en "\033[0m"
        fi

        # Check if everything is peachy keen
        if [ $mod -eq 0 ]
        then
            echo "Nothing to commit"
        fi

        cd ../
    fi

done
