# Analyzer
Runs all of our anaylsis code. 

## Using the tensor-flow based top tagger

To have easy access to TensorFlow, we need to work in a CMSSW93 release:
```
cmsrel CMSSW_9_3_3
cd CMSSW_9_3_3/src/
cmsenv
```

Then, check out the latest tagged version of the top tagger repository. 

```
git clone git@github.com:susy2015/TopTagger.git
```

Then configure and compile the tagger:
```
cd TopTagger/TopTagger/test
./configure 
make -j4
```

Now also check out our repository if not done already:
```
cd $CMSSW_BASE/src
git clone git@github.com:StealthStop/Framework.git
git clone git@github.com:susy2015/TopTaggerTools.git
git clone git@github.com:susy2015/SusyAnaTools.git
git clone git@github.com:StealthStop/Analyzer.git
cd Analyzer/Analyzer/test
source setup.csh
make -j4
```

Last step is to get the cfg and model files for the top tagger and deepESM.
```
cmsenv
getTaggerCfg.sh -t Tensorflow_Medium_Example_v1.0.2 -o
getDeepESMCfg.sh -t Keras_Tensorflow_v1.2.6 -o -s 2016
getDeepESMCfg.sh -t Keras_Tensorflow_v3.0.2 -o -s 2017
```

No changes to the analysis code should be needed. 


## Condor submission

The condor directory contains some scripts to help submit jobs via condor on the cmslpc cluster. 
The requirements for condor submission are: 
 - A script to run on the worker node. This script should set up the area, copy any needed files, call your executable with the right options, and make sure the output gets copied to where you want. The example included here is [run_Analyzer_condor.tcsh](Analyzer/test/condor/run_Analyzer_condor.tcsh)
 - One or more tarballs to unpack on the worker node, these usually contain a slimmed down CMSSW area, and your executable with any needed libraries
 - A so-called jdl file that contains the condor setup and specifies the jobs to be submitted
The last two items are produced by a python script called [condorSubmit.py](Analyzer/test/condor/condorSubmit.py). 

```
[condor]$ python condorSubmit.py -h
Usage: condorSubmit.py [options]


Options:
  -h, --help         show this help message and exit
  -n NUMFILE         number of files per job
  -d DATASETS        List of datasets, comma separated
  -l                 List all datacollections
  -L                 List all datacollections and sub collections
  -c                 Do not submit jobs.  Only create condor_submit.txt.
  --output=OUTPATH   Name of directory where output of each condor job goes
  --analyze=ANALYZE  AnalyzeBackground (b), AnalyzeEventSelection (s),
                     Analyze0Lep (z), Analyze1Lep (o), MakeNJetDists (n)
```
As you can see from the above help menu, there are a number of options. 
With the `-n` option you can specify how many files to run over per job. The `--analyze` option lets you pick which analyzer to use. 
The MyAnalysis program has been updated to have these same switches. 
MyAnalysis now also uses the samples code to keep track of datasets, their cross sections, 6nd their names. 
To see a list of available datasets, you can call the submission script with the `-l` or `-L` options. Pass the list of datasets you want to run over to the script with the option `-d`. 
Before submitting jobs, make sure to have called `voms-proxy-init`. 

## Making inputs for the fit

Running the condor jobs to produce the input histograms for the fit.

```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/condor
python condorSubmit.py -d 2016_Data_SingleElectron,2016_Data_SingleMuon,2016_TT,2016_TT_fsrUp,2016_TT_fsrDown,2016_TT_isrUp,2016_TT_isrDown,2016_WJetsToLNu,2016_DYJetsToLL_M-50,2016_QCD,2016_ST,2016_Diboson,2016_Rare,2016_AllSignal -n 10 --analyze n --output CondorOutput_Keras1.2.6_Final
python condorSubmit.py -d 2017_Data_SingleElectron,2017_Data_SingleMuon,2017_TT,2017_WJetsToLNu,2017_DYJetsToLL_M-50,2017_QCD,2017_ST,2017_Diboson,2017_Rare,2017_AllSignal -n 10 --analyze n --output CondorOutput_Keras3.0.2_Final
```

Now hadd the outputs when the jobs are done.

```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/condor
python hadder.py -d 2016_Data_SingleElectron,2016_Data_SingleMuon,2016_TT,2016_TT_fsrUp,2016_TT_fsrDown,2016_TT_isrUp,2016_TT_isrDown,2016_WJetsToLNu,2016_DYJetsToLL_M-50,2016_QCD,2016_ST,2016_Diboson,2016_Rare,2016_AllSignal -H MakeNJetsDists_Kerasv1.2.6_Final -p CondorOutput_Keras1.2.6_Final/output-files -y 2016 --doHack
python hadder.py -d  2017_Data_SingleElectron,2017_Data_SingleMuon,2017_TT,2017_WJetsToLNu,2017_DYJetsToLL_M-50,2017_QCD,2017_ST,2017_Diboson,2017_Rare,2017_AllSignal -H MakeNJetsDists_Kerasv3.0.2_Final -p CondorOutput_Keras3.0.2_Final/output-files -y 2017 --doHack
```

If there are missing jobs you should see a message in red after each sample is hadded.
After hadding now put all the samples together into one file for the fit input.

```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/
python write_fit_input.py -d condor/MakeNJetsDists_Kerasv1.2.6_Final -H FitInput/Keras_V1.2.6_Final -y 2016
python write_fit_input.py -d condor/MakeNJetsDists_Kerasv3.0.2_Final -H FitInput/Keras_V3.0.2_Final -y 2017
```
