# Analyzer
Runs all of our anaylsis code. 

## Using the tensor-flow based top tagger

To have easy access to TensorFlow and UpRoot, we need to work in a CMSSW10_2_9 release:
```
cmsrel CMSSW_10_2_9
cd CMSSW_10_2_9/src/
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
git clone -b Stealth git@github.com:susy2015/TopTaggerTools.git
git clone git@github.com:susy2015/SusyAnaTools.git
git clone git@github.com:StealthStop/Analyzer.git
cd Analyzer/Analyzer/test
source setup.sh #.csh if in tcsh
./configure
make -j4
```

We set up the top tagger cfg files for per year, because per year has different b-tagger working points (WPs).
Last step is to get the cfg and model files for the top tagger and deepESM.
```
cmsenv
getTaggerCfg.sh -t StealthStop_DeepCSV_DeepResolved_DeepAK8_wp0.98_2016_v1 -f TopTaggerCfg_2016.cfg -o
getTaggerCfg.sh -t StealthStop_DeepCSV_DeepResolved_DeepAK8_wp0.98_2017_v1 -f TopTaggerCfg_2017.cfg -o
getTaggerCfg.sh -t StealthStop_DeepCSV_DeepResolved_DeepAK8_wp0.98_2018_v1 -f TopTaggerCfg_2018.cfg -o
getDeepESMCfg.sh -t Keras_Tensorflow_2016_v1.1 -o -s 2016
getDeepESMCfg.sh -t Keras_Tensorflow_2017_v1.1 -o -s 2017
getDeepESMCfg.sh -t Keras_Tensorflow_2018pre_v1.0 -o -s 2018pre
getDeepESMCfg.sh -t Keras_Tensorflow_2018post_v1.0 -o -s 2018post
```

Example of running MyAnalysis interactively
```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/
./MyAnalysis -A MakeNJetDists -H myoutputfile.root -D 2016_TT -E 1001
```


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
python condorSubmit.py --analyze MakeNJetDists -d 2016_Data_SingleElectron,2016_Data_SingleMuon,2016_TT,2016_TT_fsrUp,2016_TT_fsrDown,2016_TT_isrUp,2016_TT_isrDown,2016_WJets,2016_DYJetsToLL_M-50,2016_QCD,2016_ST,2016_Diboson,2016_Triboson,2016_TTX,2016_AllSignal -n 15 --output CondorOutput_2016_v1.1
python condorSubmit.py --analyze MakeNJetDists -d 2017_Data_SingleElectron,2017_Data_SingleMuon,2017_TT,2017_WJets,2017_DYJetsToLL_M-50,2017_QCD,2017_ST,2017_Diboson,2017_Triboson,2017_TTX,2017_AllSignal                                                             -n 15 --output CondorOutput_2017_v1.1
python condorSubmit.py --analyze MakeNJetDists -d 2018pre_Data_SingleElectron,2018pre_Data_SingleMuon,2018pre_TT,2018pre_WJets,2018pre_DYJetsToLL_M-50,2018pre_QCD,2018pre_ST,2018pre_Diboson,2018pre_Triboson,2018pre_TTX,2018pre_AllSignal                            -n 15 --output CondorOutput_2018pre_v1.0
python condorSubmit.py --analyze MakeNJetDists -d 2018post_Data_SingleElectron,2018post_Data_SingleMuon,2018post_TT,2018post_WJets,2018post_DYJetsToLL_M-50,2018post_QCD,2018post_ST,2018post_Diboson,2018post_Triboson,2018post_TTX,2018post_AllSignal                 -n 15 --output CondorOutput_2018post_v1.0
```

Now hadd the outputs when the jobs are done.

```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/condor
python hadder.py -d  2016_Data_SingleElectron,2016_Data_SingleMuon,2016_TT,2016_TT_fsrUp,2016_TT_fsrDown,2016_TT_isrUp,2016_TT_isrDown,2016_WJets,2016_DYJetsToLL_M-50,2016_QCD,2016_ST,2016_Diboson,2016_Triboson,2016_TTX,2016_AllSignal  -H MakeNJetsDists_2016_v1.1     -p CondorOutput_2016_v1.1/output-files     -y 2016     --haddOther --haddData
python hadder.py -d  2017_Data_SingleElectron,2017_Data_SingleMuon,2017_TT,2017_WJets,2017_DYJetsToLL_M-50,2017_QCD,2017_ST,2017_Diboson,2017_Triboson,2017_TTX,2017_AllSignal                                                              -H MakeNJetsDists_2017_v1.1     -p CondorOutput_2017_v1.1/output-files     -y 2017     --haddOther --haddData
python hadder.py -d  2018pre_Data_SingleElectron,2018pre_Data_SingleMuon,2018pre_TT,2018pre_WJets,2018pre_DYJetsToLL_M-50,2018pre_QCD,2018pre_ST,2018pre_Diboson,2018pre_Triboson,2018pre_TTX,2018pre_AllSignal                             -H MakeNJetsDists_2018pre_v1.0  -p CondorOutput_2018pre_v1.0/output-files  -y 2018pre  --haddOther --haddData
python hadder.py -d  2018post_Data_SingleElectron,2018post_Data_SingleMuon,2018post_TT,2018post_WJets,2018post_DYJetsToLL_M-50,2018post_QCD,2018post_ST,2018post_Diboson,2018post_Triboson,2018post_TTX,2018post_AllSignal                  -H MakeNJetsDists_2018post_v1.0 -p CondorOutput_2018post_v1.0/output-files -y 2018post --haddOther --haddData
```

If there are missing jobs you should see a message in red after each sample is hadded.
After hadding now put all the samples together into one file for the fit input.

```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/
python write_fit_input.py -d condor/MakeNJetsDists_2016_v1.1     -H FitInput/Keras_2016_v1.1     -y 2016
python write_fit_input.py -d condor/MakeNJetsDists_2017_v1.1     -H FitInput/Keras_2017_v1.1     -y 2017
python write_fit_input.py -d condor/MakeNJetsDists_2018pre_v1.0  -H FitInput/Keras_2018pre_v1.0  -y 2018pre
python write_fit_input.py -d condor/MakeNJetsDists_2018post_v1.0 -H FitInput/Keras_2018post_v1.0 -y 2018post
```

## Deriving ttbar shape systematics

Before running `njets_systs_comp.py` one must run `run_fits4ttbar_systematics.sh` in the HiggsAnalysis-CombinedLimit repository. This will stash away results in directories called TAG_YEAR_SYSTEMATIC

Once these have been generated, one can run `njets_systs_comp.py` with three arguments:

```
--fittag TAG
--year YEAR
--fitdir Base directory containing TAG_YEAR_SYSTEMATIC results
```

An example workflow to run would be something like:

```
cd $HOME/../../CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/

./run_fits4ttbar_systematics.sh

# Let this run for several hours...
# After it is finished successfully

cd $HOME/../../CMSSW_9_3_3/src/Analyzer/Analyzer/test/

# For example running on 2018pre
python njets_systs_comp.py --year 2018pre --fitdir $HOME/../../src/HiggsAnalysis/CombinedLimit/ --fittag Approval_StatErrPlusFullDev_12JetFix
```
