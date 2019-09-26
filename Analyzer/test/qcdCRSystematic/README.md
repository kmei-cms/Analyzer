README for QCD CR Systematic Folder

##Last updated by Kelvin Mei on September 23rd, 2019

The purpose of the scripts in this folder is to create a systematic from looking at data in the QCD control region and the Powheg TT MC. The general idea is that the trends in high njets in the MVA shape of the neural network trainings are similar, and so an additional multiplicative weight can be derived from data and from MC and the ratio of the two can be the backbone of a new systematic, quantifying how different the data behaves than the MC.

To run these scripts, you will first need:

* The mini tree that has the totalEventWeight (all SFs applied) for the TT MC that pass the signal region baseline
* The mini tree that has the data (Single Muon data set) in the QCD control region that pass the control region baseline
* The njets\_for\_Aron.root file that has the final tt bar njets shape (this can be made from the mini trees if you do not have it.
* (Optional) The mini tree that has the QCD MC that pass the signal region baseline and a separate one for events that pass the control region baseline. The code currently is set so that it runs with the latter, and it produces nice plots with them, but it is not necessary.

The order of the scripts are as follows:

1. `python makeMVAHistograms.py --signal`
2. `python makeMVAHistograms.py --signal --ttbar`
3. `python plotMVAComp.py`
4. `hadd ratios_[2016,2017,2018].root mvaRatioPlots/[2016,2017,2018]*.root`
5. `python makeRatioOfRatios.py`
6. `python throwToysForSystematic.py`
7. `python makeQCDCRSystematic.py`
