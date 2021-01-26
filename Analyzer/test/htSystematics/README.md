# README for the H<sub>T</sub> systematics for the SUS-19-004 Analysis

## Last Updated by Kelvin Mei on January 26th, 2021

The purpose of the scripts in this folder is to derive the H<sub>T</sub> scale factor for the SUS-19-004 analysis.
This involves having the H<sub>T</sub> distributions for the ttbar MC after the signal region selection for events 
with 5 jets, 6 jets, and 7 jets, as well as the standard mini trees for events with 7 or more jets.

This set of scripts was derived from work done by Chris Madrid back in 2018, and utilizes his DataSetInfo python class.
The updated scritps are located in the sub-directory mainHtSystematic.

The order of scripts are as follows:
1. If you do not have the 5 jets, 6 jets, and 7 jets H<sub>T</sub> distributions, you can run makeHtStudyHistos.py, but this requires having the mini trees for ttbar in the signal region with events with 5 jets and 6 jets.
2. Afterwards, just run fitRatio.py


The second batch of scripts in this folder derives the two additional H<sub>T</sub> systematics that utilize
two other modified auxiliary scale factors:

1. The scale factor where the value is kept constant after H<sub>T</sub> > 2000 GeV in order to address low statistics at high H<sub>T</sub>
2. The scale factor where the explicit N<sub>jets</sub> dependence is removed from the scale factor (using the scale factor at N<sub>jets</sub> >= 7.