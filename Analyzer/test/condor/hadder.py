import sys, os
from os import system, environ
import subprocess
sys.path = [environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/condor/",] + sys.path

from samples import SampleCollection
import optparse
from glob import glob
import datetime
import shutil

# Parse command line arguments
parser = optparse.OptionParser("usage: %prog [options]\n")
parser.add_option('-d', dest='datasets', type='string', default='',             help="Lists of datasets, comma separated")
parser.add_option('-H', dest='outDir',   type='string', default='rootfiles',    help="Can pass in the output directory name")
parser.add_option('-p', dest='inPath',   type='string', default='output-files', help="Can pass in the input directory name")
parser.add_option('-o', action='store_true',                                    help="Overwrite output directory")
options, args = parser.parse_args()

# Checks if user specified a dataset(s)
datasets = []
if options.datasets:
    datasets = options.datasets.split(',')
else:
    print "Failed: No dataset specified"
    exit(0)

# Check if output directory exits and makes it if not
outDir = options.outDir
overwrite = options.o
if os.path.exists(outDir):
    if overwrite: 
        print "Overwriting output directory"
        shutil.rmtree(outDir)
        os.makedirs(outDir)
    else:
        print "Failed: Output directory %s already exits" % ('"'+outDir+'"')
        exit(0)

else:
    os.makedirs(outDir) 

# Get input directory path
inPath = options.inPath

# Loop over all sample options to find files to hadd
sc = SampleCollection("../sampleSets.cfg", "../sampleCollections.cfg")
scl = sc.sampleCollectionList()
for sampleCollection in scl:
    sl = sc.sampleList(sampleCollection)
    if sampleCollection in datasets:
        directory = sampleCollection
        files = ""
        print "-----------------------------------------------------------"
        print sampleCollection
        print "-----------------------------------------------------------"
        
        # copy signal root files
        if sampleCollection == "AllSignal":
            for sample in sl:
                command = "cp %s/%s/MyAnalysis_%s_0.root %s/%s.root" % (inPath, directory, sample[1], outDir, sample[1])
                print command
                system(command)

        # hadd other condor jobs
        else:
            for sample in sl:
                files += " " + " ".join(glob("%s/%s/MyAnalysis_%s_*.root" % (inPath, directory, sample[1])))
            command = "hadd %s/%s.root %s" % (outDir, sampleCollection, files)
            try:
                process = subprocess.Popen(command, shell=True)
                process.wait()
            except:
                print "\033[91m Too many files to hadd: using the exception setup \033[0m"
                command = "hadd %s/%s.root %s/%s/*" % (outDir, sampleCollection, inPath, sampleCollection)
                system(command)
                pass

# Hack to make the BG_noTT.root file
sigNttbar = ["AllSignal", "TT", "TTJets", "Data_SingleMuon", "Data_SingleElectron"]
files = ""
for sampleCollection in scl:
    sl = sc.sampleList(sampleCollection)
    if sampleCollection in datasets:
        if sampleCollection not in sigNttbar: 
            directory = sampleCollection
            files += " %s/%s.root " % (outDir, directory)
command = "hadd %s/BG_noTT.root %s" % (outDir, files)
print "-----------------------------------------------------------"
print command
print "-----------------------------------------------------------"
system(command)

# Hack to make the Data.root file (hadd all the data together)
dataFiles = ["Data_SingleMuon.root", "Data_SingleElectron.root"]
command = "hadd %s/Data.root " % outDir
for f in dataFiles:
    command += " %s/%s" % (outDir, f)
print "-----------------------------------------------------------"
print command
print "-----------------------------------------------------------"
system(command)
