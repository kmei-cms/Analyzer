import sys, os
from os import system, environ
import subprocess
sys.path = [environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/condor/",] + sys.path

from samples import SampleCollection
import optparse
from glob import glob
import datetime
import shutil
import ROOT

def red(string):
     CRED = "\033[91m"
     CEND = "\033[0m"
     return CRED + str(string) + CEND

def checkNumEvents(nEvents, rootFile):
    try:
         f = ROOT.TFile.Open(rootFile)
         try:
              h = f.Get("EventCounter")
              nNeg = h.GetBinContent(1)
              nPos = h.GetBinContent(2)
              diff = float(nEvents)-(nPos-nNeg)
              if abs(diff) > 10.0:
                   print red("------------------------------------------------------------------------------------------------")
                   print red("Num events in \"EventCounter\" doesn't match the number in \"sampleSet.cfg\"")
                   print "SampleSet nEvents: ", red(nEvents), "EventCounter nEvents: ", red(nPos-nNeg), "=", red(nPos), red(-nNeg)
                   print red("------------------------------------------------------------------------------------------------")
         except:
              print red("Problem opening and reading from histogram \"EventCounter\"")
              pass
         f.Close()
    except:
         print red("Can't open rootFile: %s" % rootFile)
         pass

def main():
    # Parse command line arguments
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-d', dest='datasets', type='string', default='',             help="Lists of datasets, comma separated")
    parser.add_option('-H', dest='outDir',   type='string', default='rootfiles',    help="Can pass in the output directory name")
    parser.add_option('-p', dest='inPath',   type='string', default='output-files', help="Can pass in the input directory name")
    parser.add_option('-y', dest='year',     type='string', default='',             help="Can pass in the year for this data")
    parser.add_option('-o', action='store_true',                                    help="Overwrite output directory")
    parser.add_option('--doHack', action='store_true',                              help="Do the hack to make Data.root and BG_noTT.root")
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
            
            # hadd signal root files
            if sampleCollection == "AllSignal" or sampleCollection == "2016_AllSignal" or sampleCollection == "2017_AllSignal":
                for sample in sl:
                    files = " " + " ".join(glob("%s/%s/MyAnalysis_%s_*.root" % (inPath, directory, sample[1])))
                    outfile = "%s/%s.root" % (outDir,sample[1])
                    command = "hadd %s/%s.root %s" % (outDir, sample[1], files)
                    system(command)
                    checkNumEvents(nEvents=sample[2], rootFile=outfile)
    
            # hadd other condor jobs
            else:
                nEvents=0
                for sample in sl:
                    files += " " + " ".join(glob("%s/%s/MyAnalysis_%s_*.root" % (inPath, directory, sample[1])))
                    nEvents+=sample[2]

                outfile = "%s/%s.root" % (outDir,sampleCollection)
                command = "hadd %s %s" % (outfile, files)
                try:
                    process = subprocess.Popen(command, shell=True)
                    process.wait()
                except:
                    print "\033[91m Too many files to hadd: using the exception setup \033[0m"
                    command = "hadd %s/%s.root %s/%s/*" % (outDir, sampleCollection, inPath, sampleCollection)
                    system(command)
                    pass

                checkNumEvents(nEvents=nEvents, rootFile=outfile)
    
    if options.doHack:
        # Hack to make the BG_noTT.root file
        sigNttbar_old = ["AllSignal", "TT", "TTJets", "Data_SingleMuon", "Data_SingleElectron"]
        sigNttbar_2016 = ["2016_AllSignal", "2016_TT", "2016_TTJets", "2016_Data_SingleMuon", "2016_Data_SingleElectron"]
        sigNttbar_2017 = ["2017_AllSignal", "2017_TT", "2017_TTJets", "2017_Data_SingleMuon", "2017_Data_SingleElectron"]
        sigNttbar = sigNttbar_old+sigNttbar_2016+sigNttbar_2017
        files = ""
        for sampleCollection in scl:
            sl = sc.sampleList(sampleCollection)
            if sampleCollection in datasets:
                if sampleCollection not in sigNttbar: 
                    directory = sampleCollection
                    files += " %s/%s.root " % (outDir, directory)
        if options.year:
            command = "hadd %s/%s_BG_noTT.root %s" % (outDir, options.year, files)
        else:
            command = "hadd %s/BG_noTT.root %s" % (outDir, files)
        print "-----------------------------------------------------------"
        print command
        print "-----------------------------------------------------------"
        system(command)
        
        # Hack to make the Data.root file (hadd all the data together)
        dataFiles = ["Data_SingleMuon.root", "Data_SingleElectron.root", 
                     "2016_Data_SingleMuon.root", "2016_Data_SingleElectron.root", 
                     "2017_Data_SingleMuon.root", "2017_Data_SingleElectron.root"]
        if options.year:
            command = "hadd %s/%s_Data.root " % (outDir,options.year)
        else:
            command = "hadd %s/Data.root " % outDir
        for f in dataFiles:
            if os.path.exists(outDir+"/"+f):
                command += " %s/%s" % (outDir, f)
        print "-----------------------------------------------------------"
        print command
        print "-----------------------------------------------------------"
        system(command)

if __name__ == "__main__":
    main()
