import os
from collections import defaultdict

privateNtuples = False

if privateNtuples:
    ## Making filelists for our own ntuples
    samplesets = [
        "WJetsToLNu",
        "diboson",
        "qcd",
        "qcd_bgenfilter",
        "signalR2",
        "tX",
        "ttWJetsToLNu",
        "ttX",
        "ttXX",
        "ttbar",
        "ttbar2",
        "tttX",
        "tttt"
        ]
    basedir = "/store/user/lpcsusyhad/StealthStop/TreeMaker_ntuples/"
    tempfilename = "tmp.txt"

    samples = defaultdict(list)

    for sampleset in samplesets:
        command = os.system("eos root://cmseos.fnal.gov ls %s/%s > %s" % (basedir, sampleset, tempfilename))
        with open(tempfilename, 'r') as tempfile:
            for line in tempfile:
                # extract sample name, Summer16_private.TTTT_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8_6_RA2AnalysisTree.root
                shortline = line.split(".")[1].rpartition("_")[0].rpartition("_")[0]
                shortline = shortline.replace("_ext1","").replace("_ext2","").replace("_backup","")
                newline = "root://cmseos.fnal.gov/" + basedir + sampleset + "/" + line
                samples[shortline].append(newline)
    
    for k,v in samples.iteritems():
        newfile = open("filelists/" + k + ".txt", 'w')
        for l in v:
            newfile.write(l)
        newfile.close()


else:

    prod = "V17"
    basedir = "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2Production%s/"%prod
    tempfilename = "tmp.txt"
    fdir = "filelists_Kevin_%s/"%prod
    if not os.path.isdir(fdir):
        os.mkdir(fdir)

    samples = defaultdict(list)

    command = os.system("eos root://cmseos.fnal.gov ls %s > %s" % (basedir, tempfilename))
    with open(tempfilename, 'r') as tempfile:
        for line in tempfile:
            if not ".root" in line: continue
            if "Fast" in line: continue
            print line
            era = ""
            if "Run2016" in line or "Summer16v3" in line: 
                era = "2016_"
            elif "Run2017" in line or "Fall17" in line: 
                era = "2017_"
            elif "Run2018" in line or "Autumn18" in line: 
                era = "2018_"
            else:
                continue
            # extract sample name, Summer16_private.TTTT_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8_6_RA2AnalysisTree.root
            # for data: Run2016B-17Jul2018_ver2-v1.SingleMuon_106_RA2AnalysisTree.root
            # for 2017 MC: RunIIFall17MiniAODv2.TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8_123_RA2AnalysisTree.root
            shortline = line.split(".")[1].rpartition("_")[0].rpartition("_")[0]
            shortline = era+shortline.replace("_ext1","").replace("_ext2","").replace("_ext3","").replace("_backup","")
            print shortline
            newline = "root://cmseos.fnal.gov/" + basedir + line
            samples[shortline].append(newline)

    for k,v in samples.iteritems():
        newfile = open("filelists_Kevin_%s/"%prod + k + ".txt", 'w')
        for l in v:
            newfile.write(l)
        newfile.close()
