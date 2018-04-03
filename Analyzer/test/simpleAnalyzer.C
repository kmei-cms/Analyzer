#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/SATException.h"

//Hack--Joes update to samples.cc/h should fix this
#include "Framework/Framework/include/samples.h"
//#include "SusyAnaTools/Tools/samples.h"

#include "TopTaggerTools/Tools/include/HistoContainer.h"
#include "TopTagger/CfgParser/include/TTException.h"

//#include "derivedTupleVariables.h"
//#include "baselineDef.h"
//#include "BTagCorrector.h"
//#include "TTbarCorrector.h"
//#include "ISRCorrector.h"
//#include "PileupWeights.h"
//#include "customize.h"

#include "TopTagger/TopTagger/include/TopTaggerResults.h"
//#include "Constituent.h"

#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

#include "math.h"

#include "Math/VectorUtil.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TChain.h"

void stripRoot(std::string &path)
{
    int dot = path.rfind(".root");
    if (dot != std::string::npos)
    {
        path.resize(dot);
    }
}

int main(int argc, char* argv[])
{

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                       //
    //                         Defining and setting up the options                           //
    //                                                                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////
    bool runOnCondor = false, enableTTbar = false, doWgt = true, runStealth = false;
    int nFiles = -1, startFile = 0, nEvts = -1;
    std::string dataSets = "TT", filename = "example.root";
    std::string sampleloc = AnaSamples::fileDir;    
    int opt;
    int option_index = 0;
    
    static struct option long_options[] = {
        {"condor",              no_argument, 0, 'c'},
        {"TTbar weight",        no_argument, 0, 't'},
        {"no event weighting",  no_argument, 0, 'd'},
        {"run stealth version", no_argument, 0, 's'},
        {"dataSets",      required_argument, 0, 'D'},
        {"numFiles",      required_argument, 0, 'N'},
        {"startFile",     required_argument, 0, 'M'},
        {"numEvts",       required_argument, 0, 'E'},
        {"output",        required_argument, 0, 'O'}
    };
    
    while((opt = getopt_long(argc, argv, "ctdsD:N:M:E:O:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'c':
            runOnCondor = true;
            std::cout << "Configured for condor compatibility." << std::endl;
            break;
    
        case 't':
            enableTTbar = true;
            std::cout << "Enabled TTbar event weighting." << std::endl;
            break;
    
        case 'd':
            doWgt = false;
            std::cout << "No Event weighting." << std::endl;
            break;
    
        case 's':
            runStealth = true;
            std::cout << "Running stealth verison" << std::endl;
            break;
    
        case 'D':
            dataSets = optarg;
            std::cout << "Running over the " << dataSets << " data sets." << std::endl;
            break;
    
        case 'N':
            nFiles = int(atoi(optarg));
            std::cout << "Running over " << nFiles << " files." << std::endl;
            break;
    
        case 'M':
            startFile = int(atoi(optarg));
            std::cout << "Starting on file #" << startFile << std::endl;
            break;
    
        case 'E':
            nEvts = int(atoi(optarg));
            std::cout << "Events: " << nEvts << std::endl;
            break;
    
        case 'O':
            filename = optarg;
            std::cout << "Filename: " << filename << std::endl;
        }
    }
    
    //if running on condor override all optional settings
    if(runOnCondor)
    {
        char thistFile[128];
        stripRoot(filename);
        sprintf(thistFile, "%s_%s_%d.root", filename.c_str(), dataSets.c_str(), startFile);
        filename = thistFile;
        std::cout << "Filename modified for use with condor: " << filename << std::endl;
        sampleloc = "condor";
    }
    
    TH1::AddDirectory(false);
    
    bool savefile = true;
    if(filename == "-"){
        savefile = false;
        std::cout << "Histogram file will not be saved." << std::endl;
    }

    std::cout << "Sample location: " << sampleloc << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                       //
    //                   Getting/defining everything for the sample set                      //
    //                                                                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////
    AnaSamples::SampleSet        ss(sampleloc, AnaSamples::luminosity);
    AnaSamples::SampleCollection sc(ss);

    //if(dataSets.find("Data") != std::string::npos){
    //   std::cout << "This looks like a data n-tuple. No weighting will be applied." << std::endl;
    //   doWgt = false;
    //}
    //
    //if(dataSets.find("TT") != std::string::npos){
    //   std::cout << "This looks like a TTbar sample. Applying TTbar weighting" << std::endl;
    //   enableTTbar = true;
    //}

    std::cout << "Dataset: " << dataSets << std::endl;

    int events = 0, peventsLep0 = 0, peventsLep1 = 0;
    HistoContainer<NTupleReader> hists0Lep("Lep0"), hists1Lep("Lep1");
    TRandom* trand = new TRandom3();

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                       //
    //                   Running over all the events in the sample set                       //
    //                                                                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////
    try
    {
        auto& fs = ss[dataSets];
        {
            TChain *t = new TChain(fs.treePath.c_str());
            fs.addFilesToChain(t, startFile, nFiles);
            double fileWgt = fs.getWeight();
            NTupleReader tr(t);            
            
            std::cout << "File: " << fs.filePath << std::endl;
            std::cout << "Tree: " << fs.treePath << std::endl;

            const int printInterval = 1000;
            int printNumber = 0;
            while(tr.getNextEvent())
            {
                events++;

                if(nEvts > 0 && tr.getEvtNum() > nEvts) break;
                if(tr.getEvtNum() / printInterval > printNumber)
                {
                    printNumber = tr.getEvtNum() / printInterval;
                    std::cout << "Event #: " << printNumber * printInterval << std::endl;
                }

                const double isData = false;//!tr.checkBranch("genDecayLVec");

                double eWeight = fileWgt;

                //Stealth Event Selection - 0 Lepton
                if( !isData  //lets not accidently unblind the stealth SR
                    //&& passNoiseEventFilter 
                    //&& passLeptonVeto
                    //&& cntNJetsPt45Eta24 >= 6 
                    //&& (ht > 500)
                    //&& nbCSV >= 1                //Atleast 1 medium B-Jet
                    )
                {
                    peventsLep0++;
                    hists0Lep.fill(tr, eWeight, trand);
                }

                //Stealth Event Selection - 1 Lepton
                if( !isData  //lets not accidently unblind the stealth SR
                    //&& passNoiseEventFilter 
                    //&& passSingleLep30
                    //&& cntNJetsPt30Eta24 >= 6
                    //&& nbCSV >= 1               //Atleast 1 medium B-Jet
                    //&& passLepTtag
                    )
                {
                    peventsLep1++;
                    hists1Lep.fill(tr, eWeight, trand);
                }
            }
        }
    }
    catch(const std::string e)
    {
        std::cout << e << std::endl;
        return 0;
    }
    catch(const TTException e)
    {
        std::cout << e << std::endl;
        return 0;
    }
    catch(const SATException e)
    {
        std::cout << e << std::endl;
        return 0;
    }

    std::cout << "Processed " << events << " events. " << peventsLep0 << " passed Lep0 selection. " << peventsLep1 << " passed Lep1 selection." << std::endl;

    if(savefile)
    {
        std::cout << "Saving root file..." << std::endl;

        TFile *f = new TFile(filename.c_str(),"RECREATE");
        if(f->IsZombie())
        {
            std::cout << "Cannot create " << filename << std::endl;
            throw "File is zombie";
        }

        hists0Lep.save(f);
        hists1Lep.save(f);

        f->Write();
        f->Close();
    }
}
