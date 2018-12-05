#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/BTagCalibrationStandalone.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "SusyAnaTools/Tools/PileupWeights.h"
#include "SusyAnaTools/Tools/MiniTupleMaker.h"

#include "TopTagger/CfgParser/include/TTException.h"

#include "Analyzer/Analyzer/include/MakeSmallChecks.h"
#include "Analyzer/Analyzer/include/MakeSlimChecks.h"
#include "Analyzer/Analyzer/include/MakeDataMCPlots.h"
#include "Analyzer/Analyzer/include/AnalyzeHighNJetQCD.h"
#include "Analyzer/Analyzer/include/AnalyzeLepTrigger.h"
#include "Analyzer/Analyzer/include/AnalyzeStealthTopTagger.h"

#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Photon.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/RunFisher.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/DeepEventShape.h"
#include "Framework/Framework/include/ScaleFactors.h"

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"

#include <iostream>
#include <getopt.h>
#include <string>
#include <functional>

std::string color(const std::string& text, const std::string& color)
{
    std::string c;
    if(color=="red") c = "31";
    else if(color=="green") c = "32";
    else if(color=="yellow") c = "33";
    else if(color=="blue") c = "34";
    else if(color=="white") c = "37";
    
    return "\033[1;"+c+"m"+ text +"\033[0m";
}

template<typename Analyze> void run(std::set<AnaSamples::FileSummary> vvf, 
                                    int startFile, int nFiles, int maxEvts, 
                                    bool isSkim, TFile* outfile)
{
    std::cout << "Initializing..." << std::endl;
    const int jecOn = 0, jerOn = 0; 
    
    if( jecOn * jerOn != 0 || std::fabs(jecOn) > 1 || std::fabs(jerOn) > 1) 
    {
        std::cerr<<color("Error: ", "red")
                 <<"Invalid values of jecOn and jerOn. "
                 <<"They should be either -1, 0, or 1. "
                 <<"you cannot have both jet energy corrections and jet energy resolutions on at the same time"<<std::endl;
        return;
    }

    std::string                 myVarSuffix = "";
    if      ( jecOn == 1 )      myVarSuffix = "JECup";
    else if ( jecOn == -1 )     myVarSuffix = "JECdown";
    else if ( jerOn == 1 )      myVarSuffix = "JERup";
    else if ( jerOn == -1 )     myVarSuffix = "JERdown";
    
    Analyze a;
    for(const AnaSamples::FileSummary& file : vvf)
    {
        // Define what is needed per sample set
        std::cout << "Running over sample " << file.tag << std::endl;
        TChain* ch = new TChain( (file.treePath).c_str() );
        file.addFilesToChain(ch, startFile, nFiles);
        NTupleReader tr(ch);
        double weight = file.getWeight(); // not used currently
        std::string runtype = (file.tag.find("Data") != std::string::npos) ? "Data" : "MC";
        std::cout << "Starting loop (in run)" << std::endl;
        printf( "runtype: %s fileWeight: %f nFiles: %i startFile: %i maxEvts: %i \n",runtype.c_str(),weight,nFiles,startFile,maxEvts ); fflush( stdout );
        tr.registerDerivedVar<std::string>("runtype",runtype);
        tr.registerDerivedVar<std::string>("filetag",file.tag);
        tr.registerDerivedVar<double>("etaCut",2.4);
        tr.registerDerivedVar<bool>("blind",true);

        // Define classes/functions that add variables on the fly
        if ( !isSkim ) 
        {
            RunTopTagger rtt;
            tr.registerFunction(rtt);
        }
        
        Muon muon;
        Electron electron;
        Photon photon;
        Jet jet(myVarSuffix);
        BJet bjet(myVarSuffix);
        CommonVariables commonVariables( myVarSuffix );
        MakeMVAVariables makeMVAVariables(false, myVarSuffix,false);
        Baseline baseline;
        DeepEventShape deepEventShape;

        // Register classes/functions that add variables on the fly
        tr.registerFunction( muon );
        tr.registerFunction( electron );
        tr.registerFunction( jet );
        tr.registerFunction( bjet );
        tr.registerFunction( photon );
        tr.registerFunction( commonVariables );
        tr.registerFunction( makeMVAVariables );
        tr.registerFunction( baseline );
        tr.registerFunction( deepEventShape );

        if( runtype == "MC" ) 
        {
            //std::string eosPath =   "root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/";
            BTagCorrectorTemplate<double> bTagCorrector("allInOne_BTagEff.root","", false, file.tag);
            bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets", "Jets_bDiscriminatorCSV", "Jets_partonFlavor");
            Pileup_SysTemplate<double> pileup("PileupHistograms_0121_69p2mb_pm4p6.root");
            ScaleFactors scaleFactors("allInOne_leptonSF_Moriond17.root");
            tr.registerFunction(bTagCorrector);
            tr.registerFunction(pileup);
            tr.registerFunction(scaleFactors);
        }
        // Loop over all of the events and fill histos
        a.Loop(tr, weight, maxEvts);

        // Cleaning up dynamic memory
        delete ch;
            
    }
    std::cout << "Writing histograms..." << std::endl;
    a.WriteHistos(outfile);
}

std::set<AnaSamples::FileSummary> setFS(const std::string& dataSets, const bool& isCondor)
{
    AnaSamples::SampleSet        ss("sampleSets.cfg", isCondor);
    AnaSamples::SampleCollection sc("sampleCollections.cfg", ss);

    std::map<std::string, std::vector<AnaSamples::FileSummary>> fileMap;
    if(ss[dataSets] != ss.null())
    {
        fileMap[dataSets] = {ss[dataSets]};
        for(const auto& colls : ss[dataSets].getCollections())
        {
            fileMap[colls] = {ss[dataSets]};
        }
    }
    else if(sc[dataSets] != sc.null())
    {
        fileMap[dataSets] = {sc[dataSets]};
        int i = 0;
        for(const auto& fs : sc[dataSets])
        {
            fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
        }
    }
    std::set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    return vvf;
}

int main(int argc, char *argv[])
{
    int opt, option_index = 0;
    bool doStealthTT = false, makeSmallChecks = false, doLepTrig = false, doDataMC = false;
    bool runOnCondor = false, doHighNJetQCD = false, makeSlimChecks  = false;
    bool isSkim = false;
    std::string histFile = "", dataSets = "";
    int nFiles = -1, startFile = 0, maxEvts = -1;

    static struct option long_options[] = {
        {"doStealthTT",        no_argument, 0, 'x'},
        {"doLepTrig",          no_argument, 0, 't'},
        {"doDataMC",           no_argument, 0, 'd'},
        {"doHighNJetQCD",      no_argument, 0, 'q'},
        {"makeSmallChecks",    no_argument, 0, 'a'},
        {"makeSlimChecks",     no_argument, 0, 'b'},
        {"condor",             no_argument, 0, 'c'},
        {"histFile",     required_argument, 0, 'H'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
    };

    while((opt = getopt_long(argc, argv, "xtdqabcH:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 'x': doStealthTT      = true;              break;
            case 't': doLepTrig        = true;              break;
            case 'd': doDataMC         = true;              break;
            case 'a': makeSmallChecks  = true;              break;
            case 'q': doHighNJetQCD    = true;              break;
            case 'b': makeSlimChecks   = true;              break;
            case 'c': runOnCondor      = true;              break;
            case 'H': histFile         = optarg;            break;
            case 'D': dataSets         = optarg;            break;
            case 'N': nFiles           = int(atoi(optarg)); break;
            case 'M': startFile        = int(atoi(optarg)); break;
            case 'E': maxEvts          = int(atoi(optarg)); break;
        }
    }

    if(dataSets.find("skim") != std::string::npos) isSkim = true;  

    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "MySlimAnalysis_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
    }

    std::set<AnaSamples::FileSummary> vvf = setFS(dataSets, runOnCondor); 
    TFile* outfile = TFile::Open(histFile.c_str(), "RECREATE");

    std::vector<std::pair<bool, std::function<void(std::set<AnaSamples::FileSummary>,int,int,int,bool,TFile*)>>> AnalyzerPairVec = {
        {doStealthTT,      run<AnalyzeStealthTopTagger>},
        {doDataMC,         run<MakeDataMCPlots>},
        {doLepTrig,        run<AnalyzeLepTrigger>},
        {makeSmallChecks,  run<MakeSmallChecks>},
        {makeSlimChecks,   run<MakeSlimChecks>},
        {doHighNJetQCD,    run<AnalyzeHighNJetQCD>},
    }; 
    
    try
    {
        for(auto& pair : AnalyzerPairVec)
        {
            if(pair.first) pair.second(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        outfile->Close();
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

    return 0;
}
