#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "TopTagger/CfgParser/include/TTException.h"

#include "SusyAnaTools/Tools/BTagCalibrationStandalone.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "SusyAnaTools/Tools/PileupWeights.h"

#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/RunFisher.h"
#include "Framework/Framework/include/DeepEventShape.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/Photon.h"
#include "Framework/Framework/include/ScaleFactors.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/Baseline.h"

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
    bool runOnCondor = false;
    bool isSkim = false;
    std::string histFile = "", dataSets = "";
    int nFiles = -1, startFile = 0, maxEvts = -1;

    static struct option long_options [] = {
        {"condor",            no_argument, 0, 'c'},
        {"histFile",    required_argument, 0, 'H'},
        {"dataSets",    required_argument, 0, 'D'},
        {"numFiles",    required_argument, 0, 'N'},
        {"startFile",   required_argument, 0, 'M'},
        {"numEvts",     required_argument, 0, 'E'}
    };

    while((opt = getopt_long(argc, argv, "cH:D:N:M:E:", long_options, &option_index)) != -1 )
    {
        switch(opt)
        {
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
        sprintf(thistFile, "MySkimAnalysis_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
    }

    std::set<AnaSamples::FileSummary> vvf = setFS(dataSets, runOnCondor); 
    TFile* outfile = TFile::Open(histFile.c_str(), "RECREATE");

    std::cout << "Initializing..." << std::endl;
    const int jecOn = 0, jerOn = 0; 
    
    if( jecOn * jerOn != 0 || std::fabs(jecOn) > 1 || std::fabs(jerOn) > 1) 
    {
        std::cerr<<color("Error: ", "red")
                 <<"Invalid values of jecOn and jerOn. "
                 <<"They should be either -1, 0, or 1. "
                 <<"you cannot have both jet energy corrections and jet energy resolutions on at the same time"<<std::endl;
        return 1;
    }

    std::string                 myVarSuffix = "";
    if      ( jecOn == 1 )      myVarSuffix = "JECup";
    else if ( jecOn == -1 )     myVarSuffix = "JECdown";
    else if ( jerOn == 1 )      myVarSuffix = "JERup";
    else if ( jerOn == -1 )     myVarSuffix = "JERdown";
    
    if( vvf.size() < 1 ) {
        std::cerr<<"No such file exists"<<std::endl;
        return 0;
    }
    
    const auto& file =  (*vvf.begin());
    // Define what is needed per sample set
    std::cout << "Running over sample " << file.tag << std::endl;
    TChain* ch = new TChain( (file.treePath).c_str() );
    file.addFilesToChain(ch, startFile, nFiles);
    TTree* mySkimTree = ch->CloneTree(0);
    
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
        tr.registerFunction( std::move(rtt) );
    }
    
    Muon muon;
    Electron electron;
    Jet jet(myVarSuffix);
    BJet bjet(myVarSuffix);
    Photon photon;
    RunFisher runFisher("v3",myVarSuffix);
    CommonVariables commonVariables;
    MakeMVAVariables makeMVAVariables(false, myVarSuffix);
    Baseline baseline;
    DeepEventShape deepEventShape;
    
    // Register classes/functions that add variables on the fly
    tr.registerFunction( std::move(muon) );
    tr.registerFunction( std::move(electron) );
    tr.registerFunction( std::move(jet) );
    tr.registerFunction( std::move(bjet) );
    tr.registerFunction( std::move(photon) );
    tr.registerFunction( std::move(runFisher) );
    tr.registerFunction( std::move(commonVariables) );
    tr.registerFunction( std::move(makeMVAVariables) );
    tr.registerFunction( std::move(baseline) );
    tr.registerFunction( std::move(deepEventShape) );
    
    if( runtype == "MC" ) 
    {
        //std::string eosPath =   "root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/";
        BTagCorrectorTemplate<double> bTagCorrector("allInOne_BTagEff.root","", false, file.tag);
        Pileup_SysTemplate<double> pileup("PileupHistograms_0121_69p2mb_pm4p6.root");
        ScaleFactors scaleFactors("allInOne_leptonSF_Moriond17.root");
        tr.registerFunction( std::move(bTagCorrector) );
        tr.registerFunction( std::move(pileup) );
        tr.registerFunction( std::move(scaleFactors) );
    }

    double                  deepESM_val;
    TBranch *b_deepESM_val = mySkimTree->Branch( "deepESM_val", &deepESM_val, "deepESM_val/D");
    

    while( tr.getNextEvent() ) 
    {
        const auto& NGoodJets_pt30  = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodLeptons    = tr.getVar<int>("NGoodLeptons");
        const auto& passTrigger     = tr.getVar<bool>("passTrigger");
        const auto& deepESMValue    = tr.getVar<double>("deepESM_val");

        if( maxEvts != -1 && tr.getEvtNum() >= maxEvts ) break;

        if( runtype == "MC" ) {
            const auto& passMadHT   = tr.getVar<bool>("passMadHT");
            if( !passMadHT ) continue;
        }

        deepESM_val = deepESMValue;

        if( NGoodJets_pt30 > 4 && NGoodLeptons > 0 && passTrigger ) mySkimTree->Fill();
    }
    
    outfile->cd();
    mySkimTree->Write();        
    outfile->Write();
    std::cout << "Writing histograms..." << std::endl;
    outfile->Close();
    // Cleaning up dynamic memory
    delete ch;
   
    return 0;
}
