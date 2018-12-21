#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/BTagCalibrationStandalone.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "SusyAnaTools/Tools/PileupWeights.h"
#include "SusyAnaTools/Tools/MiniTupleMaker.h"

#include "TopTagger/CfgParser/include/TTException.h"

#include "Analyzer/Analyzer/include/AnalyzeBackground.h"
#include "Analyzer/Analyzer/include/AnalyzeWControlRegion.h"
#include "Analyzer/Analyzer/include/AnalyzeTopTagger.h"
#include "Analyzer/Analyzer/include/AnalyzeEventSelection.h"
#include "Analyzer/Analyzer/include/AnalyzeEventShape.h"
#include "Analyzer/Analyzer/include/Analyze0Lep.h"
#include "Analyzer/Analyzer/include/Analyze1Lep.h"
#include "Analyzer/Analyzer/include/AnalyzeNjetsMinusOneCSFillDijetHists.h"
#include "Analyzer/Analyzer/include/AnalyzeNjetsMinusOneCSJetReplacement.h"
#include "Analyzer/Analyzer/include/AnalyzeStealthTopTagger.h"
#include "Analyzer/Analyzer/include/AnalyzeBTagSF.h"
#include "Analyzer/Analyzer/include/MakeNJetDists.h"
#include "Analyzer/Analyzer/include/MakeMiniTree.h"
#include "Analyzer/Analyzer/include/CalculateBTagSF.h"
#include "Analyzer/Analyzer/include/CalculateHtSF.h"

#include "Framework/Framework/include/PrepNTupleVars.h"
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
                                    TFile* outfile, bool isQuiet)
{
    std::cout << "Initializing..." << std::endl;
    Analyze a;
    for(const AnaSamples::FileSummary& file : vvf)
    {
        // Define what is needed per sample set
        std::cout << "Running over sample " << file.tag << std::endl;
        TChain* ch = new TChain( (file.treePath).c_str() );
        file.addFilesToChain(ch, startFile, nFiles);
        NTupleReader tr(ch);
        double weight = file.getWeight(); // not used currently
        const std::string runtype = (file.tag.find("Data") != std::string::npos) ? "Data" : "MC";
        const std::string runYear = (file.tag.find("2017") != std::string::npos) ? "2017" : "2016";
        const double Lumi = (runYear == "2016") ? 35900.0 : 41525.0;
        std::cout << "Starting loop (in run)" << std::endl;
        printf( "runtype: %s fileWeight: %f nFiles: %i startFile: %i maxEvts: %i \n",runtype.c_str(),weight,nFiles,startFile,maxEvts ); fflush( stdout );
        tr.registerDerivedVar("runtype",runtype);
        tr.registerDerivedVar("runYear",runYear);
        tr.registerDerivedVar("filetag",file.tag);
        tr.registerDerivedVar("etaCut",2.4); 
        tr.registerDerivedVar("Lumi",Lumi); 
        tr.registerDerivedVar("blind",true);

        // Define classes/functions that add variables on the fly
        PrepNTupleVars prep;
        RunTopTagger rtt;
        Muon muon;
        Electron electron;
        Photon photon;
        Jet jet;
        BJet bjet;
        RunFisher runFisher("v3");
        CommonVariables commonVariables;
        MakeMVAVariables makeMVAVariables;
        Baseline baseline;
        DeepEventShape deepEventShape;
        
        // Register classes/functions that add variables on the fly
        tr.registerFunction(prep);
        tr.registerFunction(rtt);
        tr.registerFunction(muon);
        tr.registerFunction(electron);
        tr.registerFunction(photon);
        tr.registerFunction(jet);
        tr.registerFunction(bjet);
        tr.registerFunction(runFisher);
        tr.registerFunction(commonVariables);
        tr.registerFunction(makeMVAVariables);
        tr.registerFunction(baseline);
        tr.registerFunction(deepEventShape);
        
        if( runtype == "MC" ) 
        {
            //std::string eosPath =   "root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/";
            BTagCorrectorTemplate<double> bTagCorrector("allInOne_BTagEff.root","", false, file.tag);
            bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets", "Jets_bDiscriminatorCSV", "Jets_partonFlavor");
            Pileup_SysTemplate<double> pileup("PileupHistograms_0121_69p2mb_pm4p6.root");
            std::string scaleFactorHistoFileName = (file.tag.find("2017") != std::string::npos ) ? "allInOne_leptonSF_2017.root" : "allInOne_leptonSF_Moriond17.root";
            ScaleFactors scaleFactors( scaleFactorHistoFileName );

            tr.registerFunction(bTagCorrector);
            tr.registerFunction(pileup);
            tr.registerFunction(scaleFactors);
        }
        // Loop over all of the events and fill histos
        a.Loop(tr, weight, maxEvts, isQuiet);

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
    bool doBackground = false, doTopTagger = false, doEventSelection = false, 
         doEventShape = false, do0Lep = false, do1Lep = false, doStealthTT = false,
         doBTagSF = false, calcBTagSF = false, calcHtSF = false, doWControlRegion = false, 
         makeMiniTree = false, makeNJetDists = false,
         doNjetsMinusOneCSFillDijetHists = false, doNjetsMinusOneCSJetReplacement = false, isQuiet = true;

    bool runOnCondor = false;
    std::string histFile = "", dataSets = "";
    int nFiles = -1, startFile = 0, maxEvts = -1;

    static struct option long_options[] = {
        {"doBackground",       no_argument, 0, 'b'},
        {"doWControlRegion",   no_argument, 0, 'w'},
        {"doTopTagger",        no_argument, 0, 't'},
        {"doEventSelection",   no_argument, 0, 's'},
        {"doEventShape",       no_argument, 0, 'p'},
        {"do0Lep",             no_argument, 0, 'z'},
        {"do1Lep",             no_argument, 0, 'o'},
        {"doNjetsMinusOneCSFillDijetHists",  no_argument, 0, 'q'},
        {"doNjetsMinusOneCSJetReplacement",  no_argument, 0, 'r'},
        {"doStealthTT",        no_argument, 0, 'x'},
        {"calcBTagSF",         no_argument, 0, 'f'},
        {"calcHtSF",           no_argument, 0, 'F'},
        {"doBTagSF",           no_argument, 0, 'g'},
        {"makeMiniTree",       no_argument, 0, 'm'},
        {"makeNJetDists",      no_argument, 0, 'n'},
        {"condor",             no_argument, 0, 'c'},
        {"verbose",            no_argument, 0, 'v'},
        {"histFile",     required_argument, 0, 'H'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
    };

    while((opt = getopt_long(argc, argv, "bwtspzoqrxfFgmncvH:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 'b': doBackground      = true;              break;
            case 'w': doWControlRegion  = true;              break;
            case 't': doTopTagger       = true;              break;
            case 's': doEventSelection  = true;              break;
            case 'p': doEventShape      = true;              break;
            case 'z': do0Lep            = true;              break;
            case 'o': do1Lep            = true;              break;
            case 'q': doNjetsMinusOneCSFillDijetHists = true;              break;
            case 'r': doNjetsMinusOneCSJetReplacement = true;              break;
            case 'x': doStealthTT       = true;              break;
            case 'f': calcBTagSF        = true;              break;
            case 'F': calcHtSF          = true;              break;
            case 'g': doBTagSF          = true;              break;
            case 'm': makeMiniTree      = true;              break;
            case 'n': makeNJetDists     = true;              break;
            case 'c': runOnCondor       = true;              break;
            case 'v': isQuiet           = false;             break;
            case 'H': histFile          = optarg;            break;
            case 'D': dataSets          = optarg;            break;
            case 'N': nFiles            = int(atoi(optarg)); break;
            case 'M': startFile         = int(atoi(optarg)); break;
            case 'E': maxEvts           = int(atoi(optarg)); break;
        }
    }

    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "MyAnalysis_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
    }

    std::set<AnaSamples::FileSummary> vvf = setFS(dataSets, runOnCondor); 
    TFile* outfile = TFile::Open(histFile.c_str(), "RECREATE");

    std::vector<std::pair<bool, std::function<void(std::set<AnaSamples::FileSummary>,int,int,int,TFile*,bool)>>> AnalyzerPairVec = {
        {doBackground,      run<AnalyzeBackground>},
        {doWControlRegion,  run<AnalyzeWControlRegion>},
        {doTopTagger,       run<AnalyzeTopTagger>},
        {doEventSelection,  run<AnalyzeEventSelection>},
        {doEventShape,      run<AnalyzeEventShape>},
        {do0Lep,            run<Analyze0Lep>},
        {do1Lep,            run<Analyze1Lep>},
        {doNjetsMinusOneCSFillDijetHists, run<AnalyzeNjetsMinusOneCSFillDijetHists>},
        {doNjetsMinusOneCSJetReplacement, run<AnalyzeNjetsMinusOneCSJetReplacement>},
        {doStealthTT,       run<AnalyzeStealthTopTagger>},
        {doBTagSF,          run<AnalyzeBTagSF>},
        {calcBTagSF,        run<CalculateBTagSF>},
        {calcHtSF,          run<CalculateHtSF>},
        {makeMiniTree,      run<MakeMiniTree>},
        {makeNJetDists,     run<MakeNJetDists>},
    }; 

    try
    {
        for(auto& pair : AnalyzerPairVec)
        {
            if(pair.first) pair.second(vvf,startFile,nFiles,maxEvts,outfile,isQuiet); 
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
