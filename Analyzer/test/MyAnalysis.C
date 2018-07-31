#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "TopTagger/CfgParser/include/TTException.h"

#include "Analyzer/Analyzer/include/AnalyzeStealthTopTagger.h"
#include "Analyzer/Analyzer/include/AnalyzeBTagSF.h"
#include "Analyzer/Analyzer/include/AnalyzeHadTrigger.h"
#include "Analyzer/Analyzer/include/CalculateBTagSF.h"

#include "SusyAnaTools/Tools/BTagCalibrationStandalone.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "SusyAnaTools/Tools/PileupWeights.h"

#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/RunFisher.h"
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/ScaleFactors.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/Baseline.h"

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"

#include <iostream>
#include <getopt.h>
#include <string>

template<typename Analyze> void run(std::set<AnaSamples::FileSummary> vvf, 
                                    int startFile, int nFiles, int maxEvts, 
                                    bool isSkim, TFile* outfile)
{
    std::cout << "Initializing..." << std::endl;
    const int jecOn = 0; 
    const int jerOn = 1;
    
    if( jecOn * jerOn != 0 || std::fabs(jecOn) > 1 || std::fabs(jerOn) > 1) {
            std::cerr<<"Invalid values of jecOn and jerOn. They should be either -1, 0, or 1, and you cannot have both jet energy corrections and jet energy resolutions on at the same time"<<std::endl;
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
        std::shared_ptr<RunTopTagger> rtt;
        if ( !isSkim ) rtt = std::make_shared<RunTopTagger>();
        RunFisher runFisher("v3",myVarSuffix);
        //RunFisher runFisher("test");
        if( runtype == "MC" ) {
            BTagCorrector bTagCorrector("allInOne_BTagEff.root","", false, file.tag);
            Pileup_Sys pileup("PileupHistograms_0121_69p2mb_pm4p6.root");
            tr.registerFunction( std::move(bTagCorrector) );
            tr.registerFunction( std::move(pileup) );
        }
        Muon muon;
        Electron electron;
        Jet jet(myVarSuffix);
        BJet bjet(myVarSuffix);
        ScaleFactors scaleFactors;
        CommonVariables commonVariables;
        Baseline baseline;

        // Register classes/functions that add variables on the fly
        if ( !isSkim ) tr.registerFunction( std::move(*rtt) );
        tr.registerFunction( std::move(muon) );
        tr.registerFunction( std::move(electron) );
        tr.registerFunction( std::move(jet) );
        tr.registerFunction( std::move(runFisher) );
        tr.registerFunction( std::move(bjet) );
        tr.registerFunction( std::move(commonVariables) );
        tr.registerFunction( std::move(scaleFactors) );
        tr.registerFunction( std::move(baseline) );

        // Loop over all of the events and fill histos
        a.InitHistos(tr);
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
    bool doBTagSF = false; bool calcBTagSF = false; bool doTrigAnalysis = false;
    bool runOnCondor = false;
    bool isSkim = false;
    std::string histFile = "", dataSets = "";
    int nFiles = -1, startFile = 0, maxEvts = -1;

    static struct option long_options[] = {
        {"calcBTagSF",         no_argument, 0, 'f'},
        {"doBTagSF",           no_argument, 0, 'g'},
        {"doTrigAnalysis",     no_argument, 0, 'a'},
        {"condor",             no_argument, 0, 'c'},
        {"histFile",     required_argument, 0, 'H'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
    };

    while((opt = getopt_long(argc, argv, "fgacH:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 'f': calcBTagSF       = true;              break;
            case 'g': doBTagSF         = true;              break;
            case 'a': doTrigAnalysis   = true;              break;
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
        sprintf(thistFile, "MyHadTrigAnalysis_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
    }

    std::set<AnaSamples::FileSummary> vvf = setFS(dataSets, runOnCondor); 
    TFile* outfile = TFile::Open(histFile.c_str(), "RECREATE");
    
    try
    {
        if(doBTagSF)
        {
            printf("\n\n running AnalyzeBTagSF\n\n") ;
            run<AnalyzeBTagSF>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        if(doTrigAnalysis)
        {
            printf("\n\n running AnalyzeHadTrigger\n\n") ;
            run<AnalyzeHadTrigger>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(calcBTagSF)
        {
            printf("\n\n running CalculateBTagSF\n\n") ;
            run<CalculateBTagSF>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
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
