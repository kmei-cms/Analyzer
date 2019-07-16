#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
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
#include "Analyzer/Analyzer/include/AnalyzeTest.h"
#include "Analyzer/Analyzer/include/AnalyzeBTagSF.h"
#include "Analyzer/Analyzer/include/MakeNJetDists.h"
#include "Analyzer/Analyzer/include/MakeMiniTree.h"
#include "Analyzer/Analyzer/include/CalculateBTagSF.h"
#include "Analyzer/Analyzer/include/CalculateSFMean.h"
#include "Analyzer/Analyzer/include/Config.h"
#include "Analyzer/Analyzer/include/Semra_Analyzer.h" // semra

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"

#include <iostream>
#include <getopt.h>
#include <string>
#include <functional>
#include <unistd.h>

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

const std::string getFullPath(const std::string& file)
{
    char buf[512];
    int count = readlink(file.c_str(), buf, sizeof(buf));
    if(count >= 0)
    {
        buf[count] = '\0';
        return std::string(buf);
    }
    else
    {
        std::cout<<"Could not get full path of "<<file<<std::endl;
        return std::string();
    }
}

template<typename Analyze> void run(const std::set<AnaSamples::FileSummary>& vvf, 
                                    const int startFile, const int nFiles, const int maxEvts, 
                                    TFile* const outfile, const bool isQuiet, const std::string& analyzer)
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
        const std::string runYear;
        if (file.tag.find("2016") != std::string::npos) {
            runYear = "2016";
        } else if (file.tag.find("2017") != std::string::npos) {
            runYear = "2017";
        } else if (file.tag.find("2018") != std::string::npos) {
            runYear = "2018";
        }
        const double Lumi;
        if (runYear == "2016") {
            Lumi = 35900.0;
        } else if (runYear == "2017") {
            Lumi = 41525.0;
        } else if (runYear == "2018") {
            Lumi = 59740.0;
        }

        const bool isSignal = (file.tag.find("_stop") != std::string::npos || file.tag.find("_mStop") != std::string::npos) ? true : false;
        const std::string DeepESMCfg = (runYear == "2016") ? "DeepEventShape_2016.cfg" : "DeepEventShape_2017.cfg";
        const std::string ModelFile = (runYear == "2016") ? "keras_frozen_2016.pb" : "keras_frozen_2017.pb";
        std::cout << "Starting loop (in run)" << std::endl;
        printf( "runtype: %s fileWeight: %f nFiles: %i startFile: %i maxEvts: %i \n",runtype.c_str(),weight,nFiles,startFile,maxEvts ); fflush( stdout );
        tr.registerDerivedVar("runtype",runtype);
        tr.registerDerivedVar("runYear",runYear);
        tr.registerDerivedVar("filetag",file.tag);
        tr.registerDerivedVar("etaCut",2.4); 
        tr.registerDerivedVar("Lumi",Lumi);
        tr.registerDerivedVar("isSignal",isSignal);
        tr.registerDerivedVar("DeepESMCfg",DeepESMCfg);
        tr.registerDerivedVar("ModelFile",ModelFile);        
        tr.registerDerivedVar("blind",false);
        tr.registerDerivedVar("analyzer",analyzer);

        // Define classes/functions that add variables on the fly        
        Config c;
        c.registerModules(tr);

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
    bool runOnCondor = false, isQuiet = true;
    std::string histFile = "", dataSets = "", analyzer = "";
    int nFiles = -1, startFile = 0, maxEvts = -1;

    static struct option long_options[] = {
        {"condor",             no_argument, 0, 'c'},
        {"verbose",            no_argument, 0, 'v'},
        {"analyzer",     required_argument, 0, 'A'},
        {"histFile",     required_argument, 0, 'H'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
    };

/// here is the options to run the codes / can add options
    while((opt = getopt_long(argc, argv, "cvA:H:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 'c': runOnCondor       = true;              break;
            case 'v': isQuiet           = false;             break;
            case 'A': analyzer          = optarg;            break;
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

    std::vector<std::pair<std::string, std::function<void(const std::set<AnaSamples::FileSummary>&,const int,const int,const int,TFile* const,const bool,const std::string&)>>> AnalyzerPairVec = {
        {"AnalyzeBackground",       run<AnalyzeBackground>},
        {"AnalyzeWControlRegion",   run<AnalyzeWControlRegion>},
        {"AnalyzeTopTagger",        run<AnalyzeTopTagger>},
        {"AnalyzeEventSelection",   run<AnalyzeEventSelection>},
        {"AnalyzeEventShape",       run<AnalyzeEventShape>},
        {"Analyze0Lep",             run<Analyze0Lep>},
        {"Analyze1Lep",             run<Analyze1Lep>},
        {"AnalyzeStealthTopTagger", run<AnalyzeStealthTopTagger>},
        {"AnalyzeBTagSF",           run<AnalyzeBTagSF>},
        {"AnalyzeTest",             run<AnalyzeTest>},
        {"CalculateBTagSF",         run<CalculateBTagSF>},
        {"CalculateSFMean",         run<CalculateSFMean>},
        {"MakeMiniTree",            run<MakeMiniTree>},
        {"MakeNJetDists",           run<MakeNJetDists>},
        {"AnalyzeNjetsMinusOneCSFillDijetHists", run<AnalyzeNjetsMinusOneCSFillDijetHists>},
        {"AnalyzeNjetsMinusOneCSJetReplacement", run<AnalyzeNjetsMinusOneCSJetReplacement>},
	{"Semra_Analyzer",          run<Semra_Analyzer>}, // SEMRA / to run my analyzer
    }; 

    try
    {
        for(auto& pair : AnalyzerPairVec)
        {
            if(pair.first==analyzer) pair.second(vvf,startFile,nFiles,maxEvts,outfile,isQuiet,analyzer); 
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
