#include "Analyzer/Analyzer/include/AnalyzeBackground.h"
#include "Analyzer/Analyzer/include/AnalyzeTopTagger.h"
#include "Analyzer/Analyzer/include/AnalyzeEventSelection.h"
#include "Analyzer/Analyzer/include/AnalyzeEventShape.h"
#include "Analyzer/Analyzer/include/Analyze0Lep.h"
#include "Framework/Framework/include/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/RunFisher.h"
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"

#include <iostream>
#include <getopt.h>

template<typename Analyze> void run(std::set<AnaSamples::FileSummary> vvf, 
                                    int startFile, int nFiles, int maxEvts, bool isSkim)
{
    std::cout << "Initializing..." << std::endl;
    Analyze t = Analyze();
    for(const AnaSamples::FileSummary& file : vvf)
    {
        // Define what is needed per sample set
        std::cout << "Running over sample " << file.tag << std::endl;
        TChain* ch = new TChain( (file.treePath).c_str() );
        file.addFilesToChain(ch, startFile, nFiles);
        NTupleReader tr(ch);
        double weight = file.getWeight();
        std::string runtype = "";
        if(file.tag.find("Data") != std::string::npos)
        {
            runtype = "Data";
        }
        else
        {
            runtype = "MC";
        }
        std::cout << "Starting loop (in run)" << std::endl;
        printf( "runtype: %s weight: %f nFiles: %i startFile: %i maxEvts: %i \n",runtype.c_str(),weight,nFiles,startFile,maxEvts ); fflush( stdout );
        tr.registerDerivedVar<std::string>("runtype",runtype);
        tr.registerDerivedVar<double>("etaCut",2.4);

        // Define classes/functions that add variables on the fly
        std::shared_ptr<RunTopTagger> rtt;
        if ( !isSkim ) rtt = std::make_shared<RunTopTagger>();
        RunFisher runFisher;
        Muon muon;
        Electron electron;
        Jet jet;
        BJet bjet;
        CommonVariables commonVariables;

        // Register classes/functions that add variables on the fly
        if ( !isSkim ) tr.registerFunction( std::move(*rtt) );
        tr.registerFunction( std::move(jet) );
        tr.registerFunction( std::move(runFisher) );
        tr.registerFunction( std::move(muon) );
        tr.registerFunction( std::move(electron) );
        tr.registerFunction( std::move(bjet) );
        tr.registerFunction( std::move(commonVariables) );

        // Loop over all of the events and fill histos
        t.Loop(tr, weight, maxEvts, file.tag);

        // Cleaning up dynamic memory
        delete ch;
    }
    std::cout << "Writing histograms..." << std::endl;
    t.WriteHistos();
}

std::set<AnaSamples::FileSummary> setFS(std::string sampleloc, std::string dataSets)
{
    AnaSamples::SampleSet        ss(sampleloc);
    AnaSamples::SampleCollection sc(ss);

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
    bool doBackground = false, doTopTagger = false, doEventSelection = false, doEventShape = false, do0Lep = false;
    bool runOnCondor = false;
    bool isSkim = false;
    std::string histFile = "", dataSets = "", sampleloc = AnaSamples::fileDir;
    int nFiles = -1, startFile = 0, maxEvts = -1;

    static struct option long_options[] = {
        {"doBackground",       no_argument, 0, 'b'},
        {"doTopTagger",        no_argument, 0, 't'},
        {"doEventSelection",   no_argument, 0, 's'},
        {"doEventShape",       no_argument, 0, 'p'},
        {"do0Lep",             no_argument, 0, 'z'},
        {"condor",             no_argument, 0, 'c'},
        {"histFile",     required_argument, 0, 'H'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
    };

    while((opt = getopt_long(argc, argv, "btspzcH:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 'b': doBackground     = true;              break;
            case 't': doTopTagger      = true;              break;
            case 's': doEventSelection = true;              break;
            case 'p': doEventShape     = true;              break;
            case 'z': do0Lep           = true;              break;
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
        sprintf(thistFile, "MyAnalysis_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
        sampleloc = "condor";
    }

    std::set<AnaSamples::FileSummary> vvf = setFS(sampleloc, dataSets); 
    TFile* myfile = TFile::Open(histFile.c_str(), "RECREATE");

    if(doBackground)
    {
        printf("\n\n running AnalyzeBackground\n\n" ) ;
        run<AnalyzeBackground>(vvf,startFile,nFiles,maxEvts,isSkim);
    }
    else if(doTopTagger)
    {
        printf("\n\n running AnalyzeTopTagger\n\n" ) ;
        run<AnalyzeTopTagger>(vvf,startFile,nFiles,maxEvts,isSkim);
    }
    else if(doEventSelection)
    {
        printf("\n\n running AnalyzeEventSelection\n\n") ;
        run<AnalyzeEventSelection>(vvf,startFile,nFiles,maxEvts,isSkim);
    }
    else if(doEventShape)
    {
        printf("\n\n running AnalyzeEventShape\n\n") ;
        run<AnalyzeEventShape>(vvf,startFile,nFiles,maxEvts,isSkim);
    }
    else if(do0Lep)
    {
        printf("\n\n running Analyze0Lep\n\n") ;
        run<Analyze0Lep>(vvf,startFile,nFiles,maxEvts,isSkim);
    }

    myfile->Close();

    return 0;
}
