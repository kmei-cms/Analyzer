#include "Analyzer/Analyzer/include/ExploreBackground.h"
#include "Analyzer/Analyzer/include/ExploreTopTagger.h"
#include "Analyzer/Analyzer/include/ExploreEventSelection.h"
#include "Framework/Framework/include/samples.h"
#include "TopTaggerTools/Tools/include/HistoContainer.h"

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"

#include <iostream>
#include <getopt.h>

template<typename Explore> void run(TChain* ch, std::set<AnaSamples::FileSummary> vvf, 
                                    std::string runType, int startFile, int nFiles, int maxEvts)
{
    Explore t = Explore(ch);
    std::cout << "Initializing..." << std::endl;
    t.InitHistos();
    for(const AnaSamples::FileSummary& file : vvf)
    {
        std::cout << "Running over sample " << file.tag << std::endl;
        TChain* new_ch = new TChain( (AnaSamples::treeName).c_str());
        t.Init(new_ch);
        file.addFilesToChain(new_ch, startFile, nFiles);
        double weight = file.getWeight();
        std::string runtype = "";
        if(file.tag.find(runType) != std::string::npos)
            runtype = runType;
        std::cout << "Starting loop" << std::endl;
        printf( "weight: %f nFiles: %i startFile: %i maxEvts: %i \n",weight,nFiles,startFile,maxEvts );
        t.Loop(weight, maxEvts, runtype, file.tag);
    }
    std::cout << "Writing histograms..." << std::endl;
    t.WriteHistos();
}

int main(int argc, char *argv[])
{

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"doBackground",       no_argument, 0, 'b'},
        {"doTopTagger",        no_argument, 0, 't'},
        {"doEventSelection",   no_argument, 0, 's'},
        {"condor",             no_argument, 0, 'c'},
        {"histFile",     required_argument, 0, 'H'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
    };

    bool doBackground = false, doTopTagger = false, doEventSelection = false;
    bool runOnCondor = false;
    std::string histFile = "", dataSets = "", sampleloc = AnaSamples::fileDir;
    int nFiles = -1, startFile = 0, maxEvts = -1;

    while((opt = getopt_long(argc, argv, "btscH:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'b':
            doBackground = true;
            break;

        case 't':
            doTopTagger = true;
            break;

        case 's':
            doEventSelection = true;
            break;

        case 'c':
            runOnCondor = true;
            break;

        case 'H':
            histFile = optarg;
            break;

        case 'D':
            dataSets = optarg;
            break;

        case 'N':
            nFiles = int(atoi(optarg));
            break;

        case 'M':
            startFile = int(atoi(optarg));
            break;

        case 'E':
            maxEvts = int(atoi(optarg));
            break;
        }
    }

    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "MyAnalysis_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
        sampleloc = "condor";
    }

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

    TFile* myfile = TFile::Open(histFile.c_str(), "RECREATE");
    TChain* ch = new TChain( (AnaSamples::treeName).c_str() ) ;

    if(doBackground)
    {
        run<ExploreBackground>(ch,vvf,"Data",startFile,nFiles,maxEvts);
    }
    else if(doTopTagger)
    {
        run<ExploreTopTagger>(ch,vvf,"qcd",startFile,nFiles,maxEvts);
    }
    else if(doEventSelection)
    {
        run<ExploreEventSelection>(ch,vvf,"Data",startFile,nFiles,maxEvts);
    }

    myfile->Close();

    return 0;
}
