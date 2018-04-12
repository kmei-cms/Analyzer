#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/SATException.h"

//Hack--Joes update to samples.cc/h should fix this
#include "Framework/Framework/include/samples.h"
//#include "SusyAnaTools/Tools/samples.h"

#include "TopTaggerTools/Tools/include/HistoContainer.h"
#include "TopTagger/CfgParser/include/TTException.h"

#include "TopTagger/TopTagger/include/TopTaggerResults.h"

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

class AnalyzeStealthTopTagger
{
private:

public:
    void Loop(NTupleReader& tr, double weight, int maxevents, std::string filetag, bool savefile = true)
    {
        
        int events = 0, peventsLep0 = 0, peventsLep1 = 0;
        HistoContainer<NTupleReader> hists0Lep("Lep0");
        TRandom* trand = new TRandom3();        
        std::string filename = "example.root";

        while(tr.getNextEvent())
        {
            events++;

            if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
            if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

            const auto& MET          = tr.getVar<double>("MET");
            const auto& HT           = tr.getVar<double>("HT");
            const auto& ntops        = tr.getVar<int>("ntops");
            const auto& ntops_3jet   = tr.getVar<int>("ntops_3jet");
            const auto& ntops_2jet   = tr.getVar<int>("ntops_2jet");
            const auto& ntops_1jet   = tr.getVar<int>("ntops_1jet");
            const auto& runtype      = tr.getVar<std::string>("runtype");     
            const auto& TriggerNames = tr.getVec<std::string>("TriggerNames");
            const auto& TriggerPass  = tr.getVec<int>("TriggerPass");
            const auto& NJets_pt30   = tr.getVar<int>("NJets_pt30");
            const auto& NJets_pt45   = tr.getVar<int>("NJets_pt45");
            const auto& NBJets       = tr.getVar<int>("NBJets");
            const auto& NBJets_pt45  = tr.getVar<int>("NBJets_pt45");
            const auto& NGoodLeptons = tr.getVar<int>("NGoodLeptons");
            const auto& HT_trigger   = tr.getVar<double>("HT_trigger");

            double eventweight = 1.0;               
            bool isData = false;
            if(runtype == "MC") 
            {
                const auto& madHT   = tr.getVar<double>("madHT");
                const auto& Weight  = tr.getVar<double>("Weight");
                double lumi = 35900; // Lumi for 2016
                    
                // Exclude events with MadGraph HT > 100 from the DY inclusive sample
                if(filetag == "DYJetsToLL_M-50_Incl" && madHT > 100) continue;
                    
                // Weight from NTuples            
                eventweight = lumi*Weight;                    
            }
            else if (runtype == "Data")
            {
                isData = true;
            }

            //Stealth Event Selection - 0 Lepton
            if( !isData
                && NGoodLeptons==0
                && NJets_pt45>=6 
                && HT_trigger > 500
                && NBJets_pt45 >= 2
                )
            {
                peventsLep0++;
                hists0Lep.fill(tr, eventweight, trand);
            }
        }

        std::cout << "Processed " << events << " events. " << peventsLep0 << " passed Lep0 selection. " << peventsLep1 << " passed Lep1 selection"<<std::endl;

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
        
            f->Write();
            f->Close();
        }
    }

    void WriteHistos()
    {            
        std::cout<<"Fix Me Please"<<std::endl;
    }

};
