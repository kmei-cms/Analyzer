#include "Analyzer/Analyzer/include/StealthHemispheres.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "Framework/Framework/include/SetUpTopTagger.h"


#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"

StealthHemispheres::StealthHemispheres() : inithisto(false)
{
}

void StealthHemispheres::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {

        const auto& eventCounter    = tr.getVar<int>("eventCounter");
        
        // ------------------------
        // -- Print event number   
        // ------------------------     
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype         = tr.getVar<std::string>("runtype");
        const auto& filetag         = tr.getVar<std::string>("filetag");
        const auto& JetID           = tr.getVar<bool>("JetID");
        const auto& NGoodLeptons    = tr.getVar<int>("NGoodLeptons");
        const auto& passTriggerMC   = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodBJets_pt45 = tr.getVar<int>("NGoodBJets_pt45");
        const auto& HT_trigger_pt45 = tr.getVar<double>("HT_trigger_pt45");
        const auto& NGoodJets_pt45  = tr.getVar<int>("NGoodJets_pt45");
        const auto& passMadHT       = tr.getVar<bool>("passMadHT");
        const auto& MET             = tr.getVar<double>("MET");
        const auto& Jets            = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt45   = tr.getVec<bool>("GoodJets_pt45");
        const auto& GoodBJets_pt45  = tr.getVec<bool>("GoodBJets_pt45");
        const auto& stopMass        = tr.getVar<double>("stopMass");
        std::cout<<stopMass<<std::endl;

        // -------------------
        // -- Define weight
        // -------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double bTagScaleFactor      = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight = tr.getVar<double>("Weight");
            const auto& lumi   = tr.getVar<double>("Lumi");
            eventweight        = lumi*Weight;

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");

            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }
    } 
} 

void StealthHemispheres::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }

    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }

    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }
}
