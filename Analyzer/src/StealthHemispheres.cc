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
    Reset();
}

void StealthHemispheres::Reset()
{
  seed_method    = -1;
  assoc_method   = -1;
  MT2            = -99999.99;
  MCT            = -99999.99;
  AlphaT         = -99999.99;
  minDHT         = -99999.99;
  maxDR          = -99999.99;
  dPhi           = -99999.99;
  //for(int i=0; i<m_jetSize; ++i){
  //  jindices1[i]  =-1;
  //  jindices2[i]  =-1;
  //}
  //for(int i=0; i<m_eleSize; ++i){
  //  eleindices1[i]=-1;
  //  eleindices2[i]=-1;
  //}
  //for(int i=0; i<m_muoSize; ++i){
  //  muoindices1[i]=-1;
  //  muoindices2[i]=-1;
  //}

  lv1. SetPxPyPzE(0, 0, 0, 0);
  lv2. SetPxPyPzE(0, 0, 0, 0);
  UTM. SetPxPyPzE(0, 0, 0, 0);
}

// ------------------------------
// -- Call GetMT2Hemi function 
// ------------------------------
Double_t StealthHemispheres::GetMT2Hemi(double testmass, bool massive, int PFJID, double minJPt, 
                                        double maxJEta, int hemi_association, int met)
{    
    //TLorentzVector MET(0., 0., 0., 0.);
    //if(met==1)      MET = pfmet[0];
    //else if(met==2) MET = MHTloose[0];
    //else if(met==3) MET = pfmet[0]; //plus OS dileptons
    //else if(met==4) MET.SetPxPyPzE(0,0,0,0);
    //else return -999;
    //
    //if( met ==3 )
    //{
    //    // adding OS dilepton LV to MET               
    //    double dilept_invmass= GetDiLeptonInvMass(0,1,0,5.0,true);
    //    if(dilept_invmass < 111 && dilept_invmass > 71)
    //    {
    //        for(int i=0; i<NEles; ++i)
    //        {
    //            MET = MET + ele[i].lv;
    //        }
    //        for(int i=0; i<NMuons; ++i)
    //        {
    //            MET = MET + muo[i].lv;
    //        }
    //    } else
    //    {return -111;}                
    //}
    //        
    //vector<float> px, py, pz, E;
    //for(int i=0; i<NJets; ++i)
    //{
    //    if(jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID) == false) continue;
    //    px.push_back(jet[i].lv.Px());
    //    py.push_back(jet[i].lv.Py());
    //    pz.push_back(jet[i].lv.Pz());
    //    E .push_back(jet[i].lv.E ());
    //}
    //
    //if (px.size()<2) return -999;
    //
    //// Get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)           
    //Hemisphere* hemi = new Hemisphere(px, py, pz, E, 2, hemi_association);
    //vector<int> grouping = hemi->getGrouping();
    //
    //TLorentzVector pseudojet1(0.,0.,0.,0.);
    //TLorentzVector pseudojet2(0.,0.,0.,0.);
    //
    //for(int i=0; i<px.size(); ++i)
    //{
    //    if(grouping[i]==1){
    //        pseudojet1.SetPx(pseudojet1.Px() + px[i]);
    //        pseudojet1.SetPy(pseudojet1.Py() + py[i]);
    //        pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
    //        pseudojet1.SetE( pseudojet1.E()  + E[i]);
    //    }else if(grouping[i] == 2)
    //    {
    //        pseudojet2.SetPx(pseudojet2.Px() + px[i]);
    //        pseudojet2.SetPy(pseudojet2.Py() + py[i]);
    //        pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
    //        pseudojet2.SetE( pseudojet2.E()  + E[i]);
    //    }
    //}
    //
    //delete hemi;
    //
    //if(met==4) {MET = -pseudojet1 - pseudojet2;}
    //if(MET.Pt()<30) {return -222;}
    //return CalcMT2(testmass, massive, pseudojet1, pseudojet2, MET);        
} 

// -------------------------------
// -- Call HemiMassTop function
// -------------------------------
Double_t StealthHemispheres::HemiMassTop()
{
    //Double_t M1 = hemi[1].lv1.M();
    //Double_t M2 = hemi[1].lv2.M();
    //if(fabs(M1-172) < fabs(M2-172)) return M1;
    //else                            return M2;
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
