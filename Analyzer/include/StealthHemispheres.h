#ifndef StealthHemispheres_H
#define StealthHemispheres_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include "TH2D.h"

#include <map>

class NTupleReader;

class StealthHemispheres 
{
private:
    bool inithisto;
    
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;

    int          seed_method;
    int          assoc_method;
    double       MT2;
    double       MCT;
    double       AlphaT;
    double       minDHT;
    double       maxDR;
    double       dPhi;

    std::vector<int> jindices1;
    std::vector<int> jindices2;
    std::vector<int> eleindices1;
    std::vector<int> eleindices2;
    std::vector<int> muoindices1;
    std::vector<int> muoindices2;
    TLorentzVector lv1;
    TLorentzVector lv2;
    TLorentzVector UTM;

    Double_t GetMT2Hemi(double testmass=0, bool massive=false, int PFJID=0, double minJPt=20, double maxJEta=2.4, int hemi_association=3, int met=1);
    Double_t CalcMT2(double testmass, bool massive, TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET); 
    Double_t HemiMassTop();

    void Reset();
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void WriteHistos(TFile* outfile);
    StealthHemispheres();

};

#endif
