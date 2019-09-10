#ifndef StealthHemispheres_H
#define StealthHemispheres_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class StealthHemispheres : public TObject {
public:
  StealthHemispheres();
  virtual ~StealthHemispheres();

  void Reset();
  Int_t          seed_method;
  Int_t          assoc_method;
  Double_t       MT2;
  Double_t       MCT;
  Double_t       AlphaT;
  Double_t       minDHT;
  Double_t       maxDR;
  Double_t       dPhi;

  Int_t          jindices1  [m_jetSize];
  Int_t          jindices2  [m_jetSize];
  Int_t          eleindices1[m_eleSize];
  Int_t          eleindices2[m_eleSize];
  Int_t          muoindices1[m_muoSize];
  Int_t          muoindices2[m_muoSize];
  TLorentzVector lv1;
  TLorentzVector lv2;
  TLorentzVector UTM;

  ClassDef(StealthHemispheres, 4)
};

class StealthHemispheres : public TObject {

    Double_t GetMT2Hemi(double testmass=0, bool massive=false, int PFJID=0, double minJPt=20, double maxJEta=2.4, int hemi_association=3, int met=1);
    Double_t CalcMT2(double testmass, bool massive, TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET); 
    Double_t HemiMassTop();

};

#endif
