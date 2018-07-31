#ifndef CalculateBTagSF_h
#define CalculateBTagSF_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <map>
#include <string>

class NTupleReader;

class CalculateBTagSF{

public :
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
   std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;

   CalculateBTagSF();
   ~CalculateBTagSF(){};

   void     Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   virtual void     InitHistos(NTupleReader &tr);
   virtual void     WriteHistos(TFile* outfile);
   bool     JetPassCuts( const TLorentzVector& jet);
};

#endif
