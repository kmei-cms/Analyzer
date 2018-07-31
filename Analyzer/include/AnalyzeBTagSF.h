#ifndef AnalyzeBTagSF_h
#define AnalyzeBTagSF_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeBTagSF{

public :
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
   std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;

   AnalyzeBTagSF();
   ~AnalyzeBTagSF(){};

   void     Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   virtual void     InitHistos(NTupleReader &tr);
   virtual void     WriteHistos(TFile* outfile);
};

#endif
