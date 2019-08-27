#ifndef AnalyzeTopTagger_h
#define AnalyzeTopTagger_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>
#include "TopTaggerTools/Tools/include/HistoContainer.h"

class NTupleReader;

class AnalyzeTopTagger
{
private:
   HistoContainer<NTupleReader> hists;
public:
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
   std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    
   AnalyzeTopTagger();
   ~AnalyzeTopTagger(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void InitHistos();
   void WriteHistos(TFile* outfile);
};

#endif
