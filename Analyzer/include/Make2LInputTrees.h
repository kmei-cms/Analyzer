#ifndef Make2LInputTrees_h
#define Make2LInputTrees_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class MiniTupleMaker;

class Make2LInputTrees{

public :
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
   std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;

   Make2LInputTrees();
   ~Make2LInputTrees(){};

   void     Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void     InitHistos();
   void     WriteHistos(TFile* outfile); 

   MiniTupleMaker *myMiniTuple;
   TTree          *myTree;
};

#endif
