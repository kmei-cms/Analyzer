#ifndef MakeNNVariables_h
#define MakeNNVariables_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;
class MiniTupleMaker;

class MakeNNVariables{

public :
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
   std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;

   MakeNNVariables();
   ~MakeNNVariables(){};

   void     Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void     InitHistos();
   void     WriteHistos(TFile* outfile); 

   MiniTupleMaker *myMiniTupleTrain_0l;
   MiniTupleMaker *myMiniTupleTrain_1l; 
   MiniTupleMaker *myMiniTupleTrain_2l;
   MiniTupleMaker *myMiniTupleTest_0l;
   MiniTupleMaker *myMiniTupleTest_1l;
   MiniTupleMaker *myMiniTupleTest_2l;
   MiniTupleMaker *myMiniTupleVal_0l;
   MiniTupleMaker *myMiniTupleVal_1l;
   MiniTupleMaker *myMiniTupleVal_2l;

   TTree          *myTreeTrain_0l;
   TTree          *myTreeTrain_1l;
   TTree          *myTreeTrain_2l;
   TTree          *myTreeTest_0l;
   TTree          *myTreeTest_1l;
   TTree          *myTreeTest_2l;
   TTree          *myTreeVal_0l;
   TTree          *myTreeVal_1l;
   TTree          *myTreeVal_2l;
};

#endif
