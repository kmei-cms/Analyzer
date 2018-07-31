#ifndef AnalyzeHadTrigger_h
#define AnalyzeHadTrigger_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeHadTrigger{

public :
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
   std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;

   AnalyzeHadTrigger();
   ~AnalyzeHadTrigger(){};

   void     Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   virtual void     InitHistos(NTupleReader &tr);
   virtual void     WriteHistos(TFile* outfile);
   bool     passTriggerGeneral(std::vector<std::string> &mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass );
   bool     passLooseSelection ( TLorentzVector &lv );
   bool     passTightSelection ( TLorentzVector &lv );
   void     printTriggerList   ( const std::vector<std::string>& TriggerNames );
   void     doesTriggerExist   ( std::vector<std::string>& myTrigTest, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass );

};

#endif
