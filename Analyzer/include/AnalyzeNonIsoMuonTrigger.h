#ifndef AnalyzeNonIsoMuonTrigger_h
#define AnalyzeNonIsoMuonTrigger_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeNonIsoMuonTrigger 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    
    AnalyzeNonIsoMuonTrigger();
    ~AnalyzeNonIsoMuonTrigger(){};
    
    void Loop(NTupleReader& tr, double weight = 1.0, int maxevents = -1, bool isQuiet = true );
    void InitHistos();
    void WriteHistos(TFile* outfile);
    bool passTriggerGeneral( std::vector<std::string>& myTriggerVector, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass );
    void printTriggerList( const std::vector<std::string>& TriggerNames );
    void doesTriggerExist( std::vector<std::string>& myTrigTestVector, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass );
    
};

#endif
