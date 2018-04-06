#ifndef AnalyzeEventShape_h
#define AnalyzeEventShape_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeEventShape 
{
public:
    std::map<std::string, TH1D*>  my_histos;
    std::map<std::string, TH2D*>  my_2d_histos;
    std::map<std::string, TEfficiency*>  my_efficiencies;
    
    AnalyzeEventShape(){ printf("\n\n In AnalyzeEventShape constructor.\n\n") ; fflush( stdout ) ; };
    ~AnalyzeEventShape(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, std::string filetag = "", bool isQuiet = false);
    void InitHistos();
    void WriteHistos();
    bool PassTriggerGeneral(std::vector<std::string> &mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass);
    bool PassTriggerAllHad(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass);
    bool PassTriggerMuon(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass);
    bool PassTriggerElectron(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass);
    
};

#endif
