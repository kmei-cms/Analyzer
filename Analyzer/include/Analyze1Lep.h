#ifndef Analyze1Lep_h
#define Analyze1Lep_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

#include "SusyAnaTools/Tools/NTupleReader.h"

class Analyze1Lep 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TProfile>>  my_tp_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    bool initHistos;
    
    Analyze1Lep();
    ~Analyze1Lep(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(const std::map<std::string, bool>& cutMap);    
    void WriteHistos(TFile* outfile);

    template<typename T> T getnJetDeepESMVar(const NTupleReader& tr, const std::string var, const std::vector<int> nJetTypes = {}, const int& nJets = -1, const bool& doNJet = false)
    {
        std::string name = "deepESM_"+var;
        T v = tr.getVar<T>(name);
        if(doNJet)        
        {
            int maxNJetTraining = *std::max_element(nJetTypes.begin(), nJetTypes.end());
            for(const auto& j : nJetTypes)
            {
                if(nJets == j || (nJets >= maxNJetTraining && j == maxNJetTraining))
                {
                    name = "deepESM_nJet"+std::to_string(j)+"_"+var;
                    v = tr.getVar<T>(name);
                }
            }        
        }
        return v;
    }    
};

#endif
