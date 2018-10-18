#ifndef Analyze1Lep_h
#define Analyze1Lep_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

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

};

#endif
