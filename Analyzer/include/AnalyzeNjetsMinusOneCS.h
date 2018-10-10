#ifndef AnalyzeNjetsMinusOneCS_h
#define AnalyzeNjetsMinusOneCS_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include "Framework/Framework/include/DeepEventShape.h"
#include <map>
#include <string>

class NTupleReader;

class AnalyzeNjetsMinusOneCS 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TProfile>>  my_tp_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    bool initHistos;
    
    AnalyzeNjetsMinusOneCS();
    ~AnalyzeNjetsMinusOneCS(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(const std::map<std::string, bool>& cutMap);
    void WriteHistos(TFile* outfile);
private:
    DeepEventShape* des_ ;
    
};

#endif
