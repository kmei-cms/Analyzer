#ifndef CalculateHtSF_h
#define CalculateHtSF_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

#include "Analyzer/Analyzer/include/Histo.h"
class NTupleReader;

class CalculateHtSF
{
public:
    std::vector<std::unique_ptr<Histo_Base>> my_Histos;
    
    CalculateHtSF();
    ~CalculateHtSF(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(std::string filetag);
    void WriteHistos(TFile* outfile); 
};

#endif
