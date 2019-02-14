#ifndef CalculateSFMean_h
#define CalculateSFMean_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TFile.h>

#include <map>
#include <string>

#include "Analyzer/Analyzer/include/Histo.h"
class NTupleReader;

class CalculateSFMean
{
public:
    std::vector<std::unique_ptr<Histo_Base>> my_Histos;
    
    CalculateSFMean();
    ~CalculateSFMean(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(std::string filetag);
    void WriteHistos(TFile* outfile); 
};

#endif
