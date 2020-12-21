#ifndef AnalyzeDoubleDisCo_h
#define AnalyzeDoubleDisCo_h

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeDoubleDisCo 
{
private:
    class TH1DInfo
    {
    public:
        std::string name;
        int nBins;
        double low;
        double high;        
    };

    class TH2DInfo
    {
    public:
        std::string name;
        int nBinsX;
        double lowX;
        double highX;        
        int nBinsY;
        double lowY;
        double highY;        
    };

public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    bool initHistos;
    
    AnalyzeDoubleDisCo();
    ~AnalyzeDoubleDisCo(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(const std::map<std::string, bool>& cutMap, const std::map<int, bool>& njetsMap, const std::map<std::string, bool>& ABCDmap, const std::vector<TH1DInfo>& histInfos, const std::vector<TH2DInfo>& hist2DInfos);
    void WriteHistos(TFile* outfile);
};

#endif
