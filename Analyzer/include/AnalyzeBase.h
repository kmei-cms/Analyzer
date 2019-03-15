#ifndef AnalyzerBase_h
#define AnalyzerBase_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TFile.h>
#include <map>
#include <string>

class NTupleReader;

class AnalyzeBase
{
protected:
    std::vector<std::unique_ptr<Histo_Base>> my_Histos;

public:
    virtual void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false) = 0;
    virtual void InitHistos(std::string filetag) = 0;
    void WriteHistos(TFile* outfile)
    {
        outfile->cd();

        for(const auto& h : my_Histos)
        {
            h->Write();
        }
    }

    void printEventNum(const int maxEvents, const int evtNum, const int divisor = 1000)
    {
        if( maxevents != -1 && evtNum >= maxevents ) break;
        if( tr.getEvtNum() % divisor == 0 ) printf( " Event %i\n", evtNum );
    }

    void Fill(NTupleReader& tr)
    {
        for(auto& h : my_Histos)
        {
            h->Fill(tr);
        }
    }
};

#endif
