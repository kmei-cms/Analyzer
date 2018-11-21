#ifndef MakeNJetDists_h
#define MakeNJetDists_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

//class NTupleReader;
#include "SusyAnaTools/Tools/NTupleReader.h"

class Histo 
{
private:
    std::string name_;
    int nBins_;
    double low_;
    double high_;
    std::string var_;
    std::vector<std::string> cuts_;
    std::vector<std::string> weights_;

    std::unique_ptr<TH1> histo_;

public:
    void Fill(const NTupleReader& tr)
    {
        bool pass = true;
        for(const auto& cutName : cuts_)
        {
            const auto& cut = tr.getVar<bool>(cutName);
            pass = pass && cut;
            if(!pass) break;
        }

        if(pass)
        {
            int variable = tr.getVar<int>(var_);
            double weight = 1.0;
            for(const auto& wName : weights_)
            {
                const auto& w = tr.getVar<double>(wName);
                weight *= w;
            }
            
            histo_->Fill(variable, weight);
        }
    }

    void Write() const
    {
        histo_->Write();
    }
    
    Histo(const std::string& name, const int nBins, const double low, const double high, const std::string& var, const std::vector<std::string>& cuts, const std::vector<std::string>& weights)
        : name_(name)
        , nBins_(nBins)
        , low_(low)
        , high_(high)
        , var_(var)
        , cuts_(cuts)
        , weights_(weights)
    {
        histo_ = std::make_unique<TH1D>( name_.c_str(), name_.c_str(), nBins_, low_, high_ );
    }
};

class MakeNJetDists
{
public:
   std::vector<Histo> my_Histos;
   std::vector<std::string> myVarSuffixs;

   MakeNJetDists();
   ~MakeNJetDists(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void InitHistos();
   void WriteHistos(TFile* outfile); 
};

#endif
