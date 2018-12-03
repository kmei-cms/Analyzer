#ifndef Histo_h
#define Histo_h

#include "SusyAnaTools/Tools/NTupleReader.h"

class Histo_Base
{
public:
    virtual void Fill(const NTupleReader&) = 0;
    virtual void Write() const = 0;
};

template<typename Type>
class Histo_FirstChild : public Histo_Base
{
protected:
    std::string name_;
    int nBinsX_;
    double lowX_;
    double highX_;
    std::string varX_;
    std::vector<std::string> cuts_;
    std::vector<std::string> weights_;
    std::unique_ptr<Type> histo_;
    
    void Write() const
    {
        histo_->Write();
    }

    bool passCuts(const NTupleReader& tr) const
    {
        bool pass = true;
        for(const auto& cutName : cuts_)
        {
            const auto& cut = tr.getVar<bool>(cutName);
            pass = pass && cut;
            if(!pass) break;
        }
        return pass;
    }

    double getWeight(const NTupleReader& tr) const
    {
        double weight = 1.0;
        for(const auto& wName : weights_)
        {
            const auto& w = tr.getVar<double>(wName);
            weight *= w;
        }        
        return weight;
    }
};

class Histo1D : public Histo_FirstChild<TH1>
{
private:
    void fillHisto(const std::string& varX, std::unique_ptr<TH1>& histo, const double weight, const NTupleReader& tr) const
    {
        std::string type;
        tr.getType(varX, type);
        
        if(type.find("double")       != std::string::npos) histo->Fill( tr.getVar<double>(varX),       weight );
        if(type.find("unsigned int") != std::string::npos) histo->Fill( tr.getVar<unsigned int>(varX), weight );
        if(type.find("int")          != std::string::npos) histo->Fill( tr.getVar<int>(varX),          weight );
        if(type.find("bool")         != std::string::npos) histo->Fill( tr.getVar<bool>(varX),         weight );
        if(type.find("float")        != std::string::npos) histo->Fill( tr.getVar<float>(varX),        weight );
        if(type.find("char")         != std::string::npos) histo->Fill( tr.getVar<char>(varX),         weight );
        if(type.find("short")        != std::string::npos) histo->Fill( tr.getVar<short>(varX),        weight );
        if(type.find("long")         != std::string::npos) histo->Fill( tr.getVar<long>(varX),         weight );        
    } 

public:
    void Fill(const NTupleReader& tr)
    {
        const bool pass = passCuts(tr);
        if(pass)
        {
            const double weight = getWeight(tr);
            fillHisto(varX_, histo_, weight, tr);
        }
    }

    Histo1D(const std::string& name, const int nBinsX, const double lowX, const double highX, const std::string& varX, const std::vector<std::string>& cuts = {}, const std::vector<std::string>& weights = {})
    {
        name_ = name;
        nBinsX_ = nBinsX;
        lowX_ = lowX;
        highX_ = highX;
        varX_ = varX;
        cuts_ = cuts;
        weights_ = weights;
        histo_ = std::make_unique<TH1D>( name_.c_str(), name_.c_str(), nBinsX_, lowX_, highX_ );
    }
};

template<typename TX, typename TY>
class Histo2D : public Histo_FirstChild<TH2>
{
private:
    int nBinsY_;
    double lowY_;
    double highY_;
    std::string varY_;

public:
    void Fill(const NTupleReader& tr)
    {
        const bool pass = passCuts(tr);
        if(pass)
        {
            TX variableX = tr.getVar<TX>(varX_);
            TY variableY = tr.getVar<TY>(varY_);
            const double weight = getWeight(tr);
            histo_->Fill(variableX, variableY, weight);
        }
    }

    Histo2D(const std::string& name, 
            const int nBinsX, const double lowX, const double highX, const std::string& varX, 
            const int nBinsY, const double lowY, const double highY, const std::string& varY, 
            const std::vector<std::string>& cuts = {}, const std::vector<std::string>& weights = {})
    {
        name_ = name;
        nBinsX_ = nBinsX;
        lowX_ = lowX;
        highX_ = highX;
        varX_ = varX;
        nBinsY_ = nBinsY;
        lowY_ = lowY;
        highY_ = highY;
        varY_ = varY;
        cuts_ = cuts;
        weights_ = weights;
        histo_ = std::make_unique<TH2D>( name_.c_str(), name_.c_str(), nBinsX_, lowX_, highX_, nBinsY_, lowY_, highY_ );
    }
};

#endif
