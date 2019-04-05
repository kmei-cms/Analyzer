#ifndef CalculateSFMean_h
#define CalculateSFMean_h

#include "Analyzer/Analyzer/include/AnalyzeBase.h"
#include "Analyzer/Analyzer/include/Histo.h"
class NTupleReader;

class CalculateSFMean : public AnalyzeBase
{
public:
    void InitHistos(std::string filetag)
    {
        TH1::SetDefaultSumw2();
        TH2::SetDefaultSumw2();

        my_Histos.emplace_back(new Histo1D(filetag+"_ht",               100, 0, 5, "htDerivedweightUncor",          {"passMadHT"}, {"Lumi", "Weight"}));
        my_Histos.emplace_back(new Histo1D(filetag+"_ht_flat2000",      100, 0, 5, "htDerivedweightFlat200Uncor",   {"passMadHT"}, {"Lumi", "Weight"}));
        my_Histos.emplace_back(new Histo1D(filetag+"_ht_njet7",         100, 0, 5, "htDerivedweightNJet7Uncor",     {"passMadHT"}, {"Lumi", "Weight"}));
        my_Histos.emplace_back(new Histo1D(filetag+"_sclUp",            100, 0, 5, "scaleWeightUpUncor",            {"passMadHT"}, {"Lumi", "Weight"}));
        my_Histos.emplace_back(new Histo1D(filetag+"_sclDown",          100, 0, 5, "scaleWeightDownUncor",          {"passMadHT"}, {"Lumi", "Weight"}));
        my_Histos.emplace_back(new Histo1D(filetag+"_pdf_Up",           100, 0, 5, "PDFweightUpUncor",              {"passMadHT"}, {"Lumi", "Weight"}));
        my_Histos.emplace_back(new Histo1D(filetag+"_pdf_Down",         100, 0, 5, "PDFweightDownUncor",            {"passMadHT"}, {"Lumi", "Weight"}));    
        my_Histos.emplace_back(new Histo1D(filetag+"_pu",               100, 0, 5, "puWeight",                      {"passMadHT"}, {"Lumi", "Weight"}));    
        my_Histos.emplace_back(new Histo1D(filetag+"_pu_Up",            100, 0, 5, "puSysUp",                       {"passMadHT"}, {"Lumi", "Weight"}));    
        my_Histos.emplace_back(new Histo1D(filetag+"_pu_Down",          100, 0, 5, "puSysDown",                     {"passMadHT"}, {"Lumi", "Weight"}));    
    }

    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false)
    {
        const auto& runtype = tr.getVar<std::string>("runtype");
        const auto& filetag = tr.getVar<std::string>("filetag");

        //-------------------------------------
        //-- Initialize histograms to be filled
        //-------------------------------------
        InitHistos(filetag);

        while( tr.getNextEvent() )
        {
            //------------------------------------
            //-- Print Event Number
            //------------------------------------
            const bool breakLoop = printEventNum(maxevents, tr.getEvtNum());
            if(breakLoop) break;

            //-----------------------------------
            //-- Code unique to this Analyzer
            //-----------------------------------
            if( runtype != "MC" ) 
            {
                std::cerr<<"Please run over an MC file since these scale factors should not be applied to data!!"<<std::endl;
                break;
            }
               
            //-----------------------------------
            //-- Fill Histograms Below
            //-----------------------------------
            Fill(tr);

        }//END of while tr.getNextEvent loop   
    }//END of function
};

#endif
