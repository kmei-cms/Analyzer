#include "Analyzer/Analyzer/test/plotter.h"
#include <iostream>
#include <getopt.h>

class FisherHolder
{
public:
    std::vector<std::string> cutNames_;
    std::string plotName_;
};

void setHistInfo(const std::string& path, std::vector<histInfo>& data, std::vector<histInfo>& bg, std::vector<histInfo>& sig, const std::string& year, int color = 0, bool doQCD = false)
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry root file                 draw options  draw color
    data = {
        {"Data_SingleLepton" , path + "/"+year+"_Data.root", "PEX0", kBlack},
        ////{"MadGraph T#bar{T}",        path + "/"+year+"_TT.root",              "hist", kBlack   },
        //{"RPV 350: 2016", "condor/Analyze1Lep_Kerasv1.2.8/2016_RPV_2t6j_mStop-350.root",               "hist", kBlack, 1.0/35.9},
        //{"T#bar{T}: 2016", "condor/Analyze1Lep_Kerasv1.2.8/2016_TT.root",               "hist", kBlack, 1.0/35.9},
    };
    
    std::vector<int> bgColor  = {kRed, kBlack, kBlue};
    std::vector<int> sigColor = {kBlack, kBlue, kRed};

    bg = {
        //{"Powheg T#bar{T}",        path + "/"+year+"_TT.root",              "hist", kRed   },
        //{"RPV 350: 2017", "condor/Analyze1Lep_Kerasv3.0.4/2017_RPV_2t6j_mStop-350.root",               "hist", kRed, 1.0/41.525},
        //{"T#bar{T}: 2017", "condor/Analyze1Lep_Kerasv3.0.4/2017_TT.root",               "hist", kRed, 1.0/41.525},

        {"Triboson",        path + "/"+year+"_Triboson.root",        "hist", kGray       },
        {"Diboson",         path + "/"+year+"_Diboson.root",         "hist", kMagenta + 1},
        {"DYJetsToLL_M-50", path + "/"+year+"_DYJetsToLL_M-50.root", "hist", kOrange + 2 },        
        {"TTX",             path + "/"+year+"_TTX.root",             "hist", kCyan + 1   },
        {"WJets",           path + "/"+year+"_WJets.root",           "hist", kYellow + 1 },
        {"ST",              path + "/"+year+"_ST.root",              "hist", kRed + 1    },
        //{"RPV 350", path + "/"+year+"_RPV_2t6j_mStop-350.root",      "hist", kCyan  + 2*color, 0.3},
        {"QCD",             path + "/"+year+"_QCD.root",             "hist", kGreen + 1  },
        {"T#bar{T}",        path + "/"+year+"_TT.root",              "hist", kBlue - 6   },

        //{"T#bar{T}",        path + "/"+year+"_TT.root",              "hist", kBlue - 6   },
        //{"QCD",             path + "/"+year+"_QCD.root",             "hist", kGreen + 1  },
    };
    if(doQCD)
    {
        bg = {
            {"Triboson",        path + "/"+year+"_Triboson.root",        "hist", kGray       },
            {"Diboson",         path + "/"+year+"_Diboson.root",         "hist", kMagenta + 1},
            {"DYJetsToLL_M-50", path + "/"+year+"_DYJetsToLL_M-50.root", "hist", kOrange + 2 },        
            {"TTX",             path + "/"+year+"_TTX.root",             "hist", kCyan + 1   },
            {"WJets",           path + "/"+year+"_WJets.root",           "hist", kYellow + 1 },
            {"ST",              path + "/"+year+"_ST.root",              "hist", kRed + 1    },
            {"T#bar{T}",        path + "/"+year+"_TT.root",              "hist", kBlue - 6   },
            {"QCD",             path + "/"+year+"_QCD.root",             "hist", kGreen + 1  },
        };
    }

    sig = {        
        {"RPV 300", path + "/"+year+"_RPV_2t6j_mStop-300.root",      "hist", kCyan    + 2*color  },
        {"RPV 550", path + "/"+year+"_RPV_2t6j_mStop-550.root",      "hist", kMagenta + 2*color  },        
        {"RPV 850", path + "/"+year+"_RPV_2t6j_mStop-850.root",      "hist", kRed     + 2*color  },
        //{"SYY 550", path + "/"+year+"_StealthSYY_2t6j_mStop-550.root",       "hist", kMagenta  + 2*color  },        
    };
}

int main(int argc, char *argv[])
{
    TH1::AddDirectory(false);

    int opt, option_index = 0;    
    std::string year = "2018post";

    static struct option long_options[] = {
        {"year",     required_argument, 0, 'y'},
    };

    while((opt = getopt_long(argc, argv, "y:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'y': year = optarg; break;
        }
    }

    std::string path;
    if     (year=="2016") path= "condor/Analyze1Lep_2016_v1.0";
    else if(year=="2017") path= "condor/Analyze1Lep_2017_v1.0";
    else if(year=="2018") path= "condor/Analyze1Lep_2018_v1.0";
    else if(year=="2018pre") path= "condor/Analyze1Lep_2018pre_v1.0";
    else if(year=="2018post") path= "condor/Analyze1Lep_2018post_v1.0";

    std::vector<histInfo> data, bg, sig;
    std::vector<histInfo> dataQCD, bgQCD, sigQCD;
    setHistInfo(path, data, bg, sig, year, 0);
    setHistInfo(path, dataQCD, bgQCD, sigQCD, year, 0, true);
    HistInfoCollection histInfoCollection(data, bg, sig);
    HistInfoCollection histInfoCollectionQCD(dataQCD, bgQCD, sigQCD);

    // vector of histInfoCollection for Roc Curves
    std::map< std::string, HistInfoCollection > rocMap = { 
                                                           {year, histInfoCollection},                                                           
    };

    //make plotter object with the required sources for histograms specified
    Plotter pltRocCompare( std::move(rocMap) ,"outputPlots/FullRun2_"+year);
    Plotter plt( std::move(histInfoCollection) ,"outputPlots/FullRun2_"+year);
    Plotter pltQCD( std::move(histInfoCollectionQCD) ,"outputPlots/FullRun2_"+year);

    // --------------------
    // - Make stack plots
    // --------------------

    std::vector<std::string> mycuts_1l = {
        "",
        "_1l_HT300",
        "_1l_HT300_ge7j",
        "_1l_HT300_ge1b",
        "_1l_HT300_ge2b",
        "_1l_HT300_ge7j_ge1b",
        "_1l_HT300_ge1b_Mbl",
        "_1l_HT300_ge4j_ge1b_Mbl",
        "_1e_HT300_ge4j_ge1b_Mbl",
        "_1m_HT300_ge4j_ge1b_Mbl",
        "_1l_HT300_ge7j_ge1b_Mbl",
        "_1l_HT300_ge7j_ge1b_Mbl_noLepWeight",
        "_1l_HT300_ge7j_ge2b_Mbl",
        "_1l_HT300_ge7j_ge1b_Mbl_lBarrel",
        "_1e_HT300_ge7j_ge1b_Mbl_lBarrel",
        "_1m_HT300_ge7j_ge1b_Mbl_lBarrel",
        "_1l_HT300_ge7j_ge1b_Mbl_lEndCap",
        "_1e_HT300_ge7j_ge1b_Mbl_lEndCap",
        "_1m_HT300_ge7j_ge1b_Mbl_lEndCap",
        "_1l_HT300_ge7j_ge1b_Mbl_noHTWeight",
        "_1e_HT300_ge7j_ge1b_Mbl",
        "_1m_HT300_ge7j_ge1b_Mbl",
        "_1l_HT300_ge7j_ge1b_Mbl_d1",
        "_1l_HT300_ge7j_ge1b_Mbl_d2",
        "_1l_HT300_ge7j_ge1b_Mbl_d3",
        "_1l_HT300_ge7j_ge1b_Mbl_d4",
        //"_1l_HT300_1j_ge1b_Mbl",
        //"_1l_HT300_2j_ge1b_Mbl",
        //"_1l_HT300_3j_ge1b_Mbl",
        //"_1l_HT300_4j_ge1b_Mbl",
        //"_1l_HT300_5j_ge1b_Mbl",
        //"_1l_HT300_6j_ge1b_Mbl",
        //"_1l_HT300_7j_ge1b_Mbl",
        //"_1l_HT300_8j_ge1b_Mbl",
        //"_1l_HT300_9j_ge1b_Mbl",
        //"_1l_HT300_10j_ge1b_Mbl",
        //"_1l_HT300_11j_ge1b_Mbl",
        //"_1l_HT300_12j_ge1b_Mbl",
        //"_1l_HT300_13j_ge1b_Mbl",
        //"_1l_HT300_14j_ge1b_Mbl",
        //"_1l_HT300_15j_ge1b_Mbl",
        //"_1l_HT300_5j_ge1b_Mbl_htCorr",
        //"_1l_HT300_6j_ge1b_Mbl_htCorr",
        //"_1l_HT300_7j_ge1b_Mbl_htCorr",
        //"_1l_HT300_8j_ge1b_Mbl_htCorr",
        
        "_passQCDCR",
    };

    //plt.plotStack("h_njets_1l_ge7j_ge1b",    "N_{J}", "Events", true, -1, false, false);
    //plt.plotStack("h_njets_1l_ge7j_ge1b_d1", "N_{J}", "Events", true, -1, false, false);
    //plt.plotStack("h_njets_1l_ge7j_ge1b_d2", "N_{J}", "Events", true, -1, false, false);
    //plt.plotStack("h_njets_1l_ge7j_ge1b_d3", "N_{J}", "Events", true, -1, false, false);
    //plt.plotStack("h_njets_1l_ge7j_ge1b_d4", "N_{J}", "Events", true, -1, false, false);
    for(std::string mycut : mycuts_1l)
    {
        plt.plotStack(    "h_njets"+mycut,         "N_{J}" ,             "Events",      true, -1, false, true);
        pltQCD.plotStack( "h_njetsQCDCR"+mycut,    "N_{J}" ,             "Events",      true, -1, false, true);
        plt.plotStack(    "h_deepESM"+mycut,       "DeepESM" ,           "Events",      true, 20, false, true);
        pltQCD.plotStack( "h_deepESMQCDCR"+mycut,  "DeepESM" ,           "Events",      false,20, false, true);
        plt.plotStack(    "h_deepESMMerged"+mycut, "DeepESM Bin",        "Events",      false,-1, false, true);
        plt.plotStack(    "h_mbl"+mycut,           "M(l,b) [GeV]",       "Events",      true, 10, false, true);
        plt.plotStack(    "h_ht"+mycut,            "H_{T} [GeV]",        "Events",      true, 10, false, true);
        pltQCD.plotStack( "h_htQCDCR"+mycut,       "H_{T} [GeV]",        "Events",      true, 10, false, true);
        plt.plotStack(    "h_lPt"+mycut,           "Lepton P_{T} [GeV]", "Num Leptons", true,  2, false, true, 0, 1000);
        plt.plotStack(    "h_lEta"+mycut,          "Lepton Eta",         "Num Leptons", true,  2, false, true);
        plt.plotStack(    "h_lPhi"+mycut,          "Lepton Phi",         "Num Leptons", false, 2, false, true);
        plt.plotStack(    "h_jPt"+mycut,           "Jet P_{T} [GeV]",    "Num Jets",    true,  2, false, true, 0, 1000);
        plt.plotStack(    "h_jEta"+mycut,          "Jet Eta",            "Num Jets",    true,  2, false, true);
        plt.plotStack(    "h_jPhi"+mycut,          "Jet Phi",            "Num Jets",   false,  2, false, true);
        plt.plotStack(    "h_ntops"+mycut,         "N_{T}" ,             "Events",      true, -1, false, true);
        plt.plotStack(    "h_nb"   +mycut,         "N_{B}" ,             "Events",      true, -1, false, true);        
        //Need these until we un blind
        plt.plotStack( "blind_njets"+mycut,         "N_{J}" ,             "Events",      true, -1, false, true);
        plt.plotStack( "blind_deepESM"+mycut,       "DeepESM" ,           "Events",      true, 10, false, true);
        plt.plotStack( "blind_deepESMMerged"+mycut, "DeepESM Bin",        "Events",      true, -1, false, true);
        plt.plotStack( "blind_mbl"+mycut,           "M(l,b) [GeV]",       "Events",      true, 10, false, true);
        plt.plotStack( "blind_ht"+mycut,            "H_{T} [GeV]",        "Events",      true, 10, false, true);
        plt.plotStack( "blind_lPt"+mycut,           "Lepton P_{T} [GeV]", "Num Leptons", true,  2, false, true, 0, 1000);
        plt.plotStack( "blind_lEta"+mycut,          "Lepton Eta",         "Num Leptons", true,  2, false, true);
        plt.plotStack( "blind_lPhi"+mycut,          "Lepton Phi",         "Num Leptons", false, 2, false, true);
        plt.plotStack( "blind_jPt"+mycut,           "Jet P_{T} [GeV]",    "Num Jets",    true,  2, false, true, 0, 1000);
        plt.plotStack( "blind_jEta"+mycut,          "Jet Eta",            "Num Jets",    true,  2, false, true);
        plt.plotStack( "blind_jPhi"+mycut,          "Jet Phi",            "Num Jets",   false,  2, false, true);
        plt.plotStack( "blind_ntops"+mycut,         "N_{T}" ,             "Events",      true, -1, false, true);
        plt.plotStack( "blind_nb"   +mycut,         "N_{B}" ,             "Events",      true, -1, false, true);        
        //
        plt.plotNormFisher("h_deepESM"+mycut, "DeepESM" , "Norm", false, 10);
        plt.plotNormFisher("h_njets"+mycut, "DeepESM" , "Events", false);
        // - Make  Roc Curve
        pltRocCompare.plotRocFisher("h_deepESM"+mycut,"Background","Signal", true, false);        
    }
    //plt.plotStack( "fwm2_top6_1l_ge7j_ge1b",  "FWM2", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "fwm3_top6_1l_ge7j_ge1b",  "FWM3", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "fwm4_top6_1l_ge7j_ge1b",  "FWM4", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "fwm5_top6_1l_ge7j_ge1b",  "FWM5", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "fwm6_top6_1l_ge7j_ge1b",  "FWM6", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "fwm7_top6_1l_ge7j_ge1b",  "FWM7", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "fwm8_top6_1l_ge7j_ge1b",  "FWM8", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "fwm9_top6_1l_ge7j_ge1b",  "FWM9", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "fwm10_top6_1l_ge7j_ge1b", "FWM10","Events", true, -1, false, 0, 2);
    //plt.plotStack( "jmt_ev0_top6_1l_ge7j_ge1b", "JMT0", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "jmt_ev1_top6_1l_ge7j_ge1b", "JMT1", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "jmt_ev2_top6_1l_ge7j_ge1b", "JMT2", "Events", true, -1, false, 0, 2);
    //for(unsigned int i = 0; i < 7; i++)
    //{
    //    plt.plotStack("Jet_cm_pt_"+std::to_string(i+1)+"_1l_ge7j_ge1b",  "Jet_"+std::to_string(i+1)+" P_{T}", "Events", true, -1, false, 0, 1500);
    //    plt.plotStack("Jet_cm_eta_"+std::to_string(i+1)+"_1l_ge7j_ge1b", "Jet_"+std::to_string(i+1)+" Eta",   "Events", true, -1, false, -6, 6);
    //    plt.plotStack("Jet_cm_phi_"+std::to_string(i+1)+"_1l_ge7j_ge1b", "Jet_"+std::to_string(i+1)+" Phi",   "Events", true, -1, false, -4, 4);
    //    plt.plotStack("Jet_cm_m_"+std::to_string(i+1)  +"_1l_ge7j_ge1b", "Jet_"+std::to_string(i+1)+" Mass",  "Events", true, -1, false,  0, 200);
    //}
    //plt.plotStack( "h_leptonweight_1l_ge7j_ge1b",    "Lepton W",         "Events", true, -1, false, 0, 2);
    //plt.plotStack( "h_weight_1l_ge7j_ge1b",          "Total Weight",     "Events", true, -1, false, 0, 2);
    //plt.plotStack( "h_htDerivedweight_1l_ge7j_ge1b", "HtDerived Weight", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "h_bTagWeight_1l_ge7j_ge1b",      "BTagWeight",       "Events", true, -1, false, 0, 2);
    
    // --------------------
    // - Make fisher plots
    // --------------------

    std::vector< FisherHolder > fisherHolder
    {
        { {"h_njets_1l_HT300_ge7j_ge1b_Mbl_d1", "h_njets_1l_HT300_ge7j_ge1b_Mbl_d2", "h_njets_1l_HT300_ge7j_ge1b_Mbl_d3"  , "h_njets_1l_HT300_ge7j_ge1b_Mbl_d4"  } , "njets_1l_HT300_ge7j_ge1b_Mbl"  },
    };
    
    for (auto& f : fisherHolder)
    {
        plt.plotFisher(f.cutNames_,      f.plotName_, "N_{J}", "Events",           true, 7, "DeepESM");
        plt.plotRatioFisher(f.cutNames_, f.plotName_, "N_{J}", "N_{J+1} / N_{J}", false, 6, "DeepESM");
    }
    
    // --------------------
    // - Compute Yields
    // --------------------

    //for(std::string mycut : mycuts_1l)
    //{
    //    const auto& yieldMap = histInfoCollection.computeYields("h_njets_1l"+mycut,"njets",12,20);
    //    //const auto& yieldMap = histInfoCollection.computeYields("h_njets_1l"+mycut,"njets",6,8);
    //    auto allbg  = yieldMap.find("AllBG"); auto data_JetHT = yieldMap.find("Data_JetHT"); auto ttbar = yieldMap.find("T#bar{T}"); auto QCD = yieldMap.find("QCD");
    //    printf("%45.45s:  %s: %12.1lf   %s: %12.1lf   %s: %12.1lf   %s: %12.1lf\n",("h_njets_1l"+mycut).c_str(), (data_JetHT->first).c_str(), data_JetHT->second, (allbg->first).c_str(), allbg->second, "ttbar", ttbar->second, (QCD->first).c_str(), QCD->second );
    //}       
}
