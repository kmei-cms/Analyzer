#include "Analyzer/Analyzer/test/plotter.h"
#include <iostream>

std::string year = "2016";

class FisherHolder
{
public:
    std::vector<std::string> cutNames_;
    std::string plotName_;
};

void setHistInfo(const std::string& path, std::vector<histInfo>& data, std::vector<histInfo>& bg, std::vector<histInfo>& sig, int color = 0)
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry root file                 draw options  draw color
    data = {
        {"Data_SingleLepton" , path + "/"+year+"_Data.root", "PEX0", kBlack},
    };
    
    std::vector<int> bgColor = {kRed, kBlack, kBlue};
    std::vector<int> sigColor = {kBlack, kBlue, kRed};

    bg = {
        {"DYJetsToLL_M-50", path + "/"+year+"_DYJetsToLL_M-50.root", "hist", kOrange + 2 },        
        {"Rare",            path + "/"+year+"_Rare.root",            "hist", kCyan + 1   },
        {"Diboson",         path + "/"+year+"_Diboson.root",         "hist", kMagenta + 1},
        {"WJetsToLNu",      path + "/"+year+"_WJetsToLNu.root",      "hist", kYellow + 1 },
        {"ST",              path + "/"+year+"_ST.root",              "hist", kRed + 1    },
        {"QCD",             path + "/"+year+"_QCD.root",             "hist", kGreen + 1  },
        {"T#bar{T}",        path + "/"+year+"_TT.root",              "hist", kBlue       },
    };
    if(year == "2016")
    {        
        sig = {
            {"RPV 850", path + "/"+year+"_rpv_stop_850.root",               "hist", kRed  + 2*color      },
            //{"RPV 650", path + "/"+year+"_rpv_stop_650.root",               "hist", kGreen  + 2*color    },        
            //{"RPV 550", path + "/"+year+"_rpv_stop_550.root",               "hist", kMagenta  + 2*color  },        
            //{"RPV 450", path + "/"+year+"_rpv_stop_450.root",               "hist", kOrange  + 2*color   },        
            {"RPV 350", path + "/"+year+"_rpv_stop_350.root",               "hist", kCyan  + 2*color     },
            //
            //{"SYY 850", path + "/"+year+"_stealth_stop_850_SYY.root",       "hist", kRed  + 2*color      },
            //{"SYY 650", path + "/"+year+"_stealth_stop_650_SYY.root",       "hist", kGreen  + 2*color    },        
            {"SYY 550", path + "/"+year+"_stealth_stop_550_SYY.root",       "hist", kMagenta  + 2*color  },        
            //{"SYY 450", path + "/"+year+"_stealth_stop_450_SYY.root",       "hist", kOrange  + 2*color   },        
            //{"SYY 350", path + "/"+year+"_stealth_stop_350_SYY.root",       "hist", kCyan  + 2*color     },
            
            //{"SHuHd 850", path + "/"+year+"_stealth_stop_850_SHuHd.root",     "hist", kRed  + 2*color      },
            //{"SHuHd 650", path + "/"+year+"_stealth_stop_650_SHuHd.root",     "hist", kGreen  + 2*color    },        
            //{"SHuHd 550", path + "/"+year+"_stealth_stop_550_SHuHd.root",     "hist", kMagenta  + 2*color  },        
            //{"SHuHd 450", path + "/"+year+"_stealth_stop_450_SHuHd.root",     "hist", kOrange  + 2*color   },        
            //{"SHuHd 350", path + "/"+year+"_stealth_stop_350_SHuHd.root",     "hist", kCyan  + 2*color     },
        };
    }
    else
    {
        sig = {        
            {"RPV 850", path + "/"+year+"_RPV_2t6j_mStop-850.root",               "hist", kRed  + 2*color      },
            //{"RPV 650", path + "/"+year+"_RPV_2t6j_mStop-650.root",               "hist", kGreen  + 2*color    },        
            //{"RPV 550", path + "/"+year+"_RPV_2t6j_mStop-550.root",               "hist", kMagenta  + 2*color  },        
            //{"RPV 450", path + "/"+year+"_RPV_2t6j_mStop-450.root",               "hist", kOrange  + 2*color   },        
            {"RPV 350", path + "/"+year+"_RPV_2t6j_mStop-350.root",               "hist", kCyan  + 2*color     },
            
            //{"SYY 850", path + "/"+year+"_StealthSYY_2t6j_mStop-850.root",       "hist", kRed  + 2*color      },
            //{"SYY 650", path + "/"+year+"_StealthSYY_2t6j_mStop-650.root",       "hist", kGreen  + 2*color    },        
            {"SYY 550", path + "/"+year+"_StealthSYY_2t6j_mStop-550.root",       "hist", kMagenta  + 2*color  },        
            //{"SYY 450", path + "/"+year+"_StealthSYY_2t6j_mStop-450.root",       "hist", kOrange  + 2*color   },        
            //{"SYY 350", path + "/"+year+"_StealthSYY_2t6j_mStop-350.root",       "hist", kCyan  + 2*color     },
            //
            //{"SHuHd 850", path + "/"+year+"_StealthSHH_2t4b_mStop-850.root",     "hist", kRed  + 2*color      },
            //{"SHuHd 650", path + "/"+year+"_StealthSHH_2t4b_mStop-650.root",     "hist", kGreen  + 2*color    },        
            //{"SHuHd 550", path + "/"+year+"_StealthSHH_2t4b_mStop-550.root",     "hist", kMagenta  + 2*color  },        
            //{"SHuHd 450", path + "/"+year+"_StealthSHH_2t4b_mStop-450.root",     "hist", kOrange  + 2*color   },        
            //{"SHuHd 350", path + "/"+year+"_StealthSHH_2t4b_mStop-350.root",     "hist", kCyan  + 2*color     },
        };
    }
}

int main()
{
    TH1::AddDirectory(false);

    //std::string pathGRfalse = "deepESM_GRfalse_1Layer_8Vars";
    //std::string pathGRfalse = "deepESM_GRfalse_1Layer_24Vars";
    //std::string pathGRfalse = "deepESM_GRfalse_1Layer_9Vars_nJets";
    //std::string pathGRfalse = "condor/Analyze1Lep_Kerasv1.3.0";

    //std::string pathGRtrue = "condor/Analyze1Lep_Kerasv1.2.0_MCTrigger_bTag_leptonWeight_ht300";
    //std::string pathGRtrue = "condor/Analyze1Lep_Kerasv1.2.3";
    //std::string pathGRtrue = "condor/Analyze1Lep_Kerasv1.2.4";
    //std::string pathGRtrue = "condor/Analyze1Lep_Kerasv1.2.4_HTSF";
    //std::string pathGRtrue = "condor/Analyze1Lep_Kerasv1.2.5_noHT";
    std::string pathGRtrue = "condor/oldTest/Analyze1Lep_Kerasv1.2.5_Freezing";
    //std::string pathGRtrue = "condor/Analyze1Lep_Kerasv3.0.0";
    //std::string pathGRtrue = "condor/Analyze1Lep_Kerasv3.0.0_v2";
    //std::string pathGRtrue = "condor/Analyze1Lep_Kerasv3.0.1";
    //std::string pathGRtrue = "condor/Analyze1Lep_Kerasv3.0.1_Freezing";
    //std::string pathFisher = "oldTest/DeepESMTests_JoesCode/deepESM_v1";
    //std::string pathPhoton = "photonCR";
    //std::string pathPhoton = "photonCR_Barrel";

    std::vector<histInfo> data_GRfalse, bg_GRfalse, sig_GRfalse;
    std::vector<histInfo> data_GRtrue, bg_GRtrue, sig_GRtrue;
    std::vector<histInfo> data_fisher, bg_fisher, sig_fisher;
    std::vector<histInfo> data_photon, bg_photon, sig_photon;

    //setHistInfo(pathGRfalse, data_GRfalse, bg_GRfalse, sig_GRfalse, 1);
    setHistInfo(pathGRtrue, data_GRtrue, bg_GRtrue, sig_GRtrue, 0);
    //setHistInfo(pathFisher, data_fisher, bg_fisher, sig_fisher, 2);
    //setHistInfo(pathPhoton, data_photon, bg_photon, sig_photon, 0);

    //make histInfoCollection
    //HistInfoCollection histInfoCollection_GRfalse(data_GRfalse, bg_GRfalse, sig_GRfalse);
    HistInfoCollection histInfoCollection_GRtrue(data_GRtrue, bg_GRtrue, sig_GRtrue);
    //HistInfoCollection histInfoCollection_fisher(data_fisher, bg_fisher, sig_fisher);
    //HistInfoCollection histInfoCollection_photon(data_photon, bg_photon, sig_photon);

    // vector of histInfoCollection for Roc Curves
    std::map< std::string, HistInfoCollection > rocMap = { //{"No 350", histInfoCollection_GRfalse},
                                                           {year, histInfoCollection_GRtrue},
                                                           //{"fisher", histInfoCollection_fisher},
    };

    //make plotter object with the required sources for histograms specified
    //Plotter pltRoc( std::move(rocMapTest) );
    Plotter pltRocCompare( std::move(rocMap) );
    Plotter plt( std::move(histInfoCollection_GRtrue) );
    //Plotter plt( std::move(histInfoCollection_GRfalse) );
    //Plotter plt( std::move(histInfoCollection_fisher) );
    //Plotter plt( std::move(histInfoCollection_photon) );

    // --------------------
    // - Make stack plots
    // --------------------

    std::vector<std::string> mycuts_1l = {
        //"",
        //"_1l",
        //"_1l_ge7j",                        
        //"_1l_ge1b",                       
        //"_1l_ge2b",                       
        //"_1e_1m_ge2b_le5j",
        //"_1l_1t",                          
        //"_1l_2t",                          
        //"_1l_ge1t",                        
        //"_1l_ge2t",                        
        //"_1l_4j_ge1b",
        //"_1l_5j_ge1b",
        //"_1l_6j_ge1b",
        //"_1l_7j_ge1b", 
        //"_1l_8j_ge1b", 
        //"_1l_9j_ge1b", 
        //"_1l_10j_ge1b", 
        //"_1l_11j_ge1b", 
        //"_1l_12j_ge1b", 
        //"_1l_13j_ge1b", 
        //"_1l_14j_ge1b", 
        //"_1l_15j_ge1b",
        //"_1l_ge7j_ge1b_noMbl",
        "_1l_ge7j_ge1b",
        //"_1e_ge7j_ge1b",
        //"_1m_ge7j_ge1b",
        //"_1l_0b_ge300ht_50to110mt_ge30MET",
        //"_1l_0b_ge300ht_50to110mt_ge30MET_even",
        //"_1l_0b_ge300ht_50to110mt_ge30MET_odd",
        //"_1l_ge7j_ge1b_d1",
        //"_1l_ge7j_ge1b_d2",
        //"_1l_ge7j_ge1b_d3",
        //"_1l_ge7j_ge1b_d4",
        //"_1l_ge7j_ge2b",                 
        //"_1l_ge7j_ge1b_1t",
        //"_1l_ge7j_ge1b_ge1t",                                                    
        //"_1l_ge7j_ge1b_1t1",               
        //"_1l_ge7j_ge1b_1t2",               
        //"_1l_ge7j_ge1b_1t3",               
        //"_1l_ge7j_ge1b_1t2or3",            
        //"_1l_ge7j_ge1b_ge1t1",             
        //"_1l_ge7j_ge1b_ge1t2",             
        //"_1l_ge7j_ge1b_ge1t3",             
    };

    for(std::string mycut : mycuts_1l)
    {
        plt.plotStack( "h_njets"+mycut, "N_{J}" ,              "Events", true, -1, true);
        plt.plotStack( "h_deepESM"+mycut, "DeepESM" ,          "Events", true, 10, true);
        plt.plotStack( "h_deepESMMerged"+mycut, "DeepESM Bin", "Events", true, -1, true);
        //plt.plotStack( "h_mbl"+mycut,    "M(l,b) [GeV]",       "Events", true, 10);
        //plt.plotStack( "h_allMbl"+mycut, "M(l,b) [GeV]",       "Events", true, 10, false);
        //plt.plotStack( "h_ht"+mycut,     "H_{T} [GeV]",        "Events", true, 10, false);
        //plt.plotStack( "h_lPt"+mycut,    "Lepton P_{T} [GeV]", "Events", true, 2, false, 0, 1000);
        //plt.plotStack( "h_lEta"+mycut,   "Lepton Eta",         "Events", true, 2);
        //plt.plotStack( "h_jPt"+mycut,     "Jet P_{T} [GeV]",     "Events", true,  2, false, 0, 1000);
        //plt.plotStack( "h_jEta"+mycut,    "Jet Eta",             "Events", true,  2, false);
        //plt.plotStack( "h_ntops"+mycut, "N_{T}" , "Events", true);
        //plt.plotStack( "h_nb"   +mycut, "N_{B}" , "Events", true);        
        //plt.plotStack( "h_fisher"+mycut, "fisher value" , "Events", true, 4);        
        //plt.plotStack( "h_photonPt"+mycut, "P_{T}^{#gamma}" , "Events", false, 10);
        //
        //// Make Normalized fisher
        //pltSkim.plotNormFisher("h_fisher_1l"+mycut, "fisher value" , "Events", false, 4);

        //plt.plotNormFisher("h_deepESM"+mycut, "DeepESM" , "Norm", false, 10);
        
        //plt.plotNormFisher("h_njets"+mycut, "DeepESM" , "Events", true);
        //plt.plotNormFisher("h_BestComboMass_1l"+mycut, "Average BestCombo Mass [GeV]" , "Events", false, 4);
        //plt.plotNormFisher("h_BestComboPt_1l"+mycut, "Average BestCombo P_{T} [GeV]" , "Events", false, 4);
        //plt.plotNormFisher("h_BestComboMassDiff_1l"+mycut, "BestCombo Mass Diff [GeV]" , "Events", false, 2);
        //plt.plotNormFisher("h_BestComboMassDiffAbs_1l"+mycut, "BestCombo Mass Abs(Diff) [GeV]" , "Events", false, 2);
        //plt.plotNormFisher("h_BestComboRelDiff_1l"+mycut, "BestCombo Rel Diff" , "Events", false, 2);
        //plt.plotNormFisher("h_BestComboRelDiffAbs_1l"+mycut, "BestCombo Rel Abs(Diff)" , "Events", false, 2);
        //
        // - Make  Roc Curve
        //pltRoc.plotRocFisher("h_deepESM_1l"+mycut,"Background","Signal", false);
        
        //pltRocCompare.plotRocFisher("h_deepESM"+mycut,"Background","Signal", true, false);
        
        //////Need these until we un blind
        //plt.plotStack( "blind_njets"+mycut,   "N_{J}",               "Events", true, -1, false);
        //plt.plotStack( "blind_deepESM"+mycut, "DeepESM",             "Events", true, 10, false);
        //plt.plotStack( "blind_deepESMMerged"+mycut, "DeepESM Bin",   "Events", true, -1, false);
        //plt.plotStack( "blind_mbl"+mycut,     "M(l,b) [GeV]",        "Events", true, 10, false);
        //plt.plotStack( "blind_allMbl"+mycut,  "M(l,b) [GeV]",        "Events", true, 10, false);
        //plt.plotStack( "blind_ht"+mycut,      "H_{T} [GeV]",         "Events", true, 10, false);
        //plt.plotStack( "blind_ntops"+mycut,   "N_{T}",               "Events", true, -1, false);
        //plt.plotStack( "blind_nb"   +mycut,   "N_{B}",               "Events", true, -1, false);        
        //plt.plotStack( "blind_lPt"+mycut,     "Lepton P_{T} [GeV]",  "Events", true,  2, false, 0, 1000);
        //plt.plotStack( "blind_lEta"+mycut,    "Lepton Eta",          "Events", true,  2, false);
        //plt.plotStack( "blind_jPt"+mycut,     "Jet P_{T} [GeV]",     "Events", true,  2, false, 0, 1000);
        //plt.plotStack( "blind_jEta"+mycut,    "Jet Eta",             "Events", true,  2, false);
        //////plt.plotStack( "blind_fisher"+mycut, "fisher value" , "Events", true, 10);        
    }
    //plt.plotStack( "h_leptonweight_1l_ge7j_ge1b",    "Lepton W",         "Events", true, -1, false, 0, 2);
    //plt.plotStack( "h_weight_1l_ge7j_ge1b",          "Total Weight",     "Events", true, -1, false, 0, 2);
    //plt.plotStack( "h_htDerivedweight_1l_ge7j_ge1b", "HtDerived Weight", "Events", true, -1, false, 0, 2);
    //plt.plotStack( "h_bTagWeight_1l_ge7j_ge1b",      "BTagWeight",       "Events", true, -1, false, 0, 2);
    
    // --------------------
    // - Make fisher plots
    // --------------------

    //std::vector< FisherHolder > fisherHolder
    //{
    //    //{ {"h_njets_1l_ge7j_ge1b_d1"  , "h_njets_1l_ge7j_ge1b_d2"  , "h_njets_1l_ge7j_ge1b_d3"  , "h_njets_1l_ge7j_ge1b_d4"  } , "njets_1l_ge7j_ge1b"  },
    //    { {"h_njets_1l_ge7j_ge1b_d1"  } , "njets_1l_ge7j_ge1b"  },
    //    //{ {"h_njets_1l_ge7j_ge1b_ge6-7esm", "h_njets_1l_ge7j_ge1b_ge7-8esm", "h_njets_1l_ge7j_ge1b_ge8-95esm", "h_njets_1l_ge7j_ge1b_ge95esm" }, "njets_1l_ge7j_ge1b_esmBins"},
    //    //{ {"h_njets_1l_ge7j_ge1b_ge5-6esm", "h_njets_1l_ge7j_ge1b_ge6-7esm", "h_njets_1l_ge7j_ge1b_ge7-8esm", "h_njets_1l_ge7j_ge1b_ge8-95esm", "h_njets_1l_ge7j_ge1b_ge95esm" }, "njets_1l_ge7j_ge1b_esmBins"},
    //    //{ {"h_njets_1l_ge7j_ge1b_f1"  , "h_njets_1l_ge7j_ge1b_f2"  , "h_njets_1l_ge7j_ge1b_f3"  , "h_njets_1l_ge7j_ge1b_f4"  } , "njets_1l_ge7j_ge1b"  },
    //};
    //
    //for (auto& f : fisherHolder)
    //{
    //    //plt.plotFisher(f.cutNames_,  f.plotName_, "N_{J}", "Events", true, 9, "DeepESM");
    //    plt.plotRatioFisher(f.cutNames_,  f.plotName_, "N_{J}", "N_{J+1} / N_{J}", false, 6, "DeepESM");
    //    //plt.plotFisher(f.cutNames_,  f.plotName_, "N_{J}", "Events", true, 9, "Fisher");
    //    //plt.plotRatioFisher(f.cutNames_,  f.plotName_, "N_{J}", "N_{J+1} / N_{J}", false, 6, "Fisher");
    //}
    
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
