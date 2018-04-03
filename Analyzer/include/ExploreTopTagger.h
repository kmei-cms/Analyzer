#ifndef ExploreTopTagger_h
#define ExploreTopTagger_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class ExploreTopTagger
{
public:
   std::map<std::string, TH1D*>  my_histos;
   std::map<std::string, TH2D*>  my_2d_histos;
   std::map<std::string, TEfficiency*>  my_efficiencies;

   ExploreTopTagger(){};
   ~ExploreTopTagger(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, std::string type = "", std::string filetag = "", bool isQuiet = false);
   void     InitHistos();
   void     WriteHistos();

};

#endif
