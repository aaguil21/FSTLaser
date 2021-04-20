#include <iostream>
#include <string>
#include "TSystem.h"
#include "TString.h"

R__LOAD_LIBRARY(../../lib/libFstTracking.dylib)

int doFstTracking_timesequence(string hv = "HV140V", bool isSavePed = true, bool isApplyCMNCorr = false, float nFstHitsCut = 4.5, int numOfUsedTimeBins = 3, float nFstThresholdCut = 3.5, string date = "04022020")
{
  std::cout << "gSystem::Load <- libFstTracking.dylib" << endl;

  std::string pedMode = "withPed";
  if(!isSavePed) pedMode = "woPed";
  std::string cmnMode = "withCMNCorr";
  if(!isApplyCMNCorr) cmnMode = "woCMNCorr";

  std::string inputlist = Form("../../list/FST/timesequence/FstCluster_%s_Th%1.1fTb%dPed%1.1f_%s_%s_%s.list",hv.c_str(),nFstHitsCut,numOfUsedTimeBins,nFstThresholdCut,pedMode.c_str(),cmnMode.c_str(),date.c_str());
  cout << "input list set to: " << inputlist.c_str() << endl;

  std::string outputfile = Form("../../output/timesequence/FstTracking_%s_Th%1.1fTb%dPed%1.1f_%s_%s_%s.root",hv.c_str(),nFstHitsCut,numOfUsedTimeBins,nFstThresholdCut,pedMode.c_str(),cmnMode.c_str(),date.c_str());
  cout << "output file set to: " << outputfile.c_str() << endl;

  FstTracking *fst = new FstTracking();
  fst->set_list(inputlist.c_str());
  fst->set_outputfile(outputfile.c_str());

  fst->Init();
  fst->Make();
  fst->Finish();

  std::cout << "All done" << endl;

  return 1;
}
