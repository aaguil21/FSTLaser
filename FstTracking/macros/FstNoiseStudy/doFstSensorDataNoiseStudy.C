#include <iostream>
#include <string>
#include "TSystem.h"

R__LOAD_LIBRARY(../../lib/libFstNoiseStudy.dylib)

int doFstSensorDataNoiseStudy(string module = "Mod03", string hv = "HV70V", string date = "07122020")
{
  std::cout << "gSystem::Load <- libFstNoiseStudy.dylib" << endl;

  FstNoiseStudy *fst = new FstNoiseStudy();

  std::string inputlist = Form("../../list/FST/data/%s/FstData_%s_%s_%s.list",module.c_str(),module.c_str(),hv.c_str(),date.c_str());
  cout << "input list set to: " << inputlist.c_str() << endl;

  std::string outputfile = Form("../../output/noise/%s/FstDataNoise_%s_%s_%s.root",module.c_str(),module.c_str(),hv.c_str(),date.c_str());
  cout << "output file set to: " << outputfile.c_str() << endl;

  fst->set_list(inputlist.c_str());
  fst->set_outputfile(outputfile.c_str());

  int stats = fst->Init();
  if(stats < 0)
  {
    return -1;
  }

  fst->Make();
  fst->Finish();
  std::cout << "All done" << endl;

  return 1;
}
