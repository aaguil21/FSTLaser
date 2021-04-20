#include <iostream>
#include <string>
#include "TSystem.h"

R__LOAD_LIBRARY(../../lib/libFstNoiseStudy.dylib)

int doIstChipNoiseStudy(string module = "Mod36")
{
  std::cout << "gSystem::Load <- libFstNoiseStudy.dylib" << endl;

  FstNoiseStudy *fst = new FstNoiseStudy();

  std::string inputlist = Form("../../list/IST/IstChipPed_%s.list",module.c_str());
  cout << "input list set to: " << inputlist.c_str() << endl;

  std::string outputfile = Form("../../output/noise/Ist%s/IstChipNoise_%s.root",module.c_str(),module.c_str());
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
