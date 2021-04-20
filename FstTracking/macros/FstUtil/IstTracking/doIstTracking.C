#include <iostream>
#include <string>
#include "TSystem.h"

R__LOAD_LIBRARY(../../lib/libIstTracking.dylib)

int doIstTracking()
{
  std::cout << "gSystem::Load <- libIstTracking.dylib" << endl;

  IstTracking *ist = new IstTracking();

  ist->Init();
  ist->Make();
  ist->Finish();


  std::cout << "All done" << endl;

  return 1;
}
