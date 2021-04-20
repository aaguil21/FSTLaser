#include <TMath.h>

#include "../../src/FstUtil/FstCons.h"

void calProjResolution()
{
  const double z0_fst = FST::pitchLayer03;
  const double z1_ist = FST::pitchLayer12 + FST::pitchLayer23;
  const double z2_ist = FST::pitchLayer23;
  const double z3_ist = 0.0;

  const double sigmaX = FST::pitchColumn/TMath::Sqrt(12.0);
  const double sigmaY = FST::pitchRow/TMath::Sqrt(12.0);

  double sigmaProjX = TMath::Sqrt(z0_fst*z0_fst*sigmaX*sigmaX/(z1_ist*z1_ist)+(1-z0_fst/z1_ist)*(1-z0_fst/z1_ist)*sigmaX*sigmaX);
  double sigmaProjY = TMath::Sqrt(z0_fst*z0_fst*sigmaY*sigmaY/(z1_ist*z1_ist)+(1-z0_fst/z1_ist)*(1-z0_fst/z1_ist)*sigmaY*sigmaY);

  cout << "sigmaX = " << sigmaX << ", sigmaProjX = " << sigmaProjX << endl;
  cout << "sigmaY = " << sigmaY << ", sigmaProjY = " << sigmaProjY << endl;
}
