#include <iostream>
#include <string>
#include <vector>
#include "TSystem.h"

R__LOAD_LIBRARY(../../lib/libFstRawHit.dylib)

int genRawHit()
{
  std::cout << "gSystem::Load <- libFstRawHit.dylib" << endl;

  FstRawHit *fstRawHit1 = new FstRawHit();

  fstRawHit1->Print();
  fstRawHit1->setLayer(0);
  fstRawHit1->setSensor(1);
  fstRawHit1->setColumn(14);
  fstRawHit1->setRow(3);
  for(int i_tb = 0; i_tb < 7; ++i_tb)
  {
    fstRawHit1->setCharge(i_tb*100.0, i_tb);
  }
  fstRawHit1->setMaxTb(4);
  fstRawHit1->setHitId(1024);
  fstRawHit1->Print();

  /*
  cout << "layer = " << fstRawHit1->getLayer() << endl;
  cout << "sensor = " << fstRawHit1->getSensor() << endl;
  cout << "column = " << fstRawHit1->getColumn() << endl;
  cout << "row = " << fstRawHit1->getRow() << endl;
  for(int i_tb = 0; i_tb < 9; ++i_tb)
  {
    cout << "i_tb = " << i_tb << ", charge = " << fstRawHit1->getCharge(i_tb) << endl;
  }
  cout << "tb = -1 & charge = " << fstRawHit1->getCharge(-1) << endl;
  cout << "tb = 9 & charge = " << fstRawHit1->getCharge(9) << endl;
  cout << "maxTb = " << fstRawHit1->getMaxTb() << endl;
  cout << "hitId = " << fstRawHit1->getHitId() << endl;
  */

  FstRawHit *fstRawHit2 = new FstRawHit();

  fstRawHit2->setLayer(3);
  fstRawHit2->setSensor(0);
  fstRawHit2->setColumn(4);
  fstRawHit2->setRow(43);
  for(int i_tb = 0; i_tb < 7; ++i_tb)
  {
    fstRawHit2->setCharge(i_tb*100.0, i_tb);
  }
  fstRawHit2->setMaxTb(5);
  fstRawHit2->setHitId(1);
  fstRawHit2->Print();


  std::vector<FstRawHit *> mRawHitVec;
  mRawHitVec.push_back(fstRawHit1);
  mRawHitVec.push_back(fstRawHit2);

  mRawHitVec[0]->Print();
  mRawHitVec[1]->Print();

  mRawHitVec[1]->Clear();
  mRawHitVec[1]->Print();

  std::cout << "All done" << std::endl;

  return 1;
}
