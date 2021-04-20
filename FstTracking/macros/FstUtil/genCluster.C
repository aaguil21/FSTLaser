#include <iostream>
#include <string>
#include <vector>
#include "TSystem.h"

R__LOAD_LIBRARY(../../lib/libFstCluster.dylib)

int genCluster()
{
  std::cout << "gSystem::Load <- libFstCluster.dylib" << endl;

  // 1st hit
  FstRawHit *fstRawHit1 = new FstRawHit();
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
  cout << "1st hit => " << endl;
  fstRawHit1->Print();

  // 2nd hit
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
  cout << "2nd hit => " << endl;
  fstRawHit2->Print();

  // 1st cluster
  FstCluster *fstCluster1 = new FstCluster();
  fstCluster1->setLayer(1);
  fstCluster1->setSensor(0);
  fstCluster1->setMeanColumn(1.5);
  fstCluster1->setMeanRow(35.5);
  fstCluster1->setTotCharge(1024.3);
  fstCluster1->setMaxTb(3);
  fstCluster1->setClusterType(1);
  fstCluster1->setNRawHits(2);
  fstCluster1->setNRawHitsR(2);
  fstCluster1->setNRawHitsPhi(1);
  fstCluster1->addRawHit(fstRawHit1);
  fstCluster1->addRawHit(fstRawHit2);
  cout << "1st cluster =>" << endl;
  fstCluster1->Print();

  // 1st cluster
  FstCluster *fstCluster2 = new FstCluster();
  fstCluster2->setLayer(10);
  fstCluster2->setSensor(5);
  fstCluster2->setMeanColumn(2.5);
  fstCluster2->setMeanRow(3.5);
  fstCluster2->setTotCharge(124.3);
  fstCluster2->setMaxTb(4);
  fstCluster2->setClusterType(1);
  fstCluster2->setNRawHits(3);
  fstCluster2->setNRawHitsR(2);
  fstCluster2->setNRawHitsPhi(2);
  fstCluster2->addRawHit(fstRawHit2);
  cout << "2nd cluster =>" << endl;
  fstCluster2->Print();

  std::vector<FstCluster *> mClusterVec;
  mClusterVec.push_back(fstCluster1);
  mClusterVec.push_back(fstCluster2);
  mClusterVec[0]->Print();
  mClusterVec[1]->Print();

  mClusterVec[0]->Clear();
  mClusterVec[0]->Print();

  std::cout << "All done" << std::endl;


  return 1;
}
