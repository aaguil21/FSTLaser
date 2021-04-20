#include <iostream>
#include <string>
#include <vector>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

R__LOAD_LIBRARY(../../lib/libFstRawHit.dylib)
R__LOAD_LIBRARY(../../lib/libFstCluster.dylib)

int testCluster()
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


  TFile *File_OutPut = new TFile("output.root","RECREATE");
  TTree *Tree_FstCluster = new TTree("Tree_FstCluster","Cluster info");
  FstCluster *fstCluster1 = new FstCluster();
  Tree_FstCluster->Branch("FstCluster","FstCluster",&fstCluster1);

  for(int i_cluster = 0; i_cluster < 10; ++i_cluster)
  {
    fstCluster1->setLayer(1);
    fstCluster1->setSensor(0);
    fstCluster1->setMeanColumn(1.5+i_cluster);
    fstCluster1->setMeanRow(35.5+i_cluster);
    fstCluster1->setTotCharge(1024.3+i_cluster*10);
    fstCluster1->setMaxTb(3);
    fstCluster1->setClusterType(1);
    fstCluster1->setNRawHits(2);
    fstCluster1->setNRawHitsR(2);
    fstCluster1->setNRawHitsPhi(1);
    fstCluster1->addRawHit(fstRawHit1);
    fstCluster1->addRawHit(fstRawHit2);
    cout << "1st cluster =>" << endl;
    fstCluster1->Print();

    // set hit array
    fstCluster1->clearRawHitsList();
    FstRawHit *fstRawHit1_cluster = fstCluster1->createRawHit();
    fstRawHit1_cluster->setLayer(fstRawHit1->getLayer());

    FstRawHit *fstRawHit2_cluster = fstCluster1->createRawHit();
    fstRawHit2_cluster->setLayer(fstRawHit2->getLayer());
    Tree_FstCluster->Fill();
  }
  Tree_FstCluster->Write();
  File_OutPut->Close();


  return 1;
}
