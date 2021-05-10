#include "./FstQAStudy.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include <TFile.h>
#include <TChain.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TVector3.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TGraph.h>
#include <vector>

using namespace std;

const int eventDivisions = 31;
const int eventRangeSize = 1000;

ClassImp(FstQAStudy)

//------------------------------------------
FstQAStudy::FstQAStudy() : mList("../../list/FST/FstPed_HV140.list"), mOutPutFile("./FstPed_HV140.root"), mApplyCMNCorr(true)
{
  cout << "FstQAStudy::FstQAStudy() -------- Constructor!  --------" << endl;
  mHome = getenv("HOME");
}

FstQAStudy::~FstQAStudy()
{
  cout << "FstQAStudy::~FstQAStudy() -------- Release Memory!  --------" << endl;
}
//------------------------------------------
int FstQAStudy::Init()
{
  cout << "FstQAStudy::Init => " << endl;
  File_mOutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

  bool isInPut = initChain(); // initialize input data/ped TChain;
  initSignalQA();
  initEventDisplay_TrackClusters();

  if(!isInPut) 
  {
    cout << "Failed to initialize input data!" << endl;
    return -1;
  }

  return 1;
}

int FstQAStudy::Make()
{
  cout << "FstQAStudy::Make => " << endl;

  long NumOfEvents = (long)mChainInPut->GetEntries();
  // if(NumOfEvents > 1000) NumOfEvents = 1000;
  // NumOfEvents = 1000;
  mChainInPut->GetEntry(0);

  int numOfUsedEvent = 0;

//Creating Tree with the total events which have hits. Used in laser scan method
  std::vector<Double_t> signalEv;
  std::vector<Double_t> noiseEv;
  std::vector<Int_t> maxTbEv;
  std::vector<Double_t> cmnNoiseEv;
  std::vector<Int_t> treeEventId;
  std::vector<Int_t> treeSensorId;
  std::vector<Int_t> treeRowNum;
  std::vector<Int_t> treeColNum;
  

  Event_Signal = new TTree("Event_Signal", "Event_Signal");
  Event_Signal->Branch("Branch1",  &(signalEv));
  Event_Signal->Branch("Branch2",  &(noiseEv));
  Event_Signal->Branch("Branch3",  &(maxTbEv)); 
  Event_Signal->Branch("Branch4",  &(cmnNoiseEv));
  Event_Signal->Branch("Branch5",  &(treeEventId));
  Event_Signal->Branch("Branch6",  &(treeSensorId));
  Event_Signal->Branch("Branch7",  &(treeRowNum));
  Event_Signal->Branch("Branch8",  &(treeColNum));
  
  Event_Signal->SetAutoSave(50000000);

 
  for(int i_event = 0; i_event < NumOfEvents; ++i_event)
  {
    if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << NumOfEvents << endl;
    mChainInPut->GetEntry(i_event);

    fillSignalQA(mFstEvent_InPut, i_event, treeEventId, maxTbEv, treeSensorId, treeColNum, treeRowNum, signalEv, noiseEv, cmnNoiseEv);
    fillEventDisplay_TrackClusters(mFstEvent_InPut);

    }
 
   
  cout << "processed events:  " << NumOfEvents << "/" << NumOfEvents << endl;

  return 1;
}

int FstQAStudy::Finish()
{
  cout << "FstQAStudy::Finish => " << endl;
  Event_Signal->Write();
  writeSignalQA();
  writeEventDisplay_TrackClusters();

  return 1;
}

// init Input TChain
bool FstQAStudy::initChain()
{
  cout << "FstQAStudy::initChain -> " << endl;

  mChainInPut = new TChain("mTree_FstEvent");

  if (!mList.empty())   // if input file is ok
  {
    cout << "Open input probability file list" << endl;
    ifstream in(mList.c_str());  // input stream
    if(in)
    {
      cout << "input file probability list is ok" << endl;
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  string addfile;
	  addfile = str;
	  mChainInPut->AddFile(addfile.c_str(),-1,"mTree_FstEvent");
	  Long64_t file_entries = mChainInPut->GetEntries();
	  cout << "File added to data chain: " << addfile.c_str() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      cout << "WARNING: input probability file input is problemtic" << endl;
      return false;
    }
  }

  mFstEvent_InPut = new FstEvent();
  mChainInPut->SetBranchAddress("FstEvent",&mFstEvent_InPut);

  long NumOfEvents = (long)mChainInPut->GetEntries();
  cout << "total number of events: " << NumOfEvents << endl;

  return true;
}
// init Input TChain

//--------------Signal QA---------------------
void FstQAStudy::initSignalQA()
{
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    string HistName;
    HistName = Form("p_mFstPedMap_Sensor%d",i_sensor);
    p_mFstPedMap[i_sensor] = new TProfile2D(HistName.c_str(),HistName.c_str(),FST::numRStrip,-0.5,FST::numRStrip-0.5,FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);
    HistName = Form("p_mFstSigMap_Sensor%d",i_sensor);
    p_mFstSigMap[i_sensor] = new TProfile2D(HistName.c_str(),HistName.c_str(),FST::numRStrip,-0.5,FST::numRStrip-0.5,FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);

    HistName = Form("h_mFstHitsMaxTb_Sensor%d",i_sensor);
    h_mFstHitsMaxTb[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),10,-0.5,9.5);
    HistName = Form("h_mFstHitsSignal_Sensor%d",i_sensor);
    h_mFstHitsSignal[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,1999.5);
    HistName = Form("h_mFstHitsNoise_Sensor%d",i_sensor);
    h_mFstHitsNoise[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),100,-0.5,99.5);
    HistName = Form("h_mFstHitsSNRatio_Sensor%d",i_sensor);
    h_mFstHitsSNRatio[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,99.5);

    HistName = Form("h_mFstSimpleClustersSignal_Sensor%d",i_sensor);
    h_mFstSimpleClustersSignal[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,1999.5);
    HistName = Form("h_mFstSimpleClustersNoise_Sensor%d",i_sensor);
    h_mFstSimpleClustersNoise[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),100,-0.5,99.5);
    HistName = Form("h_mFstSimpleClustersSNRatio_Sensor%d",i_sensor);
    h_mFstSimpleClustersSNRatio[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,99.5);

    HistName = Form("h_mFstScanClustersSignal_Sensor%d",i_sensor);
    h_mFstScanClustersSignal[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,1999.5);
    HistName = Form("h_mFstScanClustersNoise_Sensor%d",i_sensor);
    h_mFstScanClustersNoise[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),100,-0.5,99.5);
    HistName = Form("h_mFstScanClustersSNRatio_Sensor%d",i_sensor);
    h_mFstScanClustersSNRatio[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,99.5);

    for(int i_event = 0; i_event < eventDivisions; i_event++) {
      HistName = Form("h_mFstLaser_HitsSignal_Sensor%d_Event%d", i_sensor, i_event);
      h_mFstLaser_HitsSignal[i_sensor][i_event] = new TH3F(HistName.c_str(),HistName.c_str(),8, 0, FST::numRStrip, 128, 0, FST::numPhiSeg, 200,-0.5,4800.5);
      HistName = Form("h_mFstLaser_NoiseSignal_Sensor%d_Event%d",i_sensor, i_event);
      h_mFstLaser_NoiseSignal[i_sensor][i_event] = new TH3F(HistName.c_str(),HistName.c_str(),8, 0, FST::numRStrip, 128, 0, FST::numPhiSeg, 200,-0.5,40.5);
      HistName = Form("h_mFstLaser_SignalRatio_Sensor%d_Event%d", i_sensor, i_event);
      h_mFstLaser_SignalRatio[i_sensor][i_event] = new TH3F(HistName.c_str(),HistName.c_str(),8, 0, FST::numRStrip, 128, 0, FST::numPhiSeg, 200,-0.5,400.5);
      HistName = Form("h_mFstLaser_CmnNoiseSignal_Sensor%d_Event%d",i_sensor, i_event);
      h_mFstLaser_cmnNoiseSignal[i_sensor][i_event] = new TH3F(HistName.c_str(),HistName.c_str(), 8, 0, FST::numRStrip, 128, 0, FST::numPhiSeg, 200,-100.5,100.5);
      HistName = Form("h_mFstLaser_Timebin_Sensor%d_Event%d", i_sensor, i_event);
      h_mFstLaser_Timebin[i_sensor][i_event] = new TH3F(HistName.c_str(),HistName.c_str(),8, 0, FST::numRStrip, 128, 0, FST::numPhiSeg, FST::numTBins,0,FST::numTBins);

      HistName = Form("h_mFstLaser_SigTime_Sensor%d_Event%d", i_sensor, i_event);
      h_mFstLaser_SigTime[i_sensor][i_event] = new TH2F(HistName.c_str(),HistName.c_str(), FST::numTBins,0,FST::numTBins, 200, -0.5,4800.5);
      HistName = Form("h_mFstLaser_NoiseTime_Sensor%d_Event%d", i_sensor, i_event);
      h_mFstLaser_NoiseTime[i_sensor][i_event] = new TH2F(HistName.c_str(),HistName.c_str(), FST::numTBins,0,FST::numTBins, 200, -0.5,40.5);

    }

    for(int i_rstrip = 0; i_rstrip < FST::mFstNumRstripPerSensor; ++i_rstrip)
    {
      HistName = Form("h_mFstHitsMaxTb_Sensor%d_Rstrip%d",i_sensor,i_rstrip);
      if(i_sensor > 0) HistName = Form("h_mFstHitsMaxTb_Sensor%d_Rstrip%d",i_sensor,i_rstrip+4);
      h_mFstHitsMaxTb_Rstrip[i_sensor][i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),10,-0.5,9.5);
      HistName = Form("h_mFstSimpleClustersMaxTb_Sensor%d_Rstrip%d",i_sensor,i_rstrip);
      if(i_sensor > 0) HistName = Form("h_mFstSimpleClustersMaxTb_Sensor%d_Rstrip%d",i_sensor,i_rstrip+4);
      h_mFstSimpleClustersMaxTb_Rstrip[i_sensor][i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),10,-0.5,9.5);
      HistName = Form("h_mFstScanClustersMaxTb_Sensor%d_Rstrip%d",i_sensor,i_rstrip);
      if(i_sensor > 0) HistName = Form("h_mFstScanClustersMaxTb_Sensor%d_Rstrip%d",i_sensor,i_rstrip+4);
      h_mFstScanClustersMaxTb_Rstrip[i_sensor][i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),10,-0.5,9.5);
      HistName = Form("h_mFstHitsSignal_Sensor%d_Rstrip%d",i_sensor,i_rstrip);
      if(i_sensor > 0) HistName = Form("h_mFstHitsSignal_Sensor%d_Rstrip%d",i_sensor,i_rstrip+4);
      h_mFstHitsSignal_Rstrip[i_sensor][i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,1999.5);
      HistName = Form("h_mFstHitsNoise_Sensor%d_Rstrip%d",i_sensor,i_rstrip);
      if(i_sensor > 0) HistName = Form("h_mFstHitsNoise_Sensor%d_Rstrip%d",i_sensor,i_rstrip+4);
      h_mFstHitsNoise_Rstrip[i_sensor][i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),100,-0.5,99.5);
      HistName = Form("h_mFstHitsSNRatio_Sensor%d_Rstrip%d",i_sensor,i_rstrip);
      if(i_sensor > 0) HistName = Form("h_mFstHitsSNRatio_Sensor%d_Rstrip%d",i_sensor,i_rstrip+4);
      h_mFstHitsSNRatio_Rstrip[i_sensor][i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,99.5);
      HistName = Form("h_mFstSimpleClustersSignalSensor%d__Rstrip%d",i_sensor,i_rstrip);
      if(i_sensor > 0) HistName = Form("h_mFstSimpleClustersSignalSensor%d__Rstrip%d",i_sensor,i_rstrip+4);
      h_mFstSimpleClustersSignal_Rstrip[i_sensor][i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,1999.5);
      HistName = Form("h_mFstScanClustersSignal_Sensor%d_Rstrip%d",i_sensor,i_rstrip);
      if(i_sensor > 0) HistName = Form("h_mFstScanClustersSignal_Sensor%d_Rstrip%d",i_sensor,i_rstrip+4);
      h_mFstScanClustersSignal_Rstrip[i_sensor][i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,1999.5);
      for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
      {
	HistName = Form("h_mFstHitsSignal_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip,i_tb);
	if(i_sensor > 0) HistName = Form("h_mFstHitsSignal_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip+4,i_tb);
	h_mFstHitsSignal_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,1999.5);
	HistName = Form("h_mFstHitsNoise_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip,i_tb);
	if(i_sensor > 0) HistName = Form("h_mFstHitsNoise_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip+4,i_tb);
	h_mFstHitsNoise_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),100,-0.5,99.5);
	HistName = Form("h_mFstHitsSNRatio_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip,i_tb);
	if(i_sensor > 0) HistName = Form("h_mFstHitsSNRatio_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip+4,i_tb);
	h_mFstHitsSNRatio_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,99.5);
	HistName = Form("h_mFstSimpleClustersSignal_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip,i_tb);
	if(i_sensor > 0) HistName = Form("h_mFstSimpleClustersSignal_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip+4,i_tb);
	h_mFstSimpleClustersSignal_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,1999.5);
	HistName = Form("h_mFstScanClustersSignal_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip,i_tb);
	if(i_sensor > 0) HistName = Form("h_mFstScanClustersSignal_Sensor%d_Rstrip%d_TimeBin%d",i_sensor,i_rstrip+4,i_tb);
	h_mFstScanClustersSignal_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,1999.5);
      }
    }
  }
  
}

void FstQAStudy::fillSignalQA(FstEvent *fstEvent,int event, std::vector<Int_t>& eventIdVec, std::vector<Int_t>& maxTbVec, std::vector<Int_t>& sensorIdVec, std::vector<Int_t>& colNumVec, std::vector<Int_t>& rowNumVec,std::vector<Double_t>& signalVec, std::vector<Double_t>& noiseVec, std::vector<Double_t>& cmnNoiseVec)
{
  std::vector<FstRawHit *> fstRawHitVec;
  fstRawHitVec.clear(); // clear the container for hits

    
  for(int i_hit = 0; i_hit < fstEvent->getNumRawHits(); ++i_hit)
  { // get Hits info
    FstRawHit *fstRawHit = fstEvent->getRawHit(i_hit);
    if(fstRawHit->getLayer() == 3)
    {
      fstRawHitVec.push_back(fstRawHit);
    }
  }
  int numOfFstHits = fstEvent->getNumFstRawHits();
    
 
  if(numOfFstHits > 0)
  {
    for(int i_hit = 0; i_hit < numOfFstHits; ++i_hit)
    {
      const int sensorId  = fstRawHitVec[i_hit]->getSensor();
      const int apv       = fstRawHitVec[i_hit]->getAPV();
      const int channel   = fstRawHitVec[i_hit]->getChannel();
      int column          = fstRawHitVec[i_hit]->getColumn(); // RStrip
      const int row       = fstRawHitVec[i_hit]->getRow(); // phi segmentation
      const int maxTb     = fstRawHitVec[i_hit]->getMaxTb();
      const double ped    = fstRawHitVec[i_hit]->getPedMean(maxTb); // pedMean
      const double signal = fstRawHitVec[i_hit]->getCharge(maxTb); // adc - pedMean - cmn
      double noise        = fstRawHitVec[i_hit]->getPedStdDev(maxTb); // pedStdDev
      double cmnNoise     = fstRawHitVec[i_hit]->getCMNStdDev(maxTb);
      int eventId	      = fstEvent->getEventId();
 

     if(mApplyCMNCorr) noise  = fstRawHitVec[i_hit]->getRanStdDev(maxTb); // pedStdDev
      p_mFstPedMap[sensorId]->Fill(column,row,ped);
      p_mFstSigMap[sensorId]->Fill(column,row,signal);
      if(fstRawHitVec[i_hit]->getIsHit())
      {
	h_mFstHitsMaxTb[sensorId]->Fill(maxTb);
	h_mFstHitsSignal[sensorId]->Fill(signal);
	h_mFstHitsNoise[sensorId]->Fill(noise);
	h_mFstHitsSNRatio[sensorId]->Fill(signal/noise);

          
	signalVec.push_back(signal);
	noiseVec.push_back(noise);
    cmnNoiseVec.push_back(cmnNoise);
	maxTbVec.push_back(maxTb);
    eventIdVec.push_back(eventId);
    sensorIdVec.push_back(sensorId);
    rowNumVec.push_back(row);
	colNumVec.push_back(column);

    //Fill Histograms used for Laser data taking with manual method
	for (int j_event = 0; j_event < eventDivisions; j_event++) {
	  if ( eventId >= (j_event*eventRangeSize + 1000) && eventId <= ((j_event+1) * eventRangeSize + 1000) ) {

	    h_mFstLaser_cmnNoiseSignal[sensorId][j_event]->Fill(column,row,cmnNoise);
            h_mFstLaser_HitsSignal[sensorId][j_event]->Fill(column, row, sigal);
            h_mFstLaser_SignalRatio[sensorId][j_event]->Fill(column, row, signal/noise);
	    h_mFstLaser_NoiseSignal[sensorId][j_event]->Fill(column, row, noise);
            h_mFstLaser_Timebin[sensorId][j_event]->Fill(column, row, maxTb);
	    h_mFstLaser_SigTime[sensorId][j_event]->Fill(maxTb, signal);
	    h_mFstLaser_NoiseTime[sensorId][j_event]->Fill(maxTb, noise);

  	  }
	}
          
	if(sensorId > 0) column = column-4;

	h_mFstHitsMaxTb_Rstrip[sensorId][column]->Fill(maxTb);
	h_mFstHitsSignal_Rstrip[sensorId][column]->Fill(signal);
	h_mFstHitsNoise_Rstrip[sensorId][column]->Fill(noise);
	h_mFstHitsSNRatio_Rstrip[sensorId][column]->Fill(signal/noise);

	h_mFstHitsSignal_Rstrip_TimeBin[sensorId][column][maxTb]->Fill(signal);
	h_mFstHitsNoise_Rstrip_TimeBin[sensorId][column][maxTb]->Fill(noise);
	h_mFstHitsSNRatio_Rstrip_TimeBin[sensorId][column][maxTb]->Fill(signal/noise);
      }
    }
  }
  Event_Signal->Fill();

  std::vector<FstCluster *> fstClusterVec;
  fstClusterVec.clear(); // clear the container for clusters
  for(int i_cluster = 0; i_cluster < fstEvent->getNumClusters(); ++i_cluster)
  { // get Clusters info
    FstCluster *fstCluster = fstEvent->getCluster(i_cluster);
    if(fstCluster->getLayer() == 3)
    {
      fstClusterVec.push_back(fstCluster);
    }
  }
  int numOfFstClusters = fstClusterVec.size();
  if(numOfFstClusters > 0)
  {
    for(int i_cluster = 0; i_cluster < numOfFstClusters; ++i_cluster)
    {
      int sensorId = fstClusterVec[i_cluster]->getSensor();
      double meanColumn = fstClusterVec[i_cluster]->getMeanColumn();
      double meanRow = fstClusterVec[i_cluster]->getMeanRow();
      double signal = fstClusterVec[i_cluster]->getTotCharge();
      double maxTb = fstClusterVec[i_cluster]->getMaxTb();
      if(fstClusterVec[i_cluster]->getClusterType() == 1 && fstClusterVec[i_cluster]->getIsSeed())
      { // Simple Clusters
	double noise = 0.0;
	for(int i_hit = 0; i_hit < fstClusterVec[i_cluster]->getNumRawHits(); ++i_hit)
	{ // fill hits display for clusters with hit
	  FstRawHit *fstRawHit = fstClusterVec[i_cluster]->getRawHit(i_hit);
	  int maxTb_hit = fstRawHit->getMaxTb();
	  double ped_hit = fstRawHit->getPedStdDev(maxTb_hit); // pedStdDev
	  if(mApplyCMNCorr) ped_hit = fstRawHit->getRanStdDev(maxTb_hit); // ranStdDev
	  noise += ped_hit*ped_hit;
	}
 
	h_mFstSimpleClustersSignal[sensorId]->Fill(signal);
	h_mFstSimpleClustersNoise[sensorId]->Fill(TMath::Sqrt(noise));
	h_mFstSimpleClustersSNRatio[sensorId]->Fill(signal/TMath::Sqrt(noise));


	int column = -1;
	if(meanColumn > -0.5 && meanColumn <= 0.5) column = 0;
	if(meanColumn >  0.5 && meanColumn <= 1.5) column = 1;
	if(meanColumn >  1.5 && meanColumn <= 2.5) column = 2;
	if(meanColumn >  2.5 && meanColumn <= 3.5) column = 3;
	if(meanColumn >  3.5 && meanColumn <= 4.5) column = 0;
	if(meanColumn >  4.5 && meanColumn <= 5.5) column = 1;
	if(meanColumn >  5.5 && meanColumn <= 6.5) column = 2;
	if(meanColumn >  6.5 && meanColumn <= 7.5) column = 3;
	if(column > -1)
	{
	  h_mFstSimpleClustersMaxTb_Rstrip[sensorId][column]->Fill(maxTb);
	  h_mFstSimpleClustersSignal_Rstrip[sensorId][column]->Fill(signal);
	  h_mFstSimpleClustersSignal_Rstrip_TimeBin[sensorId][column][(int)maxTb]->Fill(signal);
	}
      }
      if(fstClusterVec[i_cluster]->getClusterType() == 2 && fstClusterVec[i_cluster]->getIsSeed())
      { // Scan Clusters
	double noise = 0.0;
	for(int i_hit = 0; i_hit < fstClusterVec[i_cluster]->getNumRawHits(); ++i_hit)
	{ // fill hits display for clusters with hit
	  FstRawHit *fstRawHit = fstClusterVec[i_cluster]->getRawHit(i_hit);
	  int maxTb_hit = fstRawHit->getMaxTb();
	  double ped_hit = fstRawHit->getPedStdDev(maxTb_hit); // pedStdDev
	  if(mApplyCMNCorr) ped_hit = fstRawHit->getRanStdDev(maxTb_hit); // ranStdDev
	  noise += ped_hit*ped_hit;
	}
	h_mFstScanClustersSignal[sensorId]->Fill(signal);
	h_mFstScanClustersNoise[sensorId]->Fill(TMath::Sqrt(noise));
	h_mFstScanClustersSNRatio[sensorId]->Fill(signal/TMath::Sqrt(noise));

	int column = -1;
	if(meanColumn > -0.5 && meanColumn <= 0.5) column = 0;
	if(meanColumn >  0.5 && meanColumn <= 1.5) column = 1;
	if(meanColumn >  1.5 && meanColumn <= 2.5) column = 2;
	if(meanColumn >  2.5 && meanColumn <= 3.5) column = 3;
	if(meanColumn >  3.5 && meanColumn <= 4.5) column = 0;
	if(meanColumn >  4.5 && meanColumn <= 5.5) column = 1;
	if(meanColumn >  5.5 && meanColumn <= 6.5) column = 2;
	if(meanColumn >  6.5 && meanColumn <= 7.5) column = 3;
	if(column > -1)
	{
	  h_mFstScanClustersMaxTb_Rstrip[sensorId][column]->Fill(maxTb);
	  h_mFstScanClustersSignal_Rstrip[sensorId][column]->Fill(signal);
	  h_mFstScanClustersSignal_Rstrip_TimeBin[sensorId][column][(int)maxTb]->Fill(signal);
	}
      }
    }
  }
}

void FstQAStudy::writeSignalQA()
{
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    p_mFstPedMap[i_sensor]->Write();
    p_mFstSigMap[i_sensor]->Write();

    h_mFstHitsMaxTb[i_sensor]->Write();
    h_mFstHitsSignal[i_sensor]->Write();
    h_mFstHitsNoise[i_sensor]->Write();
    h_mFstHitsSNRatio[i_sensor]->Write();
    h_mFstSimpleClustersSignal[i_sensor]->Write();
    h_mFstSimpleClustersNoise[i_sensor]->Write();
    h_mFstSimpleClustersSNRatio[i_sensor]->Write();
    h_mFstScanClustersSignal[i_sensor]->Write();
    h_mFstScanClustersNoise[i_sensor]->Write();
    h_mFstScanClustersSNRatio[i_sensor]->Write();

    for(int j_event = 0; j_event < eventDivisions; j_event++) {
      h_mFstLaser_HitsSignal[i_sensor][j_event]->Write();
      h_mFstLaser_NoiseSignal[i_sensor][j_event]->Write();
      h_mFstLaser_SignalRatio[i_sensor][j_event]->Write();
      h_mFstLaser_cmnNoiseSignal[i_sensor][j_event]->Write();
      h_mFstLaser_Timebin[i_sensor][j_event]->Write();
      h_mFstLaser_SigTime[i_sensor][j_event]->Write();
      h_mFstLaser_NoiseTime[i_sensor][j_event]->Write();

    }

    for(int i_rstrip = 0; i_rstrip < FST::mFstNumRstripPerSensor; ++i_rstrip)
    {
      h_mFstHitsMaxTb_Rstrip[i_sensor][i_rstrip]->Write();
      h_mFstSimpleClustersMaxTb_Rstrip[i_sensor][i_rstrip]->Write();
      h_mFstScanClustersMaxTb_Rstrip[i_sensor][i_rstrip]->Write();
      h_mFstHitsSignal_Rstrip[i_sensor][i_rstrip]->Write();
      h_mFstHitsNoise_Rstrip[i_sensor][i_rstrip]->Write();
      h_mFstHitsSNRatio_Rstrip[i_sensor][i_rstrip]->Write();
      h_mFstSimpleClustersSignal_Rstrip[i_sensor][i_rstrip]->Write();
      h_mFstScanClustersSignal_Rstrip[i_sensor][i_rstrip]->Write();
      for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
      {
	h_mFstHitsSignal_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb]->Write();
	h_mFstHitsNoise_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb]->Write();
	h_mFstHitsSNRatio_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb]->Write();
	h_mFstSimpleClustersSignal_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb]->Write();
	h_mFstScanClustersSignal_Rstrip_TimeBin[i_sensor][i_rstrip][i_tb]->Write();
      }
    }
  }
}
//--------------Signal QA---------------------

//--------------Event Display---------------------
void FstQAStudy::initEventDisplay_TrackClusters()
{
  const int numRDisplay = 10;
  const double rMinDisplay = FST::rInner - 1.0*FST::pitchR;
  const double rMaxDisplay = FST::rOuter + 5.0*FST::pitchR;
  const int numPhiDisplay = 2*FST::numPhiSeg;
  const double phiMinDisplay = -TMath::Pi()*30.0/180.0 - FST::gapPhi;
  const double phiMaxDisplay =  TMath::Pi()*30.0/180.0 + FST::gapPhi;
  h_mFstRawHitsDisplay              = new TH2F("h_mFstRawHitsDisplay","h_mFstRawHitsDisplay",numRDisplay,rMinDisplay,rMaxDisplay,numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstRawPedsDisplay              = new TH2F("h_mFstRawPedsDisplay","h_mFstRawPedsDisplay",numRDisplay,rMinDisplay,rMaxDisplay,numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstMaxTbDisplay                = new TH2F("h_mFstMaxTbDisplay","h_mFstMaxTbDisplay",numRDisplay,rMinDisplay,rMaxDisplay,numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstRawHitsDisplay_bTh          = new TH2F("h_mFstRawHitsDisplay_bTh","h_mFstRawHitsDisplay_bTh",numRDisplay,rMinDisplay,rMaxDisplay,numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstRawPedsDisplay_bTh          = new TH2F("h_mFstRawPedsDisplay_bTh","h_mFstRawPedsDisplay_bTh",numRDisplay,rMinDisplay,rMaxDisplay,numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstMaxTbDisplay_bTh            = new TH2F("h_mFstMaxTbDisplay_bTh","h_mFstMaxTbDisplay_bTh",numRDisplay,rMinDisplay,rMaxDisplay,numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstSimpleClustersDisplay       = new TH2F("h_mFstSimpleClustersDisplay","h_mFstSimpleClustersDisplay",10*numRDisplay,rMinDisplay,rMaxDisplay,2*numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstScanClustersDisplay         = new TH2F("h_mFstScanClustersDisplay","h_mFstScanClustersDisplay",10*numRDisplay,rMinDisplay,rMaxDisplay,2*numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstScanClustersDisplay_RawHits = new TH2F("h_mFstScanClustersDisplay_RawHits","h_mFstScanClustersDisplay_RawHits",numRDisplay,rMinDisplay,rMaxDisplay,numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstScanClustersDisplay_RawPeds = new TH2F("h_mFstScanClustersDisplay_RawPeds","h_mFstScanClustersDisplay_RawPeds",numRDisplay,rMinDisplay,rMaxDisplay,numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mFstScanClustersDisplay_MaxTb   = new TH2F("h_mFstScanClustersDisplay_MaxTb","h_mFstScanClustersDisplay_MaxTb",numRDisplay,rMinDisplay,rMaxDisplay,numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mHitTracksDisplay               = new TH2F("h_mHitTracksDisplay","h_mHitTracksDisplay",10*numRDisplay,rMinDisplay,rMaxDisplay,2*numPhiDisplay,phiMinDisplay,phiMaxDisplay);
  h_mClusterTracksDisplay           = new TH2F("h_mClusterTracksDisplay","h_mClusterTracksDisplay",10*numRDisplay,rMinDisplay,rMaxDisplay,2*numPhiDisplay,phiMinDisplay,phiMaxDisplay);

  g_mFstSimpleClustersDisplay = new TGraph();
  g_mFstSimpleClustersDisplay->SetName("g_mFstSimpleClustersDisplay");
  g_mFstScanClustersDisplay = new TGraph();
  g_mFstScanClustersDisplay->SetName("g_mFstScanClustersDisplay");
  g_mHitTracksDisplay = new TGraph();
  g_mHitTracksDisplay->SetName("g_mHitTracksDisplay");
  g_mClusterTracksDisplay = new TGraph();
  g_mClusterTracksDisplay->SetName("g_mClusterTracksDisplay");

  mTree_EventDisplay = new TTree("mTree_EventDisplay","Fst Hits Clusters Tracks Display");
  mTree_EventDisplay->Branch("mEventId",&mEventId,"mEventId/I");
  mTree_EventDisplay->Branch("mNumOfFstRawHits",&mNumOfFstRawHits,"mNumOfFstRawHits/I");
  mTree_EventDisplay->Branch("mNumOfIst1RawHits",&mNumOfIst1RawHits,"mNumOfIst1RawHits/I");
  mTree_EventDisplay->Branch("mNumOfIst2RawHits",&mNumOfIst2RawHits,"mNumOfIst2RawHits/I");
  mTree_EventDisplay->Branch("mNumOfIst3RawHits",&mNumOfIst3RawHits,"mNumOfIst3RawHits/I");
  mTree_EventDisplay->Branch("h_mFstRawHitsDisplay","TH2F",&h_mFstRawHitsDisplay);
  mTree_EventDisplay->Branch("h_mFstRawPedsDisplay","TH2F",&h_mFstRawPedsDisplay);
  mTree_EventDisplay->Branch("h_mFstMaxTbDisplay","TH2F",&h_mFstMaxTbDisplay);
  mTree_EventDisplay->Branch("h_mFstRawHitsDisplay_bTh","TH2F",&h_mFstRawHitsDisplay_bTh); // hit map below threshold
  mTree_EventDisplay->Branch("h_mFstRawPedsDisplay_bTh","TH2F",&h_mFstRawPedsDisplay_bTh);
  mTree_EventDisplay->Branch("h_mFstMaxTbDisplay_bTh","TH2F",&h_mFstMaxTbDisplay_bTh);

  mTree_EventDisplay->Branch("mNumOfFstSimpleClusters",&mNumOfFstSimpleClusters,"mNumOfFstSimpleClusters/I");
  mTree_EventDisplay->Branch("mNumOfFstScanClusters",&mNumOfFstScanClusters,"mNumOfFstScanClusters/I");
  mTree_EventDisplay->Branch("mNumOfIst1Clusters",&mNumOfIst1Clusters,"mNumOfIst1Clusters/I");
  mTree_EventDisplay->Branch("mNumOfIst2Clusters",&mNumOfIst2Clusters,"mNumOfIst2Clusters/I");
  mTree_EventDisplay->Branch("mNumOfIst3Clusters",&mNumOfIst3Clusters,"mNumOfIst3Clusters/I");
  mTree_EventDisplay->Branch("h_mFstSimpleClustersDisplay","TH2F",&h_mFstSimpleClustersDisplay);
  mTree_EventDisplay->Branch("h_mFstScanClustersDisplay","TH2F",&h_mFstScanClustersDisplay);
  mTree_EventDisplay->Branch("h_mFstScanClustersDisplay_RawHits","TH2F",&h_mFstScanClustersDisplay_RawHits);
  mTree_EventDisplay->Branch("h_mFstScanClustersDisplay_RawPeds","TH2F",&h_mFstScanClustersDisplay_RawPeds);
  mTree_EventDisplay->Branch("h_mFstScanClustersDisplay_MaxTb","TH2F",&h_mFstScanClustersDisplay_MaxTb);

  mTree_EventDisplay->Branch("mNumOfHitTracks",&mNumOfHitTracks,"mNumOfHitTracks/I");
  mTree_EventDisplay->Branch("h_mHitTracksDisplay","TH2F",&h_mHitTracksDisplay);
  mTree_EventDisplay->Branch("mNumOfClusterTracks",&mNumOfClusterTracks,"mNumOfClusterTracks/I");
  mTree_EventDisplay->Branch("mNumOfClusterTracks_2Layer",&mNumOfClusterTracks_2Layer,"mNumOfClusterTracks_2Layer/I");
  mTree_EventDisplay->Branch("mNumOfClusterTracks_3Layer",&mNumOfClusterTracks_3Layer,"mNumOfClusterTracks_3Layer/I");
  mTree_EventDisplay->Branch("h_mClusterTracksDisplay","TH2F",&h_mClusterTracksDisplay);

  mTree_EventDisplay->Branch("g_mFstSimpleClustersDisplay","TGraph",&g_mFstSimpleClustersDisplay);
  mTree_EventDisplay->Branch("g_mFstScanClustersDisplay","TGraph",&g_mFstScanClustersDisplay);
  mTree_EventDisplay->Branch("g_mHitTracksDisplay","TGraph",&g_mHitTracksDisplay);
  mTree_EventDisplay->Branch("g_mClusterTracksDisplay","TGraph",&g_mClusterTracksDisplay);

  mTree_EventDisplay->SetAutoSave(50000000);

  h_mNumFstSimpleClusters = new TH1F("h_mNumFstSimpleClusters","h_mNumFstSimpleClusters",100,-0.5,99.5);
  h_mNumFstRawHitsSimpleCluster = new TH1F("h_mNumFstRawHitsSimpleCluster","h_mNumFstRawHitsSimpleCluster",50,-0.5,49.5);
  h_mNumFstScanClusters = new TH1F("h_mNumFstScanClusters","h_mNumFstScanClusters",100,-0.5,99.5);
  h_mNumFstRawHitsScanCluster = new TH1F("h_mNumFstRawHitsScanCluster","h_mNumFstRawHitsScanCluster",50,-0.5,49.5);
}

void FstQAStudy::clearEventDisplay_TrackClusters()
{
  mEventId = -1;
  mNumOfFstRawHits = 0;
  mNumOfIst1RawHits = 0;
  mNumOfIst2RawHits = 0;
  mNumOfIst3RawHits = 0;

  mNumOfFstSimpleClusters = 0;
  mNumOfFstScanClusters = 0;
  mNumOfIst1Clusters = 0;
  mNumOfIst2Clusters = 0;
  mNumOfIst3Clusters = 0;

  mNumOfHitTracks = 0;
  mNumOfClusterTracks = 0;
  mNumOfClusterTracks_2Layer = 0;
  mNumOfClusterTracks_3Layer = 0;

  h_mFstRawHitsDisplay->Reset();
  h_mFstRawPedsDisplay->Reset();
  h_mFstMaxTbDisplay->Reset();
  h_mFstRawHitsDisplay_bTh->Reset();
  h_mFstRawPedsDisplay_bTh->Reset();
  h_mFstMaxTbDisplay_bTh->Reset();

  h_mFstSimpleClustersDisplay->Reset();
  h_mFstScanClustersDisplay->Reset();
  h_mFstScanClustersDisplay_RawHits->Reset();
  h_mFstScanClustersDisplay_RawPeds->Reset();
  h_mFstScanClustersDisplay_MaxTb->Reset();

  h_mHitTracksDisplay->Reset();
  h_mClusterTracksDisplay->Reset();

  g_mFstSimpleClustersDisplay->Set(-1);
  g_mFstScanClustersDisplay->Set(-1);
  g_mHitTracksDisplay->Set(-1);
  g_mClusterTracksDisplay->Set(-1);
}

void FstQAStudy::fillEventDisplay_TrackClusters(FstEvent *fstEvent)
{
  this->clearEventDisplay_TrackClusters();

  mEventId = fstEvent->getEventId();
  // cout << "mEventId = " << mEventId << endl;

  for(int i_hit = 0; i_hit < fstEvent->getNumRawHits(); ++i_hit)
  { // fill Hits Display
    FstRawHit *fstRawHit = fstEvent->getRawHit(i_hit);
    if(fstRawHit->getLayer() == 3) // FST
    {
      int sensorId    = fstRawHit->getSensor();
      double r_orig   = fstRawHit->getPosX(); // r for fst in sensorId frame
      double phi_orig = fstRawHit->getPosY(); // phi for fst in sensorId frame
      TVector2 vPosFst(r_orig,phi_orig); // set origninal pos in sensorId frame
      TVector2 vAlignedFST = this->getFstAlignedPos(vPosFst,sensorId); // convert to mDefSenorId frame
      double r_fst   = vAlignedFST.X(); // r for fst in mDefSenorId frame
      double phi_fst = vAlignedFST.Y(); // phi for fst in mDefSenorId frame
      int maxTb      = fstRawHit->getMaxTb();
      double adc     = fstRawHit->getCharge(maxTb); // adc - pedMean
      double ped     = fstRawHit->getPedStdDev(maxTb); // pedStdDev
      if(mApplyCMNCorr) ped = fstRawHit->getRanStdDev(maxTb); // ranStdDev

      if( fstRawHit->getIsHit() )
      { // above threshold
	mNumOfFstRawHits++;
	h_mFstRawHitsDisplay->Fill(r_fst,phi_fst,adc);
	h_mFstRawPedsDisplay->Fill(r_fst,phi_fst,ped);
	h_mFstMaxTbDisplay->Fill(r_fst,phi_fst,maxTb);
      }
      h_mFstRawHitsDisplay_bTh->Fill(r_fst,phi_fst,adc);
      h_mFstRawPedsDisplay_bTh->Fill(r_fst,phi_fst,ped);
      h_mFstMaxTbDisplay_bTh->Fill(r_fst,phi_fst,maxTb);
    }
    if(fstRawHit->getLayer() == 0) mNumOfIst1RawHits++;
    if(fstRawHit->getLayer() == 1) mNumOfIst2RawHits++;
    if(fstRawHit->getLayer() == 2) mNumOfIst3RawHits++;
  }

  for(int i_cluster = 0; i_cluster < fstEvent->getNumClusters(); ++i_cluster)
  { // fill Clusters Display
    FstCluster *fstCluster = fstEvent->getCluster(i_cluster);
    if(fstCluster->getLayer() == 3 && fstCluster->getClusterType() == 1) // FST Simple Cluster
    {
      if( fstCluster->getIsSeed() ) // select cluste with seed
      {
	mNumOfFstSimpleClusters++;
	h_mNumFstRawHitsSimpleCluster->Fill(fstCluster->getNumRawHits());

	int sensorId    = fstCluster->getSensor();
	double r_orig   = fstCluster->getMeanX(); // r for fst in sensorId frame
	double phi_orig = fstCluster->getMeanY(); // phi for fst in sensorId frame
	TVector2 vPosFst(r_orig,phi_orig); // set origninal pos in sensorId frame
	TVector2 vAlignedFST = this->getFstAlignedPos(vPosFst,sensorId); // convert to mDefSenorId frame
	double r_fst   = vAlignedFST.X(); // r for fst in mDefSenorId frame
	double phi_fst = vAlignedFST.Y(); // phi for fst in mDefSenorId frame
	double adc     = fstCluster->getTotCharge();
	h_mFstSimpleClustersDisplay->Fill(r_fst,phi_fst,adc);
	g_mFstSimpleClustersDisplay->SetPoint(mNumOfFstSimpleClusters-1,r_fst,phi_fst);
      }
    }
    if(fstCluster->getLayer() == 0 && fstCluster->getClusterType() == 1) mNumOfIst1Clusters++;
    if(fstCluster->getLayer() == 1 && fstCluster->getClusterType() == 1) mNumOfIst2Clusters++;
    if(fstCluster->getLayer() == 2 && fstCluster->getClusterType() == 1) mNumOfIst3Clusters++;

    if(fstCluster->getLayer() == 3 && fstCluster->getClusterType() == 2) // FST Scan Cluster
    {
      if( fstCluster->getIsSeed() ) // select cluste with seed
      {
	mNumOfFstScanClusters++;
	h_mNumFstRawHitsScanCluster->Fill(fstCluster->getNumRawHits());

	int sensorId    = fstCluster->getSensor();
	double r_orig   = fstCluster->getMeanX(); // r for fst in sensorId frame
	double phi_orig = fstCluster->getMeanY(); // phi for fst in sensorId frame
	TVector2 vPosFst(r_orig,phi_orig); // set origninal pos in sensorId frame
	TVector2 vAlignedFST = this->getFstAlignedPos(vPosFst,sensorId); // convert to mDefSenorId frame
	double r_fst   = vAlignedFST.X(); // r for fst in mDefSenorId frame
	double phi_fst = vAlignedFST.Y(); // phi for fst in mDefSenorId frame
	double adc     = fstCluster->getTotCharge();
	h_mFstScanClustersDisplay->Fill(r_fst,phi_fst,adc);
	g_mFstScanClustersDisplay->SetPoint(mNumOfFstScanClusters-1,r_fst,phi_fst);

	for(int i_hit = 0; i_hit < fstCluster->getNumRawHits(); ++i_hit)
	{ // fill hits display for clusters with hit
	  FstRawHit *fstRawHit = fstCluster->getRawHit(i_hit);
	  int sensorId        = fstRawHit->getSensor();
	  double r_hit_orig   = fstRawHit->getPosX(); // r for fst in sensorId frame
	  double phi_hit_orig = fstRawHit->getPosY(); // phi for fst in sensorId frame
	  TVector2 vPosHitFst(r_hit_orig,phi_hit_orig); // set origninal pos in sensorId frame
	  TVector2 vAlignedHitFST = this->getFstAlignedPos(vPosHitFst,sensorId); // convert to mDefSenorId frame
	  double r_hit   = vAlignedHitFST.X(); // r for fst in mDefSenorId frame
	  double phi_hit = vAlignedHitFST.Y(); // phi for fst in mDefSenorId frame
	  int maxTb_hit  = fstRawHit->getMaxTb();
	  double adc_hit = fstRawHit->getCharge(maxTb_hit); // adc - pedMean
	  double ped_hit = fstRawHit->getPedStdDev(maxTb_hit); // pedStdDev
	  if(mApplyCMNCorr) ped_hit = fstRawHit->getRanStdDev(maxTb_hit); // ranStdDev

	  h_mFstScanClustersDisplay_RawHits->Fill(r_hit,phi_hit,adc_hit);
	  h_mFstScanClustersDisplay_RawPeds->Fill(r_hit,phi_hit,ped_hit);
	  h_mFstScanClustersDisplay_MaxTb->Fill(r_hit,phi_hit,maxTb_hit);
	}
      }
    }
  }
  h_mNumFstSimpleClusters->Fill(mNumOfFstSimpleClusters);
  h_mNumFstScanClusters->Fill(mNumOfFstScanClusters);

  for(int i_track = 0; i_track < mFstEvent_InPut->getNumTracks(); ++i_track)
  { // fill Tracks Display
    FstTrack *fstTrack = mFstEvent_InPut->getTrack(i_track);

    TVector3 pos_ist1 = fstTrack->getOrigPosIst1();
    double y1_ist = pos_ist1.Y(); // original hit postion on IST1
    TVector3 pos_ist3 = fstTrack->getOrigPosIst3();
    double y3_ist = pos_ist3.Y(); // original hit postion on IST3

    TVector3 proj_fst = fstTrack->getProjFst(FST::mDefSenorId); // use mDefSenorId for now
    double r_proj = proj_fst.X();
    double phi_proj = proj_fst.Y();

    TVector3 proj_ist2 = fstTrack->getProjIst2();
    double x2_proj = proj_ist2.X(); // get aligned projected position w.r.t. IST2
    double y2_proj = proj_ist2.Y(); // x & y for ist2

    if(fstTrack->getTrackType() == 0) // track reconstructed with hits
    {
      mNumOfHitTracks++;
      h_mHitTracksDisplay->Fill(r_proj,phi_proj);
      g_mHitTracksDisplay->SetPoint(mNumOfHitTracks-1,r_proj,phi_proj);
    }
    if(fstTrack->getTrackType() == 1) // track reconstructed with clusters
    {
      mNumOfClusterTracks++;
      h_mClusterTracksDisplay->Fill(r_proj,phi_proj);
      g_mClusterTracksDisplay->SetPoint(mNumOfClusterTracks-1,r_proj,phi_proj);

      if( abs(y1_ist-y3_ist) < 17.0*FST::pitchRow )
      {
	if(r_proj >= FST::rInner-1.0*FST::pitchR && r_proj <= FST::rOuter+5.0*FST::pitchR && phi_proj >= -TMath::Pi()*15.0/180.0-10*FST::pitchPhi && phi_proj <= TMath::Pi()*15.0/180.0+10*FST::pitchPhi)
	{ // save only when in the area of interest
	  mNumOfClusterTracks_2Layer++; // satisfied 2-Layer tracking

	  for(int i_cluster = 0; i_cluster < fstEvent->getNumClusters(); ++i_cluster)
	  { // find IST2 matching
	    FstCluster *fstCluster = fstEvent->getCluster(i_cluster);
	    if(fstCluster->getLayer() == 2) // IST2
	    {
	      double x2_ist = fstCluster->getMeanX(); // x for ist2
	      double y2_ist = fstCluster->getMeanY(); // y for ist2
	      if( x2_proj >= 20.0*FST::pitchColumn && x2_proj <= 24.0*FST::pitchColumn )
	      {
		if( abs(x2_ist-x2_proj) < 6.0 && abs(y2_ist-y2_proj) < 0.6 )
		{ // IST2 matching cut
		  mNumOfClusterTracks_3Layer++;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  mTree_EventDisplay->Fill();
}

void FstQAStudy::writeEventDisplay_TrackClusters()
{
  mTree_EventDisplay->Write();

  h_mNumFstSimpleClusters->Write();
  h_mNumFstRawHitsSimpleCluster->Write();
  h_mNumFstScanClusters->Write();
  h_mNumFstRawHitsScanCluster->Write();
}
//--------------Event Display---------------------

//--------------Utility for tracking-------------------
TVector2 FstQAStudy::getFstAlignedPos(TVector2 vPosFst, int sensorId)
{
  double r_FstLocal   = vPosFst.X();
  double phi_FstLocal = vPosFst.Y();
  double x_FstLocal   = r_FstLocal*TMath::Cos(phi_FstLocal);
  double y_FstLocal   = r_FstLocal*TMath::Sin(phi_FstLocal);
  // transfer to IST2 frame
  double x_IST = (x_FstLocal-FST::x_shift[sensorId])*TMath::Cos(FST::phi_rot[sensorId]) - (y_FstLocal-FST::y_shift[sensorId])*TMath::Sin(FST::phi_rot[sensorId]);
  double y_IST = (y_FstLocal-FST::y_shift[sensorId])*TMath::Cos(FST::phi_rot[sensorId]) + (x_FstLocal-FST::x_shift[sensorId])*TMath::Sin(FST::phi_rot[sensorId]);

  // transfer to FST mDefSenorId frame
  double x_FST = x_IST*TMath::Cos(FST::phi_rot[FST::mDefSenorId]) + y_IST*TMath::Sin(FST::phi_rot[FST::mDefSenorId]) + FST::x_shift[FST::mDefSenorId];
  double y_FST = y_IST*TMath::Cos(FST::phi_rot[FST::mDefSenorId]) - x_IST*TMath::Sin(FST::phi_rot[FST::mDefSenorId]) + FST::y_shift[FST::mDefSenorId];
  double r_FST = TMath::Sqrt(x_FST*x_FST+y_FST*y_FST);
  double phi_FST = TMath::ATan2(y_FST,x_FST);

  TVector2 vAlignedFST;
  vAlignedFST.Set(r_FST,phi_FST);

  return vAlignedFST;
}
//--------------Utility for tracking-------------------
