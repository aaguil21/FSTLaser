#include "./FstNoiseStudy.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TVector3.h>

using namespace std;

ClassImp(FstNoiseStudy)

//------------------------------------------

FstNoiseStudy::FstNoiseStudy() : mList("../../list/FST/FstPed_HV70.list"), mOutPutFile("./FstPed_HV70.root")
{
  cout << "FstNoiseStudy::FstNoiseStudy() -------- Constructor!  --------" << endl;
  mHome = getenv("HOME");
}

FstNoiseStudy::~FstNoiseStudy()
{
  cout << "FstNoiseStudy::~FstNoiseStudy() -------- Release Memory!  --------" << endl;
}

//------------------------------------------
int FstNoiseStudy::Init()
{
  cout << "FstNoiseStudy::Init => " << endl;
  File_mOutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

  bool isInPut = initChain(); // initialize input data/ped TChain;
  bool isPed = initPedestal(); // initialize pedestal array;
  bool isChannelMap = getChannelMap(); // read in channel map

  if(!isInPut) 
  {
    cout << "Init=> Failed to initialize input data!" << endl;
    return -1;
  }
  if(!isPed) 
  {
    cout << "Init=> Failed to initialize pedestals!" << endl;
    return -1;
  }
  if(!isChannelMap)
  {
    cout << "Init=> Failed to read in channel map!" << endl;
    return -1;
  }

  return 1;
}

bool FstNoiseStudy::initChain()
{
  cout << "FstNoiseStudy::initChain -> " << endl;

  mChainInPut = new TChain("tree_Hits");

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
	  mChainInPut->AddFile(addfile.c_str(),-1,"tree_Hits");
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

  mChainInPut->SetBranchAddress("evt_rdo",evt_rdo);
  mChainInPut->SetBranchAddress("evt_arm",evt_arm);
  mChainInPut->SetBranchAddress("evt_port",evt_port);
  mChainInPut->SetBranchAddress("evt_apv",evt_apv);
  mChainInPut->SetBranchAddress("hit_ch",hit_ch);
  mChainInPut->SetBranchAddress("hit_tb",hit_tb);
  mChainInPut->SetBranchAddress("hit_adc",hit_adc);

  long NumOfEvents = (long)mChainInPut->GetEntries();
  cout << "total number of events: " << NumOfEvents << endl;

  return true;
}

//------------------------------------------
int FstNoiseStudy::Make()
{
  cout << "FstNoiseStudy::Make => " << endl;

  bool isPed = calPedestal(); // calculate pedestals with first 1000 events
  if( !isPed ) 
  {
    cout << "No Pedestal!!!" << endl;
    return -1;
  }

  return 1;
}

int FstNoiseStudy::Finish()
{
  cout << "FstNoiseStudy::Finish => " << endl;
  writePedestal();

  return 1;
}

//--------------pedestal---------------------
bool FstNoiseStudy::clearPedestal()
{
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	{
	  for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	  {
	    mPed[i_arm][i_port][i_apv][i_ch][i_tb] = -1.0;
	    mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb] = -1.0;
	    mPedRMS[i_arm][i_port][i_apv][i_ch][i_tb] = -1.0;

	    mCMNStdDev[i_arm][i_port][i_apv][i_ch][i_tb] = -1.0;
	    mRanStdDev[i_arm][i_port][i_apv][i_ch][i_tb] = -1.0;
	  }
	  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
	  {
	    mCMNMean[i_arm][i_port][i_apv][i_rstrip][i_tb] = -1.0;
	  }
	}
      }
    }
  }

  clearCMN();

  return true;
}

bool FstNoiseStudy::clearCMN()
{
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	{
	  for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	  {
	    mCMNStdDev_Evt[i_arm][i_port][i_apv][i_ch][i_tb] = -1.0;
	  }
	  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
	  {
	    mCMNMean_Evt[i_arm][i_port][i_apv][i_rstrip][i_tb] = -1.0;
	  }
	}
      }
    }
  }

  return true;
}

bool FstNoiseStudy::initPedestal()
{
  cout << "FstNoiseStudy::initPedestal -> " << endl;

  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string gName = Form("g_mPedMean_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mPedMean[i_layer][i_tb] = new TGraph();
      g_mPedMean[i_layer][i_tb]->SetName(gName.c_str());

      gName = Form("g_mPedSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mPedSigma[i_layer][i_tb] = new TGraph();
      g_mPedSigma[i_layer][i_tb]->SetName(gName.c_str());

      gName = Form("g_mCMNSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mCMNSigma[i_layer][i_tb] = new TGraph();
      g_mCMNSigma[i_layer][i_tb]->SetName(gName.c_str());

      gName = Form("g_mRanSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRanSigma[i_layer][i_tb] = new TGraph();
      g_mRanSigma[i_layer][i_tb]->SetName(gName.c_str());

      std::string HistName = Form("h_mPedDisplay_Layer%d_TimeBin%d",i_layer,i_tb);
      if(i_layer == 0) h_mPedDisplay[i_layer][i_tb] = new TH2F(HistName.c_str(),HistName.c_str(),FST::numRStrip,-0.5,FST::numRStrip-0.5,FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);
      else h_mPedDisplay[i_layer][i_tb] = new TH2F(HistName.c_str(),HistName.c_str(),FST::noColumns,-0.5,FST::noColumns-0.5,FST::noRows,-0.5,FST::noRows-0.5);
    }
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string HistName;
      HistName = Form("h_mPedMean_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedMean_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);

      HistName = Form("h_mPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedSigma_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);

      HistName = Form("h_mCMNSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mCMNSigma_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);

      HistName = Form("h_mRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mRanSigma_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);

      HistName = Form("h_mDisplayPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mDisplayPedSigma_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),200,-99.5,100.5);

      HistName = Form("h_mDisplayCMNSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mDisplayCMNSigma_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),200,-99.5,100.5);

      HistName = Form("h_mDisplayRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mDisplayRanSigma_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),200,-99.5,100.5);
    }
  }

  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string gName = Form("g_mRoPedMean_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRoPedMean[i_layer][i_tb] = new TGraph();
      g_mRoPedMean[i_layer][i_tb]->SetName(gName.c_str());

      gName = Form("g_mRoPedSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRoPedSigma[i_layer][i_tb] = new TGraph();
      g_mRoPedSigma[i_layer][i_tb]->SetName(gName.c_str());

      gName = Form("g_mRoCMNSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRoCMNSigma[i_layer][i_tb] = new TGraph();
      g_mRoCMNSigma[i_layer][i_tb]->SetName(gName.c_str());

      gName = Form("g_mRoRanSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRoRanSigma[i_layer][i_tb] = new TGraph();
      g_mRoRanSigma[i_layer][i_tb]->SetName(gName.c_str());
    }
  }

  return clearPedestal();
}

bool FstNoiseStudy::calPedestal()
{
  cout << "FstNoiseStudy::calPedestal -> " << endl;
  cout << " Only use first 1000 event for Pedstal Calculation!" << endl;

  int NumOfEvents = (int)mChainInPut->GetEntries();
  if(NumOfEvents > 1000) NumOfEvents = 1000;
  // const int NumOfEvents = 1000;
  mChainInPut->GetEntry(0);

  //  Calculate a rolling average and standard deviation
  int counters[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numChannels][FST::numTBins];
  double sumValues[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numChannels][FST::numTBins];
  double sumValuesSquared[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numChannels][FST::numTBins];

  //--------------------------------------------------------
  // 1st loop to get a rough estimate of ped (including Hits)
  cout << "=====>First Pedestal Pass...." << endl;
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	{
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    counters[i_arm][i_port][i_apv][i_ch][i_tb] = 0;
	    sumValues[i_arm][i_port][i_apv][i_ch][i_tb] = 0.0;
	    sumValuesSquared[i_arm][i_port][i_apv][i_ch][i_tb] = 0.0;
	  }
	}
      }
    }
  }

  for(int i_event = 0; i_event < NumOfEvents; ++i_event)
  { 
    // if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << NumOfEvents << endl;
    mChainInPut->GetEntry(i_event);

    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  int rdo  = evt_rdo[i_arm][i_port][i_apv];
	  int arm  = evt_arm[i_arm][i_port][i_apv];
	  int port = evt_port[i_arm][i_port][i_apv];
	  int apv  = evt_apv[i_arm][i_port][i_apv];

	  // NOTE THAT RDO/ARM/PORT ARE HARDWIRED HERE!!!
	  bool pass = ( ( rdo == 1 ) && ( (arm == 0) || (arm == 1 ) ) && ( (port == 0) || (port == 1) ) &&  ( apv >= 0 ) && ( apv < FST::numAPVs ) ) ;
	  bool bAPV = isBadAPV(arm,port,apv);
	  if(pass && !bAPV)
	  {
	    for(int i_ro = 0; i_ro < FST::numROChannels; ++i_ro)
	    {
	      int tb = hit_tb[i_arm][i_port][i_apv][i_ro]; // time bin
	      int ch = hit_ch[i_arm][i_port][i_apv][i_ro]; // real channel number 
	      int adc = hit_adc[i_arm][i_port][i_apv][i_ro];
	      sumValues[arm][port][apv][ch][tb] += adc;
	      sumValuesSquared[arm][port][apv][ch][tb] += adc * adc;
	      counters[arm][port][apv][ch][tb]++;
	    }
	  }
	}
      }
    }
  }

  // rough estimate of ped (including Hits)
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	{
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    if(counters[i_arm][i_port][i_apv][i_ch][i_tb] > 0) // eject bad channels
	    {
	      mPed[i_arm][i_port][i_apv][i_ch][i_tb] = sumValues[i_arm][i_port][i_apv][i_ch][i_tb]/counters[i_arm][i_port][i_apv][i_ch][i_tb];
	      mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb] = sqrt((sumValuesSquared[i_arm][i_port][i_apv][i_ch][i_tb]-(double)counters[i_arm][i_port][i_apv][i_ch][i_tb]*mPed[i_arm][i_port][i_apv][i_ch][i_tb]*mPed[i_arm][i_port][i_apv][i_ch][i_tb])/(double)(counters[i_arm][i_port][i_apv][i_ch][i_tb]-1));
	      mPedRMS[i_arm][i_port][i_apv][i_ch][i_tb] = sqrt(mPed[i_arm][i_port][i_apv][i_ch][i_tb]*mPed[i_arm][i_port][i_apv][i_ch][i_tb]+mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]*mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]); // rms^2 = mean^2+std^2
	    }
	  }
	}
      }
    }
  }
  //--------------------------------------------------------

  //--------------------------------------------------------
  // 2nd loop to find ped excluding Hits
  cout << "=====>Second Pedestal Pass...." << endl;
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	{
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    counters[i_arm][i_port][i_apv][i_ch][i_tb] = 0;
	    sumValues[i_arm][i_port][i_apv][i_ch][i_tb] = 0.0;
	    sumValuesSquared[i_arm][i_port][i_apv][i_ch][i_tb] = 0.0;
	  }
	}
      }
    }
  }
  for(int i_event = 0; i_event < NumOfEvents; ++i_event)
  { 
    // if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << NumOfEvents << endl;
    mChainInPut->GetEntry(i_event);

    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  int rdo  = evt_rdo[i_arm][i_port][i_apv];
	  int arm  = evt_arm[i_arm][i_port][i_apv];
	  int port = evt_port[i_arm][i_port][i_apv];
	  int apv  = evt_apv[i_arm][i_port][i_apv];

	  // NOTE THAT RDO/ARM/PORT ARE HARDWIRED HERE!!!
	  bool pass = ( ( rdo == 1 ) && ( (arm == 0) || (arm == 1 ) ) && ( (port == 0) || (port == 1) ) &&  ( apv >= 0 ) && ( apv < FST::numAPVs ) ) ;
	  bool bAPV = isBadAPV(arm,port,apv);
	  if(pass && !bAPV)
	  {
	    for(int i_ro = 0; i_ro < FST::numROChannels; ++i_ro)
	    {
	      int tb = hit_tb[i_arm][i_port][i_apv][i_ro]; // time bin
	      int ch = hit_ch[i_arm][i_port][i_apv][i_ro]; // real channel number 
	      int adc = hit_adc[i_arm][i_port][i_apv][i_ro];
	      float nPedsCut = FST::nIstPedsCut; // 3.0 for IST
	      if(i_arm == 1 && i_port == 1) nPedsCut = FST::nFstPedsCut; // 3.0 for FST
	      if ( (adc < mPed[arm][port][apv][ch][tb]+nPedsCut*mPedStdDev[arm][port][apv][ch][tb]) && ( adc >= 0 && adc < 4096) )
	      { // only adc belew ped+3sigma are considered for 2nd loop
		sumValues[arm][port][apv][ch][tb] += adc;
		sumValuesSquared[arm][port][apv][ch][tb] += adc * adc;
		counters[arm][port][apv][ch][tb]++;
	      }
	    }
	  }
	}
      }
    }
  }
  // calculate ped excluding Hits
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	{
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    if(counters[i_arm][i_port][i_apv][i_ch][i_tb] > 0) // eject bad channels
	    {
	      mPed[i_arm][i_port][i_apv][i_ch][i_tb] = sumValues[i_arm][i_port][i_apv][i_ch][i_tb]/counters[i_arm][i_port][i_apv][i_ch][i_tb];
	      mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb] = sqrt((sumValuesSquared[i_arm][i_port][i_apv][i_ch][i_tb]-(double)counters[i_arm][i_port][i_apv][i_ch][i_tb]*mPed[i_arm][i_port][i_apv][i_ch][i_tb]*mPed[i_arm][i_port][i_apv][i_ch][i_tb])/(double)(counters[i_arm][i_port][i_apv][i_ch][i_tb]-1));
	      mPedRMS[i_arm][i_port][i_apv][i_ch][i_tb] = sqrt(mPed[i_arm][i_port][i_apv][i_ch][i_tb]*mPed[i_arm][i_port][i_apv][i_ch][i_tb]+mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]*mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]); // rms^2 = mean^2+std^2
	    }

	    int layer = this->getLayer(i_arm,i_port);
	    int col = this->getColumn(i_arm,i_port,i_apv,i_ch);
	    int row = this->getRow(i_arm,i_port,i_apv,i_ch);

	    g_mPedMean[layer][i_tb]->SetPoint(i_apv*FST::numChannels+i_ch,i_apv*FST::numChannels+i_ch,mPed[i_arm][i_port][i_apv][i_ch][i_tb]);
	    g_mPedSigma[layer][i_tb]->SetPoint(i_apv*FST::numChannels+i_ch,i_apv*FST::numChannels+i_ch,mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	    h_mPedDisplay[layer][i_tb]->Fill(col,row,mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	    if(layer == 0 && col > -1) 
	    { // FST
	      h_mPedMean_FST[col][i_tb]->SetBinContent(row+1,mPed[i_arm][i_port][i_apv][i_ch][i_tb]);
	      h_mPedSigma_FST[col][i_tb]->SetBinContent(row+1,mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	    }

	    // fill noise display w.r.t. readout order
	    int roChannel = mRoChannelMap[i_ch];
	    g_mRoPedMean[layer][i_tb]->SetPoint(i_apv*FST::numChannels+roChannel,i_apv*FST::numROChannels+roChannel*FST::numTBins+i_tb,mPed[i_arm][i_port][i_apv][i_ch][i_tb]);
	    g_mRoPedSigma[layer][i_tb]->SetPoint(i_apv*FST::numChannels+roChannel,i_apv*FST::numROChannels+roChannel*FST::numTBins+i_tb,mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	    // cout << "RO Channel: " << i_apv*FST::numROChannels+roChannel*FST::numTBins+i_tb << " for PH Channel: " << i_ch << " and TB: " << i_tb << endl;
	  }
	}
      }
    }
  }
  //--------------------------------------------------------

  //--------------------------------------------------------
  // 3rd loop to find Common Mode Noise

  //  Calculate a rolling average and standard deviation
  int counters_CMN[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numRStrip][FST::numTBins];
  double sumValues_CMN[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numRStrip][FST::numTBins];
  double sumValuesSquared_CMN[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numRStrip][FST::numTBins];

  cout << "=====>Third Pass for Common Mode Noise...." << endl;
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
	{
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    counters_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb] = 0;
	    sumValues_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb] = 0.0;
	    sumValuesSquared_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb] = 0.0;
	  }
	}
      }
    }
  }

  for(int i_event = 0; i_event < NumOfEvents; ++i_event)
  { 
    // if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << NumOfEvents << endl;
    mChainInPut->GetEntry(i_event);

    // calculate a rolling average for each event => CMN distribution
    int counters_evt[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numRStrip][FST::numTBins];
    double sumValues_evt[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numRStrip][FST::numTBins];
    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
	  {
	    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	    {
	      counters_evt[i_arm][i_port][i_apv][i_rstrip][i_tb] = 0;
	      sumValues_evt[i_arm][i_port][i_apv][i_rstrip][i_tb] = 0.0;
	    }
	  }
	}
      }
    }
    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  int rdo  = evt_rdo[i_arm][i_port][i_apv];
	  int arm  = evt_arm[i_arm][i_port][i_apv];
	  int port = evt_port[i_arm][i_port][i_apv];
	  int apv  = evt_apv[i_arm][i_port][i_apv];

	  // NOTE THAT RDO/ARM/PORT ARE HARDWIRED HERE!!!
	  bool pass = ( ( rdo == 1 ) && ( (arm == 0) || (arm == 1 ) ) && ( (port == 0) || (port == 1) ) &&  ( apv >= 0 ) && ( apv < FST::numAPVs ) ) ;
	  bool bAPV = isBadAPV(arm,port,apv);
	  if(pass && !bAPV)
	  {
	    for(int i_ro = 0; i_ro < FST::numROChannels; ++i_ro)
	    {
	      int tb = hit_tb[i_arm][i_port][i_apv][i_ro]; // time bin
	      int ch = hit_ch[i_arm][i_port][i_apv][i_ro]; // real channel number 
	      int adc = hit_adc[i_arm][i_port][i_apv][i_ro];
	      float nPedsCut = FST::nIstPedsCut; // 3.0 for IST
	      if(i_arm == 1 && i_port == 1) nPedsCut = FST::nFstPedsCut; // 3.0 for FST

	      if ( (adc < mPed[arm][port][apv][ch][tb]+nPedsCut*mPedStdDev[arm][port][apv][ch][tb]) && ( adc >= 0 && adc < 4096) )
	      { // only adc below ped+3sigma on FST are considered for 3rd loop
		int layer = this->getLayer(arm,port);
		int col = this->getColumn(arm,port,apv,ch); // 0-7 for FST & 0-23 for IST => convert to 0 (even) & 1 (odd)
		if(layer > 0) col = col%2;
		sumValues_evt[arm][port][apv][col][tb] += adc-mPed[arm][port][apv][ch][tb];
		counters_evt[arm][port][apv][col][tb]++;
		if(layer == 0 && col > -1) h_mDisplayPedSigma_FST[col][tb]->Fill(adc-mPed[arm][port][apv][ch][tb]);
	      }
	    }
	  }
	}
      }
    }

    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
	    {
	      if(counters_evt[i_arm][i_port][i_apv][i_rstrip][i_tb] > 0) // eject bad channels
	      { // calculate CMN for each rstrip
		double CMN_evt = sumValues_evt[i_arm][i_port][i_apv][i_rstrip][i_tb]/counters_evt[i_arm][i_port][i_apv][i_rstrip][i_tb];
		sumValues_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb] += CMN_evt;
		sumValuesSquared_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb] += CMN_evt*CMN_evt;
		counters_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb]++;
		if(i_arm == 1 && i_port == 1) h_mDisplayCMNSigma_FST[i_rstrip][i_tb]->Fill(CMN_evt);
	      }
	    }
	  }
	}
      }
    }
  }
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	{
	  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
	  {
	    if(counters_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb] > 0) // eject bad channels
	    {
	      double meanCMN = sumValues_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb]/counters_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb];
	      mCMNMean[i_arm][i_port][i_apv][i_rstrip][i_tb] = sqrt((sumValuesSquared_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb]-(double)counters_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb]*meanCMN*meanCMN)/(double)(counters_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb]-1));
	    }
	  }
	  for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	  {
	    int layer = this->getLayer(i_arm,i_port);
	    int col = this->getColumn(i_arm,i_port,i_apv,i_ch);
	    if(layer > 0) col = col%2;
	    int row = this->getRow(i_arm,i_port,i_apv,i_ch);
	    if(counters_CMN[i_arm][i_port][i_apv][col][i_tb] > 0 && col > -1)
	    { // set CMN for each channel
	      mCMNStdDev[i_arm][i_port][i_apv][i_ch][i_tb] = mCMNMean[i_arm][i_port][i_apv][col][i_tb];
	    }
	    g_mCMNSigma[layer][i_tb]->SetPoint(i_apv*FST::numChannels+i_ch,i_apv*FST::numChannels+i_ch,mCMNStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	    if(layer == 0 && col > -1) h_mCMNSigma_FST[col][i_tb]->SetBinContent(row+1,mCMNStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);

	    // fill noise display w.r.t. readout order
	    int roChannel = mRoChannelMap[i_ch];
	    g_mRoCMNSigma[layer][i_tb]->SetPoint(i_apv*FST::numChannels+roChannel,i_apv*FST::numROChannels+roChannel*FST::numTBins+i_tb,mCMNStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	  }
	}
      }
    }
  }
  //--------------------------------------------------------

  //--------------------------------------------------------
  // 4th loop to find Differential Noise
  cout << "=====>Fourth Pass for Differential Noise...." << endl;
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	{
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    counters[i_arm][i_port][i_apv][i_ch][i_tb] = 0;
	    sumValues[i_arm][i_port][i_apv][i_ch][i_tb] = 0.0;
	    sumValuesSquared[i_arm][i_port][i_apv][i_ch][i_tb] = 0.0;
	  }
	}
      }
    }
  }

  for(int i_event = 0; i_event < NumOfEvents; ++i_event)
  { 
    // if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << NumOfEvents << endl;
    mChainInPut->GetEntry(i_event);

    clearCMN(); // clear event-by-event CMN
    // calculate a rolling average for each event => CMN distribution
    int counters_evt[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numRStrip][FST::numTBins];
    double sumValues_evt[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numRStrip][FST::numTBins];
    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
	  {
	    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	    {
	      counters_evt[i_arm][i_port][i_apv][i_rstrip][i_tb] = 0;
	      sumValues_evt[i_arm][i_port][i_apv][i_rstrip][i_tb] = 0.0;
	    }
	  }
	}
      }
    }
    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  int rdo  = evt_rdo[i_arm][i_port][i_apv];
	  int arm  = evt_arm[i_arm][i_port][i_apv];
	  int port = evt_port[i_arm][i_port][i_apv];
	  int apv  = evt_apv[i_arm][i_port][i_apv];

	  // NOTE THAT RDO/ARM/PORT ARE HARDWIRED HERE!!!
	  bool pass = ( ( rdo == 1 ) && ( (arm == 0) || (arm == 1 ) ) && ( (port == 0) || (port == 1) ) &&  ( apv >= 0 ) && ( apv < FST::numAPVs ) ) ;
	  bool bAPV = isBadAPV(arm,port,apv);
	  if(pass && !bAPV)
	  {
	    for(int i_ro = 0; i_ro < FST::numROChannels; ++i_ro)
	    {
	      int tb = hit_tb[i_arm][i_port][i_apv][i_ro]; // time bin
	      int ch = hit_ch[i_arm][i_port][i_apv][i_ro]; // real channel number 
	      int adc = hit_adc[i_arm][i_port][i_apv][i_ro];
	      float nPedsCut = FST::nIstPedsCut; // 3.0 for IST
	      if(i_arm == 1 && i_port == 1) nPedsCut = FST::nFstPedsCut; // 3.0 for FST

	      if ( (adc < mPed[arm][port][apv][ch][tb]+nPedsCut*mPedStdDev[arm][port][apv][ch][tb]) && ( adc >= 0 && adc < 4096) )
	      { // only adc below ped+3sigma on FST are considered for 4th loop
		int layer = this->getLayer(arm,port);
		int col = this->getColumn(arm,port,apv,ch); // 0-7 for FST & 0-23 for IST => convert to 0 (even) & 1 (odd)
		if(layer > 0) col = col%2;
		sumValues_evt[arm][port][apv][col][tb] += adc-mPed[arm][port][apv][ch][tb];
		counters_evt[arm][port][apv][col][tb]++;
	      }
	    }
	  }
	}
      }
    }

    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
	    {
	      if(counters_evt[i_arm][i_port][i_apv][i_rstrip][i_tb] > 0) // eject bad channels
	      { // calculate CMN for each group
		mCMNMean_Evt[i_arm][i_port][i_apv][i_rstrip][i_tb] = sumValues_evt[i_arm][i_port][i_apv][i_rstrip][i_tb]/counters_evt[i_arm][i_port][i_apv][i_rstrip][i_tb];
	      }
	    }
	    for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	    {
	      int layer = this->getLayer(i_arm,i_port);
	      int col = this->getColumn(i_arm,i_port,i_apv,i_ch);
	      if(layer > 0) col = col%2;
	      int row = this->getRow(i_arm,i_port,i_apv,i_ch);
	      if(counters_evt[i_arm][i_port][i_apv][col][i_tb] > 0 && col > -1)
	      { // set CMN for each channel
		mCMNStdDev_Evt[i_arm][i_port][i_apv][i_ch][i_tb] = mCMNMean_Evt[i_arm][i_port][i_apv][col][i_tb];
	      }
	    }
	  }
	}
      }
    }

    // calculate differential noise
    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  int rdo  = evt_rdo[i_arm][i_port][i_apv];
	  int arm  = evt_arm[i_arm][i_port][i_apv];
	  int port = evt_port[i_arm][i_port][i_apv];
	  int apv  = evt_apv[i_arm][i_port][i_apv];

	  // NOTE THAT RDO/ARM/PORT ARE HARDWIRED HERE!!!
	  bool pass = ( ( rdo == 1 ) && ( (arm == 0) || (arm == 1 ) ) && ( (port == 0) || (port == 1) ) &&  ( apv >= 0 ) && ( apv < FST::numAPVs ) ) ;
	  bool bAPV = isBadAPV(arm,port,apv);
	  if(pass && !bAPV)
	  {
	    for(int i_ro = 0; i_ro < FST::numROChannels; ++i_ro)
	    {
	      int tb = hit_tb[i_arm][i_port][i_apv][i_ro]; // time bin
	      int ch = hit_ch[i_arm][i_port][i_apv][i_ro]; // real channel number 
	      int adc = hit_adc[i_arm][i_port][i_apv][i_ro];
	      float nPedsCut = FST::nIstPedsCut; // 3.0 for IST
	      if(i_arm == 1 && i_port == 1) nPedsCut = FST::nFstPedsCut; // 3.0 for FST
	      if ( (adc < mPed[arm][port][apv][ch][tb]+nPedsCut*mPedStdDev[arm][port][apv][ch][tb]) && ( adc >= 0 && adc < 4096) )
	      { // only adc belew ped+3sigma are considered
		sumValues[arm][port][apv][ch][tb] += adc-mPed[arm][port][apv][ch][tb]-mCMNStdDev_Evt[arm][port][apv][ch][tb];
		sumValuesSquared[arm][port][apv][ch][tb] += (adc-mPed[arm][port][apv][ch][tb]-mCMNStdDev_Evt[arm][port][apv][ch][tb])*(adc-mPed[arm][port][apv][ch][tb]-mCMNStdDev_Evt[arm][port][apv][ch][tb]);
		counters[arm][port][apv][ch][tb]++;
		int col = this->getColumn(arm,port,apv,ch);
		if(arm == 1 && port == 1 && col > -1) h_mDisplayRanSigma_FST[col][tb]->Fill(adc-mPed[arm][port][apv][ch][tb]-mCMNStdDev_Evt[arm][port][apv][ch][tb]);
	      }
	    }
	  }
	}
      }
    }
  }
  for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
  {
    for(int i_port = 0; i_port < FST::numPorts; ++i_port)
    {
      for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
      {
	for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	{
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    if(counters[i_arm][i_port][i_apv][i_ch][i_tb] > 0) // eject bad channels
	    {
	      double meanDiff = sumValues[i_arm][i_port][i_apv][i_ch][i_tb]/counters[i_arm][i_port][i_apv][i_ch][i_tb];
	      mRanStdDev[i_arm][i_port][i_apv][i_ch][i_tb] = sqrt((sumValuesSquared[i_arm][i_port][i_apv][i_ch][i_tb]-(double)counters[i_arm][i_port][i_apv][i_ch][i_tb]*meanDiff*meanDiff)/(double)(counters[i_arm][i_port][i_apv][i_ch][i_tb]-1));
	    }

	    int layer = this->getLayer(i_arm,i_port);
	    int col = this->getColumn(i_arm,i_port,i_apv,i_ch);
	    int row = this->getRow(i_arm,i_port,i_apv,i_ch);
	    g_mRanSigma[layer][i_tb]->SetPoint(i_apv*FST::numChannels+i_ch,i_apv*FST::numChannels+i_ch,mRanStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	    if(layer == 0 && col > -1) h_mRanSigma_FST[col][i_tb]->SetBinContent(row+1,mRanStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);

	    // fill noise display w.r.t. readout order
	    int roChannel = mRoChannelMap[i_ch];
	    g_mRoRanSigma[layer][i_tb]->SetPoint(i_apv*FST::numChannels+roChannel,i_apv*FST::numROChannels+roChannel*FST::numTBins+i_tb,mRanStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	  }
	}
      }
    }
  }
  //--------------------------------------------------------

  return true;
}

void FstNoiseStudy::writePedestal()
{
  cout << "FstNoiseStudy::writePedestal => save ped mean & sigma!" << endl;

  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mPedMean[i_layer][i_tb]->Write();
      g_mPedSigma[i_layer][i_tb]->Write();
      g_mCMNSigma[i_layer][i_tb]->Write();
      g_mRanSigma[i_layer][i_tb]->Write();
      h_mPedDisplay[i_layer][i_tb]->Write();
    }
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      h_mPedMean_FST[i_rstrip][i_tb]->Write();
      h_mPedSigma_FST[i_rstrip][i_tb]->Write();
      h_mCMNSigma_FST[i_rstrip][i_tb]->Write();
      h_mRanSigma_FST[i_rstrip][i_tb]->Write();

      h_mDisplayPedSigma_FST[i_rstrip][i_tb]->Write();
      h_mDisplayCMNSigma_FST[i_rstrip][i_tb]->Write();
      h_mDisplayRanSigma_FST[i_rstrip][i_tb]->Write();
    }
  }

  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mRoPedMean[i_layer][i_tb]->Write();
      g_mRoPedSigma[i_layer][i_tb]->Write();
      g_mRoCMNSigma[i_layer][i_tb]->Write();
      g_mRoRanSigma[i_layer][i_tb]->Write();
    }
  }
}
//--------------pedestal---------------------

//--------------channel map---------------------
bool FstNoiseStudy::getChannelMap()
{
  cout << "FstNoiseStudy::getChannelMap->" << endl;
  string inputfile = "../../src/FstNoiseStudy/channelMap.txt";
  cout << "Open channel map inputfile: " << inputfile.c_str() << endl;
  ifstream file_channelMap (inputfile.c_str());
  if( !file_channelMap.is_open() )
  {
    cout << "Abort. Fail to read in channel map: " << inputfile.c_str() << endl;
    return false;
  }

  int readout_channel = -1;
  int physics_channel = -1;

  cout << "reading in channel map:" << endl;
  while (file_channelMap >> physics_channel >> readout_channel)
  {
    cout << "physics_channel = " << physics_channel << ", readout_channel = " << readout_channel << endl;
    mRoChannelMap[physics_channel] = readout_channel;
  }
  file_channelMap.close();

  return true;
}

//--------------channel map---------------------

//--------------Utility---------------------
int FstNoiseStudy::getLayer(int arm, int port)
{
  if(arm == 0 && port == 0) return 1; // top IST
  if(arm == 0 && port == 1) return 2; // middle IST
  if(arm == 1 && port == 0) return 3; // bottom IST
  if(arm == 1 && port == 1) return 0; // FST module

  return -1;
}

int FstNoiseStudy::getSensor(int arm, int port, int apv)
{
  int layer = this->getLayer(arm,port);
  if(layer == 0)
  { // layer = 0 for FST
    if(apv >= 0 && apv <= 3) return 0; // inner sector
    if(apv >= 4 && apv <= 5) return 1; // 1st half of outer sector
    if(apv >= 6 && apv <= 7) return 2; // 2nd half of outer sector
  }
  else
  { // layer = 1-3 for IST
    if(apv >= 0 && apv <= 5) return 0; // sensor 0 for first half group of APVs
    if(apv >= 6 && apv <= 11) return 1; // sensor 1 for second half group of APVs 
  }

  return -1;
}

int FstNoiseStudy::getColumn(int arm, int port, int apv, int ch)
{
  int col = -1;
  int layer = this->getLayer(arm,port);
  int sensor = this->getSensor(arm,port,apv);

  if(layer == 0)
  { // layer = 0 for FST
    if(this->getRStrip(apv,ch) >= 0)
    {
      if(sensor == 0) col = this->getRStrip(apv,ch);
      if(sensor > 0) col = 4 + this->getRStrip(apv,ch);
    }
  }
  else
  { // layer = 1-3 for IST
    if(ch >= 0 && ch <= 63) col = apv*2;
    else col = apv*2 + 1;
  }

  return col;
}

int FstNoiseStudy::getRow(int arm, int port, int apv, int ch)
{
  int row = -1;
  int layer = this->getLayer(arm,port);
  if(layer == 0)
  { // layer = 0 for FST
    row = this->getPhiSeg(apv,ch);
  }
  else
  { // layer = 1-3 for IST
    if(ch >= 0 && ch <= 63) row = ch;
    else row = 127 - ch;
  }

  return row;
}

int FstNoiseStudy::getRStrip(int apv, int ch)
{
  int r_strip = -1;

  // inner sector
  if(apv == 0 || apv == 1) r_strip = ch%4;
  if(apv == 2 || apv == 3) r_strip = 3-ch%4;

  // outer sector
  if(apv == 4 || apv == 6) r_strip = ch%4;
  if(apv == 5 || apv == 7) r_strip = 3-ch%4;

  return r_strip;
}

int FstNoiseStudy::getPhiSeg(int apv, int ch)
{
  int phi_seg = -1;

  if(apv >= 0 && apv <= 3) phi_seg = apv*32 + ch/4; // inner
  if(apv >= 4 && apv <= 7) phi_seg = (apv-4)*32 + ch/4; // outer

  return phi_seg;
}

double FstNoiseStudy::getPosX(int arm, int port, int apv, int ch)
{
  double posX = -999.9;
  int layer = this->getLayer(arm,port);
  if(layer == 0) // FST
  {
    // double r_fst = FST::rOuter + (this->getColumn(arm,port,apv,ch)-4)*FST::pitchR + 0.5*FST::pitchR;
    // double phi_fst = (63-this->getRow(arm,port,apv,ch))*FST::pitchPhi + 0.5*FST::pitchPhi;
    // posX = r_fst*TMath::Cos(phi_fst); // x = r*cos(phi)
    posX = FST::rOuter + (this->getColumn(arm,port,apv,ch)-4)*FST::pitchR + 0.5*FST::pitchR;
  }
  else // IST1-3
  {
    posX = this->getColumn(arm,port,apv,ch)*FST::pitchColumn + 0.5*FST::pitchColumn;
  }

  return posX;
}

double FstNoiseStudy::getPosY(int arm, int port, int apv, int ch)
{
  double posY = -999.9;
  int layer = this->getLayer(arm,port);
  if(layer == 0) // FST
  {
    // double r_fst = FST::rOuter + (this->getColumn(arm,port,apv,ch)-4)*FST::pitchR + 0.5*FST::pitchR;
    // double phi_fst = (63-this->getRow(arm,port,apv,ch))*FST::pitchPhi + 0.5*FST::pitchPhi;
    // posY = r_fst*TMath::Sin(phi_fst); // y = r*sin(phi)
    posY = (63-this->getRow(arm,port,apv,ch))*FST::pitchPhi + 0.5*FST::pitchPhi;
  }
  else // IST1-3
  {
    posY = (63-this->getRow(arm,port,apv,ch))*FST::pitchRow + 0.5*FST::pitchRow;
  }

  return posY;
}

double FstNoiseStudy::getMeanX(int layer, double meanColumn)
{
  double meanX = -999.9;
  if(layer == 0) // FST
  {
    // double r_fst = FST::rOuter + (this->getColumn(arm,port,apv,ch)-4)*FST::pitchR + 0.5*FST::pitchR;
    // double phi_fst = (63-this->getRow(arm,port,apv,ch))*FST::pitchPhi + 0.5*FST::pitchPhi;
    // posX = r_fst*TMath::Cos(phi_fst); // x = r*cos(phi)
    meanX = FST::rOuter + (meanColumn-4)*FST::pitchR + 0.5*FST::pitchR;
  }
  else // IST1-3
  {
    meanX = meanColumn*FST::pitchColumn + 0.5*FST::pitchColumn;
  }

  return meanX;
}

double FstNoiseStudy::getMeanY(int layer, double meanRow)
{
  double meanY = -999.9;
  if(layer == 0) // FST
  {
    meanY = (63-meanRow)*FST::pitchPhi + 0.5*FST::pitchPhi;
  }
  else // IST1-3
  {
    meanY = (63-meanRow)*FST::pitchRow + 0.5*FST::pitchRow;
  }

  return meanY;
}

bool FstNoiseStudy::isBadAPV(int arm, int port, int apv)
{
  bool bAPV = false;

  if(arm == 0 && port==0 && (apv == 1 || apv == 5||apv == 6)) bAPV = true;
  if(arm == 0 && port==1 && (apv <= 9)) bAPV = true;
  // if(arm == 1 && port==0 && (apv == 0||apv == 1)) bAPV = true;
  if(arm == 1 && port==0 && (apv == 0)) bAPV = true;
  if(arm == 1 && port==1 && (apv > 7)) bAPV = true;

  return bAPV;
}
//--------------Utility---------------------
