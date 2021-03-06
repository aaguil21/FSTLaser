#include "./FstClusterMaker.h"
#include "../FstUtil/FstRawHit.h"
#include "../FstUtil/FstCluster.h"
#include "../FstUtil/FstTrack.h"
#include "../FstUtil/FstEvent.h"

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

ClassImp(FstClusterMaker)

//------------------------------------------

FstClusterMaker::FstClusterMaker() : mList("../../list/FST/FstData_HV140.list"), mOutPutFile("./FstData_HV140.root"), mSavePed(true), mApplyCMNCorr(true), mFstHitsCut(4.5), mNumOfUsedTimeBins(3), mFstThresholdCut2(2.5), mFstThresholdCut1(3.5)
{
  cout << "FstClusterMaker::FstClusterMaker() -------- Constructor!  --------" << endl;
  mHome = getenv("HOME");
}

FstClusterMaker::~FstClusterMaker()
{
  cout << "FstClusterMaker::~FstClusterMaker() -------- Release Memory!  --------" << endl;
}

//------------------------------------------
int FstClusterMaker::Init()
{
  cout << "FstClusterMaker::Init => " << endl;
  cout << "Configurations: mSavePed: " <<  mSavePed << ", mApplyCMNCorr: " << mApplyCMNCorr << endl;
  cout << "mFstHitsCut: " << mFstHitsCut << ", mNumOfUsedTimeBins: " << mNumOfUsedTimeBins << ", mFstThresholdCut2: " << mFstThresholdCut2 << ", mFstThresholdCut1: " << mFstThresholdCut1 << endl;
  File_mOutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

  bool isInPut = initChain(); // initialize input data/ped TChain;
  bool isPed = initPedestal(); // initialize pedestal array;
  bool isSig = initSignal(); // initialize signal array;
  bool isTree = initTree();

  if(!isInPut) 
  {
    cout << "Failed to initialize input data!" << endl;
    return -1;
  }
  if(!isPed) 
  {
    cout << "Failed to initialize pedestals!" << endl;
    return -1;
  }
  if(!isSig) 
  {
    cout << "Failed to initialize signal!" << endl;
    return -1;
  }
  if(!isTree) 
  {
    cout << "Failed to initialize TTree for FST & IST Hits & Clusters!" << endl;
    return -1;
  }

  initHitDisplay(); // initialize hit display

  return 1;
}

bool FstClusterMaker::initChain()
{
  cout << "FstClusterMaker::initChain -> " << endl;

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

bool FstClusterMaker::initSignal()
{
  cout << "FstClusterMaker::initSignal -> " << endl;
  return clearSignal();
}

//------------------------------------------
bool FstClusterMaker::clearSignal()
{
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
	    mRawSig[i_arm][i_port][i_apv][i_ch][i_tb] = -1.0;
	    mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb] = -1.0;
	    mSigCMNCorr[i_arm][i_port][i_apv][i_ch][i_tb] = -1.0;
	  }
	}
      }
    }
  }

  return true;
}

//------------------------------------------
int FstClusterMaker::Make()
{
  cout << "FstClusterMaker::Make => " << endl;

  bool isPed = calPedestal(); // calculate pedestals with first 1000 events
  if( !isPed ) 
  {
    cout << "No Pedestal!!!" << endl;
    return -1;
  }

  long NumOfEvents = (long)mChainInPut->GetEntries();
  // if(NumOfEvents > 1000) NumOfEvents = 1000;
  // NumOfEvents = 1000;
  mChainInPut->GetEntry(0);

  for(int i_event = 0; i_event < NumOfEvents; ++i_event)
  {
    if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << NumOfEvents << endl;
    mChainInPut->GetEntry(i_event);
    clearSignal();

    //--------------------------------
    // Calculate event-by-event CMN
    if(mApplyCMNCorr)
    {
      clearCMN(); // clear event-by-event CMN
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
		  if(layer < 3) col = col%2;
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
		if(layer < 3) col = col%2;
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
    }
    //--------------------------------

    // calculate ped & cmn corrected signal for each ch & time bin
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

	      mRawSig[i_arm][i_port][i_apv][ch][tb] = adc;
	      mSigPedCorr[i_arm][i_port][i_apv][ch][tb] = adc-mPed[i_arm][i_port][i_apv][ch][tb]; // adc - ped
	      mSigCMNCorr[i_arm][i_port][i_apv][ch][tb] = adc-mPed[i_arm][i_port][i_apv][ch][tb]-mCMNStdDev_Evt[i_arm][i_port][i_apv][ch][tb]; // adc - ped - cmn
	    }
	  }
	}
      }
    }

    // set signal and noise for each event based on mApplyCMNCorr flag
    // mApplyCMNCorr = true => Ped & CMN corrected signal and random noise | mApplyCMNCorr = false => Ped corrected signal and total noise
    double signalEvt[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numChannels][FST::numTBins]; // signal for each event
    double noiseEvt[FST::numARMs][FST::numPorts][FST::numAPVs][FST::numChannels][FST::numTBins]; // noise for each event
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
	      if(!mApplyCMNCorr)
	      {
		signalEvt[i_arm][i_port][i_apv][i_ch][i_tb] = mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb];
		noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb] = mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb];
	      }
	      if(mApplyCMNCorr)
	      {
		signalEvt[i_arm][i_port][i_apv][i_ch][i_tb] = mSigCMNCorr[i_arm][i_port][i_apv][i_ch][i_tb];
		noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb] = mRanStdDev[i_arm][i_port][i_apv][i_ch][i_tb];
	      }
	    }
	  }
	}
      }
    }

    // find Hits
    std::vector<FstRawHit *> rawHitVec_orig;
    rawHitVec_orig.clear();
    std::vector<FstRawHit *> rawHitVec_Used;
    rawHitVec_Used.clear();
    int numOfHits = 0; // Hits & Peds on FST & IST
    int numOfHits_Used = 0; // Hits on FST & IST
    int numOfHits_IST = 0; // Hits on IST

    for(int i_arm = 0; i_arm < FST::numARMs; ++i_arm)
    {
      for(int i_port = 0; i_port < FST::numPorts; ++i_port)
      {
	for(int i_apv = 0; i_apv < FST::numAPVs; ++i_apv)
	{
	  for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	  {
	    if( // some hit quality cuts => 1st one is questionable | 2nd made sure ch shows reasonable noise
		// ( mSigPedCorr[i_arm][i_port][i_apv][i_ch][0] < mSigPedCorr[i_arm][i_port][i_apv][i_ch][3] ) &&
		( mPedStdDev[i_arm][i_port][i_apv][i_ch][FST::pedTimeBin] > FST::MinNoise)  
	      )
	    {
	      bool isHit = false;
	      int maxTB = 0;

	      if( !(i_arm == 1 && i_port == 1) ) // IST Hits always required 3 time bin
	      {
		double maxADC = mSigPedCorr[i_arm][i_port][i_apv][i_ch][0]; // init with 1st tb
		double preADC = maxADC;
		float nHitsCut = FST::nIstHitsCut; // 6.0 for IST
		for(int i_tb = 1; i_tb < FST::numTBins-1; ++i_tb)
		{ // only if 3 consequetive timebins of a ch exceed the threshold cut is considered as a hit
		  if( 
		      ( mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb-1] > nHitsCut*mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb-1]) &&
		      ( mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb] > nHitsCut*mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]) &&
		      ( mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb+1] > nHitsCut*mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb+1])
		    ) 
		  {
		    isHit = true; // set isHit to true if 3 consequetive time bins exceed the threshold
		    if(mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb] > preADC)
		    { // find time bin with max adc for 0-FST::numTBins-2
		      maxADC = mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb];
		      maxTB = i_tb;
		      preADC = maxADC;
		    }
		    if(i_tb == FST::numTBins-2 && mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb+1] > preADC)
		    { // check if last time bin has the max ADC
		      maxADC = mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb+1];
		      maxTB = i_tb+1;
		    }
		  }
		}
	      }

	      if(i_arm == 1 && i_port == 0) // FST Hits
	      {
		double maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][0]; // init with 1st tb
		double preADC = maxADC;
		// float nHitsCut = FST::nFstHitsCut; // 4.5 for FST
		float nHitsCut = mFstHitsCut; // 4.5 for FST
		if(mNumOfUsedTimeBins == 3) // FST Hits with 3 Time Bins
		{
		  for(int i_tb = 1; i_tb < FST::numTBins-1; ++i_tb)
		  { // only if 3 consequetive timebins of a ch exceed the threshold cut is considered as a hit
		    if( 
			( signalEvt[i_arm][i_port][i_apv][i_ch][i_tb-1] > nHitsCut*noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb-1]) &&
			( signalEvt[i_arm][i_port][i_apv][i_ch][i_tb]   > nHitsCut*noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb]) &&
			( signalEvt[i_arm][i_port][i_apv][i_ch][i_tb+1] > nHitsCut*noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb+1])
		      ) 
		    {
		      isHit = true; // set isHit to true if 3 consequetive time bins exceed the threshold
		      if(signalEvt[i_arm][i_port][i_apv][i_ch][i_tb] > preADC)
		      { // find time bin with max adc for 0-FST::numTBins-2
			maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][i_tb];
			maxTB = i_tb;
			preADC = maxADC;
		      }
		      if(i_tb == FST::numTBins-2 && signalEvt[i_arm][i_port][i_apv][i_ch][i_tb+1] > preADC)
		      { // check if last time bin has the max ADC
			maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][i_tb+1];
			maxTB = i_tb+1;
		      }
		    }
		  }
		}
		if(mNumOfUsedTimeBins == 2) // FST Hits with 2 Time Bins
		{
		  for(int i_tb = 1; i_tb < FST::numTBins; ++i_tb)
		  { // only if 2 consequetive timebins of a ch exceed the threshold cut is considered as a hit
		    if( 
			( signalEvt[i_arm][i_port][i_apv][i_ch][i_tb-1] > nHitsCut*noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb-1]) &&
			( signalEvt[i_arm][i_port][i_apv][i_ch][i_tb]   > nHitsCut*noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb])
		      ) 
		    {
		      isHit = true; // set isHit to true if 2 consequetive time bins exceed the threshold
		      if(signalEvt[i_arm][i_port][i_apv][i_ch][i_tb] > preADC)
		      { // find time bin with max adc for 0-FST::numTBins-1
			maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][i_tb];
			maxTB = i_tb;
			preADC = maxADC;
		      }
		      if(signalEvt[i_arm][i_port][i_apv][i_ch][i_tb-1] > preADC)
		      { // check if i_tb-1 has the max ADC
			maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][i_tb-1];
			maxTB = i_tb-1;
			preADC = maxADC;
		      }
		    }
		  }
		}
		if(mNumOfUsedTimeBins == 1) // FST Hits with 1 Time Bin
		{
		  for(int i_tb = 1; i_tb < FST::numTBins; ++i_tb)
		  { // only if 1 timebin of a ch exceed the threshold cut is considered as a hit
		    if( 
			( signalEvt[i_arm][i_port][i_apv][i_ch][i_tb] > nHitsCut*noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb])
		      ) 
		    {
		      isHit = true; // set isHit to true if 1 consequetive time bins exceed the threshold
		      if(signalEvt[i_arm][i_port][i_apv][i_ch][i_tb] > preADC)
		      { // find time bin with max adc for 0-FST::numTBins-1
			maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][i_tb];
			maxTB = i_tb;
			preADC = maxADC;
		      }
		    }
		  }
		}
	      }

	      bool isPed = false;
	      // if( !isHit && mSavePed)
	      if( !isHit )
	      {
		maxTB = 0;
		if( !(i_arm == 1 && i_port == 1) ) // IST Peds always required 3 time bin
		{
		  double maxADC = mSigPedCorr[i_arm][i_port][i_apv][i_ch][0]; // init with 1st tb
		  double preADC = maxADC;
		  float nPedsCut = FST::nIstThresholdCut; // 3.0 for IST
		  for(int i_tb = 1; i_tb < FST::numTBins-1; ++i_tb)
		  { // only if 3 consequetive timebins of a ch exceed the threshold cut is considered as a hit
		    if( 
			( mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb-1] > nPedsCut*mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb-1]) &&
			( mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb]   > nPedsCut*mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]) &&
			( mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb+1] > nPedsCut*mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb+1])
		      ) 
		    {
		      isPed = true; // set isHit to true if 3 consequetive time bins exceed the threshold
		      if(mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb] > preADC)
		      { // find time bin with max adc for 0-FST::numTBins-2
			maxADC = mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb];
			maxTB = i_tb;
			preADC = maxADC;
		      }
		      if(i_tb == FST::numTBins-2 && mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb+1] > preADC)
		      { // check if last time bin has the max ADC
			maxADC = mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb+1];
			maxTB = i_tb+1;
		      }
		    }
		  }
		}

		if(i_arm == 1 && i_port == 0) // FST Peds 
		{
		  double maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][0]; // init with 1st tb
		  double preADC = maxADC;
		  // float nPedsCut = FST::nFstThresholdCut; // 2.0 for FST
		  float nPedsCut2 = mFstThresholdCut2; // 2.5 for FST & TimeBin 2
		  for(int i_tb = 1; i_tb < FST::numTBins; ++i_tb)
		  { // only if 2 consequetive timebins of a ch exceed the threshold cut is considered as a hit
		    if( 
			( signalEvt[i_arm][i_port][i_apv][i_ch][i_tb-1] > nPedsCut2*noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb-1]) &&
			( signalEvt[i_arm][i_port][i_apv][i_ch][i_tb]   > nPedsCut2*noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb])
		      ) 
		    {
		      isPed = true; // set isPed to true if 2 consequetive time bins exceed the threshold
		      if(signalEvt[i_arm][i_port][i_apv][i_ch][i_tb] > preADC)
		      { // find time bin with max adc for 0-FST::numTBins-1
			maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][i_tb];
			maxTB = i_tb;
			preADC = maxADC;
		      }
		      if(signalEvt[i_arm][i_port][i_apv][i_ch][i_tb-1] > preADC)
		      { // check if i_tb-1 has the max ADC
			maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][i_tb-1];
			maxTB = i_tb-1;
			preADC = maxADC;
		      }
		    }
		  }
		  if( !isPed )
		  {
		    float nPedsCut1 = mFstThresholdCut1; // 3.5 for FST & TimeBin 1
		    for(int i_tb = 1; i_tb < FST::numTBins; ++i_tb)
		    { // only if 1 timebin of a ch exceed the threshold cut is considered as a hit
		      if( 
			  ( signalEvt[i_arm][i_port][i_apv][i_ch][i_tb] > nPedsCut1*noiseEvt[i_arm][i_port][i_apv][i_ch][i_tb])
			) 
		      {
			isPed = true; // set isPed to true if 1 time bins exceed the threshold
			if(signalEvt[i_arm][i_port][i_apv][i_ch][i_tb] > preADC)
			{ // find time bin with max adc for 0-FST::numTBins-1
			  maxADC = signalEvt[i_arm][i_port][i_apv][i_ch][i_tb];
			  maxTB = i_tb;
			  preADC = maxADC;
			}
		      }
		    }
		  }
		}
	      }

	      // if( (isHit || isPed) && numOfHits < FST::maxNHits)
	      // if( (isHit || isPed) && numOfHits_Used < FST::maxNHits)
	      if( (isHit || isPed) && numOfHits_IST < FST::maxNHits)
	      { // set Hit info
		FstRawHit *fstRawHit = new FstRawHit();
		fstRawHit->Clear();
		fstRawHit->setLayer(getLayer(i_arm,i_port));
		fstRawHit->setSensor(getSensor(i_arm,i_port,i_apv));
		fstRawHit->setAPV(i_apv);
		fstRawHit->setChannel(i_ch);
		fstRawHit->setColumn(getColumn(i_arm,i_port,i_apv,i_ch));
		fstRawHit->setRow(getRow(i_arm,i_port,i_apv,i_ch));
		fstRawHit->setPosX(getPosX(i_arm,i_port,i_apv,i_ch));
		fstRawHit->setPosY(getPosY(i_arm,i_port,i_apv,i_ch));
		for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
		{
		  fstRawHit->setPedMean(mPed[i_arm][i_port][i_apv][i_ch][i_tb], i_tb);
		  fstRawHit->setPedStdDev(mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb], i_tb);
		  fstRawHit->setPedRMS(mPedRMS[i_arm][i_port][i_apv][i_ch][i_tb], i_tb);
		  fstRawHit->setCMNStdDev(mCMNStdDev_Evt[i_arm][i_port][i_apv][i_ch][i_tb], i_tb);
		  fstRawHit->setRanStdDev(mRanStdDev[i_arm][i_port][i_apv][i_ch][i_tb], i_tb);
		  fstRawHit->setRawCharge(mRawSig[i_arm][i_port][i_apv][i_ch][i_tb], i_tb);
		  if(getLayer(i_arm,i_port) < 3) 
		  {
		    fstRawHit->setCharge(mSigPedCorr[i_arm][i_port][i_apv][i_ch][i_tb], i_tb); // IST
		    fstRawHit->setWeight(FST::mSigWeight_IST); // weight for ISt is always 1.0
		  }
		  if(getLayer(i_arm,i_port) == 3) 
		  {
		    int rstrip = getColumn(i_arm,i_port,i_apv,i_ch);
		    fstRawHit->setCharge(signalEvt[i_arm][i_port][i_apv][i_ch][i_tb], i_tb); // FST
		    fstRawHit->setWeight(FST::mSigWeight_FST[rstrip]);
		  }
		}
		fstRawHit->setMaxTb(maxTB);
		fstRawHit->setHitId(numOfHits);
		fstRawHit->setDefaultTb(FST::defaultTimeBin);
		fstRawHit->setIsHit(isHit && !isPed);

		// hit container for IST RawHits and FST RawHits & Peds
		if(fstRawHit->getLayer() == 3)
		{
		  rawHitVec_orig.push_back(fstRawHit);
		  numOfHits++;
		}
		if(fstRawHit->getLayer() != 0 && isHit && !isPed) 
		{
		  rawHitVec_orig.push_back(fstRawHit); 
		  numOfHits++;
		}

		// hit container for IST & FST RawHits only
		if(isHit && !isPed) 
		{ 
		  rawHitVec_Used.push_back(fstRawHit);
		  numOfHits_Used++;
		  if(fstRawHit->getLayer() != 0)
		  {
		    numOfHits_IST++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    // if(numOfHits > 0 && numOfHits <= FST::maxNHitsPerEvent) // maximum hits to expect per event is 10
    // if(numOfHits_Used > 0 && numOfHits_Used <= FST::maxNHitsPerEvent) // maximum hits to expect per event is 10
    if(numOfHits_IST > 0 && numOfHits_IST <= FST::maxNHitsPerEvent) // maximum hits to expect per event is 20
    {
      fillHitDisplay(rawHitVec_orig); // fill hit display

      // set up FstRawHit with Hits & Peds
      mFstEvent->clearRawHitsList();
      int nHits = rawHitVec_orig.size();
      int numOfFstHits = 0;
      for(int i_hit = 0; i_hit < nHits; ++i_hit)
      {
	mFstRawHit = mFstEvent->createRawHit();
	mFstRawHit->setLayer(rawHitVec_orig[i_hit]->getLayer());
	mFstRawHit->setSensor(rawHitVec_orig[i_hit]->getSensor());
	mFstRawHit->setAPV(rawHitVec_orig[i_hit]->getAPV());
	mFstRawHit->setChannel(rawHitVec_orig[i_hit]->getChannel());
	mFstRawHit->setColumn(rawHitVec_orig[i_hit]->getColumn());
	mFstRawHit->setRow(rawHitVec_orig[i_hit]->getRow());
	mFstRawHit->setPosX(rawHitVec_orig[i_hit]->getPosX());
	mFstRawHit->setPosY(rawHitVec_orig[i_hit]->getPosY());
	for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	{
	  mFstRawHit->setPedMean(rawHitVec_orig[i_hit]->getPedMean(i_tb), i_tb);
	  mFstRawHit->setPedStdDev(rawHitVec_orig[i_hit]->getPedStdDev(i_tb), i_tb);
	  mFstRawHit->setPedRMS(rawHitVec_orig[i_hit]->getPedRMS(i_tb), i_tb);
	  mFstRawHit->setCMNStdDev(rawHitVec_orig[i_hit]->getCMNStdDev(i_tb), i_tb);
	  mFstRawHit->setRanStdDev(rawHitVec_orig[i_hit]->getRanStdDev(i_tb), i_tb);
	  mFstRawHit->setRawCharge(rawHitVec_orig[i_hit]->getRawCharge(i_tb), i_tb);
	  mFstRawHit->setCharge(rawHitVec_orig[i_hit]->getCharge(i_tb), i_tb);
	}
	mFstRawHit->setWeight(rawHitVec_orig[i_hit]->getWeight());
	mFstRawHit->setMaxTb(rawHitVec_orig[i_hit]->getMaxTb());
	mFstRawHit->setHitId(rawHitVec_orig[i_hit]->getHitId());
	mFstRawHit->setDefaultTb(rawHitVec_orig[i_hit]->getDefaultTb());
	mFstRawHit->setIsHit(rawHitVec_orig[i_hit]->getIsHit());
	if(mFstRawHit->getLayer() == 3) numOfFstHits++;
      }
      mFstEvent->setNumFstRawHits(numOfFstHits);

      // hit container for clusters
      std::vector<FstRawHit *> rawHitVec_Cluster;
      rawHitVec_Cluster.clear();
      if(mSavePed)
      { // use seed and hits for FST clusters
	for(int i_hit = 0; i_hit < rawHitVec_orig.size(); ++i_hit)
	{
	  rawHitVec_Cluster.push_back(rawHitVec_orig[i_hit]);
	}
      }
      if(!mSavePed)
      { // use seed only for FST clusters
	for(int i_hit = 0; i_hit < rawHitVec_Used.size(); ++i_hit)
	{
	  rawHitVec_Cluster.push_back(rawHitVec_Used[i_hit]);
	}
      }

      // set up FstCluster with cluster from Hits
      mFstEvent->clearClustersList();
      // std::vector<FstCluster *> cluster_simple = findCluster_Simple(rawHitVec_Cluster);
      // std::vector<FstCluster *> cluster_simple = findCluster_Simple(rawHitVec_orig);
      std::vector<FstCluster *> cluster_simple = findCluster_Simple(rawHitVec_Used);
      int nClusters_simple = cluster_simple.size();
      int numOfFstClusters_simple = 0;
      for(int i_cluster = 0; i_cluster < nClusters_simple; ++i_cluster)
      {
	mFstCluster = mFstEvent->createCluster();
	mFstCluster->setLayer(cluster_simple[i_cluster]->getLayer());
	mFstCluster->setSensor(cluster_simple[i_cluster]->getSensor());
	mFstCluster->setMeanColumn(cluster_simple[i_cluster]->getMeanColumn());
	mFstCluster->setMeanRow(cluster_simple[i_cluster]->getMeanRow());
	mFstCluster->setMeanX(cluster_simple[i_cluster]->getMeanX());
	mFstCluster->setMeanY(cluster_simple[i_cluster]->getMeanY());
	mFstCluster->setTotCharge(cluster_simple[i_cluster]->getTotCharge());
	mFstCluster->setMaxTb(cluster_simple[i_cluster]->getMaxTb());
	mFstCluster->setClusterType(cluster_simple[i_cluster]->getClusterType());
	mFstCluster->setNRawHits(cluster_simple[i_cluster]->getNRawHits());
	mFstCluster->setNRawHitsR(cluster_simple[i_cluster]->getNRawHitsR());
	mFstCluster->setNRawHitsPhi(cluster_simple[i_cluster]->getNRawHitsPhi());
	mFstCluster->setClusterId(cluster_simple[i_cluster]->getClusterId());
	
	std::vector<FstRawHit *> rawHitsVec = cluster_simple[i_cluster]->getRawHitVec(); // get raw hits vec from cluster finder
	mFstCluster->clearRawHitsList(); // set up FstRawHit in the FstCluster
	for(int i_hit = 0; i_hit < rawHitsVec.size(); ++i_hit)
	{
	  // mFstCluster->setHitId(i_hit, rawHitsVec[i_hit]->getHitId());
	  mFstClusteredRawHit = mFstCluster->createRawHit();
	  mFstClusteredRawHit->setLayer(rawHitsVec[i_hit]->getLayer());
	  mFstClusteredRawHit->setSensor(rawHitsVec[i_hit]->getSensor());
	  mFstClusteredRawHit->setAPV(rawHitsVec[i_hit]->getAPV());
	  mFstClusteredRawHit->setChannel(rawHitsVec[i_hit]->getChannel());
	  mFstClusteredRawHit->setColumn(rawHitsVec[i_hit]->getColumn());
	  mFstClusteredRawHit->setRow(rawHitsVec[i_hit]->getRow());
	  mFstClusteredRawHit->setPosX(rawHitsVec[i_hit]->getPosX());
	  mFstClusteredRawHit->setPosY(rawHitsVec[i_hit]->getPosY());
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    mFstClusteredRawHit->setPedMean(rawHitsVec[i_hit]->getPedMean(i_tb), i_tb);
	    mFstClusteredRawHit->setPedStdDev(rawHitsVec[i_hit]->getPedStdDev(i_tb), i_tb);
	    mFstClusteredRawHit->setPedRMS(rawHitsVec[i_hit]->getPedRMS(i_tb), i_tb);
	    mFstClusteredRawHit->setCMNStdDev(rawHitsVec[i_hit]->getCMNStdDev(i_tb), i_tb);
	    mFstClusteredRawHit->setRanStdDev(rawHitsVec[i_hit]->getRanStdDev(i_tb), i_tb);
	    mFstClusteredRawHit->setRawCharge(rawHitsVec[i_hit]->getRawCharge(i_tb), i_tb);
	    mFstClusteredRawHit->setCharge(rawHitsVec[i_hit]->getCharge(i_tb), i_tb);
	  }
	  mFstClusteredRawHit->setWeight(rawHitsVec[i_hit]->getWeight());
	  mFstClusteredRawHit->setMaxTb(rawHitsVec[i_hit]->getMaxTb());
	  mFstClusteredRawHit->setHitId(rawHitsVec[i_hit]->getHitId());
	  mFstClusteredRawHit->setDefaultTb(rawHitsVec[i_hit]->getDefaultTb());
	  mFstClusteredRawHit->setIsHit(rawHitsVec[i_hit]->getIsHit());
	  if( rawHitsVec[i_hit]->getIsHit() ) mFstCluster->setIsSeed(rawHitsVec[i_hit]->getIsHit()); // clusters with seed
	}
	if(mFstCluster->getLayer() == 3) numOfFstClusters_simple++;
      }
      mFstEvent->setNumFstClusters_Simple(numOfFstClusters_simple);

      // std::vector<FstCluster *> cluster_scan = findCluster_Scan(rawHitVec_Cluster);
      // std::vector<FstCluster *> cluster_scan = findCluster_ScanWeight(rawHitVec_Cluster);
      std::vector<FstCluster *> cluster_scan = findCluster_ScanRadius(rawHitVec_Cluster);
      int nClusters_scan = cluster_scan.size();
      int numOfFstClusters_scan = 0;
      for(int i_cluster = 0; i_cluster < nClusters_scan; ++i_cluster)
      {
	mFstCluster = mFstEvent->createCluster();
	mFstCluster->setLayer(cluster_scan[i_cluster]->getLayer());
	mFstCluster->setSensor(cluster_scan[i_cluster]->getSensor());
	mFstCluster->setMeanColumn(cluster_scan[i_cluster]->getMeanColumn());
	mFstCluster->setMeanRow(cluster_scan[i_cluster]->getMeanRow());
	mFstCluster->setMeanX(cluster_scan[i_cluster]->getMeanX());
	mFstCluster->setMeanY(cluster_scan[i_cluster]->getMeanY());
	mFstCluster->setTotCharge(cluster_scan[i_cluster]->getTotCharge());
	mFstCluster->setMaxTb(cluster_scan[i_cluster]->getMaxTb());
	mFstCluster->setClusterType(cluster_scan[i_cluster]->getClusterType());
	mFstCluster->setNRawHits(cluster_scan[i_cluster]->getNRawHits());
	mFstCluster->setNRawHitsR(cluster_scan[i_cluster]->getNRawHitsR());
	mFstCluster->setNRawHitsPhi(cluster_scan[i_cluster]->getNRawHitsPhi());
	mFstCluster->setClusterId(cluster_scan[i_cluster]->getClusterId());
	
	std::vector<FstRawHit *> rawHitsVec = cluster_scan[i_cluster]->getRawHitVec(); // get raw hits vec from cluster finder
	mFstCluster->clearRawHitsList(); // set up FstRawHit in the FstCluster
	for(int i_hit = 0; i_hit < rawHitsVec.size(); ++i_hit)
	{
	  // mFstCluster->setHitId(i_hit, rawHitsVec[i_hit]->getHitId());
	  mFstClusteredRawHit = mFstCluster->createRawHit();
	  mFstClusteredRawHit->setLayer(rawHitsVec[i_hit]->getLayer());
	  mFstClusteredRawHit->setSensor(rawHitsVec[i_hit]->getSensor());
	  mFstClusteredRawHit->setAPV(rawHitsVec[i_hit]->getAPV());
	  mFstClusteredRawHit->setChannel(rawHitsVec[i_hit]->getChannel());
	  mFstClusteredRawHit->setColumn(rawHitsVec[i_hit]->getColumn());
	  mFstClusteredRawHit->setRow(rawHitsVec[i_hit]->getRow());
	  mFstClusteredRawHit->setPosX(rawHitsVec[i_hit]->getPosX());
	  mFstClusteredRawHit->setPosY(rawHitsVec[i_hit]->getPosY());
	  for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
	  {
	    mFstClusteredRawHit->setPedMean(rawHitsVec[i_hit]->getPedMean(i_tb), i_tb);
	    mFstClusteredRawHit->setPedStdDev(rawHitsVec[i_hit]->getPedStdDev(i_tb), i_tb);
	    mFstClusteredRawHit->setPedRMS(rawHitsVec[i_hit]->getPedRMS(i_tb), i_tb);
	    mFstClusteredRawHit->setCMNStdDev(rawHitsVec[i_hit]->getCMNStdDev(i_tb), i_tb);
	    mFstClusteredRawHit->setRanStdDev(rawHitsVec[i_hit]->getRanStdDev(i_tb), i_tb);
	    mFstClusteredRawHit->setRawCharge(rawHitsVec[i_hit]->getRawCharge(i_tb), i_tb);
	    mFstClusteredRawHit->setCharge(rawHitsVec[i_hit]->getCharge(i_tb), i_tb);
	  }
	  mFstClusteredRawHit->setWeight(rawHitsVec[i_hit]->getWeight());
	  mFstClusteredRawHit->setMaxTb(rawHitsVec[i_hit]->getMaxTb());
	  mFstClusteredRawHit->setHitId(rawHitsVec[i_hit]->getHitId());
	  mFstClusteredRawHit->setDefaultTb(rawHitsVec[i_hit]->getDefaultTb());
	  mFstClusteredRawHit->setIsHit(rawHitsVec[i_hit]->getIsHit());
	  if( rawHitsVec[i_hit]->getIsHit() ) mFstCluster->setIsSeed(rawHitsVec[i_hit]->getIsHit()); // clusters with seed
	}
	if(mFstCluster->getLayer() == 3) numOfFstClusters_scan++;
      }
      mFstEvent->setNumFstClusters_Scan(numOfFstClusters_scan);

      // set up FstTrack 
      mFstEvent->clearTracksList(); // FstTrack
      std::vector<FstTrack *> fstTrackVec_Hits = findTrack_Hits(rawHitVec_Cluster); // find tracks with Hits
      for(int i_track = 0; i_track < fstTrackVec_Hits.size(); ++i_track)
      { // get track from hits
	// if(mFstEvent->getNumTracks() > 30) continue; // only save first 30 tracks | will drop the event with more than 10 trakcs anyway
	mFstTrack = mFstEvent->createTrack();
	mFstTrack->setTrackId(fstTrackVec_Hits[i_track]->getTrackId());
	mFstTrack->setTrackType(fstTrackVec_Hits[i_track]->getTrackType());
	mFstTrack->setIdIst1(fstTrackVec_Hits[i_track]->getIdIst1());
	mFstTrack->setOrigPosIst1(fstTrackVec_Hits[i_track]->getOrigPosIst1());
	mFstTrack->setAlignedIst1(fstTrackVec_Hits[i_track]->getAlignedIst1());
	mFstTrack->setIdIst3(fstTrackVec_Hits[i_track]->getIdIst3());
	mFstTrack->setOrigPosIst3(fstTrackVec_Hits[i_track]->getOrigPosIst3());
	mFstTrack->setAlignedIst3(fstTrackVec_Hits[i_track]->getAlignedIst3());
	mFstTrack->setProjIst2(fstTrackVec_Hits[i_track]->getProjIst2());
	for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
	{
	  mFstTrack->setProjFst(i_sensor, fstTrackVec_Hits[i_track]->getProjFst(i_sensor));
	}
      }
      mFstEvent->setNumRawHitTraks(fstTrackVec_Hits.size());

      std::vector<FstTrack *> fstTrackVec_Clusters = findTrack_Clusters(cluster_simple); // find tracks with Clusters
      for(int i_track = 0; i_track < fstTrackVec_Clusters.size(); ++i_track)
      { // get track from clusters
	mFstTrack = mFstEvent->createTrack();
	// if(mFstEvent->getNumTracks() > 30) continue; // only save first 30 tracks | will drop the event with more than 10 trakcs anyway
	mFstTrack->setTrackId(fstTrackVec_Clusters[i_track]->getTrackId());
	mFstTrack->setTrackType(fstTrackVec_Clusters[i_track]->getTrackType());
	mFstTrack->setIdIst1(fstTrackVec_Clusters[i_track]->getIdIst1());
	mFstTrack->setOrigPosIst1(fstTrackVec_Clusters[i_track]->getOrigPosIst1());
	mFstTrack->setAlignedIst1(fstTrackVec_Clusters[i_track]->getAlignedIst1());
	mFstTrack->setIdIst3(fstTrackVec_Clusters[i_track]->getIdIst3());
	mFstTrack->setOrigPosIst3(fstTrackVec_Clusters[i_track]->getOrigPosIst3());
	mFstTrack->setAlignedIst3(fstTrackVec_Clusters[i_track]->getAlignedIst3());
	mFstTrack->setProjIst2(fstTrackVec_Clusters[i_track]->getProjIst2());
	for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
	{
	  mFstTrack->setProjFst(i_sensor, fstTrackVec_Clusters[i_track]->getProjFst(i_sensor));
	}
      }
      mFstEvent->setNumClusterTraks(fstTrackVec_Clusters.size());
      mFstEvent->setEventId(i_event);
      mTree_FstEvent->Fill();
      
    }
  }
  mTree_FstEvent->Print();
  cout << "processed events:  " << NumOfEvents << "/" << NumOfEvents << endl;

  return 1;
}

int FstClusterMaker::Finish()
{
  cout << "FstClusterMaker::Finish => " << endl;
  writePedestal();
  writeHitDisplay();
  writeTree();

  return 1;
}

//--------------pedestal---------------------
bool FstClusterMaker::clearPedestal()
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

bool FstClusterMaker::clearCMN()
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

bool FstClusterMaker::initPedestal()
{
  cout << "FstClusterMaker::initPedestal -> " << endl;

  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string gName = Form("g_mPedMean_Off_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mPedMean[i_layer][i_tb] = new TGraph();
      g_mPedMean[i_layer][i_tb]->SetName(gName.c_str());

      gName = Form("g_mPedSigma_Off_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mPedSigma[i_layer][i_tb] = new TGraph();
      g_mPedSigma[i_layer][i_tb]->SetName(gName.c_str());

      gName = Form("g_mCMNSigma_Off_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mCMNSigma[i_layer][i_tb] = new TGraph();
      g_mCMNSigma[i_layer][i_tb]->SetName(gName.c_str());

      gName = Form("g_mRanSigma_Off_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRanSigma[i_layer][i_tb] = new TGraph();
      g_mRanSigma[i_layer][i_tb]->SetName(gName.c_str());

      std::string HistName = Form("h_mPedDisplay_Off_Layer%d_TimeBin%d",i_layer,i_tb);
// made change to Layer 3 for the fst when using the lser test stand
      if(i_layer == 3) h_mPedDisplay[i_layer][i_tb] = new TH2F(HistName.c_str(),HistName.c_str(),FST::numRStrip,-0.5,FST::numRStrip-0.5,FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);
      else h_mPedDisplay[i_layer][i_tb] = new TH2F(HistName.c_str(),HistName.c_str(),FST::noColumns,-0.5,FST::noColumns-0.5,FST::noRows,-0.5,FST::noRows-0.5);
    }
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string HistName;
      HistName = Form("h_mPedMean_Off_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedMean_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);

      HistName = Form("h_mPedSigma_Off_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedSigma_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);

      HistName = Form("h_mCMNSigma_Off_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mCMNSigma_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);

      HistName = Form("h_mRanSigma_Off_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mRanSigma_FST[i_rstrip][i_tb] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);
    }
  }

  return clearPedestal();
}

bool FstClusterMaker::calPedestal()
{
  cout << "FstClusterMaker::calPedestal -> " << endl;
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
	    if(layer == 3 && col > -1) 
	    { // FST
	      h_mPedMean_FST[col][i_tb]->SetBinContent(row+1,mPed[i_arm][i_port][i_apv][i_ch][i_tb]);
	      h_mPedSigma_FST[col][i_tb]->SetBinContent(row+1,mPedStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	    }
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
		if(layer < 3) col = col%2;
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
	      { // calculate CMN for each rstrip
		double CMN_evt = sumValues_evt[i_arm][i_port][i_apv][i_rstrip][i_tb]/counters_evt[i_arm][i_port][i_apv][i_rstrip][i_tb];
		sumValues_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb] += CMN_evt;
		sumValuesSquared_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb] += CMN_evt*CMN_evt;
		counters_CMN[i_arm][i_port][i_apv][i_rstrip][i_tb]++;
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
	    if(layer < 3) col = col%2;
	    int row = this->getRow(i_arm,i_port,i_apv,i_ch);
	    if(counters_CMN[i_arm][i_port][i_apv][col][i_tb] > 0 && col > -1)
	    { // set CMN for each channel
	      mCMNStdDev[i_arm][i_port][i_apv][i_ch][i_tb] = mCMNMean[i_arm][i_port][i_apv][col][i_tb];
	    }
	    g_mCMNSigma[layer][i_tb]->SetPoint(i_apv*FST::numChannels+i_ch,i_apv*FST::numChannels+i_ch,mCMNStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	    if(layer == 3 && col > -1) h_mCMNSigma_FST[col][i_tb]->SetBinContent(row+1,mCMNStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
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
		if(layer < 3) col = col%2;
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
	      if(layer < 3) col = col%2;
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
	    if(layer == 3 && col > -1) h_mRanSigma_FST[col][i_tb]->SetBinContent(row+1,mRanStdDev[i_arm][i_port][i_apv][i_ch][i_tb]);
	  }
	}
      }
    }
  }
  //--------------------------------------------------------

  return true;
}

void FstClusterMaker::writePedestal()
{
  cout << "FstClusterMaker::writePedestal => save ped mean & sigma!" << endl;

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
    }
  }
}
//--------------pedestal---------------------

//--------------hit display---------------------
bool FstClusterMaker::initHitDisplay()
{
  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    std::string HistName = Form("h_mHitDisplay_Layer%d",i_layer);
    if(i_layer == 3) h_mHitDisplay[i_layer] = new TH2F(HistName.c_str(),HistName.c_str(),FST::numRStrip,-0.5,FST::numRStrip-0.5,FST::numPhiSeg,-0.5,FST::numPhiSeg-0.5);
    else h_mHitDisplay[i_layer] = new TH2F(HistName.c_str(),HistName.c_str(),FST::noColumns,-0.5,FST::noColumns-0.5,FST::noRows,-0.5,FST::noRows-0.5);
    
    HistName = Form("h_mMaxTb_Layer%d",i_layer);
    h_mMaxTb[i_layer] = new TH1F(HistName.c_str(),HistName.c_str(),10,-0.5,9.5); 
  }

  return true;
}

void FstClusterMaker::fillHitDisplay(std::vector<FstRawHit *> rawHitsVec)
{
  for(int i_hit = 0; i_hit < rawHitsVec.size(); ++i_hit)
  {
    int layer = rawHitsVec[i_hit]->getLayer();
    h_mHitDisplay[layer]->Fill(rawHitsVec[i_hit]->getColumn(),rawHitsVec[i_hit]->getRow());
    h_mMaxTb[layer]->Fill(rawHitsVec[i_hit]->getMaxTb());
  }
}

void FstClusterMaker::writeHitDisplay()
{
  cout << "FstClusterMaker::writeHitDisplay => save Hits at each Layer!" << endl;
  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    h_mHitDisplay[i_layer]->Write();
    h_mMaxTb[i_layer]->Write();
  }
}
//--------------hit display---------------------

//--------------cluster with Simple Algorithm---------------------
std::vector<FstCluster *> FstClusterMaker::findCluster_Simple(std::vector<FstRawHit *> rawHitsVec_orig)
{
  int nRawHits = 1, nRawHitsR = 1, nRawHitsPhi = 1;

  int numOfHits = rawHitsVec_orig.size();
  std::vector<FstRawHit *> rawHitsVec;
  rawHitsVec.clear();
  rawHitsVec.reserve(rawHitsVec_orig.size());
  for(int i_hit = 0; i_hit < numOfHits; ++i_hit)
  { // set temp ist hit container
    rawHitsVec.push_back(rawHitsVec_orig[i_hit]);
  }

  // find cluster
  std::vector<FstCluster *> clustersVec_Simple;
  clustersVec_Simple.clear();
  int clusterId = 0;

  if(rawHitsVec.size() > 0)
  {
    // set 1st cluster to 1st hit
    std::vector<FstRawHit *>::iterator rawHitsIt = rawHitsVec.begin();

    FstCluster *fstCluster_tmp = new FstCluster();
    fstCluster_tmp->Clear();
    fstCluster_tmp->setLayer((*rawHitsIt)->getLayer());
    fstCluster_tmp->setSensor((*rawHitsIt)->getSensor());
    fstCluster_tmp->setMeanColumn((*rawHitsIt)->getColumn());
    fstCluster_tmp->setMeanRow((*rawHitsIt)->getRow());
    fstCluster_tmp->setMeanX((*rawHitsIt)->getPosX());
    fstCluster_tmp->setMeanY((*rawHitsIt)->getPosY());
    fstCluster_tmp->setTotCharge((*rawHitsIt)->getCharge((*rawHitsIt)->getMaxTb()));
    fstCluster_tmp->setMaxTb((*rawHitsIt)->getMaxTb());
    fstCluster_tmp->setClusterType(1);
    fstCluster_tmp->setNRawHits(nRawHits);
    fstCluster_tmp->setNRawHitsR(nRawHitsR);
    fstCluster_tmp->setNRawHitsPhi(nRawHitsPhi);
    fstCluster_tmp->addRawHit((*rawHitsIt)); // save hits for ith cluster
    fstCluster_tmp->setClusterId(clusterId);

    clustersVec_Simple.push_back(fstCluster_tmp);

    rawHitsVec.erase(rawHitsIt); // remove 1st hit from ist hit container

    if( rawHitsVec.size() !=0 )
    {
      double weight, tempSumAdc;
      std::vector<FstCluster *>::iterator clusterIt = clustersVec_Simple.begin(); // get the begin of cluster container

      while( clusterIt != clustersVec_Simple.end() && !rawHitsVec.empty() )
      { //loop the existed clusters vector constainer
	std::vector<FstRawHit *> clusteredHitsVec = (*clusterIt)->getRawHitVec(); // get hits stored in current cluster
	for(std::vector<FstRawHit *>::iterator clusteredHitsIt = clusteredHitsVec.begin(); clusteredHitsIt != clusteredHitsVec.end(); clusteredHitsIt++) 
	{ // loop over the hit belong to ith cluster
	  rawHitsIt = rawHitsVec.begin(); // get hit from the rest of rawHitsVec

	  while( rawHitsIt != rawHitsVec.end() ) 
	  { // loop over the rest of hits in rawHitsVec
	    if( ( ((*rawHitsIt)->getLayer() == (*clusteredHitsIt)->getLayer()) && ((*rawHitsIt)->getSensor() == (*clusteredHitsIt)->getSensor()) && ((*rawHitsIt)->getRow() == (*clusteredHitsIt)->getRow()) && ( ((*rawHitsIt)->getColumn() == (*clusteredHitsIt)->getColumn() + 1) || ((*rawHitsIt)->getColumn() == (*clusteredHitsIt)->getColumn() - 1) )) ||
		( ((*rawHitsIt)->getLayer() == (*clusteredHitsIt)->getLayer()) && ((*rawHitsIt)->getSensor() == (*clusteredHitsIt)->getSensor()) && ((*rawHitsIt)->getColumn() == (*clusteredHitsIt)->getColumn()) && ( ((*rawHitsIt)->getRow() == (*clusteredHitsIt)->getRow() + 1) || ((*rawHitsIt)->getRow() == (*clusteredHitsIt)->getRow() - 1) ))   ) 
	    {
	      nRawHits = (*clusterIt)->getNRawHits() + 1;
	      if( (*rawHitsIt)->getRow() == (*clusteredHitsIt)->getRow() ) nRawHitsR = (*clusterIt)->getNRawHitsR() + 1; //same phi
	      if( (*rawHitsIt)->getColumn() == (*clusteredHitsIt)->getColumn() ) nRawHitsPhi = (*clusterIt)->getNRawHitsPhi() + 1; //same R

	      int maxTb_temp = (*rawHitsIt)->getMaxTb();
	      double currentAdc = (*rawHitsIt)->getCharge(maxTb_temp);
	      tempSumAdc = (*clusterIt)->getTotCharge() + currentAdc;
	      weight = currentAdc/tempSumAdc;

	      int layer_temp         = (*clusterIt)->getLayer();
	      int sensor_temp        = (*clusterIt)->getSensor();
	      double meanColumn_temp = (1.0 - weight) * (*clusterIt)->getMeanColumn() + weight * (*rawHitsIt)->getColumn();
	      double meanRow_temp    = (1.0 - weight) * (*clusterIt)->getMeanRow()    + weight * (*rawHitsIt)->getRow();
	      double meanX_temp      = (1.0 - weight) * (*clusterIt)->getMeanX()      + weight * (*rawHitsIt)->getPosX();
	      double meanY_temp      = (1.0 - weight) * (*clusterIt)->getMeanY()      + weight * (*rawHitsIt)->getPosY();
	      double meanTb_temp     = (1.0 - weight) * (*clusterIt)->getMaxTb()      + weight * (*rawHitsIt)->getMaxTb();
	      double totAdc_temp     = tempSumAdc;
	      // double meanX_temp      = this->getMeanX(layer_temp,meanColumn_temp);
	      // double meanY_temp      = this->getMeanY(layer_temp,meanRow_temp);

	      (*clusterIt)->setLayer(layer_temp);
	      (*clusterIt)->setSensor(sensor_temp);
	      (*clusterIt)->setMeanColumn(meanColumn_temp);
	      (*clusterIt)->setMeanRow(meanRow_temp);
	      (*clusterIt)->setMeanX(meanX_temp);
	      (*clusterIt)->setMeanY(meanY_temp);
	      (*clusterIt)->setTotCharge(totAdc_temp);
	      (*clusterIt)->setMaxTb(meanTb_temp);
	      (*clusterIt)->setClusterType(1);
	      (*clusterIt)->setNRawHits(nRawHits);
	      (*clusterIt)->setNRawHitsR(nRawHitsR);
	      (*clusterIt)->setNRawHitsPhi(nRawHitsPhi);
	      (*clusterIt)->addRawHit((*rawHitsIt));
	      (*clusterIt)->setClusterId(clusterId);

	      //include the hit to the cluster's component vector
	      int itPosition = std::distance(clusteredHitsVec.begin(), clusteredHitsIt);
	      clusteredHitsVec = (*clusterIt)->getRawHitVec();
	      clusteredHitsIt = clusteredHitsVec.begin() + itPosition;

	      //remove the clustered ith raw hit from the raw hits list
	      int distance = std::distance(rawHitsVec.begin(), rawHitsIt);
	      rawHitsVec.erase(rawHitsIt);
	      rawHitsIt = rawHitsVec.begin() + distance;
	    }
	    else
	    {
	      rawHitsIt++;
	    }
	  }
	}

	//if the rawHitsIt_th hit does not belong to the existed ith clusters then create a new cluster.
	if(rawHitsVec.size() < 1) continue;

	rawHitsIt = rawHitsVec.begin();
	clusterId++;

	FstCluster *fstCluster_next = new FstCluster();
	fstCluster_next->Clear();
	fstCluster_next->setLayer((*rawHitsIt)->getLayer());
	fstCluster_next->setSensor((*rawHitsIt)->getSensor());
	fstCluster_next->setMeanColumn((*rawHitsIt)->getColumn());
	fstCluster_next->setMeanRow((*rawHitsIt)->getRow());
	fstCluster_next->setMeanX((*rawHitsIt)->getPosX());
	fstCluster_next->setMeanY((*rawHitsIt)->getPosY());
	fstCluster_next->setTotCharge((*rawHitsIt)->getCharge((*rawHitsIt)->getMaxTb()));
	fstCluster_next->setMaxTb((*rawHitsIt)->getMaxTb());
	fstCluster_next->setClusterType(1);
	fstCluster_next->setNRawHits(1);
	fstCluster_next->setNRawHitsR(1);
	fstCluster_next->setNRawHitsPhi(1);
	fstCluster_next->addRawHit((*rawHitsIt)); // save hits for ith cluster
	fstCluster_next->setClusterId(clusterId);

	rawHitsVec.erase(rawHitsIt); //remove the new clustered 1st hit from the hits list

	int distanceCluster = std::distance(clustersVec_Simple.begin(), clusterIt);
	clustersVec_Simple.push_back(fstCluster_next); // push new cluster to container
	clusterIt = clustersVec_Simple.begin() + distanceCluster;
	clusterIt++;
      }
    }
  }

  // cout << "rawHitsVec_orig.size = " << rawHitsVec_orig.size() << ", clustersVec_Simple.size = " << clustersVec_Simple.size() << endl;

  return clustersVec_Simple;
}
//--------------cluster with Simple Algorithm---------------------

//--------------cluster with Scan Algorithm---------------------
std::vector<FstCluster *> FstClusterMaker::findCluster_Scan(std::vector<FstRawHit *> rawHitsVec_orig)
{
  // temp vector for raw hits and clusters
  // Merge all hits in the same phi-seg for FST Cosmic Test Stand
  // => might separate clusters in R strip in the future
  // then form clusters in neighboring phi seg
  // applied to sensor 1 only right now => can be extended to the full module
  // need to think geometry more careful when both inner and outer sectors are used
  // => one might need to be fliped by 180 degree w.r.t. to original position 
  std::vector<FstRawHit *> rawHitsVec[FST::mFstNumSensorsPerModule][FST::mFstNumPhiSegPerSensor];
  std::vector<FstCluster *> clustersVec[FST::mFstNumSensorsPerModule][FST::mFstNumPhiSegPerSensor];
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    for(int i_phi = 0; i_phi < FST::mFstNumPhiSegPerSensor; ++i_phi)
    {
      rawHitsVec[i_sensor][i_phi].reserve(FST::mFstNumRstripPerSensor);
      clustersVec[i_sensor][i_phi].reserve(FST::mFstNumRstripPerSensor);
    }
  }

  // copy raw hits to the temp container
  int numOfHits = rawHitsVec_orig.size();
  for(int i_hit = 0; i_hit < numOfHits; ++i_hit)
  {
    int layer = rawHitsVec_orig[i_hit]->getLayer();
    int sensor = rawHitsVec_orig[i_hit]->getSensor();
    int row = rawHitsVec_orig[i_hit]->getRow(); // phi segmentation
    if(layer == 3) rawHitsVec[sensor][row].push_back(rawHitsVec_orig[i_hit]); // findCluster_Scan for FST only
  }

  // start cluster finder
  std::vector<FstCluster *> clustersVec_Scan;
  int clusterId = 100;

  FstRawHit *rawHitFirst = 0;
  FstRawHit *rawHitNext = 0;

  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    // step 1: do clustering for each phi segmentation
    for(int i_phi = 0; i_phi < FST::mFstNumPhiSegPerSensor; ++i_phi)
    {
      if( !rawHitsVec[i_sensor][i_phi].empty() )
      {
	FstCluster *fstClusterTemp = new FstCluster();
	rawHitFirst= rawHitsVec[i_sensor][i_phi].back();
	rawHitsVec[i_sensor][i_phi].pop_back();

	fstClusterTemp->Clear();
	fstClusterTemp->setLayer(rawHitFirst->getLayer());
	fstClusterTemp->setSensor(rawHitFirst->getSensor());
	fstClusterTemp->setMeanColumn(rawHitFirst->getColumn());
	fstClusterTemp->setMeanRow(rawHitFirst->getRow());
	fstClusterTemp->setMeanX(rawHitFirst->getPosX());
	fstClusterTemp->setMeanY(rawHitFirst->getPosY());
	fstClusterTemp->setTotCharge(rawHitFirst->getCharge(rawHitFirst->getMaxTb()));
	fstClusterTemp->setMaxTb(rawHitFirst->getMaxTb());
	fstClusterTemp->setClusterType(2);
	fstClusterTemp->setNRawHits(1);
	fstClusterTemp->setNRawHitsR(1);
	fstClusterTemp->setNRawHitsPhi(1);
	fstClusterTemp->addRawHit(rawHitFirst); // save hits for ith cluster
	fstClusterTemp->setClusterId(clusterId); // save hits for ith cluster

	while( !rawHitsVec[i_sensor][i_phi].empty() ) 
	{ // loop over the rest of hits in rawHitsVec
	  rawHitNext = rawHitsVec[i_sensor][i_phi].back();
	  rawHitsVec[i_sensor][i_phi].pop_back();

	  int nRawHits_temp  = fstClusterTemp->getNRawHits() + 1;
	  int nRawHitsR_temp = fstClusterTemp->getNRawHitsR() + 1; //same phi

	  int maxTb_temp    = rawHitNext->getMaxTb();
	  double currentAdc = rawHitNext->getCharge(maxTb_temp);
	  double tempSumAdc = fstClusterTemp->getTotCharge() + currentAdc;
	  double weight     = currentAdc/tempSumAdc;

	  int layer_temp         = fstClusterTemp->getLayer();
	  int sensor_temp        = fstClusterTemp->getSensor();
	  double meanColumn_temp = (1.0 - weight)*fstClusterTemp->getMeanColumn() + weight*rawHitNext->getColumn();
	  double meanRow_temp    = (1.0 - weight)*fstClusterTemp->getMeanRow()    + weight*rawHitNext->getRow();
	  double meanX_temp      = (1.0 - weight)*fstClusterTemp->getMeanX()      + weight*rawHitNext->getPosX();
	  double meanY_temp      = (1.0 - weight)*fstClusterTemp->getMeanY()      + weight*rawHitNext->getPosY();
	  double meanTb_temp     = (1.0 - weight)*fstClusterTemp->getMaxTb()      + weight*rawHitNext->getMaxTb();
	  double totAdc_temp     = tempSumAdc;
	  // double meanX_temp      = this->getMeanX(layer_temp,meanColumn_temp);
	  // double meanY_temp      = this->getMeanY(layer_temp,meanRow_temp);

	  fstClusterTemp->setLayer(layer_temp);
	  fstClusterTemp->setSensor(sensor_temp);
	  fstClusterTemp->setMeanColumn(meanColumn_temp);
	  fstClusterTemp->setMeanRow(meanRow_temp);
	  fstClusterTemp->setMeanX(meanX_temp);
	  fstClusterTemp->setMeanY(meanY_temp);
	  fstClusterTemp->setTotCharge(totAdc_temp);
	  fstClusterTemp->setMaxTb(meanTb_temp);
	  fstClusterTemp->setClusterType(2);
	  fstClusterTemp->setNRawHits(nRawHitsR_temp);
	  fstClusterTemp->setNRawHitsR(nRawHitsR_temp);
	  fstClusterTemp->addRawHit(rawHitNext);
	  fstClusterTemp->setClusterId(clusterId); // save hits for ith cluster
	}
	clustersVec[i_sensor][i_phi].push_back(fstClusterTemp);
	clusterId++;
      }
    }

    // step 2: do clustering for neighboring phi segmentation
    std::vector<FstCluster *>::iterator clusterIt1, clusterIt2;
    for(int phiIndex1 = 0; phiIndex1 < FST::mFstNumPhiSegPerSensor-1; ++phiIndex1)
    {
      int phiIndex2 = phiIndex1 + 1;

      if(clustersVec[i_sensor][phiIndex1].size() > 0 && clustersVec[i_sensor][phiIndex2].size() > 0)
      {
	for(clusterIt1 = clustersVec[i_sensor][phiIndex1].begin(); clusterIt1 != clustersVec[i_sensor][phiIndex1].end() && !clustersVec[i_sensor][phiIndex1].empty(); clusterIt1++)
	{
	  for(clusterIt2 = clustersVec[i_sensor][phiIndex2].begin(); clusterIt2 != clustersVec[i_sensor][phiIndex2].end() && !clustersVec[i_sensor][phiIndex2].empty(); clusterIt2++)
	  {
	    // distance of mean RStrip between two clusters in neighboring phi segmentation
	    // currently, there is only one cluster for each phi segmentation
	    double columnDistance = (*clusterIt1)->getMeanColumn() - (*clusterIt2)->getMeanColumn();
	    if(TMath::Abs(columnDistance) < 3.5)
	    {
	      clusterId = (*clusterIt1)->getClusterId();
	      if( (*clusterIt1)->getTotCharge() < (*clusterIt2)-> getTotCharge() )
	      {
		clusterId = (*clusterIt2)->getClusterId();
	      }

	      double totCharge_temp  = (*clusterIt1)->getTotCharge() + (*clusterIt2)->getTotCharge();
	      int layer_temp         = (*clusterIt1)->getLayer();
	      int sensor_temp        = (*clusterIt1)->getSensor();
	      double meanColumn_temp = (*clusterIt1)->getMeanColumn()*(*clusterIt1)->getTotCharge()/totCharge_temp + (*clusterIt2)->getMeanColumn()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      double meanRow_temp    = (*clusterIt1)->getMeanRow()*(*clusterIt1)->getTotCharge()/totCharge_temp    + (*clusterIt2)->getMeanRow()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      double meanX_temp      = (*clusterIt1)->getMeanX()*(*clusterIt1)->getTotCharge()/totCharge_temp      + (*clusterIt2)->getMeanX()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      double meanY_temp      = (*clusterIt1)->getMeanY()*(*clusterIt1)->getTotCharge()/totCharge_temp      + (*clusterIt2)->getMeanY()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      double meanTb_temp     = (*clusterIt1)->getMaxTb()*(*clusterIt1)->getTotCharge()/totCharge_temp      + (*clusterIt2)->getMaxTb()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      int nRawHits_temp      = (*clusterIt1)->getNRawHits() + (*clusterIt2)->getNRawHits();
	      int nRawHitsR_temp     = (*clusterIt1)->getNRawHitsR() + (*clusterIt2)->getNRawHitsR();
	      int nRawHitsPhi_temp   = (*clusterIt1)->getNRawHitsPhi() + (*clusterIt2)->getNRawHitsPhi();
	      // double meanX_temp      = this->getMeanX(layer_temp,meanColumn_temp);
	      // double meanY_temp      = this->getMeanY(layer_temp,meanRow_temp);

	      (*clusterIt2)->setLayer(layer_temp);
	      (*clusterIt2)->setSensor(sensor_temp);
	      (*clusterIt2)->setMeanColumn(meanColumn_temp);
	      (*clusterIt2)->setMeanRow(meanRow_temp);
	      (*clusterIt2)->setMeanX(meanX_temp);
	      (*clusterIt2)->setMeanY(meanY_temp);
	      (*clusterIt2)->setTotCharge(totCharge_temp);
	      (*clusterIt2)->setMaxTb(meanTb_temp);
	      (*clusterIt2)->setClusterType(2);
	      (*clusterIt2)->setNRawHits(nRawHits_temp);
	      (*clusterIt2)->setNRawHitsR(nRawHitsR_temp);
	      (*clusterIt2)->setNRawHitsPhi(nRawHitsR_temp);
	      (*clusterIt2)->setClusterId(clusterId);

	      std::vector<FstRawHit *> rawHitsVec_clusterIt1 = (*clusterIt1)->getRawHitVec();
	      for(int i_hit = 0; i_hit < rawHitsVec_clusterIt1.size(); ++i_hit)
	      {
		(*clusterIt2)->addRawHit(rawHitsVec_clusterIt1[i_hit]);
	      }

	      int distance1 = std::distance(clustersVec[i_sensor][phiIndex1].begin(), clusterIt1);
	      clustersVec[i_sensor][phiIndex1].erase(clusterIt1);

	      if(distance1 == 0)
		clusterIt1 = clustersVec[i_sensor][phiIndex1].begin();
	      else
		--clusterIt1;
	    }
	  }
	}
      }
    }
  }

  // fill output container
  std::vector<FstCluster *>::iterator clusterIt;
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    for(int i_phi = 0; i_phi < FST::mFstNumPhiSegPerSensor; ++i_phi)
    {
      if(clustersVec[i_sensor][i_phi].size() <= 0) continue;

      for(clusterIt = clustersVec[i_sensor][i_phi].begin(); clusterIt != clustersVec[i_sensor][i_phi].end(); ++clusterIt)
      {
	clustersVec_Scan.push_back(*clusterIt);
      }

      rawHitsVec[i_sensor][i_phi].clear();
      clustersVec[i_sensor][i_phi].clear();
    }
  }

  return clustersVec_Scan;
}

std::vector<FstCluster *> FstClusterMaker::findCluster_ScanWeight(std::vector<FstRawHit *> rawHitsVec_orig)
{
  // temp vector for raw hits and clusters
  // Merge all hits in the same phi-seg for FST Cosmic Test Stand
  // => might separate clusters in R strip in the future
  // then form clusters in neighboring phi seg
  // applied to sensor 1 only right now => can be extended to the full module
  // need to think geometry more careful when both inner and outer sectors are used
  // => one might need to be fliped by 180 degree w.r.t. to original position 
  std::vector<FstRawHit *> rawHitsVec[FST::mFstNumSensorsPerModule][FST::mFstNumPhiSegPerSensor];
  std::vector<FstCluster *> clustersVec[FST::mFstNumSensorsPerModule][FST::mFstNumPhiSegPerSensor];
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    for(int i_phi = 0; i_phi < FST::mFstNumPhiSegPerSensor; ++i_phi)
    {
      rawHitsVec[i_sensor][i_phi].reserve(FST::mFstNumRstripPerSensor);
      clustersVec[i_sensor][i_phi].reserve(FST::mFstNumRstripPerSensor);
    }
  }

  // copy raw hits to the temp container
  int numOfHits = rawHitsVec_orig.size();
  for(int i_hit = 0; i_hit < numOfHits; ++i_hit)
  {
    int layer = rawHitsVec_orig[i_hit]->getLayer();
    int sensor = rawHitsVec_orig[i_hit]->getSensor();
    int row = rawHitsVec_orig[i_hit]->getRow(); // phi segmentation
    if(layer == 3) rawHitsVec[sensor][row].push_back(rawHitsVec_orig[i_hit]); // findCluster_Scan for FST only
  }

  // start cluster finder
  std::vector<FstCluster *> clustersVec_Scan;
  int clusterId = 200;
  // int clusterType = 3;
  int clusterType = 2;

  FstRawHit *rawHitFirst = 0;
  FstRawHit *rawHitNext = 0;

  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    // step 1: do clustering for each phi segmentation
    for(int i_phi = 0; i_phi < FST::mFstNumPhiSegPerSensor; ++i_phi)
    {
      if( !rawHitsVec[i_sensor][i_phi].empty() )
      {
	FstCluster *fstClusterTemp = new FstCluster();
	rawHitFirst= rawHitsVec[i_sensor][i_phi].back();
	rawHitsVec[i_sensor][i_phi].pop_back();

	fstClusterTemp->Clear();
	fstClusterTemp->setLayer(rawHitFirst->getLayer());
	fstClusterTemp->setSensor(rawHitFirst->getSensor());
	fstClusterTemp->setMeanColumn(rawHitFirst->getColumn());
	fstClusterTemp->setMeanRow(rawHitFirst->getRow());
	fstClusterTemp->setMeanX(rawHitFirst->getPosX());
	fstClusterTemp->setMeanY(rawHitFirst->getPosY());
	// fstClusterTemp->setTotCharge(rawHitFirst->getCharge(rawHitFirst->getMaxTb())/rawHitFirst->getWeight());
	fstClusterTemp->setTotCharge(rawHitFirst->getCharge(rawHitFirst->getMaxTb()));
	fstClusterTemp->setMaxTb(rawHitFirst->getMaxTb());
	fstClusterTemp->setClusterType(clusterType);
	fstClusterTemp->setNRawHits(1);
	fstClusterTemp->setNRawHitsR(1);
	fstClusterTemp->setNRawHitsPhi(1);
	fstClusterTemp->addRawHit(rawHitFirst); // save hits for ith cluster
	fstClusterTemp->setClusterId(clusterId); // save hits for ith cluster

	while( !rawHitsVec[i_sensor][i_phi].empty() ) 
	{ // loop over the rest of hits in rawHitsVec
	  rawHitNext = rawHitsVec[i_sensor][i_phi].back();
	  rawHitsVec[i_sensor][i_phi].pop_back();

	  int nRawHits_temp  = fstClusterTemp->getNRawHits() + 1;
	  int nRawHitsR_temp = fstClusterTemp->getNRawHitsR() + 1; //same phi

	  int maxTb_temp    = rawHitNext->getMaxTb();
	  // double currentAdc = rawHitNext->getCharge(maxTb_temp)/rawHitNext->getWeight();
	  double currentAdc = rawHitNext->getCharge(maxTb_temp);
	  double tempSumAdc = fstClusterTemp->getTotCharge() + currentAdc;

	  double currentAdcWeight = rawHitNext->getCharge(maxTb_temp)/rawHitNext->getWeight();
	  double tempSumAdcWeight = currentAdcWeight;
	  std::vector<FstRawHit *> rawHitsTemp = fstClusterTemp->getRawHitVec(); 
	  for(int i_hit = 0; i_hit < rawHitsTemp.size(); ++i_hit)
	  {
	    tempSumAdcWeight += rawHitsTemp[i_hit]->getCharge(rawHitsTemp[i_hit]->getMaxTb())/rawHitsTemp[i_hit]->getWeight();
	  }
	  double weight = currentAdcWeight/tempSumAdcWeight;

	  int layer_temp         = fstClusterTemp->getLayer();
	  int sensor_temp        = fstClusterTemp->getSensor();
	  double meanColumn_temp = (1.0 - weight)*fstClusterTemp->getMeanColumn() + weight*rawHitNext->getColumn();
	  double meanRow_temp    = (1.0 - weight)*fstClusterTemp->getMeanRow()    + weight*rawHitNext->getRow();
	  double meanX_temp      = (1.0 - weight)*fstClusterTemp->getMeanX()      + weight*rawHitNext->getPosX();
	  double meanY_temp      = (1.0 - weight)*fstClusterTemp->getMeanY()      + weight*rawHitNext->getPosY();
	  double meanTb_temp     = (1.0 - weight)*fstClusterTemp->getMaxTb()      + weight*rawHitNext->getMaxTb();
	  double totAdc_temp     = tempSumAdc;
	  // double meanX_temp      = this->getMeanX(layer_temp,meanColumn_temp);
	  // double meanY_temp      = this->getMeanY(layer_temp,meanRow_temp);

	  fstClusterTemp->setLayer(layer_temp);
	  fstClusterTemp->setSensor(sensor_temp);
	  fstClusterTemp->setMeanColumn(meanColumn_temp);
	  fstClusterTemp->setMeanRow(meanRow_temp);
	  fstClusterTemp->setMeanX(meanX_temp);
	  fstClusterTemp->setMeanY(meanY_temp);
	  fstClusterTemp->setTotCharge(totAdc_temp);
	  fstClusterTemp->setMaxTb(meanTb_temp);
	  fstClusterTemp->setClusterType(clusterType);
	  fstClusterTemp->setNRawHits(nRawHitsR_temp);
	  fstClusterTemp->setNRawHitsR(nRawHitsR_temp);
	  fstClusterTemp->addRawHit(rawHitNext);
	  fstClusterTemp->setClusterId(clusterId); // save hits for ith cluster
	}
	clustersVec[i_sensor][i_phi].push_back(fstClusterTemp);
	clusterId++;
      }
    }

    // step 2: do clustering for neighboring phi segmentation
    std::vector<FstCluster *>::iterator clusterIt1, clusterIt2;
    for(int phiIndex1 = 0; phiIndex1 < FST::mFstNumPhiSegPerSensor-1; ++phiIndex1)
    {
      int phiIndex2 = phiIndex1 + 1;

      if(clustersVec[i_sensor][phiIndex1].size() > 0 && clustersVec[i_sensor][phiIndex2].size() > 0)
      {
	for(clusterIt1 = clustersVec[i_sensor][phiIndex1].begin(); clusterIt1 != clustersVec[i_sensor][phiIndex1].end() && !clustersVec[i_sensor][phiIndex1].empty(); clusterIt1++)
	{
	  for(clusterIt2 = clustersVec[i_sensor][phiIndex2].begin(); clusterIt2 != clustersVec[i_sensor][phiIndex2].end() && !clustersVec[i_sensor][phiIndex2].empty(); clusterIt2++)
	  {
	    // distance of mean RStrip between two clusters in neighboring phi segmentation
	    // currently, there is only one cluster for each phi segmentation
	    double columnDistance = (*clusterIt1)->getMeanColumn() - (*clusterIt2)->getMeanColumn();
	    if(TMath::Abs(columnDistance) < 3.5)
	    {
	      clusterId = (*clusterIt1)->getClusterId();
	      if( (*clusterIt1)->getTotCharge() < (*clusterIt2)-> getTotCharge() )
	      {
		clusterId = (*clusterIt2)->getClusterId();
	      }

	      // calculate weight
	      double totCharge_clusterIt1 = 0.0;
	      std::vector<FstRawHit *> rawHitsClusterIt1 = (*clusterIt1)->getRawHitVec(); 
	      for(int i_hit = 0; i_hit < rawHitsClusterIt1.size(); ++i_hit)
	      {
		int maxTb_temp = rawHitsClusterIt1[i_hit]->getMaxTb();
		totCharge_clusterIt1 += rawHitsClusterIt1[i_hit]->getCharge(maxTb_temp)/rawHitsClusterIt1[i_hit]->getWeight();
	      }
	      double totCharge_clusterIt2 = 0.0;
	      std::vector<FstRawHit *> rawHitsClusterIt2 = (*clusterIt2)->getRawHitVec(); 
	      for(int i_hit = 0; i_hit < rawHitsClusterIt2.size(); ++i_hit)
	      {
		int maxTb_temp = rawHitsClusterIt2[i_hit]->getMaxTb();
		totCharge_clusterIt2 += rawHitsClusterIt2[i_hit]->getCharge(maxTb_temp)/rawHitsClusterIt2[i_hit]->getWeight();
	      }
	      double totCharge_Weight = totCharge_clusterIt1 + totCharge_clusterIt2;

	      double totCharge_temp  = (*clusterIt1)->getTotCharge() + (*clusterIt2)->getTotCharge();
	      int layer_temp         = (*clusterIt1)->getLayer();
	      int sensor_temp        = (*clusterIt1)->getSensor();
	      // double meanColumn_temp = (*clusterIt1)->getMeanColumn()*(*clusterIt1)->getTotCharge()/totCharge_temp + (*clusterIt2)->getMeanColumn()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      // double meanRow_temp    = (*clusterIt1)->getMeanRow()*(*clusterIt1)->getTotCharge()/totCharge_temp    + (*clusterIt2)->getMeanRow()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      // double meanX_temp      = (*clusterIt1)->getMeanX()*(*clusterIt1)->getTotCharge()/totCharge_temp      + (*clusterIt2)->getMeanX()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      // double meanY_temp      = (*clusterIt1)->getMeanY()*(*clusterIt1)->getTotCharge()/totCharge_temp      + (*clusterIt2)->getMeanY()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      // double meanTb_temp     = (*clusterIt1)->getMaxTb()*(*clusterIt1)->getTotCharge()/totCharge_temp      + (*clusterIt2)->getMaxTb()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      double meanColumn_temp = (*clusterIt1)->getMeanColumn()*totCharge_clusterIt1/totCharge_Weight + (*clusterIt2)->getMeanColumn()*totCharge_clusterIt2/totCharge_Weight;
	      double meanRow_temp    = (*clusterIt1)->getMeanRow()*totCharge_clusterIt1/totCharge_Weight    + (*clusterIt2)->getMeanRow()*totCharge_clusterIt2/totCharge_Weight;
	      double meanX_temp      = (*clusterIt1)->getMeanX()*totCharge_clusterIt1/totCharge_Weight      + (*clusterIt2)->getMeanX()*totCharge_clusterIt2/totCharge_Weight;
	      double meanY_temp      = (*clusterIt1)->getMeanY()*totCharge_clusterIt1/totCharge_Weight      + (*clusterIt2)->getMeanY()*totCharge_clusterIt2/totCharge_Weight;
	      double meanTb_temp     = (*clusterIt1)->getMaxTb()*totCharge_clusterIt1/totCharge_Weight      + (*clusterIt2)->getMaxTb()*totCharge_clusterIt2/totCharge_Weight;
	      int nRawHits_temp      = (*clusterIt1)->getNRawHits() + (*clusterIt2)->getNRawHits();
	      int nRawHitsR_temp     = (*clusterIt1)->getNRawHitsR() + (*clusterIt2)->getNRawHitsR();
	      int nRawHitsPhi_temp   = (*clusterIt1)->getNRawHitsPhi() + (*clusterIt2)->getNRawHitsPhi();
	      // double meanX_temp      = this->getMeanX(layer_temp,meanColumn_temp);
	      // double meanY_temp      = this->getMeanY(layer_temp,meanRow_temp);

	      (*clusterIt2)->setLayer(layer_temp);
	      (*clusterIt2)->setSensor(sensor_temp);
	      (*clusterIt2)->setMeanColumn(meanColumn_temp);
	      (*clusterIt2)->setMeanRow(meanRow_temp);
	      (*clusterIt2)->setMeanX(meanX_temp);
	      (*clusterIt2)->setMeanY(meanY_temp);
	      (*clusterIt2)->setTotCharge(totCharge_temp);
	      (*clusterIt2)->setMaxTb(meanTb_temp);
	      (*clusterIt2)->setClusterType(clusterType);
	      (*clusterIt2)->setNRawHits(nRawHits_temp);
	      (*clusterIt2)->setNRawHitsR(nRawHitsR_temp);
	      (*clusterIt2)->setNRawHitsPhi(nRawHitsR_temp);
	      (*clusterIt2)->setClusterId(clusterId);

	      std::vector<FstRawHit *> rawHitsVec_clusterIt1 = (*clusterIt1)->getRawHitVec();
	      for(int i_hit = 0; i_hit < rawHitsVec_clusterIt1.size(); ++i_hit)
	      {
		(*clusterIt2)->addRawHit(rawHitsVec_clusterIt1[i_hit]);
	      }

	      int distance1 = std::distance(clustersVec[i_sensor][phiIndex1].begin(), clusterIt1);
	      clustersVec[i_sensor][phiIndex1].erase(clusterIt1);

	      if(distance1 == 0)
		clusterIt1 = clustersVec[i_sensor][phiIndex1].begin();
	      else
		--clusterIt1;
	    }
	  }
	}
      }
    }
  }

  // fill output container
  std::vector<FstCluster *>::iterator clusterIt;
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    for(int i_phi = 0; i_phi < FST::mFstNumPhiSegPerSensor; ++i_phi)
    {
      if(clustersVec[i_sensor][i_phi].size() <= 0) continue;

      for(clusterIt = clustersVec[i_sensor][i_phi].begin(); clusterIt != clustersVec[i_sensor][i_phi].end(); ++clusterIt)
      {
	clustersVec_Scan.push_back(*clusterIt);
      }

      rawHitsVec[i_sensor][i_phi].clear();
      clustersVec[i_sensor][i_phi].clear();
    }
  }

  return clustersVec_Scan;
}

std::vector<FstCluster *> FstClusterMaker::findCluster_ScanRadius(std::vector<FstRawHit *> rawHitsVec_orig)
{
  // temp vector for raw hits and clusters
  // Merge all hits in the same phi-seg for FST Cosmic Test Stand
  // => might separate clusters in R strip in the future
  // then form clusters in neighboring phi seg
  // applied to sensor 1 only right now => can be extended to the full module
  // need to think geometry more careful when both inner and outer sectors are used
  // => one might need to be fliped by 180 degree w.r.t. to original position 
  std::vector<FstRawHit *> rawHitsVec[FST::mFstNumSensorsPerModule][FST::mFstNumPhiSegPerSensor];
  std::vector<FstCluster *> clustersVec[FST::mFstNumSensorsPerModule][FST::mFstNumPhiSegPerSensor];
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    for(int i_phi = 0; i_phi < FST::mFstNumPhiSegPerSensor; ++i_phi)
    {
      rawHitsVec[i_sensor][i_phi].reserve(FST::mFstNumRstripPerSensor);
      clustersVec[i_sensor][i_phi].reserve(FST::mFstNumRstripPerSensor);
    }
  }

  // copy raw hits to the temp container
  int numOfHits = rawHitsVec_orig.size();
  for(int i_hit = 0; i_hit < numOfHits; ++i_hit)
  {
    int layer = rawHitsVec_orig[i_hit]->getLayer();
    int sensor = rawHitsVec_orig[i_hit]->getSensor();
    int row = rawHitsVec_orig[i_hit]->getRow(); // phi segmentation
    if(layer == 3) rawHitsVec[sensor][row].push_back(rawHitsVec_orig[i_hit]); // findCluster_Scan for FST only
  }

  // start cluster finder
  std::vector<FstCluster *> clustersVec_Scan;
  int clusterId = 300;
  // int clusterType = 4;
  int clusterType = 2;

  FstRawHit *rawHitFirst = 0;
  FstRawHit *rawHitNext = 0;

  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    // step 1: do clustering for each phi segmentation
    for(int i_phi = 0; i_phi < FST::mFstNumPhiSegPerSensor; ++i_phi)
    {
      if( !rawHitsVec[i_sensor][i_phi].empty() )
      {
	FstCluster *fstClusterTemp = new FstCluster();
	rawHitFirst= rawHitsVec[i_sensor][i_phi].back();
	rawHitsVec[i_sensor][i_phi].pop_back();

	fstClusterTemp->Clear();
	fstClusterTemp->setLayer(rawHitFirst->getLayer());
	fstClusterTemp->setSensor(rawHitFirst->getSensor());
	fstClusterTemp->setMeanColumn(rawHitFirst->getColumn());
	fstClusterTemp->setMeanRow(rawHitFirst->getRow());
	fstClusterTemp->setMeanX(rawHitFirst->getPosX());
	fstClusterTemp->setMeanY(rawHitFirst->getPosY());
	fstClusterTemp->setTotCharge(rawHitFirst->getCharge(rawHitFirst->getMaxTb()));
	fstClusterTemp->setMaxTb(rawHitFirst->getMaxTb());
	fstClusterTemp->setClusterType(clusterType);
	fstClusterTemp->setNRawHits(1);
	fstClusterTemp->setNRawHitsR(1);
	fstClusterTemp->setNRawHitsPhi(1);
	fstClusterTemp->addRawHit(rawHitFirst); // save hits for ith cluster
	fstClusterTemp->setClusterId(clusterId); // save hits for ith cluster

	while( !rawHitsVec[i_sensor][i_phi].empty() ) 
	{ // loop over the rest of hits in rawHitsVec
	  rawHitNext = rawHitsVec[i_sensor][i_phi].back();
	  rawHitsVec[i_sensor][i_phi].pop_back();

	  int nRawHits_temp  = fstClusterTemp->getNRawHits() + 1;
	  int nRawHitsR_temp = fstClusterTemp->getNRawHitsR() + 1; //same phi

	  int maxTb_temp    = rawHitNext->getMaxTb();
	  double currentAdc = rawHitNext->getCharge(maxTb_temp);
	  double tempSumAdc = fstClusterTemp->getTotCharge() + currentAdc;
	  double weight     = currentAdc/tempSumAdc;

	  int layer_temp         = fstClusterTemp->getLayer();
	  int sensor_temp        = fstClusterTemp->getSensor();
	  // double meanColumn_temp = (1.0 - weight)*fstClusterTemp->getMeanColumn() + weight*rawHitNext->getColumn();
	  double meanColumn_temp = fstClusterTemp->getMeanColumn(); 
	  if(rawHitNext->getColumn() > fstClusterTemp->getMeanColumn()) meanColumn_temp = rawHitNext->getColumn();
	  double meanRow_temp    = (1.0 - weight)*fstClusterTemp->getMeanRow()    + weight*rawHitNext->getRow();
	  // double meanX_temp      = (1.0 - weight)*fstClusterTemp->getMeanX()      + weight*rawHitNext->getPosX();
	  double meanX_temp      = fstClusterTemp->getMeanX();
	  if(rawHitNext->getPosX() > fstClusterTemp->getMeanX()) meanX_temp = rawHitNext->getPosX();
	  double meanY_temp      = (1.0 - weight)*fstClusterTemp->getMeanY()      + weight*rawHitNext->getPosY();
	  double meanTb_temp     = (1.0 - weight)*fstClusterTemp->getMaxTb()      + weight*rawHitNext->getMaxTb();
	  double totAdc_temp     = tempSumAdc;
	  // double meanX_temp      = this->getMeanX(layer_temp,meanColumn_temp);
	  // double meanY_temp      = this->getMeanY(layer_temp,meanRow_temp);

	  fstClusterTemp->setLayer(layer_temp);
	  fstClusterTemp->setSensor(sensor_temp);
	  fstClusterTemp->setMeanColumn(meanColumn_temp);
	  fstClusterTemp->setMeanRow(meanRow_temp);
	  fstClusterTemp->setMeanX(meanX_temp);
	  fstClusterTemp->setMeanY(meanY_temp);
	  fstClusterTemp->setTotCharge(totAdc_temp);
	  fstClusterTemp->setMaxTb(meanTb_temp);
	  fstClusterTemp->setClusterType(clusterType);
	  fstClusterTemp->setNRawHits(nRawHitsR_temp);
	  fstClusterTemp->setNRawHitsR(nRawHitsR_temp);
	  fstClusterTemp->addRawHit(rawHitNext);
	  fstClusterTemp->setClusterId(clusterId); // save hits for ith cluster
	}
	clustersVec[i_sensor][i_phi].push_back(fstClusterTemp);
	clusterId++;
      }
    }

    // step 2: do clustering for neighboring phi segmentation
    std::vector<FstCluster *>::iterator clusterIt1, clusterIt2;
    for(int phiIndex1 = 0; phiIndex1 < FST::mFstNumPhiSegPerSensor-1; ++phiIndex1)
    {
      int phiIndex2 = phiIndex1 + 1;

      if(clustersVec[i_sensor][phiIndex1].size() > 0 && clustersVec[i_sensor][phiIndex2].size() > 0)
      {
	for(clusterIt1 = clustersVec[i_sensor][phiIndex1].begin(); clusterIt1 != clustersVec[i_sensor][phiIndex1].end() && !clustersVec[i_sensor][phiIndex1].empty(); clusterIt1++)
	{
	  for(clusterIt2 = clustersVec[i_sensor][phiIndex2].begin(); clusterIt2 != clustersVec[i_sensor][phiIndex2].end() && !clustersVec[i_sensor][phiIndex2].empty(); clusterIt2++)
	  {
	    // distance of mean RStrip between two clusters in neighboring phi segmentation
	    // currently, there is only one cluster for each phi segmentation
	    double columnDistance = (*clusterIt1)->getMeanColumn() - (*clusterIt2)->getMeanColumn();
	    if(TMath::Abs(columnDistance) < 3.5)
	    {
	      clusterId = (*clusterIt1)->getClusterId();
	      if( (*clusterIt1)->getTotCharge() < (*clusterIt2)-> getTotCharge() )
	      {
		clusterId = (*clusterIt2)->getClusterId();
	      }

	      double totCharge_temp  = (*clusterIt1)->getTotCharge() + (*clusterIt2)->getTotCharge();
	      int layer_temp         = (*clusterIt1)->getLayer();
	      int sensor_temp        = (*clusterIt1)->getSensor();
	      // double meanColumn_temp = (*clusterIt1)->getMeanColumn()*(*clusterIt1)->getTotCharge()/totCharge_temp + (*clusterIt2)->getMeanColumn()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      double meanColumn_temp = (*clusterIt1)->getMeanColumn();
	      if((*clusterIt2)->getMeanColumn() > (*clusterIt1)->getMeanColumn()) meanColumn_temp = (*clusterIt2)->getMeanColumn();
	      double meanRow_temp    = (*clusterIt1)->getMeanRow()*(*clusterIt1)->getTotCharge()/totCharge_temp    + (*clusterIt2)->getMeanRow()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      // double meanX_temp      = (*clusterIt1)->getMeanX()*(*clusterIt1)->getTotCharge()/totCharge_temp      + (*clusterIt2)->getMeanX()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      double meanX_temp      = (*clusterIt1)->getMeanX();
	      if((*clusterIt2)->getMeanX() > (*clusterIt1)->getMeanX()) meanX_temp = (*clusterIt2)->getMeanX();
	      double meanY_temp      = (*clusterIt1)->getMeanY()*(*clusterIt1)->getTotCharge()/totCharge_temp      + (*clusterIt2)->getMeanY()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      double meanTb_temp     = (*clusterIt1)->getMaxTb()*(*clusterIt1)->getTotCharge()/totCharge_temp      + (*clusterIt2)->getMaxTb()*(*clusterIt2)->getTotCharge()/totCharge_temp;
	      int nRawHits_temp      = (*clusterIt1)->getNRawHits() + (*clusterIt2)->getNRawHits();
	      int nRawHitsR_temp     = (*clusterIt1)->getNRawHitsR() + (*clusterIt2)->getNRawHitsR();
	      int nRawHitsPhi_temp   = (*clusterIt1)->getNRawHitsPhi() + (*clusterIt2)->getNRawHitsPhi();
	      // double meanX_temp      = this->getMeanX(layer_temp,meanColumn_temp);
	      // double meanY_temp      = this->getMeanY(layer_temp,meanRow_temp);

	      (*clusterIt2)->setLayer(layer_temp);
	      (*clusterIt2)->setSensor(sensor_temp);
	      (*clusterIt2)->setMeanColumn(meanColumn_temp);
	      (*clusterIt2)->setMeanRow(meanRow_temp);
	      (*clusterIt2)->setMeanX(meanX_temp);
	      (*clusterIt2)->setMeanY(meanY_temp);
	      (*clusterIt2)->setTotCharge(totCharge_temp);
	      (*clusterIt2)->setMaxTb(meanTb_temp);
	      (*clusterIt2)->setClusterType(clusterType);
	      (*clusterIt2)->setNRawHits(nRawHits_temp);
	      (*clusterIt2)->setNRawHitsR(nRawHitsR_temp);
	      (*clusterIt2)->setNRawHitsPhi(nRawHitsR_temp);
	      (*clusterIt2)->setClusterId(clusterId);

	      std::vector<FstRawHit *> rawHitsVec_clusterIt1 = (*clusterIt1)->getRawHitVec();
	      for(int i_hit = 0; i_hit < rawHitsVec_clusterIt1.size(); ++i_hit)
	      {
		(*clusterIt2)->addRawHit(rawHitsVec_clusterIt1[i_hit]);
	      }

	      int distance1 = std::distance(clustersVec[i_sensor][phiIndex1].begin(), clusterIt1);
	      clustersVec[i_sensor][phiIndex1].erase(clusterIt1);

	      if(distance1 == 0)
		clusterIt1 = clustersVec[i_sensor][phiIndex1].begin();
	      else
		--clusterIt1;
	    }
	  }
	}
      }
    }
  }

  // fill output container
  std::vector<FstCluster *>::iterator clusterIt;
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    for(int i_phi = 0; i_phi < FST::mFstNumPhiSegPerSensor; ++i_phi)
    {
      if(clustersVec[i_sensor][i_phi].size() <= 0) continue;

      for(clusterIt = clustersVec[i_sensor][i_phi].begin(); clusterIt != clustersVec[i_sensor][i_phi].end(); ++clusterIt)
      {
	clustersVec_Scan.push_back(*clusterIt);
      }

      rawHitsVec[i_sensor][i_phi].clear();
      clustersVec[i_sensor][i_phi].clear();
    }
  }

  return clustersVec_Scan;
}
//--------------cluster with Scan Algorithm---------------------

//--------------Find Track with Hits---------------------
std::vector<FstTrack *> FstClusterMaker::findTrack_Hits(std::vector<FstRawHit *> rawHitVec_orig)
{
  int numOfHits[4]; // 0 for fst | 1-3 for ist
  std::vector<FstRawHit *> rawHitVec[4];
  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    numOfHits[i_layer] = 0;
    rawHitVec[i_layer].clear();
  }

  for(int i_hit = 0; i_hit < rawHitVec_orig.size(); ++i_hit)
  { // set temp ist hit container for each layer
    int layer = rawHitVec_orig[i_hit]->getLayer();
    numOfHits[layer]++;
    rawHitVec[layer].push_back(rawHitVec_orig[i_hit]);
  }

  const double z0_fst = FST::pitchLayer03;
  const double z1_ist = FST::pitchLayer12 + FST::pitchLayer23;
  const double z2_ist = FST::pitchLayer23;
  const double z3_ist = 0.0;

  std::vector<FstTrack *> trackVec;
  trackVec.clear();
  int numOfTracks = 0;
  if(numOfHits[1] > 0 && numOfHits[3] > 0)
  { // only do tracking when at least 1 hit is found in ist1 & ist3
    TVector3 vPosIst1_orig, vPosIst3_orig;
    TVector3 vPosIst1, vPosIst3;
    TVector3 vAlignedIst1, vAlignedIst3;
    TVector3 vProjFst, vProjIst2;

    for(int i_ist1 = 0; i_ist1 < numOfHits[1]; ++i_ist1)
    {
      for(int i_ist3 = 0; i_ist3 < numOfHits[3]; ++i_ist3)
      {
	double x1_ist = rawHitVec[1][i_ist1]->getPosX();
	double y1_ist = rawHitVec[1][i_ist1]->getPosY();
	vPosIst1_orig.SetXYZ(x1_ist,y1_ist,z1_ist); // set original pos of hits on ist1
	// aligned to IST2
	double x1_corr_IST = x1_ist*TMath::Cos(FST::phi_rot_ist1) + y1_ist*TMath::Sin(FST::phi_rot_ist1) + FST::x1_shift;
	double y1_corr_IST = y1_ist*TMath::Cos(FST::phi_rot_ist1) - x1_ist*TMath::Sin(FST::phi_rot_ist1) + FST::y1_shift;
	vPosIst1.SetXYZ(x1_corr_IST,y1_corr_IST,z1_ist);
	// aligned to FST mDefSenorId
	// double x1_corr_FST = x1_corr_IST*TMath::Cos(FST::phi_rot_ist2) + y1_corr_IST*TMath::Sin(FST::phi_rot_ist2) + FST::x2_shift;
	// double y1_corr_FST = y1_corr_IST*TMath::Cos(FST::phi_rot_ist2) - x1_corr_IST*TMath::Sin(FST::phi_rot_ist2) + FST::y2_shift;
	// vAlignedIst1.SetXYZ(x1_corr_FST,y1_corr_FST,z1_ist); // set aligned pos of hits on ist1
	vAlignedIst1 = this->getIstAlignedPos(vPosIst1);

	double x3_ist = rawHitVec[3][i_ist3]->getPosX();
	double y3_ist = rawHitVec[3][i_ist3]->getPosY();
	vPosIst3_orig.SetXYZ(x3_ist,y3_ist,z3_ist); // set original pos of hits on ist3
	// aligned to IST2
	double x3_corr_IST = x3_ist*TMath::Cos(FST::phi_rot_ist3) + y3_ist*TMath::Sin(FST::phi_rot_ist3) + FST::x3_shift;
	double y3_corr_IST = y3_ist*TMath::Cos(FST::phi_rot_ist3) - x3_ist*TMath::Sin(FST::phi_rot_ist3) + FST::y3_shift;
	vPosIst3.SetXYZ(x3_corr_IST,y3_corr_IST,z3_ist);
	// aligned to FST
	// double x3_corr_FST = x3_corr_IST*TMath::Cos(FST::phi_rot_ist2) + y3_corr_IST*TMath::Sin(FST::phi_rot_ist2) + FST::x2_shift;
	// double y3_corr_FST = y3_corr_IST*TMath::Cos(FST::phi_rot_ist2) - x3_corr_IST*TMath::Sin(FST::phi_rot_ist2) + FST::y2_shift;
	// vAlignedIst3.SetXYZ(x3_corr_FST,y3_corr_FST,z3_ist); // set aligned pos of hits on ist3
	vAlignedIst3 = this->getIstAlignedPos(vPosIst3);

	// double x2_proj = x3_corr_FST + (x1_corr_FST-x3_corr_FST)*z2_ist/z1_ist;
	// double y2_proj = y3_corr_FST + (y1_corr_FST-y3_corr_FST)*z2_ist/z1_ist;
	double x2_proj = x3_corr_IST + (x1_corr_IST-x3_corr_IST)*z2_ist/z1_ist;
	double y2_proj = y3_corr_IST + (y1_corr_IST-y3_corr_IST)*z2_ist/z1_ist;
	vProjIst2.SetXYZ(x2_proj,y2_proj,z2_ist); // x & y before aligned to FST

	FstTrack *fstTrack = new FstTrack(); // generate a FstTrack when both ist1 and ist3 have a hit
	fstTrack->Clear();
	fstTrack->setTrackId(numOfTracks);
	fstTrack->setTrackType(0);

	// set IST1 & IST3 position projection on IST2
	fstTrack->setIdIst1(rawHitVec[1][i_ist1]->getHitId());
	fstTrack->setOrigPosIst1(vPosIst1_orig);
	fstTrack->setAlignedIst1(vAlignedIst1); // aligned w.r.t. mDefSenorId
	fstTrack->setIdIst3(rawHitVec[3][i_ist3]->getHitId());
	fstTrack->setOrigPosIst3(vPosIst3_orig);
	fstTrack->setAlignedIst3(vAlignedIst3);
	fstTrack->setProjIst2(vProjIst2);

	for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
	{
	  // double x0_proj = x3_corr_FST + (x1_corr_FST-x3_corr_FST)*z0_fst/z1_ist;
	  // double y0_proj = y3_corr_FST + (y1_corr_FST-y3_corr_FST)*z0_fst/z1_ist;
	  // double r0_proj = TMath::Sqrt(x0_proj*x0_proj+y0_proj*y0_proj);
	  // double phi0_proj = TMath::ATan2(y0_proj,x0_proj);
	  // vProjFst.SetXYZ(r0_proj,phi0_proj,z0_fst); // r & phi for FST
	  vProjFst = this->getIstProjection(vPosIst1, vPosIst3, i_sensor); // get aligned with mSensorId and projected to FST i_sensor

	  // set projection on FST
	  fstTrack->setProjFst(i_sensor,vProjFst);
	}
	trackVec.push_back(fstTrack);
	// cout << "numOfTracks = " << numOfTracks << endl;
	// trackVec[numOfTracks]->Print();
	numOfTracks++;
      }
    }
  }

  return trackVec;
}
//--------------Find Track with Hits---------------------

//--------------Find Track with Clusters---------------------
std::vector<FstTrack *> FstClusterMaker::findTrack_Clusters(std::vector<FstCluster *> clusterVec_orig)
{
  int numOfClusters[4]; // 0 for fst | 1-3 for ist
  std::vector<FstCluster *> clusterVec[4];
  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    numOfClusters[i_layer] = 0;
    clusterVec[i_layer].clear();
  }

  for(int i_cluster = 0; i_cluster < clusterVec_orig.size(); ++i_cluster)
  { // set temp ist cluster container for each layer
    int layer = clusterVec_orig[i_cluster]->getLayer();
    numOfClusters[layer]++;
    clusterVec[layer].push_back(clusterVec_orig[i_cluster]);
  }

  const double z0_fst = FST::pitchLayer03;
  const double z1_ist = FST::pitchLayer12 + FST::pitchLayer23;
  const double z2_ist = FST::pitchLayer23;
  const double z3_ist = 0.0;

  std::vector<FstTrack *> trackVec;
  trackVec.clear();
  int numOfTracks = 100; // 100 => reconstructed from simple cluster 
  if(numOfClusters[1] > 0 && numOfClusters[3] > 0)
  { // only do tracking when at least 1 hit is found in ist1 & ist3
    TVector3 vPosIst1_orig, vPosIst3_orig;
    TVector3 vPosIst1, vPosIst3;
    TVector3 vAlignedIst1, vAlignedIst3;
    TVector3 vProjFst, vProjIst2;

    for(int i_ist1 = 0; i_ist1 < numOfClusters[1]; ++i_ist1)
    {
      for(int i_ist3 = 0; i_ist3 < numOfClusters[3]; ++i_ist3)
      {
	double x1_ist = clusterVec[1][i_ist1]->getMeanX();
	double y1_ist = clusterVec[1][i_ist1]->getMeanY();
	vPosIst1_orig.SetXYZ(x1_ist,y1_ist,z1_ist); // set original pos of hits on ist1
	// aligned to IST2
	double x1_corr_IST = x1_ist*TMath::Cos(FST::phi_rot_ist1) + y1_ist*TMath::Sin(FST::phi_rot_ist1) + FST::x1_shift;
	double y1_corr_IST = y1_ist*TMath::Cos(FST::phi_rot_ist1) - x1_ist*TMath::Sin(FST::phi_rot_ist1) + FST::y1_shift;
	vPosIst1.SetXYZ(x1_corr_IST,y1_corr_IST,z1_ist);
	// aligned to FST 
	// double x1_corr_FST = x1_corr_IST*TMath::Cos(FST::phi_rot_ist2) + y1_corr_IST*TMath::Sin(FST::phi_rot_ist2) + FST::x2_shift;
	// double y1_corr_FST = y1_corr_IST*TMath::Cos(FST::phi_rot_ist2) - x1_corr_IST*TMath::Sin(FST::phi_rot_ist2) + FST::y2_shift;
	// vAlignedIst1.SetXYZ(x1_corr_FST,y1_corr_FST,z1_ist); // set aligned pos of hits on ist1
	vAlignedIst1 = this->getIstAlignedPos(vPosIst1);

	double x3_ist = clusterVec[3][i_ist3]->getMeanX();
	double y3_ist = clusterVec[3][i_ist3]->getMeanY();
	vPosIst3_orig.SetXYZ(x3_ist,y3_ist,z3_ist); // set original pos of hits on ist3
	// aligned to IST2
	double x3_corr_IST = x3_ist*TMath::Cos(FST::phi_rot_ist3) + y3_ist*TMath::Sin(FST::phi_rot_ist3) + FST::x3_shift;
	double y3_corr_IST = y3_ist*TMath::Cos(FST::phi_rot_ist3) - x3_ist*TMath::Sin(FST::phi_rot_ist3) + FST::y3_shift;
	vPosIst3.SetXYZ(x3_corr_IST,y3_corr_IST,z3_ist); // set aligned pos of hits on ist3
	// aligned to FST
	// double x3_corr_FST = x3_corr_IST*TMath::Cos(FST::phi_rot_ist2) + y3_corr_IST*TMath::Sin(FST::phi_rot_ist2) + FST::x2_shift;
	// double y3_corr_FST = y3_corr_IST*TMath::Cos(FST::phi_rot_ist2) - x3_corr_IST*TMath::Sin(FST::phi_rot_ist2) + FST::y2_shift;
	// vAlignedIst3.SetXYZ(x3_corr_FST,y3_corr_FST,z3_ist); // set aligned pos of hits on ist3
	vAlignedIst3 = this->getIstAlignedPos(vPosIst3);

	// double x2_proj = x3_corr_FST + (x1_corr_FST-x3_corr_FST)*z2_ist/z1_ist;
	// double y2_proj = y3_corr_FST + (y1_corr_FST-y3_corr_FST)*z2_ist/z1_ist;
	double x2_proj = x3_corr_IST + (x1_corr_IST-x3_corr_IST)*z2_ist/z1_ist;
	double y2_proj = y3_corr_IST + (y1_corr_IST-y3_corr_IST)*z2_ist/z1_ist;
	vProjIst2.SetXYZ(x2_proj,y2_proj,z2_ist); // x & y before aligned to FST

	FstTrack *fstTrack = new FstTrack(); // generate a FstTrack when both ist1 and ist3 have a hit
	fstTrack->Clear();
	fstTrack->setTrackId(numOfTracks);
	fstTrack->setTrackType(1);

	// set IST1 & IST3 position projection on IST2
	fstTrack->setIdIst1(clusterVec[1][i_ist1]->getClusterId());
	fstTrack->setOrigPosIst1(vPosIst1_orig);
	fstTrack->setAlignedIst1(vAlignedIst1); // aligned w.r.t. mDefSenorId
	fstTrack->setIdIst3(clusterVec[3][i_ist3]->getClusterId());
	fstTrack->setOrigPosIst3(vPosIst3_orig);
	fstTrack->setAlignedIst3(vAlignedIst3); // aligned w.r.t. mDefSenorId
	fstTrack->setProjIst2(vProjIst2);

	for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
	{
	  // double x0_proj = x3_corr_FST + (x1_corr_FST-x3_corr_FST)*z0_fst/z1_ist;
	  // double y0_proj = y3_corr_FST + (y1_corr_FST-y3_corr_FST)*z0_fst/z1_ist;
	  // double r0_proj = TMath::Sqrt(x0_proj*x0_proj+y0_proj*y0_proj);
	  // double phi0_proj = TMath::ATan2(y0_proj,x0_proj);
	  // vProjFst.SetXYZ(r0_proj,phi0_proj,z0_fst); // r & phi for FST
	  vProjFst = this->getIstProjection(vPosIst1, vPosIst3, i_sensor); // get aligned with mSensorId and projected to FST i_sensor

	  // set projection on FST
	  fstTrack->setProjFst(i_sensor,vProjFst);
	}
	trackVec.push_back(fstTrack);
	// cout << "numOfTracks = " << numOfTracks << endl;
	// trackVec[numOfTracks]->Print();
	numOfTracks++;
      }
    }
  }

  return trackVec;
}
//--------------Find Track with Clusters---------------------

//--------------Output TTree---------------------
bool FstClusterMaker::initTree()
{
  mTree_FstEvent = new TTree("mTree_FstEvent","Fst Hits and Clusters Info");
  mFstEvent = new FstEvent();
  mTree_FstEvent->Branch("FstEvent","FstEvent",&mFstEvent);
  mTree_FstEvent->SetAutoSave(50000000);

  return true;
}

void FstClusterMaker::writeTree()
{
  mTree_FstEvent->Write();
}
//--------------Output TTree---------------------

//--------------Utility---------------------
int FstClusterMaker::getLayer(int arm, int port)
{
  if(arm == 0 && port == 0) return 1; // top IST
  if(arm == 0 && port == 1) return 2; // middle IST
  if(arm == 1 && port == 0) return 3; // bottom IST
  if(arm == 1 && port == 1) return 0; // FST module

  return -1;
}

int FstClusterMaker::getSensor(int arm, int port, int apv)
{
  int layer = this->getLayer(arm,port);
  if(layer == 3)
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

int FstClusterMaker::getColumn(int arm, int port, int apv, int ch)
{
  int col = -1;
  int layer = this->getLayer(arm,port);
  int sensor = this->getSensor(arm,port,apv);

  if(layer == 3)
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

int FstClusterMaker::getRow(int arm, int port, int apv, int ch)
{
  int row = -1;
  int layer = this->getLayer(arm,port);
  if(layer == 3)
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

int FstClusterMaker::getRStrip(int apv, int ch)
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

int FstClusterMaker::getPhiSeg(int apv, int ch)
{
  int phi_seg = -1;

  if(apv >= 0 && apv <= 3) phi_seg = apv*32 + ch/4; // inner
  if(apv >= 4 && apv <= 7) phi_seg = (apv-4)*32 + ch/4; // outer

  return phi_seg;
}

double FstClusterMaker::getPosX(int arm, int port, int apv, int ch)
{
  double posX = -999.9;
  int layer = this->getLayer(arm,port);
  if(layer == 3) // FST
  {
    int sensor = this->getSensor(arm,port,apv);
    if(sensor == 0) posX = FST::rInner + (this->getRStrip(apv,ch))*FST::pitchR + 0.5*FST::pitchR;
    if(sensor > 0)  posX = FST::rOuter + (this->getRStrip(apv,ch))*FST::pitchR + 0.5*FST::pitchR;
  }
  else // IST1-3
  {
    posX = this->getColumn(arm,port,apv,ch)*FST::pitchColumn + 0.5*FST::pitchColumn;
  }

  return posX;
}

double FstClusterMaker::getPosY(int arm, int port, int apv, int ch)
{
  double posY = -999.9;
  int layer = this->getLayer(arm,port);
  if(layer == 3) // FST
  {
    int sensor = this->getSensor(arm,port,apv);
    // inner sector face up
    if(sensor == 0) posY = (this->getPhiSeg(apv,ch)-63)*FST::pitchPhi - 0.5*FST::pitchPhi;
    // outer sector face down
    if(sensor == 1) posY = (63-this->getPhiSeg(apv,ch))*FST::pitchPhi + 0.5*FST::pitchPhi + 0.5*FST::gapPhi; // 1st half
    if(sensor == 2) posY = (63-this->getPhiSeg(apv,ch))*FST::pitchPhi + 0.5*FST::pitchPhi - 0.5*FST::gapPhi; // 2nd half
  }
  else // IST1-3
  {
    posY = (63-this->getRow(arm,port,apv,ch))*FST::pitchRow + 0.5*FST::pitchRow;
  }

  return posY;
}

bool FstClusterMaker::isBadAPV(int arm, int port, int apv)
{
  bool bAPV = false;

  if(arm == 0 && port==0 && (apv == 1 || apv == 5||apv == 6)) bAPV = true;
  if(arm == 0 && port==1 && (apv <= 9)) bAPV = true;
  // if(arm == 1 && port==0 && (apv == 0||apv == 1)) bAPV = true;
  if(arm == 1 && port==0 && (apv == 0)) bAPV = true;
  if(arm == 1 && port==1 && (apv > 7)) bAPV = true;

  return bAPV;
}

TVector3 FstClusterMaker::getIstProjection(TVector3 vPosIst1, TVector3 vPosIst3, int sensorId)
{
  // double x1_IST = vPosIst1.X(); // aligned to IST2
  // double y1_IST = vPosIst1.Y();
  // double z1_IST = vPosIst1.Z();
  // double x1_FST = x1_IST*TMath::Cos(FST::phi_rot[FST::mDefSenorId]) + y1_IST*TMath::Sin(FST::phi_rot[FST::mDefSenorId]) + FST::x_shift[FST::mDefSenorId]; // aligned to FST sensorId mDefSenorId
  // double y1_FST = y1_IST*TMath::Cos(FST::phi_rot[FST::mDefSenorId]) - x1_IST*TMath::Sin(FST::phi_rot[FST::mDefSenorId]) + FST::y_shift[FST::mDefSenorId];
  TVector3 vAlignedIst1 = this->getIstAlignedPos(vPosIst1); // aligned to FST sensorId mDefSenorId
  double x1_FST = vAlignedIst1.X();
  double y1_FST = vAlignedIst1.Y();

  // double x3_IST = vPosIst3.X(); // aligned to IST2
  // double y3_IST = vPosIst3.Y();
  // double z3_IST = vPosIst3.Z();
  // double x3_FST = x3_IST*TMath::Cos(FST::phi_rot[FST::mDefSenorId]) + y3_IST*TMath::Sin(FST::phi_rot[FST::mDefSenorId]) + FST::x_shift[FST::mDefSenorId]; // aligned to FST sensorId mDefSenorId
  // double y3_FST = y3_IST*TMath::Cos(FST::phi_rot[FST::mDefSenorId]) - x3_IST*TMath::Sin(FST::phi_rot[FST::mDefSenorId]) + FST::y_shift[FST::mDefSenorId];
  TVector3 vAlignedIst3 = this->getIstAlignedPos(vPosIst3); // aligned to FST sensorId mDefSenorId
  double x3_FST = vAlignedIst3.X();
  double y3_FST = vAlignedIst3.Y();

  double x0_FST = x3_FST + (x1_FST-x3_FST)*(FST::pitchLayer03+FST::z_shift[sensorId])/(FST::pitchLayer12+FST::pitchLayer23);
  double y0_FST = y3_FST + (y1_FST-y3_FST)*(FST::pitchLayer03+FST::z_shift[sensorId])/(FST::pitchLayer12+FST::pitchLayer23);
  double z0_FST = FST::pitchLayer03+FST::z_shift[sensorId];
  double r0_FST = TMath::Sqrt(x0_FST*x0_FST+y0_FST*y0_FST);
  double phi0_FST = TMath::ATan2(y0_FST,x0_FST);

  TVector3 vProjFst;
  vProjFst.SetXYZ(r0_FST,phi0_FST,z0_FST);

  return vProjFst;
}

TVector3 FstClusterMaker::getIstAlignedPos(TVector3 vPosIst)
{
  double x_IST = vPosIst.X(); // aligned to IST2
  double y_IST = vPosIst.Y();
  double z_IST = vPosIst.Z();
  // aligned to FST sensorId mDefSenorId
  double x_AlignedIST = x_IST*TMath::Cos(FST::phi_rot[FST::mDefSenorId]) + y_IST*TMath::Sin(FST::phi_rot[FST::mDefSenorId]) + FST::x_shift[FST::mDefSenorId];
  double y_AlignedIST = y_IST*TMath::Cos(FST::phi_rot[FST::mDefSenorId]) - x_IST*TMath::Sin(FST::phi_rot[FST::mDefSenorId]) + FST::y_shift[FST::mDefSenorId];

  TVector3 vAlignedIst;
  vAlignedIst.SetXYZ(x_AlignedIST,y_AlignedIST,z_IST);

  return vAlignedIst;
}
//--------------Utility---------------------
