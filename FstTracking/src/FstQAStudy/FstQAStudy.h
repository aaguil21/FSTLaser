#ifndef FstQAStudy_h
#define FstQAStudy_h

#include <TObject.h>
#include "../FstUtil/FstCons.h"
#include "../FstUtil/FstRawHit.h"
#include "../FstUtil/FstCluster.h"
#include "../FstUtil/FstTrack.h"
#include "../FstUtil/FstEvent.h"
#include <vector>
#include <fstream>

class TFile;
class TChain;
class TTree;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TProfile2D;
class TGraph;

class FstQAStudy : public TObject
{
  public:
    FstQAStudy();
    ~FstQAStudy();

    void set_list(std::string inputlist)
    {
      mList = inputlist;
    }
    void set_outputfile(std::string outputfile)
    {
      mOutPutFile = outputfile;
    }
    void set_cmncorrection(bool applyCMNCorr)
    {
      mApplyCMNCorr = applyCMNCorr;
    }

    int Init();
    int Make(); // subtract pedestal => signal & find Hits
    int Finish();

    // init input TChain
    bool initChain();

    // signal & noise
    void initSignalQA();
    //void fillSignalQA(FstEvent *fstEvent, int event, std::vector<Int_t>& treeEvent, std::vector<Double_t>& sigVec);
    void writeSignalQA();
    void fillSignalQA(FstEvent *fstEvent,int event, std::vector<Int_t>& eventIdVec, std::vector<Int_t>& maxTbVec, std::vector<Int_t>& sensorIdVec, std::vector<Int_t>& colNumVec, std::vector<Int_t>& rowNumVec,std::vector<Double_t>& signalVec, std::vector<Double_t>& noiseVec, std::vector<Double_t>& cmnNoiseVec);

    // Event Display
    void initEventDisplay_TrackClusters();
    void clearEventDisplay_TrackClusters();
    void fillEventDisplay_TrackClusters(FstEvent *fstEvent);
    void writeEventDisplay_TrackClusters();

  private:
    std::string mHome, mList;
    std::string mOutPutFile;

    TFile *File_mOutPut;

    // Signal & noise for max time bin
    // hits
    TProfile2D *p_mFstPedMap[FST::mFstNumSensorsPerModule];
    TProfile2D *p_mFstSigMap[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstHitsMaxTb[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstHitsSignal[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstHitsNoise[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstHitsSNRatio[FST::mFstNumSensorsPerModule];

    TH1F *h_mFstHitsMaxTb_Rstrip[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH1F *h_mFstHitsSignal_Rstrip[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH1F *h_mFstHitsNoise_Rstrip[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH1F *h_mFstHitsSNRatio_Rstrip[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];

    TH1F *h_mFstHitsSignal_Rstrip_TimeBin[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor][FST::numTBins]; // time bin differential
    TH1F *h_mFstHitsNoise_Rstrip_TimeBin[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor][FST::numTBins];
    TH1F *h_mFstHitsSNRatio_Rstrip_TimeBin[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor][FST::numTBins];

    TH3F *h_mFstLaser_HitsSignal[FST::mFstNumSensorsPerModule][34];
    TH3F *h_mFstLaser_NoiseSignal[FST::mFstNumSensorsPerModule][34];    
    TH3F *h_mFstLaser_SignalRatio[FST::mFstNumSensorsPerModule][34];
    TH3F *h_mFstLaser_Timebin[FST::mFstNumSensorsPerModule][34];
    TH3F *h_mFstLaser_cmnNoiseSignal[FST::mFstNumSensorsPerModule][34];
    TH2F *h_mFstLaser_SigTime[FST::mFstNumSensorsPerModule][34];
    TH2F *h_mFstLaser_NoiseTime[FST::mFstNumSensorsPerModule][34];    
    
    //Laser Scan Event Tree and Branches
    TTree *Event_Signal;
/*    std::vector<Double_t> signalEv;
    std::vector<Double_t> noiseEv;
    std::vector<Int_t> maxTbEv;
    std::vector<Double_t> cmnNoiseEv;
    std::vector<Int_t> treeEventId;
    std::vector<Int_t> treeSensorId;
    std::vector<Int_t> treeRowNum;
    std::vector<Int_t> treeColNum;
*/

    //clusters
    TH1F *h_mFstSimpleClustersSignal[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstSimpleClustersNoise[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstSimpleClustersSNRatio[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstSimpleClustersSignal_Rstrip[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH1F *h_mFstSimpleClustersMaxTb_Rstrip[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH1F *h_mFstSimpleClustersSignal_Rstrip_TimeBin[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor][FST::numTBins];

    TH1F *h_mFstScanClustersSignal[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstScanClustersNoise[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstScanClustersSNRatio[FST::mFstNumSensorsPerModule];
    TH1F *h_mFstScanClustersSignal_Rstrip[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH1F *h_mFstScanClustersMaxTb_Rstrip[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH1F *h_mFstScanClustersSignal_Rstrip_TimeBin[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor][FST::numTBins];

    // Output TTree for Event Display
    bool mApplyCMNCorr;
    TTree *mTree_EventDisplay;
    int mEventId;
    int mNumOfFstRawHits;
    int mNumOfIst1RawHits;
    int mNumOfIst2RawHits;
    int mNumOfIst3RawHits;
    int mNumOfFstSimpleClusters;
    int mNumOfFstScanClusters;
    int mNumOfIst1Clusters;
    int mNumOfIst2Clusters;
    int mNumOfIst3Clusters;
    int mNumOfHitTracks;
    int mNumOfClusterTracks; // all tracks with 2-Layer tracking
    int mNumOfClusterTracks_2Layer; // 2-Layer tracking within area of interest
    int mNumOfClusterTracks_3Layer; // 3-Layer tracking within area of interest

    // Event Display
    TH2F *h_mFstRawHitsDisplay; // save event display histogram into a TTree
    TH2F *h_mFstRawPedsDisplay;
    TH2F *h_mFstMaxTbDisplay;

    TH2F *h_mFstRawHitsDisplay_bTh; // save event display histogram into a TTree
    TH2F *h_mFstRawPedsDisplay_bTh;
    TH2F *h_mFstMaxTbDisplay_bTh;

    TH2F *h_mFstSimpleClustersDisplay;
    TH2F *h_mFstScanClustersDisplay;
    TH2F *h_mFstScanClustersDisplay_RawHits; // hits & peds in a cluster
    TH2F *h_mFstScanClustersDisplay_RawPeds;
    TH2F *h_mFstScanClustersDisplay_MaxTb;

    TH2F *h_mHitTracksDisplay;
    TH2F *h_mClusterTracksDisplay;

    TGraph *g_mFstSimpleClustersDisplay; // display only
    TGraph *g_mFstScanClustersDisplay; // display only
    TGraph *g_mHitTracksDisplay;
    TGraph *g_mClusterTracksDisplay;

    TH1F *h_mNumFstSimpleClusters;
    TH1F *h_mNumFstRawHitsSimpleCluster;
    TH1F *h_mNumFstScanClusters;
    TH1F *h_mNumFstRawHitsScanCluster;

    // Input TChain for hits and clusters
    TChain *mChainInPut; // input TTree
    FstEvent *mFstEvent_InPut;
    FstRawHit *mFstRawHit;
    FstCluster *mFstCluster;
    FstTrack *mFstTrack;

    // Utility for tracking
    TVector2 getFstAlignedPos(TVector2 vPosFst, int sensorId); // get alignmed sensorId w.r.t. IST2 then to mDefSenorId

    ClassDef(FstQAStudy,1)
};

#endif
