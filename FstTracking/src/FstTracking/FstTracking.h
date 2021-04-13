#ifndef FstTracking_h
#define FstTracking_h

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

class FstTracking : public TObject
{
  public:
    FstTracking();
    ~FstTracking();

    void set_list(std::string inputlist)
    {
      mList = inputlist;
    }
    void set_outputfile(std::string outputfile)
    {
      mOutPutFile = outputfile;
    }

    int Init();
    int Make(); // subtract pedestal => signal & find Hits
    int Finish();

    // init input TChain
    bool initChain();

    // track QA
    bool initTrackingQA_Hits();
    void fillTrackingQA_Hits(std::vector<FstRawHit *> rawHitVec_orig);
    void writeTrackingQA_Hits();

    // incident angle
    void initTrackAngle();
    void fillTrackAngle(FstEvent *fstEvent);
    void writeTrackAngle();

    // track resolution
    void initTracking_Hits();
    void calResolution_Hits(FstEvent *fstEvent);
    void writeTracking_Hits();

    void initTracking_Clusters();
    void calResolution_SimpleClusters(FstEvent *fstEvent);
    void calResolution_ScanClusters(FstEvent *fstEvent);
    void writeTracking_Clusters();

    // efficiency
    void initEfficiency_Hits();
    void calEfficiency_Hits(FstEvent *fstEvent);
    void calEfficiencyModule_Hits(FstEvent *fstEvent);
    void writeEfficiency_Hits();

    void initEfficiency_Clusters();
    void calEfficiency_SimpleClusters(FstEvent *fstEvent);
    void calEfficiency_ScanClusters(FstEvent *fstEvent);
    void calEfficiencyModule_SimpleClusters(FstEvent *fstEvent);
    void calEfficiencyModule_ScanClusters(FstEvent *fstEvent);
    void writeEfficiency_Clusters();

  private:
    std::string mHome, mList;
    std::string mOutPutFile;

    TFile *File_mOutPut;

    // Tracking QA based on Hits
    TH2F *h_mHitsCorrXR[4]; // 0: ist1x vs. ist3x | 1: ist1x vs. fstr | 2: ist3x vs. fstr | 3: ist1x+ist3x vs. fstr
    TH2F *h_mHitsCorrYPhi[4]; // 0: ist1y vs. ist3y | 1: ist1y vs. fstphi | 2: ist3y vs. fstphi | 3: ist1y+ist3y vs. fstphi

    // incident angle
    TH1F *h_mClustersTrackAngle; // cosmic incident angle on IST

    // Tracking
    // hits
    TH1F *h_mHitsTrackFstResX[FST::mFstNumSensorsPerModule];
    TH1F *h_mHitsTrackFstResY[FST::mFstNumSensorsPerModule];
    TH2F *h_mHitsTrackFstResXY[FST::mFstNumSensorsPerModule]; // 2D distribution

    TH1F *h_mHitsTrackFstResR[FST::mFstNumSensorsPerModule];
    TH1F *h_mHitsTrackFstResPhi[FST::mFstNumSensorsPerModule];
    TH2F *h_mHitsTrackFstResRPhi[FST::mFstNumSensorsPerModule]; // 2D distribution

    // IST simple clusters
    TH1F *h_mSimpleClustersTrackIstResX_2Layer;
    TH1F *h_mSimpleClustersTrackIstResY_2Layer;
    TH2F *h_mSimpleClustersTrackIstResXY_2Layer; // 2D distribution

    // FST simple clusters
    TH1F *h_mSimpleClustersTrackFstResX_2Layer[FST::mFstNumSensorsPerModule]; // 0 for inner | 1 & 2 for outer
    TH1F *h_mSimpleClustersTrackFstResY_2Layer[FST::mFstNumSensorsPerModule];
    TH2F *h_mSimpleClustersTrackFstResXY_2Layer[FST::mFstNumSensorsPerModule]; // 2D distribution

    TH1F *h_mSimpleClustersTrackFstResR_2Layer[FST::mFstNumSensorsPerModule];
    TH1F *h_mSimpleClustersTrackFstResPhi_2Layer[FST::mFstNumSensorsPerModule];
    TH2F *h_mSimpleClustersTrackFstResRPhi_2Layer[FST::mFstNumSensorsPerModule]; // 2D distribution
    
    TH1F *h_mSimpleClustersTrackFstResX_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH1F *h_mSimpleClustersTrackFstResY_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH2F *h_mSimpleClustersTrackFstResXY_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor]; // 2D distribution
    
    TH1F *h_mSimpleClustersTrackFstResR_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH1F *h_mSimpleClustersTrackFstResPhi_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
    TH2F *h_mSimpleClustersTrackFstResRPhi_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor]; // 2D distribution

    // FST scan clusters
    TH1F *h_mScanClustersTrackFstResX_2Layer[FST::mFstNumSensorsPerModule];
    TH1F *h_mScanClustersTrackFstResY_2Layer[FST::mFstNumSensorsPerModule];
    TH2F *h_mScanClustersTrackFstResXY_2Layer[FST::mFstNumSensorsPerModule]; // 2D distribution

    TH1F *h_mScanClustersTrackFstResR_2Layer[FST::mFstNumSensorsPerModule];
    TH1F *h_mScanClustersTrackFstResPhi_2Layer[FST::mFstNumSensorsPerModule];
    TH2F *h_mScanClustersTrackFstResRPhi_2Layer[FST::mFstNumSensorsPerModule]; // 2D distribution
    
    TH1F *h_mScanClustersTrackFstResX_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::numRStrip];
    TH1F *h_mScanClustersTrackFstResY_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::numRStrip];
    TH2F *h_mScanClustersTrackFstResXY_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::numRStrip]; // 2D distribution
    
    TH1F *h_mScanClustersTrackFstResR_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::numRStrip];
    TH1F *h_mScanClustersTrackFstResPhi_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::numRStrip];
    TH2F *h_mScanClustersTrackFstResRPhi_2Layer_Rstrips[FST::mFstNumSensorsPerModule][FST::numRStrip]; // 2D distribution

    // Efficiency based on Hits
    // 0: no matching | 1-7 matching within (1-7)*pitchR in r & 3.0(6.0)*pitchPhi in phi
    TH2F *h_mHitsTrackIstCounts[FST::mFstNumSensorsPerModule][FST::mFstNumMatching]; // position on FST from IST projection
    TH2F *h_mHitsTrackFstCounts[FST::mFstNumSensorsPerModule][FST::mFstNumMatching]; // position on FST from FST measurement
    // whole module
    TH2F *h_mHitsTrackIstCountsModule[FST::mFstNumMatching]; // position on FST from IST projection
    TH2F *h_mHitsTrackFstCountsModule[FST::mFstNumMatching]; // position on FST from FST measurement

    // Efficiency based on Clusters
    // simple clusters
    TH2F *h_mSimpleClustersTrackIstCounts_2Layer[FST::mFstNumSensorsPerModule][FST::mFstNumMatching]; // position on FST from IST projection
    TH2F *h_mSimpleClustersTrackFstCounts_2Layer[FST::mFstNumSensorsPerModule][FST::mFstNumMatching]; // position on FST from FST measurement
    // whole module
    TH2F *h_mSimpleClustersTrackIstCountsModule_2Layer[FST::mFstNumMatching]; // position on FST from IST projection
    TH2F *h_mSimpleClustersTrackFstCountsModule_2Layer[FST::mFstNumMatching]; // position on FST from FST measurement

    // scan clusters
    TH2F *h_mScanClustersTrackIstCounts_2Layer[FST::mFstNumSensorsPerModule][FST::mFstNumMatching]; // position on FST from IST projection
    TH2F *h_mScanClustersTrackFstCounts_2Layer[FST::mFstNumSensorsPerModule][FST::mFstNumMatching]; // position on FST from FST measurement
    // whole module
    TH2F *h_mScanClustersTrackIstCountsModule_2Layer[FST::mFstNumMatching]; // position on FST from IST projection
    TH2F *h_mScanClustersTrackFstCountsModule_2Layer[FST::mFstNumMatching]; // position on FST from FST measurement

    // Input TChain for hits and clusters
    TChain *mChainInPut; // input TTree
    FstEvent *mFstEvent_InPut;
    FstRawHit *mFstRawHit;
    FstCluster *mFstCluster;
    FstTrack *mFstTrack;

    // Utility for tracking
    TVector2 getFstAlignedPos(TVector2 vPosFst, int sensorId); // get alignmed sensorId w.r.t. IST2 then to mDefSenorId

    ClassDef(FstTracking,1)
};

#endif
