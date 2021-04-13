#ifndef IstTracking_h
#define IstTracking_h

#include <TObject.h>
#include <./IstTrackingCons.h>
#include <vector>

class TFile;
class TChain;
class TH1F;
class TH2F;
class TH3F;
class TGraph;

typedef struct
{
  int layer;
  int sensor;
  int column;
  int row;
  double maxAdc;
  int maxTb;
  bool filled;
} IstHit;

typedef struct
{
  int layer;
  int sensor;
  double meanColumn;
  double meanRow;
  double totAdc;
  int maxTb;
  int clusterSize;
  int clusterSizeX;
  int clusterSizeY;
  int clusterType; // 0 for ARMDisplay | 1 for Simple | 2 for Scan
} IstCluster;

class IstTracking : public TObject
{
  public:
    IstTracking();
    ~IstTracking();

    void set_list(std::string inputlist)
    {
      mList = inputlist;
    }
    void set_outputfile(std::string outputfile)
    {
      mOutPutFile = outputfile;
    }

    bool clearSignal();
    bool clearHit();

    int Init();
    bool initChain();
    bool initSignal();
    bool initHit();
    bool initCluster();

    int Make(); // subtract pedestal => signal & find Hits

    int Finish();

    // pedestal
    bool clearPedestal();
    bool initPedestal();
    bool calPedestal(); // extract pedestal for each ch and fill TGraphs for ped mean & sigma (noise)
    void writePedestal();

    // hit display
    bool initHitDisplay();
    void fillHitDisplay(std::vector<IstHit> isthitvec);
    void writeHitDisplay();

    // cluster with ARMDisplay
    bool clearCluster_ARMDisplay();
    bool initTracking_ARMDisplay();
    bool findCluster_ARMDisplay(std::vector<IstHit> isthit_orig);
    void fillTracking_ARMDisplay(std::vector<IstCluster> istcluster_orig);
    void writeTracking_ARMDisplay();

    // cluster with Simple Algorithm
    bool clearCluster_Simple();
    bool initTracking_Simple();
    bool findCluster_Simple(std::vector<IstHit> isthit_orig);
    // void fillTracking_Simple(std::vector<IstCluster> istcluster_orig);
    void writeTracking_Simple();

  private:
    std::string mHome, mList;
    std::string mOutPutFile;

    TFile *File_mOutPut;

    // Pedestal
    double mPed[IST::numARMs][IST::numPorts][IST::numAPVs][IST::numChannels];
    double mPedStdDev[IST::numARMs][IST::numPorts][IST::numAPVs][IST::numChannels];
    // double pedRMS[IST::numARMs][IST::numPorts][IST::numAPVs][IST::numChannels];

    TGraph *g_mPedMean_Layer1;
    TGraph *g_mPedSigma_Layer1;
    TGraph *g_mPedMean_Layer2;
    TGraph *g_mPedSigma_Layer2;
    TGraph *g_mPedMean_Layer3;
    TGraph *g_mPedSigma_Layer3;

    TH2F *h_mHitDisplay_Layer1; // col & row
    TH2F *h_mHitDisplay_Layer2;
    TH2F *h_mHitDisplay_Layer3;
    TH1F *h_mMaxTb_Layer1; 
    TH1F *h_mMaxTb_Layer2;
    TH1F *h_mMaxTb_Layer3;

    std::vector<IstHit> mIstHitVec; // store hit information after ped subtraction
    double mSigPedCorr[IST::numARMs][IST::numPorts][IST::numAPVs][IST::numChannels][IST::numTBins];
    double mRawSig[IST::numARMs][IST::numPorts][IST::numAPVs][IST::numChannels][IST::numTBins];

    std::vector<IstCluster> mIstClusterVec_ARMDisplay; // cluster with ARMDisplay
    TH1F *h_mXResidual_ARMDisplay;
    TH1F *h_mYResidual_ARMDisplay;
    TH1F *h_mAdc_Layer1_ARMDisplay;
    TH1F *h_mAdc_Layer3_ARMDisplay;
    TH1F *h_mAdcAngleCorr_Layer1_ARMDisplay;
    TH1F *h_mAdcAngleCorr_Layer3_ARMDisplay;
    TH1F *h_mTrackAngle_ARMDisplay;

    std::vector<IstCluster> mIstClusterVec_Simple; // cluster with Simple Algorithm
    TH1F *h_mXResidual_Simple;
    TH1F *h_mYResidual_Simple;
    TH1F *h_mAdc_Layer1_Simple;
    TH1F *h_mAdc_Layer3_Simple;
    TH1F *h_mAdcAngleCorr_Layer1_Simple;
    TH1F *h_mAdcAngleCorr_Layer3_Simple;
    TH1F *h_mTrackAngle_Simple;

    // Utility for tracking
    int getLayer(int arm, int port); // return layer based on arm & port
    int getSensor(int apv); // return sensor based on APV
    int getColumn(int apv, int ch); // return column based on arm & port & apv & ch
    int getRow(int apv, int ch); // return row based on arm & port & apv & ch

    TChain *mChainInPut; // input TTree
    int evt_rdo[IST::numARMs][IST::numPorts][IST::numAPVs];
    int evt_arm[IST::numARMs][IST::numPorts][IST::numAPVs];
    int evt_port[IST::numARMs][IST::numPorts][IST::numAPVs];
    int evt_apv[IST::numARMs][IST::numPorts][IST::numAPVs];
    int hit_ch[IST::numARMs][IST::numPorts][IST::numAPVs][IST::numROChannels];
    int hit_tb[IST::numARMs][IST::numPorts][IST::numAPVs][IST::numROChannels];
    int hit_adc[IST::numARMs][IST::numPorts][IST::numAPVs][IST::numROChannels];

    ClassDef(IstTracking,1)
};

#endif
