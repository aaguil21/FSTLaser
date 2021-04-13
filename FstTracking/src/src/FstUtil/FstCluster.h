#ifndef FstCluster_h
#define FstCluster_h

#include <TObject.h>
#include <TClonesArray.h>
#include <vector>
#include "./FstRawHit.h"
#include "./FstCons.h"

class FstCluster : public TObject
{
  private:
    int mLayer; // 0 for FST | 1-3 for IST
    int mSensor;
    double mMeanColumn; // r_strip for FST | x for IST
    double mMeanRow;  // phi for FST | y for IST
    double mMeanX; // extracted from column & row
    double mMeanY; // r/phi for FST | x/y for IST
    double mTotCharge; // sum of pedestal subtracted ADC for max time bin of each hit in the cluster
    double mMaxTb;
    int mClusterType; // 0 for ARMDisplay | 1 for Simple | 2 for Scan
    bool mIsSeed; // true: cluster with seed | false cluster without seed
    int mNRawHits; // number of hits in the cluster
    int mNRawHitsR; // number of hits in the cluster along the R direction
    int mNRawHitsPhi; // number of hits in the cluster along the phi direction
    std::vector<FstRawHit *> mRawHitVec; //! container to save raw hits which contribute to the cluster => don't save in TTree
    int mClusterId; // cluster Id
    // int mHitId[FST::maxNHits];

    int mNumOfRawHits; // number of raw hits in the cluster
    TClonesArray *mFstRawHits; //->

  public:
    FstCluster();
    ~FstCluster();

    void setLayer(int layer);
    void setSensor(int sensor);
    void setMeanColumn(double column);
    void setMeanRow(double row);
    void setMeanX(double meanX);
    void setMeanY(double meanY);
    void setTotCharge(double charge);
    void setMaxTb(double maxTb);
    void setClusterType(int type);
    void setIsSeed(bool isSeed);
    void setNRawHits(int nRawHits);
    void setNRawHitsR(int nRawHitsR);
    void setNRawHitsPhi(int nRawHitsPhi);
    void addRawHit(FstRawHit *fstRawHit);
    void setClusterId(int clusterId);
    // void setHitId(int i_hit, int hitId);

    int getLayer() const;
    int getSensor() const;
    double getMeanColumn() const;
    double getMeanRow() const;
    double getMeanX() const;
    double getMeanY() const;
    double getTotCharge() const;
    double getMaxTb() const;
    int getClusterType() const;
    bool getIsSeed() const;
    int getNRawHits() const;
    int getNRawHitsR() const;
    int getNRawHitsPhi() const;
    std::vector<FstRawHit *> getRawHitVec() const;
    int getClusterId() const;
    // int getHitId(int i_hit) const;

    // FstRawHit
    FstRawHit* createRawHit();
    void clearRawHitsList();
    int getNumRawHits() const;
    FstRawHit* getRawHit(int i_hit) const;

    void Print() const;
    void Clear();

    ClassDef(FstCluster,1)
};

#endif
