#ifndef FstRawHit_h
#define FstRawHit_h

#include <TObject.h>
#include "./FstCons.h"

class FstRawHit : public TObject
{
  private:
    int mLayer; // 0 for FST | 1-3 for IST
    int mSensor;
    int mAPV;
    int mChannel;
    int mColumn; // r_strip for FST | x for IST
    int mRow; // phi for FST | y for IST
    double mPosX; // extracted from column & row
    double mPosY; // r/phi for FST | x/y for IST
    double mPedMean[FST::numTBins]; // mean pedestal for each time bin
    double mPedStdDev[FST::numTBins]; // rms of pedestal for each time bin => noise level
    double mPedRMS[FST::numTBins]; // rms of pedestal for each time bin => noise level
    double mCMNStdDev[FST::numTBins]; // rms of pedestal for each time bin => noise level
    double mRanStdDev[FST::numTBins]; // rms of pedestal for each time bin => noise level
    double mRawCharge[FST::numTBins]; // FST: adc | IST: adc
    double mCharge[FST::numTBins]; // FST: adc-ped-cmn | IST: adc-ped
    double mWeight; // weight to account difference in FST RStrip
    int mMaxTb;
    int mHitId; // hit id
    int mDefaultTb;
    bool mIsHit;

  public:
    FstRawHit();
    ~FstRawHit();

    void setLayer(int layer);
    void setSensor(int sensor);
    void setAPV(int apv);
    void setChannel(int channel);
    void setColumn(int column);
    void setRow(int row);
    void setPosX(double posX);
    void setPosY(double posY);
    void setPedMean(double pedMean, int tb);
    void setPedStdDev(double pedStdDev, int tb);
    void setPedRMS(double pedRMS, int tb);
    void setCMNStdDev(double cmnStdDev, int tb);
    void setRanStdDev(double ranStdDev, int tb);
    void setRawCharge(double charge, int tb);
    void setCharge(double charge, int tb);
    void setWeight(double weight);
    void setMaxTb(int maxTb);
    void setHitId(int hitId);
    void setDefaultTb(int tb);
    void setIsHit(bool isHit);

    int getLayer() const;
    int getSensor() const;
    int getAPV() const;
    int getChannel() const;
    int getColumn() const;
    int getRow() const;
    double getPosX() const;
    double getPosY() const;
    double getPedMean(int tb = 0) const;
    double getPedStdDev(int tb = 0) const;
    double getPedRMS(int tb = 0) const;
    double getCMNStdDev(int tb = 0) const;
    double getRanStdDev(int tb = 0) const;
    double getRawCharge(int tb = 0) const;
    double getCharge(int tb = 0) const;
    double getWeight() const;
    int getMaxTb() const;
    int getHitId() const;
    int getDefaultTb() const;
    bool getIsHit() const;

    void Print() const;
    void Clear();

    ClassDef(FstRawHit,1)
};

#endif
