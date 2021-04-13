#ifndef FstTrack_h
#define FstTrack_h

#include <TObject.h>
#include <TVector3.h>
#include "./FstCons.h"

class FstTrack : public TObject
{
  private:
    int mTrackId; // 0 => reconstructed from hits | 100 => reconstructed from simple cluster
    int mTrackType; // 0: reconstructed from hits | 1: reconstructed from simple cluster

    int mIdIst1;
    TVector3 mOrigPosIst1; // measured original position in IST1 frame
    TVector3 mAlignedIst1; // IST1 position in mDefSenorId frame
    int mIdIst3;
    TVector3 mOrigPosIst3; // measured original position in IST3 frame
    TVector3 mAlignedIst3; // IST3 position in mDefSenorId frame
    TVector3 mProjIst2; // projected cosmic position on IST2 in IST2 frame
    TVector3 mProjFst[FST::mFstNumSensorsPerModule]; // projected cosmic position on FST for each sensorId in mDefSenorId frame

  public:
    FstTrack();
    ~FstTrack();

    void setTrackId(int trackId);
    void setTrackType(int trackType);
    void setIdIst1(int id);
    void setOrigPosIst1(TVector3 pos);
    void setAlignedIst1(TVector3 pos);
    void setIdIst3(int id);
    void setOrigPosIst3(TVector3 pos);
    void setAlignedIst3(TVector3 pos);
    void setProjIst2(TVector3 pos);
    void setProjFst(int sensorId, TVector3 pos);

    int getTrackId() const;
    int getTrackType() const;
    int getIdIst1() const;
    TVector3 getOrigPosIst1() const;
    TVector3 getAlignedIst1() const;
    int getIdIst3() const;
    TVector3 getOrigPosIst3() const;
    TVector3 getAlignedIst3() const;
    TVector3 getProjIst2() const;
    TVector3 getProjFst(int sensorId) const;

    void Print() const;
    void Clear();

    ClassDef(FstTrack,1)
};

#endif

