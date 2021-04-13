#include <iostream>
#include "./FstTrack.h"

using namespace std;

ClassImp(FstTrack)

//------------------------------------------

FstTrack::FstTrack()
{
  Clear();
}

FstTrack::~FstTrack()
{
}

//------------------------------------------
// modifiers
void FstTrack::setTrackId(int trackId)
{
  mTrackId = trackId;
}
void FstTrack::setTrackType(int trackType)
{
  mTrackType = trackType;
}
void FstTrack::setIdIst1(int id)
{
  mIdIst1 = id;
}
void FstTrack::setOrigPosIst1(TVector3 pos)
{
  mOrigPosIst1.SetXYZ(pos.X(),pos.Y(),pos.Z());
}
void FstTrack::setAlignedIst1(TVector3 pos)
{
  mAlignedIst1.SetXYZ(pos.X(),pos.Y(),pos.Z());
}
void FstTrack::setIdIst3(int id)
{
  mIdIst3 = id;
}
void FstTrack::setOrigPosIst3(TVector3 pos)
{
  mOrigPosIst3.SetXYZ(pos.X(),pos.Y(),pos.Z());
}
void FstTrack::setAlignedIst3(TVector3 pos)
{
  mAlignedIst3.SetXYZ(pos.X(),pos.Y(),pos.Z());
}
void FstTrack::setProjIst2(TVector3 pos)
{
  mProjIst2.SetXYZ(pos.X(),pos.Y(),pos.Z());
}
void FstTrack::setProjFst(int sensorId, TVector3 pos)
{
  mProjFst[sensorId].SetXYZ(pos.X(),pos.Y(),pos.Z());
}

// accessors
int FstTrack::getTrackId() const
{
  return mTrackId;
}
int FstTrack::getTrackType() const
{
  return mTrackType;
}
int FstTrack::getIdIst1() const
{
  return mIdIst1;
}
TVector3 FstTrack::getOrigPosIst1() const
{
  return mOrigPosIst1;
}
TVector3 FstTrack::getAlignedIst1() const
{
  return mAlignedIst1;
}
int FstTrack::getIdIst3() const
{
  return mIdIst3;
}
TVector3 FstTrack::getOrigPosIst3() const
{
  return mOrigPosIst3;
}
TVector3 FstTrack::getAlignedIst3() const
{
  return mAlignedIst3;
}
TVector3 FstTrack::getProjIst2() const
{
  return mProjIst2;
}
TVector3 FstTrack::getProjFst(int sensorId) const
{
  return mProjFst[sensorId];
}
//------------------------------------------
void FstTrack::Print() const
{
  cout << "mTrackId = " << mTrackId << endl;
  cout << "mTrackType = " << mTrackType << endl;

  cout << "Ist1: " << endl;
  cout << "mIdIst1 = " << mIdIst1 << endl;
  cout << "Original Ist1 X = " << mOrigPosIst1.X() << ", Y = " << mOrigPosIst1.Y() << ", Z = " << mOrigPosIst1.Z() << endl;;
  cout << "Aligned  Ist1 X = " << mAlignedIst1.X() << ", Y = " << mAlignedIst1.Y() << ", Z = " << mAlignedIst1.Z() << endl;;

  cout << "Ist3: " << endl;
  cout << "mIdIst3 = " << mIdIst1 << endl;
  cout << "Original Ist3 X = " << mOrigPosIst1.X() << ", Y = " << mOrigPosIst1.Y() << ", Z = " << mOrigPosIst1.Z() << endl;;
  cout << "Aligned  Ist3 X = " << mAlignedIst1.X() << ", Y = " << mAlignedIst1.Y() << ", Z = " << mAlignedIst1.Z() << endl;;

  cout << "Projection in IST2: " << endl;
  cout << "Projection X = " << mProjIst2.X() << ", Y = " << mProjIst2.Y() << ", Z = " << mProjIst2.Z() << endl;;

  cout << "Projection in FST: " << endl;
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    cout << "i_sensor = " << i_sensor << "Projection X = " << mProjFst[i_sensor].X() << ", Y = " << mProjFst[i_sensor].Y() << ", Z = " << mProjFst[i_sensor].Z() << endl;;
  }
}

void FstTrack::Clear()
{
  mTrackId = -1;
  mTrackType = -1;

  mIdIst1 = -1;
  mOrigPosIst1.SetXYZ(-999.9,-999.9,-999.9); 
  mAlignedIst1.SetXYZ(-999.9,-999.9,-999.9);
  mIdIst3 = -1;
  mOrigPosIst3.SetXYZ(-999.9,-999.9,-999.9);
  mAlignedIst3.SetXYZ(-999.9,-999.9,-999.9);
  mProjIst2.SetXYZ(-999.9,-999.9,-999.9);
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    mProjFst[i_sensor].SetXYZ(-999.9,-999.9,-999.9);
  }
}
//------------------------------------------
