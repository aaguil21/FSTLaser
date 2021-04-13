#include <iostream>
#include "./FstCluster.h"

using namespace std;

ClassImp(FstCluster)

//------------------------------------------

FstCluster::FstCluster()
{
  Clear();

  mNumOfRawHits = 0;
  mFstRawHits = new TClonesArray("FstRawHit", 10);
}

FstCluster::~FstCluster()
{
  delete mFstRawHits;
  mFstRawHits = NULL;
}

//------------------------------------------

// modifiers
void FstCluster::setLayer(int layer)
{
  mLayer = layer;
}
void FstCluster::setSensor(int sensor)
{
  mSensor = sensor;
}
void FstCluster::setMeanColumn(double column)
{
  mMeanColumn = column;
}
void FstCluster::setMeanRow(double row)
{
  mMeanRow = row;
}
void FstCluster::setMeanX(double meanX)
{
  mMeanX = meanX;
}
void FstCluster::setMeanY(double meanY)
{
  mMeanY = meanY;
}
void FstCluster::setTotCharge(double charge)
{
  mTotCharge = charge;
}
void FstCluster::setMaxTb(double maxTb)
{
  mMaxTb = maxTb;
}
void FstCluster::setClusterType(int type)
{
  mClusterType = type;
}
void FstCluster::setIsSeed(bool isSeed)
{
  mIsSeed = isSeed;
}
void FstCluster::setNRawHits(int nRawHits)
{
  mNRawHits = nRawHits;
}
void FstCluster::setNRawHitsR(int nRawHitsR)
{
  mNRawHitsR = nRawHitsR;
}
void FstCluster::setNRawHitsPhi(int nRawHitsPhi)
{
  mNRawHitsPhi = nRawHitsPhi;
}
void FstCluster::addRawHit(FstRawHit *fstRawHit)
{
  mRawHitVec.push_back(fstRawHit);
}
void FstCluster::setClusterId(int clusterId)
{
  mClusterId = clusterId;
}
// void FstCluster::setHitId(int i_hit, int hitId)
// {
//   mHitId[i_hit] = hitId;
// }

// accessors
int FstCluster::getLayer() const
{
  return mLayer;
}
int FstCluster::getSensor() const
{
  return mSensor;
}
double FstCluster::getMeanColumn() const
{
  return mMeanColumn;
}
double FstCluster::getMeanRow() const
{
  return mMeanRow;
}
double FstCluster::getMeanX() const
{
  return mMeanX;
}
double FstCluster::getMeanY() const
{
  return mMeanY;
}
double FstCluster::getTotCharge() const
{
  return mTotCharge;
}
double FstCluster::getMaxTb() const
{
  return mMaxTb;
}
int FstCluster::getClusterType() const
{
  return mClusterType;
}
bool FstCluster::getIsSeed() const
{
  return mIsSeed;
}
int FstCluster::getNRawHits() const
{
  return mNRawHits;
}
int FstCluster::getNRawHitsR() const
{
  return mNRawHitsR;
}
int FstCluster::getNRawHitsPhi() const
{
  return mNRawHitsPhi;
}
std::vector<FstRawHit *> FstCluster::getRawHitVec() const
{
  return mRawHitVec;
}
int FstCluster::getClusterId() const
{
  return mClusterId;
}
// int FstCluster::getHitId(int i_hit) const
// {
//   return mHitId[i_hit];
// }

//------------------------------------------

void FstCluster::Print() const
{
  cout << "mLayer = " << mLayer << endl;
  cout << "mSensor = " << mSensor << endl;
  cout << "mMeanColumn = " << mMeanColumn << endl;
  cout << "mMeanRow = " << mMeanRow << endl;
  cout << "mMeanX = " << mMeanX << endl;
  cout << "mMeanY = " << mMeanY << endl;
  cout << "mTotCharge = " << mTotCharge << endl;
  cout << "mMaxTb =" << mMaxTb << endl;
  cout << "mClusterType =" << mClusterType << endl;
  cout << "mIsSeed =" << mIsSeed << endl;
  cout << "mNRawHits = " << mNRawHits << endl;
  cout << "mNRawHitsR = " << mNRawHitsR << endl;
  cout << "mNRawHitsPhi = " << mNRawHitsPhi << endl;
  cout << "mClusterId = " << mClusterId << endl;
  cout << endl;
  for(int i_hit = 0; i_hit < mNRawHits; ++i_hit)
  {
    cout << "Print " << i_hit << "th hit:" << endl;
    // cout << "mHitId = " << mHitId[i_hit] << endl;
    FstRawHit *fstRawHit = this->getRawHit(i_hit);
    fstRawHit->Print();
  }
  cout << endl;
}

void FstCluster::Clear()
{
  mLayer = -1;
  mSensor = -1;
  mMeanColumn = -1;
  mMeanRow = -1;
  mMeanX = -1;
  mMeanY = -1;
  mTotCharge = -1;
  mMaxTb = -1;
  mClusterType = -1;
  mIsSeed = false;
  mNRawHits = -1;
  mNRawHitsR = -1;
  mNRawHitsPhi = -1;
  mRawHitVec.clear();
  mClusterId = -1;
  // for(int i_hit = 0; i_hit < FST::maxNHits; ++i_hit)
  // {
  //   mHitId[i_hit] = -1;
  // }

  // clearRawHitsList();
}

// FstRawHit
FstRawHit* FstCluster::createRawHit()
{
  if(mNumOfRawHits == mFstRawHits->GetSize())
    mFstRawHits->Expand( mNumOfRawHits + 10 );
  if(mNumOfRawHits >= 1000)
  {
    Fatal( "FstCluster::createRawHit()", "ERROR: Too many hits (>1000)!" );
    exit( 2 );
  }

  new((*mFstRawHits)[mNumOfRawHits++]) FstRawHit;
  return (FstRawHit*)((*mFstRawHits)[mNumOfRawHits - 1]);
}

void FstCluster::clearRawHitsList()
{
  mNumOfRawHits = 0;
  mFstRawHits->Clear();
}

int FstCluster::getNumRawHits() const
{
  return mNumOfRawHits;
}

FstRawHit* FstCluster::getRawHit(int i_hit) const
{
  return i_hit < mNumOfRawHits ? (FstRawHit*)((*mFstRawHits)[i_hit]) : NULL;
}
// FstRawHit
