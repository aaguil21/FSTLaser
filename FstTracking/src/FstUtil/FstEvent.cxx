#include <iostream>
#include "./FstEvent.h"

using namespace std;

ClassImp(FstEvent)

//------------------------------------------
FstEvent::FstEvent()
{
  mNumOfHits               = 0;
  mNumOfFstHits            = 0;
  mNumOfClusters           = 0;
  mNumOfFstClusters_Simple = 0;
  mNumOfFstClusters_Scan   = 0;
  mNumOfTracks             = 0;
  mNumOfRawHitTracks       = 0;
  mNumOfClusterTracks      = 0;
  mEventId                 = -1;

  mRawHits  = new TClonesArray("FstRawHit", 10);
  mClusters = new TClonesArray("FstCluster", 10);
  mTracks   = new TClonesArray("FstTrack", 10);
}

FstEvent::~FstEvent()
{
  delete mRawHits;
  mRawHits = NULL;
  delete mClusters;
  mClusters = NULL;
  delete mTracks;
  mTracks = NULL;
}
//------------------------------------------

// FstRawHit
FstRawHit* FstEvent::createRawHit()
{
  if(mNumOfHits == mRawHits->GetSize())
    mRawHits->Expand( mNumOfHits + 10 );
  // if(mNumOfHits >= 25)
  if(mNumOfHits >= 1000)
  {
    Fatal( "FstEvent::createRawHit()", "ERROR: Too many hits (>1000)!" );
    exit( 2 );
  }

  new((*mRawHits)[mNumOfHits++]) FstRawHit;
  return (FstRawHit*)((*mRawHits)[mNumOfHits - 1]);
}

void FstEvent::clearRawHitsList()
{
  mNumOfHits = 0;
  mNumOfFstHits = 0;
  mRawHits->Clear();
}

int FstEvent::getNumRawHits() const
{
  return mNumOfHits;
}

FstRawHit* FstEvent::getRawHit(int i_hit) const
{
  return i_hit < mNumOfHits ? (FstRawHit*)((*mRawHits)[i_hit]) : NULL;
}

void FstEvent::setNumFstRawHits(int numOfFstHits)
{
  mNumOfFstHits = numOfFstHits;
}

int FstEvent::getNumFstRawHits() const
{
  return mNumOfFstHits;
}
// FstRawHit

// FstCluster
FstCluster* FstEvent::createCluster()
{
  if(mNumOfClusters == mClusters->GetSize())
    mClusters->Expand( mNumOfClusters + 10 );
  if(mNumOfClusters >= 1000)
  {
    Fatal( "FstEvent::createCluster()", "ERROR: Too many clusters (>1000)!" );
    exit( 2 );
  }

  new((*mClusters)[mNumOfClusters++]) FstCluster;
  return (FstCluster*)((*mClusters)[mNumOfClusters - 1]);
}

void FstEvent::clearClustersList()
{
  mNumOfClusters = 0;
  mNumOfFstClusters_Simple = 0;
  mNumOfFstClusters_Scan = 0;
  mClusters->Clear();
}

int FstEvent::getNumClusters() const
{
  return mNumOfClusters;
}

FstCluster* FstEvent::getCluster(int i_cluster) const
{
  return i_cluster < mNumOfClusters ? (FstCluster*)((*mClusters)[i_cluster]) : NULL;
}

void FstEvent::setNumFstClusters_Simple(int numOfFstClusters)
{
  mNumOfFstClusters_Simple = numOfFstClusters;
}

int FstEvent::getNumFstClusters_Simple() const
{
  return mNumOfFstClusters_Simple;
}

void FstEvent::setNumFstClusters_Scan(int numOfFstClusters)
{
  mNumOfFstClusters_Scan = numOfFstClusters;
}

int FstEvent::getNumFstClusters_Scan() const
{
  return mNumOfFstClusters_Scan;
}
// FstCluster

// FstTrack
FstTrack* FstEvent::createTrack()
{
  if(mNumOfTracks == mTracks->GetSize())
    mTracks->Expand( mNumOfTracks + 10 );
  if(mNumOfTracks >= 1000)
  {
    Fatal( "FstEvent::createTrack()", "ERROR: Too many tracks (>1000)!" );
    exit( 2 );
  }

  new((*mTracks)[mNumOfTracks++]) FstTrack;
  return (FstTrack*)((*mTracks)[mNumOfTracks - 1]);
}

void FstEvent::clearTracksList()
{
  mNumOfTracks = 0;
  mTracks->Clear();
}

int FstEvent::getNumTracks() const
{
  return mNumOfTracks;
}

FstTrack* FstEvent::getTrack(int i_track) const
{
  return i_track < mNumOfTracks ? (FstTrack*)((*mTracks)[i_track]) : NULL;
}

void FstEvent::setNumRawHitTraks(int numOfRawHitTracks)
{
  mNumOfRawHitTracks = numOfRawHitTracks;
}
int FstEvent::getNumRawHitTraks() const
{
  return mNumOfRawHitTracks;
}
void FstEvent::setNumClusterTraks(int numOfClusterTracks)
{
  mNumOfClusterTracks = numOfClusterTracks;
}
int FstEvent::getNumClusterTraks() const
{
  return mNumOfClusterTracks;
}
// FstTrack

// Event Id
void FstEvent::setEventId(int eventId)
{
  mEventId = eventId;
}

int FstEvent::getEventId() const
{
  return mEventId;
}
// Event Id
