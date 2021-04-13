#ifndef FstEvent_h
#define FstEvent_h

#include <TObject.h>
#include <TClonesArray.h>
#include "./FstRawHit.h"
#include "./FstCluster.h"
#include "./FstTrack.h"

class FstEvent : public TObject
{
  private:
    int mNumOfHits;
    int mNumOfFstHits;
    int mNumOfClusters;
    int mNumOfFstClusters_Simple;
    int mNumOfFstClusters_Scan;
    int mNumOfTracks;
    int mNumOfRawHitTracks;
    int mNumOfClusterTracks;
    int mEventId;

    TClonesArray *mRawHits;      //->
    TClonesArray *mClusters;     //->
    TClonesArray *mTracks;       //->

  public:
    FstEvent();
    ~FstEvent();

    FstRawHit* createRawHit();
    void clearRawHitsList();
    int getNumRawHits() const;
    FstRawHit* getRawHit(int i_hit) const;
    void setNumFstRawHits(int numOfFstHits);
    int getNumFstRawHits() const;

    FstCluster* createCluster();
    void clearClustersList();
    int getNumClusters() const;
    FstCluster* getCluster(int i_cluster) const;
    void setNumFstClusters_Simple(int numOfFstClusters);
    int getNumFstClusters_Simple() const;
    void setNumFstClusters_Scan(int numOfFstClusters);
    int getNumFstClusters_Scan() const;

    FstTrack* createTrack();
    void clearTracksList();
    int getNumTracks() const;
    FstTrack* getTrack(int i_track) const;
    void setNumRawHitTraks(int numOfRawHitTracks);
    int getNumRawHitTraks() const;
    void setNumClusterTraks(int numOfClusterTracks);
    int getNumClusterTraks() const;

    void setEventId(int eventId);
    int getEventId() const;

    ClassDef(FstEvent,1)
};

#endif

