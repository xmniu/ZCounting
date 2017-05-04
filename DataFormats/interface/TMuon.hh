#ifndef ZCOUNTING_DATAFORMATS_TMUON_HH
#define ZCOUNTING_DATAFORMATS_TMUON_HH

#include "ZCounting/DataFormats/interface/MiniBaconDefs.hh"
#include <TObject.h>

namespace baconhep
{
  class TMuon : public TObject
  {
    public:
      TMuon():
      pt(0), eta(0), phi(0), ptErr(0),
      staPt(0), staEta(0), staPhi(0),
      pfPt(0), pfEta(0), pfPhi(0),
      trkIso(-1), ecalIso(-1), hcalIso(-1),
      chHadIso(-1), gammaIso(-1), neuHadIso(-1), puIso(-1),
      d0(-999.), dz(-999.), sip3d(-999.),
      validFraction(0),
      segmentCompatibility(0),
      tkNchi2(-999.), muNchi2(-999.),
      trkKink(0), glbKink(0), chi2LocalPosition(0),
      q(0),
      nValidHits(0),
      typeBits(0),
      nTkHits(0), nPixHits(0),
      nTkLayers(0), nPixLayers(0),
      nMatchStn(0),
      trkID(-1),
      hltMatchBits(0)
      {}
      ~TMuon(){}
    
      float          pt, eta, phi, ptErr;                   // kinematics
      float          staPt, staEta, staPhi;                 // STA track kinematics
      float          pfPt, pfEta, pfPhi;                    // matched PFCandidate
      float          trkIso, ecalIso, hcalIso;              // detector isolation
      float          chHadIso, gammaIso, neuHadIso, puIso;  // PF isolation variables
      float          d0, dz, sip3d;                         // impact parameter
      float	     validFraction;			    // fraction of valid tracker hits
      float	     segmentCompatibility;		    // segment compatibility
      float          tkNchi2, muNchi2;                      // track fit normalized chi-square
      float          trkKink, glbKink;                      // track kink
      float          chi2LocalPosition;			    // tracker-standalone position match
      int            q;                                     // charge
      int            nValidHits;                            // number of valid muon hits in global fit
      unsigned int   typeBits;                              // muon type bits
      unsigned int   nTkHits, nPixHits;                     // number of hits in tracker
      unsigned int   nTkLayers, nPixLayers;                 // number of hit layers in tracker
      unsigned int   nMatchStn;                             // number of stations with muon segments
      int            trkID;                                 // tracker track ID (unique per event)
      TriggerObjects hltMatchBits;                          // HLT matching
          
    ClassDef(TMuon,2)
  };

  enum EMuType
  {
    // following convention in DataFormats/MuonReco/interface/Muon.h
    kGlobal     = 2,
    kTracker    = 4,
    kStandalone = 8,
    kCaloMuon   = 16,
    kPFMuon     = 32,
    kRPCMuon    = 64
  };

}
#endif
