#ifndef SelectZmm_H
#define SelectZmm_H

#include "FWCore/Framework/interface/MakerMacros.h"      // definitions for declaring plug-in modules
#include "FWCore/Framework/interface/Frameworkfwd.h"     // declaration of EDM types
#include "FWCore/Framework/interface/EDAnalyzer.h"       // EDAnalyzer class
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // Parameters
#include "FWCore/Utilities/interface/InputTag.h"
#include <string>                                        // string class
#include "ZCounting/DataFormats/interface/MiniBaconDefs.hh"
#include "ZCounting/DataFormats/interface/TMuon.hh"
#include "ZCounting/Utils/interface/TTrigger.hh"

#include <TMath.h>
#include <cassert>

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//#include <memory>
//#include <fstream>
//
////Framework
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Utilities/interface/InputTag.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//
////event
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/EventSetup.h"
//
////DQM
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
//
////Candidate handling
//#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/Candidate/interface/CandidateFwd.h"
//
//// Vertex utilities
//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//
//// Trigger
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/HLTReco/interface/TriggerObject.h"
//#include "DataFormats/HLTReco/interface/TriggerEvent.h"
//#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "FWCore/Common/interface/TriggerNames.h"
//
////Muon
//#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
//
////Tracker
//#include "DataFormats/TrackReco/interface/Track.h"
class TFile;
class TH1D;
class TTree;
class TClonesArray;
namespace edm {
  class TriggerResults;
  class TriggerNames;
}
namespace baconhep {
  class TEventInfo;
  class TTrigger;
  class FillerEventInfo;
  class FillerVertex;
  class FillerMuon;
}
 
class SelectZmm: public DQMEDAnalyzer{

public:

  SelectZmm(const edm::ParameterSet& ps);
  virtual ~SelectZmm();
  
protected:

  void dqmBeginRun(edm::Run const &, edm::EventSetup const &) override;
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(edm::Event const& e, edm::EventSetup const& eSetup) override;
  void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup) override ;
  void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup) override;
  void endRun(edm::Run const& run, edm::EventSetup const& eSetup) override;

private:
  //other functions
  bool isMuonTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits);
  bool isMuonTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits);
  bool passMuonID(const baconhep::TMuon *muon, const std::string idtype);
  bool passMuonIso(const baconhep::TMuon *muon, const std::string isotype, const float isocut);

  // specify trigger paths of interest
  void setTriggers();

  // initialization from HLT menu; needs to be called on every change in HLT menu
  void initHLT(const edm::TriggerResults&, const edm::TriggerNames&);

  bool fSkipOnHLTFail;

  // variables to handle triggers
  edm::ParameterSetID fTriggerNamesID;
  std::string         fHLTFile;
  edm::InputTag       fHLTObjTag;
  edm::InputTag       fHLTTag;
  edm::EDGetTokenT<trigger::TriggerEvent> fHLTObjTag_token;
  edm::EDGetTokenT<edm::TriggerResults> fHLTTag_token;
  edm::EDGetTokenT<reco::VertexCollection> fPVName_token;

  // AOD collection names
  std::string fPVName;

  // bacon fillers
  baconhep::FillerEventInfo *fFillerEvtInfo;
  baconhep::FillerVertex    *fFillerPV;
  baconhep::FillerMuon      *fFillerMuon;
  baconhep::TTrigger        *fTrigger;

  // Objects and arrays for output file
  std::string              fOutputName;
  TFile                   *fOutputFile;
  TH1D                    *fTotalEvents;
  TTree                   *fEventTree;
  baconhep::TEventInfo    *fEvtInfo;
  TClonesArray            *fMuonArr;
  TClonesArray            *fPVArr;

  std::string IDType_;
  std::string IsoType_;
  double IsoCut_;

  double PtCutL1_;
  double PtCutL2_;
  double EtaCutL1_;
  double EtaCutL2_;

  int    MassBin_;
  double MassMin_;
  double MassMax_;

  int    LumiBin_;
  double LumiMin_;
  double LumiMax_;

  int    PVBin_;
  double PVMin_;
  double PVMax_;

  const Double_t MUON_MASS  = 0.105658369;
  const Double_t MUON_BOUND = 0.9;

  // Histograms
  MonitorElement* h_mass_HLT_pass_central;
  MonitorElement* h_mass_HLT_pass_forward;
  MonitorElement* h_mass_HLT_fail_central;
  MonitorElement* h_mass_HLT_fail_forward;

  MonitorElement* h_mass_SIT_pass_central;
  MonitorElement* h_mass_SIT_pass_forward;
  MonitorElement* h_mass_SIT_fail_central;
  MonitorElement* h_mass_SIT_fail_forward;

  MonitorElement* h_mass_Sta_pass_central;
  MonitorElement* h_mass_Sta_pass_forward;
  MonitorElement* h_mass_Sta_fail_central;
  MonitorElement* h_mass_Sta_fail_forward;

  MonitorElement* h_npv;
  MonitorElement* h_yield_Z;
};


#endif
