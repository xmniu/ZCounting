#include "ZCounting/DataFillers/interface/FillerEventInfo.hh"
#include "ZCounting/DataFormats/interface/TEventInfo.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <TLorentzVector.h>
#include <string>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerEventInfo::FillerEventInfo(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC):
  fPUInfoName (iConfig.getUntrackedParameter<std::string>("edmPileupInfoName","addPileupInfo")),
  fBSName     (iConfig.getUntrackedParameter<std::string>("edmBeamspotName","offlineBeamSpot"))
{
  fPUInfoName_token = iC.consumes< std::vector<PileupSummaryInfo> >(fPUInfoName);
  fBSName_token = iC.consumes<reco::BeamSpot>(fBSName);
}

//--------------------------------------------------------------------------------------------------
FillerEventInfo::~FillerEventInfo(){}

//--------------------------------------------------------------------------------------------------
void FillerEventInfo::fill(TEventInfo *evtInfo,
                           const edm::Event &iEvent, const reco::Vertex &pv,
                           const bool hasGoodPV,
			   const TriggerBits triggerBits)//,TSusyGen *iSusyGen)
{
  assert(evtInfo);
  
  evtInfo->runNum  = iEvent.id().run();
  evtInfo->lumiSec = iEvent.luminosityBlock();
  evtInfo->evtNum  = iEvent.id().event();

  //
  // primary vertex info
  //==============================
  evtInfo->hasGoodPV  = hasGoodPV;
 
  //
  // fired triggers
  //==============================
  evtInfo->triggerBits = triggerBits;
}
