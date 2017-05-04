#ifndef ZCOUNTING_DATAFILLERS_FILLEREVENTINFO_HH
#define ZCOUNTING_DATAFILLERS_FILLEREVENTINFO_HH

#include <string>
#include "FWCore/Framework/interface/Frameworkfwd.h"     // declaration of EDM types
#include "FWCore/Framework/interface/EDAnalyzer.h"       // EDAnalyzer class
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // Parameters
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Framework/interface/MakerMacros.h"      // definitions for declaring plug-in modules
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "ZCounting/DataFormats/interface/MiniBaconDefs.hh"
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class TEventInfo;  // foward declaration
  class FillerEventInfo
  {
    public:
    FillerEventInfo(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);
      ~FillerEventInfo();
      
      void fill(TEventInfo         *evtInfo,       // output object to be filled
                const edm::Event   &iEvent,        // EDM event info
		const reco::Vertex &pv,            // event primary vertex
		const bool          hasGoodPV,     // flag for if PV passing cuts is found
		const TriggerBits   triggerBits);//,   // bits for corresponding fired triggers
//		TSusyGen           *susyGen=0);      // output for SUSY objects
	       
    protected:
    
      // EDM object collection names
      std::string fPUInfoName;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > fPUInfoName_token;
      std::string fBSName;
      edm::EDGetTokenT<reco::BeamSpot> fBSName_token;
  };
}
#endif
