// bacon classes and constants
 #include "ZCounting/DataFormats/interface/MiniBaconDefs.hh"
 #include "ZCounting/DataFormats/interface/TEventInfo.hh"
 #include "ZCounting/DataFormats/interface/TMuon.hh"
 #include "ZCounting/DataFormats/interface/TVertex.hh"
 #include "ZCounting/Utils/interface/TTrigger.hh"

 #include "ZCounting/DataFillers/interface/FillerEventInfo.hh"
 #include "ZCounting/DataFillers/interface/FillerVertex.hh"
 #include "ZCounting/DataFillers/interface/FillerMuon.hh"

 // tools to parse HLT name patterns
 #include <boost/foreach.hpp>
 #include "FWCore/Utilities/interface/RegexMatch.h"

 #include "FWCore/Common/interface/TriggerNames.h"
 #include "DataFormats/Common/interface/TriggerResults.h"
 #include "DataFormats/HLTReco/interface/TriggerEvent.h"

 // data format classes
 #include "FWCore/Framework/interface/Event.h"
 #include "DataFormats/Common/interface/Handle.h"
 #include "DataFormats/VertexReco/interface/Vertex.h"

 // ROOT classes
 #include <TFile.h>
 #include <TH1D.h>
 #include <TTree.h>
 #include <TClonesArray.h>
 #include <TLorentzVector.h>
 #include <TMath.h>

#include "ZCounting/ZSelectors/plugins/SelectZmm.hh"
//#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Utilities/interface/Exception.h"
//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
//
//#include "FWCore/Framework/interface/Event.h"
//#include "DataFormats/Common/interface/Handle.h"
//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/HLTReco/interface/TriggerEvent.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "TrackingTools/IPTools/interface/IPTools.h"
//#include <TClonesArray.h>
//#include <TLorentzVector.h>
//#include <TMath.h>
//
//#include "TLorentzVector.h"
//#include "TFile.h"
//#include <vector>
//#include <algorithm>
//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <iomanip>
//#include <stdio.h>
//#include <string>
//#include <sstream>
//#include <math.h>

//
// -------------------------------------- Constructor --------------------------------------------
//
SelectZmm::SelectZmm(const edm::ParameterSet& iConfig):
  fSkipOnHLTFail     (iConfig.getUntrackedParameter<bool>("skipOnHLTFail",false)),
  fHLTFile           (iConfig.getUntrackedParameter<std::string>("TriggerFile","HLT")),
  fHLTObjTag         (iConfig.getParameter<edm::InputTag>("TriggerEvent")),
  fHLTTag            (iConfig.getParameter<edm::InputTag>("TriggerResults")),
  fPVName (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fFillerEvtInfo     (0),
  fFillerPV          (0),
  fFillerMuon        (0),
  fTrigger           (0),
  fOutputName        (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile        (0),
  fTotalEvents       (0),
  fEventTree         (0),
  fEvtInfo           (0),
  fMuonArr           (0),
  fPVArr             (0)
{
  edm::LogInfo("SelectZmm") <<  "Constructor  SelectZmm::SelectZmm " << std::endl;
  std::cout<<"Constructor  SelectZmm::SelectZmm"<<std::endl;
  
  //Get parameters from configuration file
  fHLTTag_token = consumes<edm::TriggerResults>(fHLTTag);
  fHLTObjTag_token = consumes<trigger::TriggerEvent>(fHLTObjTag);
  fPVName_token = consumes<reco::VertexCollection>(fPVName);

  // Don't write TObject part of the objects
  baconhep::TEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TMuon::Class()->IgnoreTObjectStreamer();
  baconhep::TVertex::Class()->IgnoreTObjectStreamer();

  // Set up bacon objects and configure fillers
  if(iConfig.existsAs<edm::ParameterSet>("Info",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Info"));

    fEvtInfo       = new baconhep::TEventInfo();                             assert(fEvtInfo);
    fFillerEvtInfo = new baconhep::FillerEventInfo(cfg,consumesCollector()); assert(fFillerEvtInfo);
  }

  if(iConfig.existsAs<edm::ParameterSet>("PV",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("PV"));

    fPVArr         = new TClonesArray("baconhep::TVertex");                  assert(fPVArr);
    fFillerPV      = new baconhep::FillerVertex(cfg);                        assert(fFillerPV);
  }

  if(iConfig.existsAs<edm::ParameterSet>("Muon",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Muon"));

    fMuonArr       = new TClonesArray("baconhep::TMuon");                    assert(fMuonArr);
    fFillerMuon    = new baconhep::FillerMuon(cfg,consumesCollector());      assert(fFillerMuon);
  }

  //Cuts
  IDType_   = iConfig.getUntrackedParameter<std::string>("IDType");
  IsoType_  = iConfig.getUntrackedParameter<std::string>("IsoType");
  IsoCut_   = iConfig.getUntrackedParameter<double>("IsoCut");

  PtCutL1_  = iConfig.getUntrackedParameter<double>("PtCutL1");
  PtCutL2_  = iConfig.getUntrackedParameter<double>("PtCutL2");
  EtaCutL1_ = iConfig.getUntrackedParameter<double>("EtaCutL1");
  EtaCutL2_ = iConfig.getUntrackedParameter<double>("EtaCutL2");

  MassBin_  = iConfig.getUntrackedParameter<int>("MassBin");
  MassMin_  = iConfig.getUntrackedParameter<double>("MassMin");
  MassMax_  = iConfig.getUntrackedParameter<double>("MassMax");

  LumiBin_  = iConfig.getUntrackedParameter<int>("LumiBin");
  LumiMin_  = iConfig.getUntrackedParameter<double>("LumiMin");
  LumiMax_  = iConfig.getUntrackedParameter<double>("LumiMax");

  PVBin_    = iConfig.getUntrackedParameter<int>("PVBin");
  PVMin_    = iConfig.getUntrackedParameter<double>("PVMin");
  PVMax_    = iConfig.getUntrackedParameter<double>("PVMax");
}

//
//  -------------------------------------- Destructor --------------------------------------------
//
SelectZmm::~SelectZmm()
{
  edm::LogInfo("SelectZmm") <<  "Destructor SelectZmm::~SelectZmm " << std::endl;
  std::cout<<"Destructor SelectZmm::~SelectZmm "<<std::endl;

  delete fFillerEvtInfo;
  delete fFillerPV;
  delete fFillerMuon;

  delete fTrigger;
  delete fEvtInfo;
  delete fMuonArr;
  delete fPVArr;
}

//
// -------------------------------------- beginRun --------------------------------------------
//
void SelectZmm::dqmBeginRun(edm::Run const &, edm::EventSetup const &)
{
  edm::LogInfo("SelectZmm") <<  "SelectZmm::beginRun" << std::endl;
  std::cout<<"SelectZmm::beginRun"<<std::endl;

  // Create output file, trees, and histograms
  fOutputFile  = new TFile(fOutputName.c_str(), "RECREATE");
  fTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);
  fEventTree   = new TTree("Events","Events");

  fEventTree->Branch("Info",fEvtInfo);
  fEventTree->Branch("Muon",     &fMuonArr);
  fEventTree->Branch("PV",       &fPVArr);

  // Triggers
  setTriggers();
}
//
// -------------------------------------- bookHistos --------------------------------------------
//
void SelectZmm::bookHistograms(DQMStore::IBooker & ibooker_, edm::Run const &, edm::EventSetup const &)
{
  edm::LogInfo("SelectZmm") <<  "SelectZmm::bookHistograms" << std::endl;
  std::cout<<"SelectZmm::bookHistograms"<<std::endl;

  ibooker_.cd();
  ibooker_.setCurrentFolder("ZCounting/Histograms");
/*
  h_mass_HLT_pass_central = ibooker_.book1D("h_mass_HLT_pass_central", "h_mass_HLT_pass_central", MassBin_, MassMin_, MassMax_);
  h_mass_HLT_pass_forward = ibooker_.book1D("h_mass_HLT_pass_forward", "h_mass_HLT_pass_forward", MassBin_, MassMin_, MassMax_);
  h_mass_HLT_fail_central = ibooker_.book1D("h_mass_HLT_fail_central", "h_mass_HLT_fail_central", MassBin_, MassMin_, MassMax_);
  h_mass_HLT_fail_forward = ibooker_.book1D("h_mass_HLT_fail_forward", "h_mass_HLT_fail_forward", MassBin_, MassMin_, MassMax_);

  h_mass_SIT_pass_central = ibooker_.book1D("h_mass_SIT_pass_central", "h_mass_SIT_pass_central", MassBin_, MassMin_, MassMax_);
  h_mass_SIT_pass_forward = ibooker_.book1D("h_mass_SIT_pass_forward", "h_mass_SIT_pass_forward", MassBin_, MassMin_, MassMax_);
  h_mass_SIT_fail_central = ibooker_.book1D("h_mass_SIT_fail_central", "h_mass_SIT_fail_central", MassBin_, MassMin_, MassMax_);
  h_mass_SIT_fail_forward = ibooker_.book1D("h_mass_SIT_fail_forward", "h_mass_SIT_fail_forward", MassBin_, MassMin_, MassMax_);

  h_mass_Sta_pass_central = ibooker_.book1D("h_mass_Sta_pass_central", "h_mass_Sta_pass_central", MassBin_, MassMin_, MassMax_);
  h_mass_Sta_pass_forward = ibooker_.book1D("h_mass_Sta_pass_forward", "h_mass_Sta_pass_forward", MassBin_, MassMin_, MassMax_);
  h_mass_Sta_fail_central = ibooker_.book1D("h_mass_Sta_fail_central", "h_mass_Sta_fail_central", MassBin_, MassMin_, MassMax_);
  h_mass_Sta_fail_forward = ibooker_.book1D("h_mass_Sta_fail_forward", "h_mass_Sta_fail_forward", MassBin_, MassMin_, MassMax_);
*/
  h_mass_HLT_pass_central = ibooker_.book2D("h_mass_HLT_pass_central", "h_mass_HLT_pass_central", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
  h_mass_HLT_pass_forward = ibooker_.book2D("h_mass_HLT_pass_forward", "h_mass_HLT_pass_forward", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
  h_mass_HLT_fail_central = ibooker_.book2D("h_mass_HLT_fail_central", "h_mass_HLT_fail_central", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
  h_mass_HLT_fail_forward = ibooker_.book2D("h_mass_HLT_fail_forward", "h_mass_HLT_fail_forward", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

  h_mass_SIT_pass_central = ibooker_.book2D("h_mass_SIT_pass_central", "h_mass_SIT_pass_central", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
  h_mass_SIT_pass_forward = ibooker_.book2D("h_mass_SIT_pass_forward", "h_mass_SIT_pass_forward", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
  h_mass_SIT_fail_central = ibooker_.book2D("h_mass_SIT_fail_central", "h_mass_SIT_fail_central", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
  h_mass_SIT_fail_forward = ibooker_.book2D("h_mass_SIT_fail_forward", "h_mass_SIT_fail_forward", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

  h_mass_Sta_pass_central = ibooker_.book2D("h_mass_Sta_pass_central", "h_mass_Sta_pass_central", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
  h_mass_Sta_pass_forward = ibooker_.book2D("h_mass_Sta_pass_forward", "h_mass_Sta_pass_forward", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
  h_mass_Sta_fail_central = ibooker_.book2D("h_mass_Sta_fail_central", "h_mass_Sta_fail_central", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
  h_mass_Sta_fail_forward = ibooker_.book2D("h_mass_Sta_fail_forward", "h_mass_Sta_fail_forward", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

  h_npv                   = ibooker_.book2D("h_npv",     "h_npv",     LumiBin_, LumiMin_, LumiMax_, PVBin_, PVMin_, PVMax_);
  h_yield_Z               = ibooker_.book1D("h_yield_Z", "h_yield_Z", LumiBin_, LumiMin_, LumiMax_);

}
//
// -------------------------------------- beginLuminosityBlock --------------------------------------------
//
void SelectZmm::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, 
                                            edm::EventSetup const& context) 
{
  edm::LogInfo("SelectZmm") <<  "SelectZmm::beginLuminosityBlock" << std::endl;
  std::cout<<"SelectZmm::beginLuminosityBlock"<<std::endl;

}


//
// -------------------------------------- Analyze --------------------------------------------
//
//--------------------------------------------------------------------------------------------------
void SelectZmm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{// Fill event tree on the fly 
  edm::LogInfo("SelectZmm") <<  "SelectZmm::analyze" << std::endl;
//  std::cout<<"SelectZmm::analyze"<<std::endl;
//  std::cout<<iEvent.luminosityBlock()<<std::endl;

  fTotalEvents->Fill(1);

  edm::Handle<edm::TriggerResults> hTrgRes;
  iEvent.getByToken(fHLTTag_token,hTrgRes);
  assert(hTrgRes.isValid());

  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*hTrgRes);
  Bool_t config_changed = false;
  if(fTriggerNamesID != triggerNames.parameterSetID()) {
    fTriggerNamesID = triggerNames.parameterSetID();
    config_changed  = true;
  }
  if(config_changed) {
    initHLT(*hTrgRes, triggerNames);
  }

  TriggerBits triggerBits; 
  for(unsigned int irec=0; irec<fTrigger->fRecords.size(); irec++) {
    if(fTrigger->fRecords[irec].hltPathIndex == (unsigned int)-1) continue;
    if(hTrgRes->accept(fTrigger->fRecords[irec].hltPathIndex)) {
      triggerBits [fTrigger->fRecords[irec].baconTrigBit] = 1;
    }
  }
  if(fSkipOnHLTFail && triggerBits == 0) return;

  // Fill vertex
  fPVArr->Clear();
  int nvertices = 0; 
  const reco::Vertex *pv = fFillerPV->fill(fPVArr, nvertices, iEvent);
  assert(pv);

  // Fill event info 
  fFillerEvtInfo->fill(fEvtInfo, iEvent, *pv, (nvertices>0), triggerBits);//,fSusyGen);

  edm::Handle<trigger::TriggerEvent> hTrgEvt;
  iEvent.getByToken(fHLTObjTag_token,hTrgEvt);

  // Fill muon info 
  fMuonArr->Clear();  
  fFillerMuon->fill(fMuonArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);

  fEventTree->Fill();

}

//
// -------------------------------------- endLuminosityBlock --------------------------------------------
//
void SelectZmm::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup) 
{
  edm::LogInfo("SelectZmm") <<  "SelectZmm::endLuminosityBlock" << std::endl;
  std::cout<<"SelectZmm::endLuminosityBlock"<<std::endl;
}

//
// -------------------------------------- endRun --------------------------------------------
//
void SelectZmm::endRun(edm::Run const& run, edm::EventSetup const& eSetup)
{
  edm::LogInfo("SelectZmm") <<  "SelectZmm::endRun" << std::endl;
  std::cout<<"SelectZmm::endRun"<<std::endl;

  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();

  std::string cmssw_base_src = getenv("CMSSW_BASE"); cmssw_base_src+="/src/";
  const baconhep::TTrigger triggerMenu(cmssw_base_src + fHLTFile);
//  const baconhep::TTrigger triggerMenu("/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_9_0_0/src/ZCounting/DataFormats/data/HLT_50nsGRun");
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  TClonesArray *muonArr    = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  TFile *infile=0;
  TTree *eventTree=0;

  infile = TFile::Open(fOutputName.c_str());
  assert(infile);
  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);
  eventTree->SetBranchAddress("Info", &info);      TBranch *infoBr = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon", &muonArr);   TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");

  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    infoBr->GetEntry(ientry);

    vertexArr->Clear();
    vertexBr->GetEntry(ientry);
    h_npv->Fill(info->lumiSec,vertexArr->GetEntries());

    // Trigger requirement
    if(!isMuonTrigger(triggerMenu, info->triggerBits)) continue;

    // Good vertex requirement
    if(!(info->hasGoodPV)) continue;

    muonArr->Clear();
    muonBr->GetEntry(ientry);

    TLorentzVector vTag(0,0,0,0);
    for(Int_t i1=0; i1<muonArr->GetEntriesFast(); i1++){
      const baconhep::TMuon *tag = (baconhep::TMuon*)((*muonArr)[i1]);

      // Tag selection
      if(tag->pt        < PtCutL1_)        continue;  // lepton pT cut
      if(fabs(tag->eta) > EtaCutL1_)       continue;  // lepton |eta| cut
      if(!(passMuonID(tag, IDType_) && passMuonIso(tag, IsoType_, IsoCut_))) continue;  // lepton selection
      if(!isMuonTriggerObj(triggerMenu, tag->hltMatchBits)) continue; // trigger matching

      vTag.SetPtEtaPhiM(tag->pt,tag->eta,tag->phi,MUON_MASS);

      TLorentzVector vProbe(0,0,0,0);
      for(Int_t i2=0; i2<muonArr->GetEntriesFast(); i2++){
        if(i1==i2) continue;
        const baconhep::TMuon *probe = (baconhep::TMuon*)((*muonArr)[i2]);

        // Probe selection
        if(probe->pt        < PtCutL2_)  continue;  // lepton pT cut
        if(fabs(probe->eta) > EtaCutL2_) continue;  // lepton |eta| cut
        if(tag->q == probe->q)           continue;    // opposite charge requirement

        vProbe.SetPtEtaPhiM(probe->pt,probe->eta,probe->phi,MUON_MASS);

        // Mass window
        TLorentzVector vDilep = vTag + vProbe;
        if((vDilep.M() < MassMin_) || (vDilep.M() > MassMax_)) continue;

        bool isProbeCentral = false;
        bool isTagCentral = false;
        if(fabs(probe->eta) < MUON_BOUND) isProbeCentral = true;
        if(fabs(tag->eta) < MUON_BOUND)   isTagCentral   = true;

        // Determine event category for efficiency calculation
        if(passMuonID(probe, IDType_) && passMuonIso(probe, IsoType_, IsoCut_)){
          if(isMuonTriggerObj(triggerMenu, probe->hltMatchBits)){
	    // category 2HLT: both muons passing trigger requirements 
            if(i1>i2) continue;  // make sure we don't double count MuMu2HLT category
            if(isProbeCentral){
              h_mass_HLT_pass_central->Fill(info->lumiSec, vDilep.M());
              h_mass_SIT_pass_central->Fill(info->lumiSec, vDilep.M());
              h_mass_Sta_pass_central->Fill(info->lumiSec, vDilep.M());
            }
            else
            {
              h_mass_HLT_pass_forward->Fill(info->lumiSec, vDilep.M());
              h_mass_SIT_pass_forward->Fill(info->lumiSec, vDilep.M());
              h_mass_Sta_pass_forward->Fill(info->lumiSec, vDilep.M());
            }

            if(isTagCentral){
              h_mass_HLT_pass_central->Fill(info->lumiSec, vDilep.M());
              h_mass_SIT_pass_central->Fill(info->lumiSec, vDilep.M());
              h_mass_Sta_pass_central->Fill(info->lumiSec, vDilep.M());
            }
            else
            {
              h_mass_HLT_pass_forward->Fill(info->lumiSec, vDilep.M());
              h_mass_SIT_pass_forward->Fill(info->lumiSec, vDilep.M());
              h_mass_Sta_pass_forward->Fill(info->lumiSec, vDilep.M());
            }

          }
          else{
            // category 1HLT: probe passing selection but not trigger
            if(isProbeCentral){
              h_mass_HLT_fail_central->Fill(info->lumiSec, vDilep.M());
              h_mass_SIT_pass_central->Fill(info->lumiSec, vDilep.M());
              h_mass_Sta_pass_central->Fill(info->lumiSec, vDilep.M());
            }
            else
            {
              h_mass_HLT_fail_forward->Fill(info->lumiSec, vDilep.M());
              h_mass_SIT_pass_forward->Fill(info->lumiSec, vDilep.M());
              h_mass_Sta_pass_forward->Fill(info->lumiSec, vDilep.M());
            }
          }

          h_yield_Z->Fill(info->lumiSec);
        }
        else if(probe->typeBits & baconhep::EMuType::kGlobal){
          // category NoSel: probe is a GLB muon but failing selection 
          if(isProbeCentral){
            h_mass_SIT_fail_central->Fill(info->lumiSec, vDilep.M());
            h_mass_Sta_pass_central->Fill(info->lumiSec, vDilep.M());
          }
          else
          {
            h_mass_SIT_fail_forward->Fill(info->lumiSec, vDilep.M());
            h_mass_Sta_pass_forward->Fill(info->lumiSec, vDilep.M());
          }
        }
        else if(probe->typeBits & baconhep::EMuType::kStandalone){
          // category STA: probe is a STA muon
          if(isProbeCentral){
            h_mass_Sta_fail_central->Fill(info->lumiSec, vDilep.M());
          }
          else
          {
            h_mass_Sta_fail_forward->Fill(info->lumiSec, vDilep.M());
          }
        }
        else if(probe->nTkLayers>=6 && probe->nPixHits>=1){
          // cateogry Trk: probe is a tracker track
          if(isProbeCentral){
            h_mass_Sta_fail_central->Fill(info->lumiSec, vDilep.M());
          }
          else
          {
            h_mass_Sta_fail_forward->Fill(info->lumiSec, vDilep.M());
          }
        }

      }// End of probe loop
    }// End of tag loop
  }// End of event loop
}

//
// -------------------------------------- functions --------------------------------------------
//
void SelectZmm::setTriggers()
{
  std::string cmssw_base_src = getenv("CMSSW_BASE"); cmssw_base_src+="/src/";
  fTrigger = new baconhep::TTrigger(cmssw_base_src + fHLTFile);
}

//--------------------------------------------------------------------------------------------------
void SelectZmm::initHLT(const edm::TriggerResults& result, const edm::TriggerNames& triggerNames)
{
  for(unsigned int irec=0; irec<fTrigger->fRecords.size(); irec++) {
    fTrigger->fRecords[irec].hltPathName  = "";
    fTrigger->fRecords[irec].hltPathIndex = (unsigned int)-1;
    const std::string pattern = fTrigger->fRecords[irec].hltPattern;
    if(edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
      std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(triggerNames.triggerNames(), pattern);
      if(matches.empty()) {
        std::cout << "requested pattern [" << pattern << "] does not match any HLT paths" << std::endl;
      } else {
        BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches) {
          fTrigger->fRecords[irec].hltPathName = *match;
        }
      }
    } else {  // take full HLT path name given
      fTrigger->fRecords[irec].hltPathName = pattern;
    }
    // Retrieve index in trigger menu corresponding to HLT path
    unsigned int index = triggerNames.triggerIndex(fTrigger->fRecords[irec].hltPathName);
    if(index < result.size()) {  // check for valid index
      fTrigger->fRecords[irec].hltPathIndex = index;
    }
  }
}

//--------------------------------------------------------------------------------------------------
bool SelectZmm::isMuonTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits)
{
  return triggerMenu.pass("HLT_IsoMu24_v*",hltBits);
}

//--------------------------------------------------------------------------------------------------
bool SelectZmm::isMuonTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits)
{
  return triggerMenu.passObj("HLT_IsoMu24_v*","hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09",hltMatchBits);
}

//--------------------------------------------------------------------------------------------------
bool SelectZmm::passMuonID(const baconhep::TMuon *muon, const std::string idtype)
{
  if(idtype == "Loose"){

    if(!(muon->typeBits & baconhep::EMuType::kGlobal) && !(muon->typeBits & baconhep::EMuType::kTracker)) return kFALSE;
    if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;

    return kTRUE;
  }
  else if(idtype == "Medium"){

    if(!(muon->typeBits & baconhep::EMuType::kGlobal) && !(muon->typeBits & baconhep::EMuType::kTracker)) return kFALSE;
    if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
    if(muon->validFraction <= 0.8) return kFALSE;

    bool goodGlob = muon->typeBits & baconhep::EMuType::kGlobal && muon->muNchi2 < 3 && muon->chi2LocalPosition < 12 && muon->trkKink < 20;
    bool isMedium = muon->segmentCompatibility > (goodGlob ? 0.303 : 0.451);

    return isMedium;
  }
  else if(idtype == "Tight"){

    if(!(muon->typeBits & baconhep::EMuType::kGlobal)) return kFALSE;
    if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
    if(muon->muNchi2    > 10)                          return kFALSE;
    if(muon->nValidHits < 1)                           return kFALSE;
    if(muon->nMatchStn  < 2)                           return kFALSE;
    if(fabs(muon->d0)   > 0.2)                         return kFALSE;
    if(fabs(muon->dz)   > 0.5)                         return kFALSE;
    if(muon->nPixHits   < 1)                           return kFALSE;
    if(muon->nTkLayers  < 6)                           return kFALSE;

    return kTRUE;
  }
  else if(idtype == "NULL"){
   return kTRUE;
  }
  else{
    std::cout<<"Wrong input for IDType..."<<std::endl;
    assert(0);
  }
}

//--------------------------------------------------------------------------------------------------
bool SelectZmm::passMuonIso(const baconhep::TMuon *muon, const std::string isotype, const float isocut)
{
  if(isotype == "Tracker-based"){
    float iso = muon->trkIso;
    return (iso < isocut);
  }
  else if(isotype == "PF-based"){
    float iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso),Double_t(0));
    return (iso < isocut);
  }
  else if(isotype == "NULL"){
    return kTRUE;
  }
  else{
    std::cout<<"Wrong input for IsoType..."<<std::endl;
    assert(0);
  }
}

DEFINE_FWK_MODULE(SelectZmm);
