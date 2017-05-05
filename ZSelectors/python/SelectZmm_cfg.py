import FWCore.ParameterSet.Config as cms

process = cms.Process('RECODQM')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')

# load DQM
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")

# my analyzer
#process.load('ZCounting.ZSelectors.SelectZmm_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v3'

# trigger filter
import os
cmssw_base = os.environ['CMSSW_BASE']
hlt_filename = "ZCounting/DataFormats/data/HLT_50nsGRun"
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(False)
process.hltHighLevel.HLTPaths = cms.vstring()


from CondCore.CondDB.CondDB_cfi import *

hlt_file = open(cmssw_base + "/src/" + hlt_filename, "r")
for line in hlt_file.readlines():
  line = line.strip()              # strip preceding and trailing whitespaces
  if (line[0:3] == 'HLT'):         # assumes typical lines begin with HLT path name (e.g. HLT_Mu15_v1)
    hlt_path = line.split()[0]
    process.hltHighLevel.HLTPaths.extend(cms.untracked.vstring(hlt_path))

    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
    process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'/store/data/Run2016H/SingleMuon/RECO/PromptReco-v2/000/281/693/00000/000146F7-1686-E611-BD78-FA163E90B8E2.root',
'/store/data/Run2016H/SingleMuon/RECO/PromptReco-v2/000/281/727/00000/000744C8-D286-E611-BBEE-02163E012677.root'
                                )
                                )
    process.source.inputCommands = cms.untracked.vstring("keep *",
                                                         "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(False),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
  )

is_data_flag = True
do_hlt_filter = False
process.zcounting = cms.EDAnalyzer('SelectZmm',
                                 #outputName     = cms.untracked.string('Output.root'),
                                 skipOnHLTFail  = cms.untracked.bool(do_hlt_filter),
                                 TriggerFile    = cms.untracked.string(hlt_filename),
                                 TriggerEvent   = cms.InputTag('hltTriggerSummaryAOD','','HLT'),
                                 TriggerResults = cms.InputTag('TriggerResults','','HLT'),

                                 IDType   = cms.untracked.string("Medium"),# Tight, Medium, Loose
                                 IsoType  = cms.untracked.string("NULL"),  # Tracker-based, PF-based
                                 IsoCut   = cms.untracked.double(0.),     # {0.05, 0.10} for Tracker-based, {0.15, 0.25} for PF-based

                                 PtCutL1  = cms.untracked.double(27.0),
                                 PtCutL2  = cms.untracked.double(27.0),
                                 EtaCutL1 = cms.untracked.double(2.4),
                                 EtaCutL2 = cms.untracked.double(2.4),

                                 MassBin  = cms.untracked.int32(50),
                                 MassMin  = cms.untracked.double(66.0),
                                 MassMax  = cms.untracked.double(116.0),

                                 LumiBin  = cms.untracked.int32(400),
                                 LumiMin  = cms.untracked.double(0.0),
                                 LumiMax  = cms.untracked.double(2000.0),

                                 PVBin    = cms.untracked.int32(60),
                                 PVMin    = cms.untracked.double(0.0),
                                 PVMax    = cms.untracked.double(60.0),

                                 Info = cms.untracked.PSet(
    edmPileupInfoName    = cms.untracked.string('addPileupInfo'),
    edmBeamspotName      = cms.untracked.string('offlineBeamSpot')
    ),

                                 PV = cms.untracked.PSet(
    edmName       = cms.untracked.string('offlinePrimaryVertices'),
    minNTracksFit = cms.untracked.uint32(0),
    minNdof       = cms.untracked.double(4),
    maxAbsZ       = cms.untracked.double(24),
    maxRho        = cms.untracked.double(2)
    ),

                                 Muon = cms.untracked.PSet(
    minPt         = cms.untracked.double(20),
    edmName       = cms.untracked.string('muons'),
    edmPFCandName = cms.untracked.string('particleFlow'),

    # save general tracker tracks in our muon collection (used in tag-and-probe for muons)
    minTrackPt   = cms.untracked.double(20),
    edmTrackName = cms.untracked.string('generalTracks')

    ),
                                 )

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
                                     fileName = cms.untracked.string("OUT_step1.root"))

# Path and EndPath definitions
process.dqmoffline_step = cms.Path(process.zcounting)
process.dqmsave_step = cms.Path(process.DQMSaver)
#process.DQMoutput_step = cms.EndPath(process.DQMoutput)


# Schedule definition
process.schedule = cms.Schedule(
    process.dqmoffline_step,
#    process.DQMoutput_step
    process.dqmsave_step
    )

process.dqmSaver.workflow = '/SingleMuon/Run2016H-PromptReco-v2/RECO'
