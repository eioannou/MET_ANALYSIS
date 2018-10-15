import FWCore.ParameterSet.Config as cms

process = cms.Process("Analyzer")

#__Options_:
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')

#__Logger_:
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuples')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(1)
)

#__MaxEvents_:
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#__InputFiles_:
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'file:pickevents.root'
        )
                            )

#__GlobalTag_:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag( process.GlobalTag, "94X_dataRun2_ReReco_EOY17_v2")

#__Adding JECs_:
from CondCore.CondDB.CondDB_cfi import *
CondDBSetup = CondDB.clone()
CondDBSetup.__delattr__('connect')
process.jec = cms.ESSource("PoolDBESSource",
                          CondDBSetup,
                          connect = cms.string('sqlite_file:/afs/cern.ch/work/a/aioannou/MET_STUDY/CMSSW_9_4_0/src/MET_Stuff/METStuffanalyzer/JECDatabase/SQLiteFiles/Fall17_17Nov2017F_V11_DATA.db'),
                          toGet = cms.VPSet(
        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag    = cms.string("JetCorrectorParametersCollection_Fall17_17Nov2017F_V11_DATA_AK4PFchs"),
            label  = cms.untracked.string("AK4PFchs")
            )
        )
                          )
process.es_prefer_jec = cms.ESPrefer("PollDBESSource", "jec")

#__CSC Halo Filter_:
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
process.CSCTightHaloFilter.taggingMode = cms.bool(True)
process.load('RecoMET.METFilters.CSCTightHalo2015Filter_cfi')
process.CSCTightHalo2015Filter.taggingMode = cms.bool(True)
process.load('RecoMET.METFilters.CSCTightHaloTrkMuUnvetoFilter_cfi')
process.CSCTightHaloTrkMuUnvetoFilter.taggingMode = cms.bool(True)

#__Global Halo Filter_:
process.load('RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi')
process.globalSuperTightHalo2016Filter.taggingMode = cms.bool(True)

#__HCAL Noise Filter_:
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS4ifJetInLowBVRegion = cms.bool(False)
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('RecoMET.METFilters.HcalStripHaloFilter_cfi')
process.HcalStripHaloFilter.taggingMode = cms.bool(True)

#__Tracking Failure FIlter_:
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.taggingMode = cms.bool(True)

#__Ecal Dead Cell Filter_:
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

#__EE Bad Sc Filter_:
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.taggingMode = cms.bool(True)

#__Bad Charged Candidate Filter_:
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

#__Bad PF Muon Filter_:
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.taggingMode = cms.bool(True)

#__Bad PF Muon Summer 16 Filter_:
process.load('RecoMET.METFilters.BadPFMuonSummer16Filter_cfi')
process.BadPFMuonSummer16Filter.taggingMode = cms.bool(True)


#__Good Vertices_:
process.goodVertices = cms.EDFilter(
    "VertexSelector",
    filter = cms.bool(False),
    src    = cms.InputTag("offlinePrimaryVertices"),
    cut    = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho <2")
)

#__Primary Vertex_:
process.primaryVertexFilter = cms.EDFilter(
    "GoodVertexFilter",
    vertexCollection = cms.InputTag( 'offlinePrimaryVertices'),
    minimumNDOF      = cms.uint32(4),
    maxAbsZ          = cms.double(24),
    maxd0            = cms.double(2)
    )

#__MET Producer_:
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(
    process,
    isData = True,
    )


#__Analyzer_:
process.Analyzer = cms.EDAnalyzer( 'MissingEtAnalyzer',
                                   metSrc  = cms.untracked.InputTag("slimmedMETs"),
                                   jetSrc  = cms.untracked.InputTag("slimmedJets"),
                                   muonSrc = cms.untracked.InputTag("slimmedMuons"),
                                   noiseFilterTag = cms.InputTag("TriggerResults"),
                                   verticesSrc = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                                   HLTSrc = cms.InputTag("TriggerResults", "", "HLT"),
                                   triggerobjects = cms.InputTag("slimmedPatTrigger"),
                                   GoodVtxNoiseFilter_Selector = cms.string("Flag_goodVertices"),
                                   GlobalTightHalo2016NoiseFilter_Selector = cms.string("Flag_globalTightHalo2016Filter"),
                                   HBHENoiseFilter_Selector = cms.string("Flag_HBHENoiseFilter"),
                                   EcalDeadCellTriggerPrimitiveNoiseFilter_Selector = cms.string("Flag_EcalDeadCellTriggerPrimitiveFilter"),
                                   EEBadScNoiseFilter_Selector = cms.string("Flag_eeBadScFilter"),
                                   BadPFMuonFilter_Selector = cms.string("Flag_BadPFMuonFilter"),
                                   BadChargedCandidateFilter_Selector = cms.string("Flag_BadChargedCandidateFilter"),
                                   ecalBadCalibFilter_Selector = cms.string("Flag_ecalBadCalibFilter")

                                   )

#__Output File_:
#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("test.root")
#                                   )

#__Paths Declaration_:
process.p = cms.EndPath(
    process.fullPatMetSequence*
    process.Analyzer
    )


    
