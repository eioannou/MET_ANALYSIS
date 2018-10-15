import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuples')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(10)
)

process.maxEvents = cms.untracked.PSet( input =  cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall17MiniAOD/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/0CF65340-0200-E811-ABB7-0025905C53F0.root'
        )
                            )

#_Global_Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag                                           
process.GlobalTag = GlobalTag( process.GlobalTag, "94X_mc2017_realistic_v10")

#_JEC:
#

#_Filters:                                                                                           
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
process.CSCTightHaloFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.CSCTightHalo2015Filter_cfi')
process.CSCTightHalo2015Filter.taggingMode = cms.bool(True)
process.load('RecoMET.METFilters.CSCTightHaloTrkMuUnvetoFilter_cfi')
process.CSCTightHaloTrkMuUnvetoFilter.taggingMode = cms.bool(True)

                                                                                    
process.load('RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi')
process.globalSuperTightHalo2016Filter.taggingMode = cms.bool(True)

                                                                                        
process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion = cms.bool(False)
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('RecoMET.METFilters.HcalStripHaloFilter_cfi')
process.HcalStripHaloFilter.taggingMode = cms.bool(True)
                                                                                             
process.goodVertices = cms.EDFilter(
    "VertexSelector",
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
    )

process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadPFMuonSummer16Filter_cfi')
process.BadPFMuonSummer16Filter.taggingMode = cms.bool(True)

process.primaryVertexFilter = cms.EDFilter(
    "GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
    )

process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB = cms.vint32(12, 13, 14)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE = cms.vint32(12, 13, 14)

#__MET:
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(
    process,
    isData = False,
    postfix = "Modified"
)

process.TEST = cms.EDAnalyzer('MissingEtAnalyzer',
                              metSrc = cms.untracked.InputTag("slimmedMETs"),
                              metmodifiedSrc = cms.untracked.InputTag("slimmedMETsModified"),
                              jetSrc = cms.untracked.InputTag("slimmedJets"),
                              muonSrc = cms.untracked.InputTag("slimmedMuons"),
                              noiseFilterTag = cms.InputTag("TriggerResults"),
                              verticesSrc = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                              HLTSrc = cms.InputTag("TriggerResults", "", "HLT"),
                              triggerobjects = cms.InputTag("slimmedPatTrigger"),
                              GoodVtxNoiseFilter_Selector = cms.string("Flag_goodVertices"),
                              GlobalTightHalo2016NoiseFilter_Selector = cms.string("Flag_globalTightHalo2016Filter"),
                              HBHENoiseFilter_Selector = cms.string("Flag_HBHENoiseFilter"),
                              HBHENoiseIsoFilter_Selector = cms.string("Flag_HBHENoiseIsoFilter"),
                              EcalDeadCellTriggerPrimitiveNoiseFilter_Selector = cms.string("Flag_EcalDedCellTriggerPrimitiveFilter"),
                              EEBadScNoiseFilter_Selector = cms.string("Flag_eeBadScFilter"),
                              BadPFMuonFilter_Selector = cms.string("Flag_BadPFMuonFilter"),
                              BadChargedCandidateFilter_Selector = cms.string("Flag_BadChargedCandidateFilter"),
                              ecalBadCalibFilter_Selector = cms.string("Flag_ecalBadCalibFilter")
                              )

#_output_file:
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("METSignificance_test.root")
                                   )

process.p = cms.EndPath(
    process.fullPatMetSequenceModified*
    process.TEST
)
