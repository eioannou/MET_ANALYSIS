import FWCore.ParameterSet.Config as cms

process = cms.Process("Analyzer")

#__Options__#
RUNONMC = False

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')

#__Logger__#
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuples')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(1)
)

#process.load("FWCore.MessageService.MessageLogger_cfi")

#__MaX_Events__#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:pickevents.root'
        '/store/mc/RunIIFall17MiniAOD/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/0CF65340-0200-E811-ABB7-0025905C53F0.root'
        #'file:/afs/cern.ch/user/w/wverbeke/public/pickevents.root'
        )
                            )

#__Global_Tag__
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#if RUNONMC:
    #process.GlobalTag = GlobalTag(process.GlobalTag, "94X_mc2017_realistic_v12")
#    process.GlobalTag  = GlobalTag(process.GlobalTag, "80X_mcRun2_asymptotic_2016_TrancheIV_v6")
#if not RUNONMC:
#process.GlobalTag = GlobalTag( process.GlobalTag, "94X_dataRun2_ReReco_EOY17_v2")
process.GlobalTag = GlobalTag( process.GlobalTag, "94X_mc2017_realistic_v10")

#__Adding_JeC__#
#from CondCore.CondDB.CondDB_cfi import *
#CondDBSetup = CondDB.clone()
#CondDBSetup.__delattr__('connect')
#process.jec = cms.ESSource("PoolDBESSource",
#                           CondDBSetup,
#                           connect = cms.string('sqlite_file:/afs/cern.ch/work/a/aioannou/MET_STUDY/CMSSW_9_4_0/src/MET_Stuff/METStuffanalyzer/JECDatabase/SQLiteFiles/Fall17_17Nov2017F_V11_DATA.db'),
#                           toGet = cms.VPSet(
#        cms.PSet(
#            record = cms.string("JetCorrectionsRecord"),
#            tag = cms.string("JetCorrectorParametersCollection_Fall17_17Nov2017F_V11_DATA_AK4PFchs"),
#            label = cms.untracked.string("AK4PFchs")
#            )
#        )
#                           )

#process.es_prefer_jec = cms.ESPrefer("PoolDBESSource", "jec")

#_New_Jet
#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

#updateJetCollection(
#    process,
#    jetSource = cms.InputTag('slimmedJets'),
#    labelName = 'UpdatedJEC',
#    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')
#)



#________________________________FILTERS_____________________________________________#           

#__CSC_Halo_Filter__#                                                                            
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
process.CSCTightHaloFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.CSCTightHalo2015Filter_cfi')
process.CSCTightHalo2015Filter.taggingMode = cms.bool(True)
process.load('RecoMET.METFilters.CSCTightHaloTrkMuUnvetoFilter_cfi')
process.CSCTightHaloTrkMuUnvetoFilter.taggingMode = cms.bool(True)

#__Global_Halo_Filter__#                                                                         
process.load('RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi')
process.globalSuperTightHalo2016Filter.taggingMode = cms.bool(True)

#__HCAL_Noise_Filter__#                                                                          
process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion = cms.bool(False)
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('RecoMET.METFilters.HcalStripHaloFilter_cfi')
process.HcalStripHaloFilter.taggingMode = cms.bool(True)

#__Good_Vertex__#                                                                                
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


#__Primary_Vertex__#                                                                             
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

#__Updating_MET_with_new_JEC__                                                                   
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(
    process,
    isData = False,
    #fixEE2017=True,
    #fixEE2017Params = {'userawPt': False, 'PtThreshold':25.0, 'MinEtaThreshold':2.50, 'MaxEtaThreshold':3.139},
    #jetSelection = "pt>15 && abs(eta)<9.9 && !(pt<75. && abs(eta)>2.65 && abs(eta)<3.139)",
    #recoMetFromPFCs = True,
    #reclusterJets = True,
    #postfix = "ModifiedMET"

    )


#if RUNONMC:
#    process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
#    process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
#    process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
#    process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
#    process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
#    process.shifterPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")


process.Analyzer = cms.EDAnalyzer('MissingEtAnalyzer',
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
                                  HBHENoiseIsoFilter_Selector = cms.string("Flag_HBHENoiseIsoFilter"),
                                  EcalDeadCellTriggerPrimitiveNoiseFilter_Selector = cms.string("Flag_EcalDeadCellTriggerPrimitiveFilter"),
                                  EEBadScNoiseFilter_Selector = cms.string("Flag_eeBadScFilter"),
                                  BadPFMuonFilter_Selector = cms.string("Flag_BadPFMuonFilter"),
                                  BadChargedCandidateFilter_Selector = cms.string("Flag_BadChargedCandidateFilter"),
                                  ecalBadCalibFilter_Selector = cms.string("Flag_ecalBadCalibFilter")

                                  )

#__Output_File__                                                                                 
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("Significance_trial_.root")
                                   )






process.p = cms.EndPath(
    process.fullPatMetSequence*
    #process.patJets *
    #pfCandidatesGoodEE2017 *
    process.Analyzer
    )
