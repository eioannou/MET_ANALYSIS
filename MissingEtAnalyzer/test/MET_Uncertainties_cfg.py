import FWCore.ParameterSet.Config as cms

process = cms.Process("METUncertainties")

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/483975CE-8142-E811-A72D-001E67E71E20.root'
        )
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag( process.GlobalTag, "94X_mc2017_realistic_v14")

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

# Puppi
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD(process, True)


runMetCorAndUncFromMiniAOD(
    process,
    isData = False,
    postfix = "Modified"
    )

runMetCorAndUncFromMiniAOD(process,
                           isData = False,
                           metType= "Puppi",
                           jetFlavor = "AK4PFPuppi",
                           postfix = "Puppi"
                           )
process.puppiNoLep.useExistingWeights = False
process.puppi.useExistingWeights = False

process.METUncertainties = cms.EDAnalyzer("METUncertainties",
                                          metSrc = cms.untracked.InputTag("slimmedMETs"),
                                          metmodifiedSrc = cms.untracked.InputTag("slimmedMETsModified"),
                                          metPuppiSrc = cms.untracked.InputTag("slimmedMETsPuppi")
                                          )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("METUncertainties_ntuples.root")
                                   )

process.p = cms.EndPath(
    process.egmPhotonIDSequence*
    process.puppiMETSequence*
    process.fullPatMetSequencePuppi*
    process.fullPatMetSequenceModified*
    process.METUncertainties
    )
