import FWCore.ParameterSet.Config as cms

process = cms.Process("SignificanceTest")

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuples')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(10)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'/store/mc/RunIISpring18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/100X_upgrade2018_realistic_v10-v1/70000/FEA6B290-5424-E811-B6F2-008CFAC93E5C.root'
        'file:patMiniAOD_standard.root'
        )
                            )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag( process.GlobalTag, "94X_mc2017_realistic_v10")

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(
    process,
    isData = False,
    postfix = "Modified"
    )


process.SignificanceTest = cms.EDAnalyzer("SignificanceTest",
                                          metSrc = cms.untracked.InputTag("slimmedMETs"),
                                          metmodifiedSrc = cms.untracked.InputTag("slimmedMETsModified")
                                          )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("Significance_test_MET.root")
                                   )

process.p = cms.EndPath(
    process.fullPatMetSequenceModified*
    process.SignificanceTest
    )
