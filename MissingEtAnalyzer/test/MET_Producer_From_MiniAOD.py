import FWCore.ParameterSet.Config as cms

process = cms.Process("RERUN")

import PhysicsTools.PatAlgos.tools.helpers as configtools
patAlgosToolsTask = configtools.getPatAlgosToolsTask(process)

process.load('Configuration.StandardSequences.Services_cff')
patAlgosToolsTask.add(process.randomEngineStateProducer)
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = [ 'cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['run2_mc']

fname = '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/483975CE-8142-E811-A72D-001E67E71E20.root'

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring([fname])
)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(process,
                           isData = False,
                           )

from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True );

runMetCorAndUncFromMiniAOD(process,
                           isData = False,
                           metType = "Puppi",
                           jetFlavor = "AK4PFPuppi",
                           postfix = "Puppi"
                           )

runMetCorAndUncFromMiniAOD(process,
                           isData = False,
                           metType = "Puppi",
                           pfCandColl = cms.InputTag("puppiForMET"),
                           recoMetFromPFCs = True,
                           reclusterJets = True,
                           jetFlavor = "AK4PFPuppi",
                           postfix = "PuppiRecorrect",
                           )

process.puppiNoLep.useExistingWeights = False
process.puppi.useExistingWeights = False

process.MINIAODSIMoutput = cms.OutputModule(
    "PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = cms.untracked.vstring(
        "keep *_slimmedMETs_*_*",
        "keep *_patPFMet_*_*",
        "keep *_patPFMetT1*_*_*",
        "keep *_slimmedMETsPuppi*_*_*",
        ),
    fileName = cms.untracked.string('corrected_MET_MINIAOD.root'),
    dataset  = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier   = cms.untracked.string('')
        ),
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning  = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
    )

process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput, patAlgosToolsTask)
