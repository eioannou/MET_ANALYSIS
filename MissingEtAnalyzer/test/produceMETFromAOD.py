from PhysicsTools.PatAlgos.patTemplate_cfg import cms, process, patAlgosToolsTask

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
patAlgosToolsTask.add( process.patCandidatesTask)

process.patTaus.skipMissingTauID = True

process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
patAlgosToolsTask.add(process.selectedPatCandidatesTask)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
patAlgosToolsTask.add( process.inclusiveVertexingTask )
patAlgosToolsTask.add( process.inclusiveCandidateVertexingTask )
patAlgosToolsTask.add( process.inclusiveCandidateVertexingCvsLTask )

process.load("PhysicsTools.PatAlgos.slimming.slimming_cff")
patAlgosToolsTask.add( process.slimmingTask )

from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeCommon, miniAOD_customizeMC
miniAOD_customizeCommon(process)
miniAOD_customizeMC(process)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'/store/mc/RunIISpring18DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/100X_upgrade2018_realistic_v10-v1/70000/FEDCEFC1-8224-E811-9B1B-A0369F8363BE.root'
        '/store/mc/RunIIFall17DRPremix/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/AODSIM/94X_mc2017_realistic_v10-v1/50000/26338CBF-BAE8-E711-AE63-FA163E022194.root'
        )
                            )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.out.outputCommands = process.MicroEventContentMC.outputCommands
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeOutput
miniAOD_customizeOutput(process.out)

process.out.fileName = 'patMiniAOD_standard.root'
