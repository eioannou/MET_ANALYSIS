import FWCore.ParameterSet.Config as cms

process = cms.Process("jectxt")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# define your favourite global tag
process.GlobalTag.globaltag = cms.string('94X_mc2017_realistic_v10')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.source = cms.Source("EmptySource")

process.readAK4PFchs = cms.EDAnalyzer('JetCorrectorDBReader',
                                      #below is the communication to the database
                                      payloadName = cms.untracked.string('AK4PFchs'),
                                      globalTag = cms.untracked.string('94X_mc2017_realistic_v10'),
                                      printScreen = cms.untracked.bool(False),
                                      createTextFile = cms.untracked.bool(True)
                                      )

process.readAK8PFchs = process.readAK4PFchs.clone( payloadName = 'AK8PFchs')
process.readAK4PF    = process.readAK4PFchs.clone( payloadName = 'AK4PF')

process.p = cms.Path(
    process.readAK4PFchs *
    process.readAK8PFchs *
    process.readkAK4PF
)
