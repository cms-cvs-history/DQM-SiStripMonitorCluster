import FWCore.ParameterSet.Config as cms

process = cms.Process("DQMOnlineRealData")
#--------------------------
# Event Source
#--------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        '/store/data/Commissioning08/Cosmics/RAW/v1/000/067/647/0000721C-35A3-DD11-9132-001D09F291D7.root'
        '/store/data/CRAFT09/Cosmics/RAW/v1/000/110/998/001404E1-0F8A-DE11-ADB3-000423D99EEE.root'
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

#-------------------------------------------------
# Message Logger
#-------------------------------------------------
process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('siStripDigis',
                                         'siStripZeroSuppression',
                                         'siStripClusters'
                                         'SiStripMonitorCluster'),
    cout = cms.untracked.PSet(threshold = cms.untracked.string('ERROR')),
    destinations = cms.untracked.vstring('cout')
)

#-------------------------------------------------
# Geometry
#-------------------------------------------------
process.load("Configuration.StandardSequences.Geometry_cff")

#-------------------------------------------------
# Calibration
#-------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR09_31X_V6P::All"
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

#-----------------------
#  Reconstruction Modules
#-----------------------
process.load("EventFilter.SiStripRawToDigi.SiStripDigis_cfi")
process.siStripDigis.ProductLabel = 'source'

process.load("RecoLocalTracker.Configuration.RecoLocalTracker_Cosmics_cff")

#--------------------------
# DQM Services
#--------------------------
process.DQMStore = cms.Service("DQMStore",
    referenceFileName = cms.untracked.string(''),
    verbose = cms.untracked.int32(0)
)

process.TkDetMap = cms.Service("TkDetMap")
process.SiStripDetInfoFileReader = cms.Service("SiStripDetInfoFileReader")

#--------------------------
# SiStrip MonitorCluster
#--------------------------
process.load("DQM.SiStripMonitorCluster.SiStripMonitorCluster_cfi")
process.SiStripMonitorCluster.CreateTrendMEs = True
process.SiStripMonitorCluster.TkHistoMap_On = True
process.SiStripMonitorCluster.SelectAllDetectors = True
process.SiStripMonitorCluster.TProfTotalNumberOfClusters.subdetswitchon = True
process.SiStripMonitorCluster.TH1TotalNumberOfClusters.subdetswitchon = True
process.SiStripMonitorCluster.TProfClustersApvCycle.subdetswitchon = True
process.SiStripMonitorCluster.TH2ClustersApvCycle.subdetswitchon = True
process.SiStripMonitorCluster.TH2ClustersApvCycle.yfactor = 0.005
process.SiStripMonitorCluster.OutputMEsInRootFile = True
process.SiStripMonitorCluster.OutputFileName = 'SiStripMonitorCluster.root'


process.outP = cms.OutputModule("AsciiOutputModule")

process.AdaptorConfig = cms.Service("AdaptorConfig")

process.RecoForDQM = cms.Sequence(process.siStripDigis*process.siStripZeroSuppression*process.siStripClusters)
process.p = cms.Path(process.RecoForDQM*process.SiStripMonitorCluster)
process.ep = cms.EndPath(process.outP)

