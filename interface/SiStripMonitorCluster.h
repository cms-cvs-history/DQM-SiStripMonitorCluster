#ifndef SiStripMonitorCluster_SiStripMonitorCluster_h
#define SiStripMonitorCluster_SiStripMonitorCluster_h
// -*- C++ -*-
// Package:     SiStripMonitorCluster
// Class  :     SiStripMonitorCluster
/**\class SiStripMonitorCluster SiStripMonitorCluster.h DQM/SiStripMonitorCluster/interface/SiStripMonitorCluster.h
   Data Quality Monitoring source of the Silicon Strip Tracker. Produces histograms related to clusters.
*/
// Original Author:  dkcira
//         Created:  Wed Feb  1 16:47:14 CET 2006
// $Id: SiStripMonitorCluster.h,v 1.16 2008/04/19 20:13:06 dutta Exp $
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/MonitorElement.h"


#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"

#include "DQM/SiStripCommon/interface/SiStripFolderOrganizer.h"
#include "AnalysisDataFormats/SiStripClusterInfo/interface/SiStripClusterInfo.h"
#include <vector>

#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/Common/interface/DetSetNew.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/DetSetVector.h"

class DQMStore;

class SiStripMonitorCluster : public edm::EDAnalyzer {
 public:
  explicit SiStripMonitorCluster(const edm::ParameterSet&);
  ~SiStripMonitorCluster();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void beginJob(edm::EventSetup const&) ;
  virtual void endJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);
  struct ModMEs{ // MEs for one single detector module
    ModMEs():    
      LayernClusters(0),
	 LayernClustersTrend(0),
	 LayerClusterStoN(0),
	 LayerClusterStoNTrend(0),
	 LayerClusterCharge(0),
	 LayerClusterChargeTrend(0),
	 LayerClusterNoise(0),
	 LayerClusterNoiseTrend(0),
	 LayerClusterWidth(0),
	 LayerClusterWidthTrend(0),
	 LayerClusterPos(0),
	 LayerClusterPGV(0){};

    // MEs for one single detector module
    // MEs at subdetector level
    MonitorElement* LayernClusters;
    MonitorElement* LayernClustersTrend;
    MonitorElement* LayerClusterStoN;
    MonitorElement* LayerClusterStoNTrend;
    MonitorElement* LayerClusterCharge;
    MonitorElement* LayerClusterChargeTrend;
    MonitorElement* LayerClusterNoise;
    MonitorElement* LayerClusterNoiseTrend;
    MonitorElement* LayerClusterWidth;
    MonitorElement* LayerClusterWidthTrend;
    MonitorElement* LayerClusterPos;
    MonitorElement* LayerClusterPGV;


    MonitorElement* NumberOfClusters;
    MonitorElement* NumberOfClustersTrend;
    MonitorElement* ClusterPosition;
    MonitorElement* ClusterPositionTrend;
    MonitorElement* ClusterWidth;
    MonitorElement* ClusterWidthTrend;
    MonitorElement* ClusterCharge;
    MonitorElement* ClusterChargeTrend;
    MonitorElement* ClusterNoise;
    MonitorElement* ClusterNoiseTrend;
    MonitorElement* ClusterSignalOverNoise;
    MonitorElement* ClusterSignalOverNoiseTrend;
    MonitorElement* ModuleLocalOccupancy;
    MonitorElement* ModuleLocalOccupancyTrend;
    MonitorElement* NrOfClusterizedStrips; // can be used at client level for occupancy calculations
    MonitorElement* NrOfClusterizedStripsTrend; // can be used at client level for occupancy calculations
  };

 private:
  void ResetModuleMEs(uint32_t idet);
  void createMEs(const edm::EventSetup& es);
  void bookSubDetMEs(TString name,TString flag);
  void AllClusters( const edm::EventSetup& es);
/*   bool clusterInfos(SiStripClusterInfo* cluster, const uint32_t& detid,std::string flag, const LocalVector LV); */
/*   void fillTrendMEs(SiStripClusterInfo* cluster,std::string name,float cos, std::string flag); */
  bool clusterInfos(SiStripClusterInfo* cluster, const uint32_t& detid);
  void fillTrendMEs(SiStripClusterInfo* cluster,std::string name);
  void fillTrend(MonitorElement* me ,float value);
  void bookTrendMEs(TString name,int32_t layer,uint32_t id,std::string flag);
  void book(); 
  inline void fillME(MonitorElement* ME,float value1){if (ME!=0)ME->Fill(value1);}
  inline void fillME(MonitorElement* ME,float value1,float value2){if (ME!=0)ME->Fill(value1,value2);}
  inline void fillME(MonitorElement* ME,float value1,float value2,float value3){if (ME!=0)ME->Fill(value1,value2,value3);}
  inline void fillME(MonitorElement* ME,float value1,float value2,float value3,float value4){if (ME!=0)ME->Fill(value1,value2,value3,value4);}
  MonitorElement * bookMETrend(const char*, const char*);
  MonitorElement* bookME1D(const char* ParameterSetLabel, const char* HistoName);

 private:
  DQMStore* dqmStore_;
  edm::ParameterSet conf_;
  std::map<uint32_t, ModMEs> ClusterMEs;
  // flags
  bool show_mechanical_structure_view, show_readout_view, show_control_view, select_all_detectors, reset_each_run, fill_signal_noise;
  unsigned long long m_cacheID_;

  TString name;
  LocalVector LV;

  edm::ESHandle<SiStripDetCabling> SiStripDetCabling_;
  std::vector<uint32_t> ModulesToBeExcluded_;
  SiStripFolderOrganizer folder_organizer;

  std::map<TString, ModMEs> ModMEsMap;
  std::map<TString, MonitorElement*> MEMap;

  edm::ParameterSet Parameters;

  edm::Handle< edm::DetSetVector<SiStripCluster> > cluster_detsetvektor;
  std::map<std::pair<std::string,int32_t>,bool> DetectedLayers;

  int runNb, eventNb;
  int firstEvent;
  int count, NClus[4];

  bool layerswitchncluson;
  bool layerswitchcluschargeon;
  bool layerswitchclusstonon;
  bool layerswitchclusposon;
  bool layerswitchclusnoiseon;
  bool layerswitchcluswidthon;

  bool moduleswitchncluson;
  bool moduleswitchcluschargeon;
  bool moduleswitchclusstonon;
  bool moduleswitchclusposon;
  bool moduleswitchclusnoiseon;
  bool moduleswitchcluswidthon;

  bool tibon;
  bool tidon;
  bool tobon;
  bool tecon;

};
#endif
