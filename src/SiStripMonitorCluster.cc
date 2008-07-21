// -*- C++ -*-
// Package:    SiStripMonitorCluster
// Class:      SiStripMonitorCluster
/**\class SiStripMonitorCluster SiStripMonitorCluster.cc DQM/SiStripMonitorCluster/src/SiStripMonitorCluster.cc
 */
// Original Author:  Dorian Kcira
//         Created:  Wed Feb  1 16:42:34 CET 2006
//<<<<<<< SiStripMonitorCluster.cc
// $Id: SiStripMonitorCluster.cc,v 1.40 2008/04/29 15:00:43 dutta Exp $
//=======
// $Id: SiStripMonitorCluster.cc,v 1.37.2.1 2008/04/29 16:26:38 dutta Exp $
//>>>>>>> 1.37.2.1
#include <vector>
#include <numeric>
#include <fstream>
#include <math.h>
#include "TNamed.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CondFormats/DataRecord/interface/SiStripNoisesRcd.h"
#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
#include "CalibTracker/Records/interface/SiStripGainRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripGain.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DQM/SiStripCommon/interface/SiStripFolderOrganizer.h"
#include "DQM/SiStripCommon/interface/SiStripHistoId.h"
#include "DQM/SiStripMonitorCluster/interface/SiStripMonitorCluster.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "TMath.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"



//--------------------------------------------------------------------------------------------
SiStripMonitorCluster::SiStripMonitorCluster(const edm::ParameterSet& iConfig) : dqmStore_(edm::Service<DQMStore>().operator->()), conf_(iConfig), show_mechanical_structure_view(true), show_readout_view(false), show_control_view(false), select_all_detectors(false), reset_each_run(false), m_cacheID_(0), fill_signal_noise (true) , folder_organizer()
{
  for(int i=0;i<4;++i) NClus[i]=0;
  firstEvent = -1;
  eventNb = 0;


  //get on/off option for every cluster from cfi
  edm::ParameterSet ParametersnClusters =  conf_.getParameter<edm::ParameterSet>("TH1nClusters");
  layerswitchncluson = ParametersnClusters.getParameter<bool>("layerswitchon");
  std::cout<<"layerswitchncluson "<<layerswitchncluson<<std::endl;
  moduleswitchncluson = ParametersnClusters.getParameter<bool>("moduleswitchon");
  std::cout<<"moduleswitchncluson "<<moduleswitchncluson<<std::endl;
  
  edm::ParameterSet ParametersClusterCharge =  conf_.getParameter<edm::ParameterSet>("TH1ClusterCharge");
  layerswitchcluschargeon = ParametersClusterCharge.getParameter<bool>("layerswitchon");
  std::cout<<"layerswitchcluschargeon "<<layerswitchcluschargeon<<std::endl;
  moduleswitchcluschargeon = ParametersClusterCharge.getParameter<bool>("moduleswitchon");
  std::cout<<"moduleswitchcluschargeon "<<moduleswitchcluschargeon<<std::endl;
  
  edm::ParameterSet ParametersClusterStoN =  conf_.getParameter<edm::ParameterSet>("TH1ClusterStoN");
  layerswitchclusstonon = ParametersClusterStoN.getParameter<bool>("layerswitchon");
  std::cout<<"layerswitchclusstonon "<<layerswitchclusstonon<<std::endl;
  moduleswitchclusstonon = ParametersClusterStoN.getParameter<bool>("moduleswitchon");
  std::cout<<"moduleswitchclusstonon "<<moduleswitchclusstonon<<std::endl;
  
  edm::ParameterSet ParametersClusterPos =  conf_.getParameter<edm::ParameterSet>("TH1ClusterPos");
  layerswitchclusposon = ParametersClusterPos.getParameter<bool>("layerswitchon");
  std::cout<<"layerswitchclusposon "<<layerswitchclusposon<<std::endl;
  moduleswitchclusposon = ParametersClusterPos.getParameter<bool>("moduleswitchon");
  std::cout<<"moduleswitchclusposon "<<moduleswitchclusposon<<std::endl;
  
  edm::ParameterSet ParametersClusterNoise =  conf_.getParameter<edm::ParameterSet>("TH1ClusterNoise");
  layerswitchclusnoiseon = ParametersClusterNoise.getParameter<bool>("layerswitchon");
  std::cout<<"layerswitchclusnoiseon "<<layerswitchclusnoiseon<<std::endl;
  moduleswitchclusnoiseon = ParametersClusterNoise.getParameter<bool>("moduleswitchon");
  std::cout<<"moduleswitchclusnoiseon "<<moduleswitchclusnoiseon<<std::endl;
  
  edm::ParameterSet ParametersClusterWidth =  conf_.getParameter<edm::ParameterSet>("TH1ClusterWidth");
  layerswitchcluswidthon = ParametersClusterWidth.getParameter<bool>("layerswitchon");
  std::cout<<"layerswitchcluswidthon "<<layerswitchcluswidthon<<std::endl;
  moduleswitchcluswidthon = ParametersClusterWidth.getParameter<bool>("moduleswitchon");
  std::cout<<"moduleswitchcluswidthon "<<moduleswitchcluswidthon<<std::endl;
  

  edm::ParameterSet ParametersDetsOn =  conf_.getParameter<edm::ParameterSet>("detectorson");
  tibon = ParametersDetsOn.getParameter<bool>("tibon");
  tidon = ParametersDetsOn.getParameter<bool>("tidon");
  tobon = ParametersDetsOn.getParameter<bool>("tobon");
  tecon = ParametersDetsOn.getParameter<bool>("tecon");

} 


SiStripMonitorCluster::~SiStripMonitorCluster() { }

//--------------------------------------------------------------------------------------------
void SiStripMonitorCluster::beginRun(const edm::Run& run, const edm::EventSetup& es){

  if (show_mechanical_structure_view) {
    unsigned long long cacheID = es.get<SiStripDetCablingRcd>().cacheIdentifier();
    if (m_cacheID_ != cacheID) {
      m_cacheID_ = cacheID;       
      edm::LogInfo("SiStripMonitorCluster") <<"SiStripMonitorCluster::beginRun: " 
					    << " Creating MEs for new Cabling ";     

      createMEs(es);
    } 
  } else if (reset_each_run) {
    edm::LogInfo("SiStripMonitorCluster") <<"SiStripMonitorCluster::beginRun: " 
					  << " Resetting MEs ";        
    for (std::map<uint32_t, ModMEs >::const_iterator idet = ClusterMEs.begin() ; idet!=ClusterMEs.end() ; idet++) {
      ResetModuleMEs(idet->first);
    }
  }

  es.get<SiStripDetCablingRcd>().get( SiStripDetCabling_ );
  book();


}

//--------------------------------------------------------------------------------------------
void SiStripMonitorCluster::endRun(const edm::Run&, const edm::EventSetup&){
}




//--------------------------------------------------------------------------------------------
void SiStripMonitorCluster::createMEs(const edm::EventSetup& es){
  if ( show_mechanical_structure_view ){
    // take from eventSetup the SiStripDetCabling object - here will use SiStripDetControl later on
    edm::ESHandle<SiStripDetCabling> tkmechstruct;
    es.get<SiStripDetCablingRcd>().get(tkmechstruct);
    
    // get list of active detectors from SiStripDetCabling - this will change and be taken from a SiStripDetControl object
    std::vector<uint32_t> activeDets;
    activeDets.clear(); // just in case
    tkmechstruct->addActiveDetectorsRawIds(activeDets);
    
    std::vector<uint32_t> SelectedDetIds;
    if(select_all_detectors){
      // select all detectors if appropriate flag is set,  for example for the mtcc
      SelectedDetIds = activeDets;
    }else{
      // use SiStripSubStructure for selecting certain regions
      SiStripSubStructure substructure;
      
      if(tibon) substructure.getTIBDetectors(activeDets, SelectedDetIds, 0, 0, 0, 0); // this adds rawDetIds to SelectedDetIds
      if(tobon) substructure.getTOBDetectors(activeDets, SelectedDetIds, 0, 0, 0);    // this adds rawDetIds to SelectedDetIds
      if(tidon) substructure.getTIDDetectors(activeDets, SelectedDetIds, 0, 0, 0, 0); // this adds rawDetIds to SelectedDetIds
      if(tecon) substructure.getTECDetectors(activeDets, SelectedDetIds, 0, 0, 0, 0, 0, 0); // this adds rawDetIds to SelectedDetIds
      
    }
    
    // remove any eventual zero elements - there should be none, but just in case
    for(std::vector<uint32_t>::iterator idets = SelectedDetIds.begin(); idets != SelectedDetIds.end(); idets++){
      if(*idets == 0) SelectedDetIds.erase(idets);
    }
    
    // use SistripHistoId for producing histogram id (and title)
    SiStripHistoId hidmanager;
    // create SiStripFolderOrganizer
    //already done in class
    //SiStripFolderOrganizer folder_organizer;
    
    folder_organizer.setSiStripFolder();

    // loop over detectors and book MEs
    edm::LogInfo("SiStripTkDQM|SiStripMonitorCluster")<<"nr. of SelectedDetIds:  "<<SelectedDetIds.size();
    for(std::vector<uint32_t>::const_iterator detid_iterator = SelectedDetIds.begin(); detid_iterator!=SelectedDetIds.end(); detid_iterator++){
      ModMEs modSingle;
      std::string hid;
      // set appropriate folder using SiStripFolderOrganizer
      folder_organizer.setDetectorFolder(*detid_iterator); // pass the detid to this method
      if (reset_each_run) ResetModuleMEs(*detid_iterator);

      //nr. of clusters per module
      if(moduleswitchncluson) {
	hid = hidmanager.createHistoId("NumberOfClusters","det",*detid_iterator);
	modSingle.NumberOfClusters = dqmStore_->book1D(hid, hid, 5,-0.5,4.5); dqmStore_->tag(modSingle.NumberOfClusters, *detid_iterator);
	modSingle.NumberOfClusters->setAxisTitle("number of clusters in one detector module");
      }

      //ClusterPosition
      if(moduleswitchclusposon) {
	hid = hidmanager.createHistoId("ClusterPosition","det",*detid_iterator);
	modSingle.ClusterPosition = dqmStore_->book1D(hid, hid, 768,-0.5,767.5); dqmStore_->tag(modSingle.ClusterPosition, *detid_iterator); // 6 APVs -> 768 strips
	modSingle.ClusterPosition->setAxisTitle("cluster position [strip number +0.5]");
      }

      //ClusterWidth
      if(moduleswitchcluswidthon) {
	hid = hidmanager.createHistoId("ClusterWidth","det",*detid_iterator);
	modSingle.ClusterWidth = dqmStore_->book1D(hid, hid, 20,-0.5,19.5); dqmStore_->tag(modSingle.ClusterWidth, *detid_iterator);
	modSingle.ClusterWidth->setAxisTitle("cluster width [nr strips]");
      }

      //ClusterCharge
      if(moduleswitchcluschargeon) {
	hid = hidmanager.createHistoId("ClusterCharge","det",*detid_iterator);
	modSingle.ClusterCharge = dqmStore_->book1D(hid, hid, 100,0.,500.); dqmStore_->tag(modSingle.ClusterCharge, *detid_iterator);
	modSingle.ClusterCharge->setAxisTitle("cluster charge [ADC]");
      }

      //ClusterNoise
      if(moduleswitchclusnoiseon) {
	hid = hidmanager.createHistoId("ClusterNoise","det",*detid_iterator);
	modSingle.ClusterNoise = dqmStore_->book1D(hid, hid, 20,0.,10.); dqmStore_->tag(modSingle.ClusterNoise, *detid_iterator);
	modSingle.ClusterNoise->setAxisTitle("cluster noise");
      }

      //ClusterSignalOverNoise
      if(moduleswitchclusstonon) {
	hid = hidmanager.createHistoId("ClusterSignalOverNoise","det",*detid_iterator);
	modSingle.ClusterSignalOverNoise = dqmStore_->book1D(hid, hid, 60,0.,200.); dqmStore_->tag(modSingle.ClusterSignalOverNoise, *detid_iterator);
	modSingle.ClusterSignalOverNoise->setAxisTitle("ratio of signal to noise for each cluster");
      }

      //ModuleLocalOccupancy
      hid = hidmanager.createHistoId("ModuleLocalOccupancy","det",*detid_iterator);
      // occupancy goes from 0 to 1, probably not over some limit value (here 0.1)
      modSingle.ModuleLocalOccupancy = dqmStore_->book1D(hid, hid, 20,-0.005,0.05); dqmStore_->tag(modSingle.ModuleLocalOccupancy, *detid_iterator);
      modSingle.ModuleLocalOccupancy->setAxisTitle("module local occupancy [% of clusterized strips]");

      //NrOfClusterizedStrips
      hid = hidmanager.createHistoId("NrOfClusterizedStrips","det",*detid_iterator);
      modSingle.NrOfClusterizedStrips = dqmStore_->book1D(hid, hid, 10,-0.5,9.5); dqmStore_->tag(modSingle.NrOfClusterizedStrips, *detid_iterator);
      modSingle.NrOfClusterizedStrips->setAxisTitle("number of clusterized strips");


      // append to ClusterMEs
      ClusterMEs.insert( std::make_pair(*detid_iterator, modSingle));

    }//end of loop over detectors

  }//end of if

}//end of method



//--------------------------------------------------------------------------------------------
void SiStripMonitorCluster::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;

  runNb   = iEvent.id().run();
  //   eventNb = iEvent.id().event();
  eventNb++;
  //std::cout << " run " << iEvent.id().run() << runNb << " event " << iEvent.id().event() << eventNb << std::endl;
  count = 0;

  edm::ESHandle<SiStripNoises> noiseHandle;
  iSetup.get<SiStripNoisesRcd>().get(noiseHandle);

  edm::ESHandle<SiStripGain> gainHandle;
  iSetup.get<SiStripGainRcd>().get(gainHandle);

  edm::ESHandle<SiStripQuality> qualityHandle;
  iSetup.get<SiStripQualityRcd>().get(qualityHandle);


  // retrieve producer name of input StripClusterCollection
  std::string clusterProducer = conf_.getParameter<std::string>("ClusterProducer");
  // get collection of DetSetVector of clusters from Event
  //   edm::Handle< edmNew::DetSetVector<SiStripCluster> > cluster_detsetvektor;
  //edm::Handle< edm::DetSetVector<SiStripCluster> > cluster_detsetvektor;
  iEvent.getByLabel(clusterProducer, cluster_detsetvektor);

  //if (!cluster_detsetvektor.isValid()) std::cout<<" collection not valid"<<std::endl;
  if (!cluster_detsetvektor.isValid()) return;

  // loop over MEs. Mechanical structure view. No need for condition here. If map is empty, nothing should happen.
  for (std::map<uint32_t, ModMEs>::const_iterator iterMEs = ClusterMEs.begin() ; iterMEs!=ClusterMEs.end() ; iterMEs++) 
    {
      uint32_t detid = iterMEs->first;  ModMEs modSingle = iterMEs->second;
      // get from DetSetVector the DetSet of clusters belonging to one detid - first make sure there exists clusters with this id

      edm::DetSetVector<SiStripCluster>::const_iterator isearch = cluster_detsetvektor->find(detid); // search  clusters of detid

      
      if(isearch==cluster_detsetvektor->end()){
	if(moduleswitchncluson) {
	  if(modSingle.NumberOfClusters != NULL){
	    (modSingle.NumberOfClusters)->Fill(0.,1.); // no clusters for this detector module, so fill histogram with 0
	  }
	}
	continue; // no clusters for this detid => jump to next step of loop
      }

      //cluster_detset is a structure, cluster_detset.data is a std::vector<SiStripCluster>, cluster_detset.id is uint32_t
      edm::DetSet<SiStripCluster> cluster_detset = (*cluster_detsetvektor)[detid]; // the statement above makes sure there exists an element with 'detid'

      if(moduleswitchncluson) {
	if(modSingle.NumberOfClusters != NULL){ // nr. of clusters per module
	  (modSingle.NumberOfClusters)->Fill(static_cast<float>(cluster_detset.size()),1.);
	}
      }

      short total_clusterized_strips = 0;
      //
      float clusterNoise = 0.;
      float clusterNoise2 = 0;
      int nrnonzeroamplitudes = 0;

      SiStripNoises::Range detNoiseRange = noiseHandle->getRange(detid);
      SiStripApvGain::Range detGainRange =  gainHandle->getRange(detid); 
      SiStripQuality::Range qualityRange = qualityHandle->getRange(detid);

      for(edm::DetSet<SiStripCluster>::const_iterator clusterIter = cluster_detset.begin(); clusterIter!= cluster_detset.end(); clusterIter++){
      
	if(moduleswitchclusposon) {
	  if(modSingle.ClusterPosition != NULL){ // position of cluster
	    (modSingle.ClusterPosition)->Fill(clusterIter->barycenter(),1.);
	  }
	}
    
	const std::vector<uint16_t>& ampls = clusterIter->amplitudes();
	short local_size = ampls.size(); // width defined as nr. of strips that belong to cluster
	total_clusterized_strips = total_clusterized_strips + local_size; // add nr of strips of this cluster to total nr. of clusterized strips
	
	if(moduleswitchcluswidthon) {
	  if(modSingle.ClusterWidth != NULL){ // width of cluster, calculate yourself, no method for getting it
	    (modSingle.ClusterWidth)->Fill(static_cast<float>(local_size),1.);
	  }
	}


	if(moduleswitchclusstonon) {
	  if( fill_signal_noise && modSingle.ClusterSignalOverNoise){
	    const std::vector<uint16_t>& ampls = clusterIter->amplitudes();
	    float clusterSignal = 0;
	    float noise;
	    for(uint iamp=0; iamp<ampls.size(); iamp++){
	      if(ampls[iamp]>0){ // nonzero amplitude
		clusterSignal += ampls[iamp];
		try{
		  if(!qualityHandle->IsStripBad(qualityRange, clusterIter->firstStrip()+iamp)){
		    noise = noiseHandle->getNoise(clusterIter->firstStrip()+iamp,detNoiseRange)/gainHandle->getStripGain(clusterIter->firstStrip()+iamp, detGainRange);
		  }
		}
		catch(cms::Exception& e){
		  edm::LogError("SiStripTkDQM|SiStripMonitorCluster|DB")<<" cms::Exception:  detid="<<detid<<" firstStrip="<<clusterIter->firstStrip()<<" iamp="<<iamp<<e.what();
		}
		clusterNoise2 += noise*noise;
		nrnonzeroamplitudes++;
	      }
	    }
	    clusterNoise = sqrt(clusterNoise2/nrnonzeroamplitudes);
	    if(modSingle.ClusterNoise) (modSingle.ClusterNoise)->Fill(clusterNoise,1.);
	    if(modSingle.ClusterSignalOverNoise) (modSingle.ClusterSignalOverNoise)->Fill(clusterSignal/clusterNoise,1.);
	  }
	}
	  
	//
	if(moduleswitchcluschargeon) {
	  if(modSingle.ClusterCharge != NULL){ // charge of cluster
	    const std::vector<uint16_t>& ampls = clusterIter->amplitudes();
	    short local_charge = 0;
	    for(std::vector<uint16_t>::const_iterator iampls = ampls.begin(); iampls<ampls.end(); iampls++){
	      local_charge += *iampls;
	    }
	    (modSingle.ClusterCharge)->Fill(static_cast<float>(local_charge),1.);
	  }
	}
	
	  
      } // end loop on clusters for the given detid
	
      if(modSingle.NrOfClusterizedStrips != NULL){ // nr of clusterized strips
	modSingle.NrOfClusterizedStrips->Fill(static_cast<float>(total_clusterized_strips),1.);
      }
	
      short total_nr_strips = 6 * 128; // assume 6 APVs per detector for the moment. later ask FedCabling object
      float local_occupancy = static_cast<float>(total_clusterized_strips)/static_cast<float>(total_nr_strips);
      if(modSingle.ModuleLocalOccupancy != NULL){ // nr of clusterized strips
	modSingle.ModuleLocalOccupancy->Fill(local_occupancy,1.);
      }

    }//end of loop over MEs
  
      
  //Now fill layer histos
  AllClusters(iSetup);
  

  //Summary Counts of clusters
  std::map<TString, MonitorElement*>::iterator iME;
  std::map<TString, ModMEs>::iterator          iModME ;
  
  int nTot=0;
  for (int i=0;i<4;++i){ // loop over TIB, TID, TOB, TEC
    
    //     std::cout<< " NClus[i] " << NClus[i] << std::endl;
    //     nTot+=NClus[i];
    //     std::cout<< " nTot " << nTot << std::endl;
    //     NClus[i]=0;
	
  }// loop over TIB, TID, TOB, TEC
      
      
  name="NumberOfClusters";
  iME = MEMap.find(name);
  if(iME!=MEMap.end() && layerswitchncluson) iME->second->Fill(nTot);
  
  iME = MEMap.find(name+"Trend");
  if(iME!=MEMap.end() && layerswitchncluson) fillTrend(iME->second,nTot);
  
  
}//end of method




//--------------------------------------------------------------------------------------------
void SiStripMonitorCluster::endJob(void){
  bool outputMEsInRootFile = conf_.getParameter<bool>("OutputMEsInRootFile");
  std::string outputFileName = conf_.getParameter<std::string>("OutputFileName");
  if(outputMEsInRootFile){

    std::ofstream monitor_summary("monitor_cluster_summary.txt");
    monitor_summary<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    monitor_summary<<"SiStripMonitorCluster::endJob ClusterMEs.size()="<<ClusterMEs.size()<<std::endl;

    for(std::map<uint32_t, ModMEs>::const_iterator idet = ClusterMEs.begin(); idet!= ClusterMEs.end(); idet++ ){

      monitor_summary<<"SiStripTkDQM|SiStripMonitorCluster"<<"      ++++++detid  "<<idet->first<<std::endl<<std::endl;

      if(moduleswitchncluson) {
	monitor_summary<<"SiStripTkDQM|SiStripMonitorCluster"<<"              +++ NumberOfClusters "<<(idet->second).NumberOfClusters->getEntries()<<" "<<(idet->second).NumberOfClusters->getMean()<<" "<<(idet->second).NumberOfClusters->getRMS()<<std::endl;
      }

      if(moduleswitchclusposon) {
	monitor_summary<<"SiStripTkDQM|SiStripMonitorCluster"<<"              +++ ClusterPosition "<<(idet->second).ClusterPosition->getEntries()<<" "<<(idet->second).ClusterPosition->getMean()<<" "<<(idet->second).ClusterPosition->getRMS()<<std::endl;
      }

      if(moduleswitchcluswidthon) {
	monitor_summary<<"SiStripTkDQM|SiStripMonitorCluster"<<"              +++ ClusterWidth "<<(idet->second).ClusterWidth->getEntries()<<" "<<(idet->second).ClusterWidth->getMean()<<" "<<(idet->second).ClusterWidth->getRMS()<<std::endl;
      }

      if(moduleswitchcluschargeon) {
	monitor_summary<<"SiStripTkDQM|SiStripMonitorCluster"<<"              +++ ClusterCharge "<<(idet->second).ClusterCharge->getEntries()<<" "<<(idet->second).ClusterCharge->getMean()<<" "<<(idet->second).ClusterCharge->getRMS()<<std::endl;
      }

      if(moduleswitchclusnoiseon) {
	monitor_summary<<"SiStripTkDQM|SiStripMonitorCluster"<<"              +++ ClusterNoise "<<(idet->second).ClusterNoise->getEntries()<<" "<<(idet->second).ClusterNoise->getMean()<<" "<<(idet->second).ClusterNoise->getRMS()<<std::endl;
      }

      if(moduleswitchclusstonon) {
	monitor_summary<<"SiStripTkDQM|SiStripMonitorCluster"<<"              +++ ClusterSignalOverNoise "<<(idet->second).ClusterSignalOverNoise->getEntries()<<" "<<(idet->second).ClusterSignalOverNoise->getMean()<<" "<<(idet->second).ClusterSignalOverNoise->getRMS()<<std::endl;
      }

      monitor_summary<<"SiStripTkDQM|SiStripMonitorCluster"<<"              +++ ModuleLocalOccupancy "<<(idet->second).ModuleLocalOccupancy->getEntries()<<" "<<(idet->second).ModuleLocalOccupancy->getMean()<<" "<<(idet->second).ModuleLocalOccupancy->getRMS()<<std::endl;

      monitor_summary<<"SiStripTkDQM|SiStripMonitorCluster"<<"              +++ NrOfClusterizedStrips "<<(idet->second).NrOfClusterizedStrips->getEntries()<<" "<<(idet->second).NrOfClusterizedStrips->getMean()<<" "<<(idet->second).NrOfClusterizedStrips->getRMS()<<std::endl;
 
    }

    monitor_summary<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  }

  // save histos in a file
  dqmStore_->save(outputFileName);
  
}




//--------------------------------------------------------------------------------------------
void SiStripMonitorCluster::ResetModuleMEs(uint32_t idet){
  std::map<uint32_t, ModMEs >::iterator pos = ClusterMEs.find(idet);
  ModMEs mod_me = pos->second;

  if(moduleswitchncluson) mod_me.NumberOfClusters->Reset();
  if(moduleswitchclusposon) mod_me.ClusterPosition->Reset();
  if(moduleswitchcluswidthon) mod_me.ClusterWidth->Reset();
  if(moduleswitchcluschargeon) mod_me.ClusterCharge->Reset();
  if(moduleswitchclusnoiseon) mod_me.ClusterNoise->Reset();
  if(moduleswitchclusstonon) mod_me.ClusterSignalOverNoise->Reset();
  mod_me.ModuleLocalOccupancy->Reset();
  mod_me.NrOfClusterizedStrips->Reset(); // can be used at client level for occupancy calculations
}






//------------------------------------------------------------------------  
void SiStripMonitorCluster::book() 
{
  
  std::vector<uint32_t> vdetId_;
  SiStripDetCabling_->addActiveDetectorsRawIds(vdetId_);

  std::vector<uint32_t> LayerSelectedDetIds;
  if(select_all_detectors){
    // select all detectors if appropriate flag is set,  for example for the mtcc
    LayerSelectedDetIds = vdetId_;
  }else{
    // use SiStripSubStructure for selecting certain regions
    SiStripSubStructure substructure;

    if(tibon) substructure.getTIBDetectors(vdetId_, LayerSelectedDetIds, 0, 0, 0, 0); // this adds rawDetIds to LayerSelectedDetIds
    if(tobon) substructure.getTOBDetectors(vdetId_, LayerSelectedDetIds, 0, 0, 0);    // this adds rawDetIds to LayerSelectedDetIds
    if(tidon) substructure.getTIDDetectors(vdetId_, LayerSelectedDetIds, 0, 0, 0, 0); // this adds rawDetIds to LayerSelectedDetIds
    if(tecon) substructure.getTECDetectors(vdetId_, LayerSelectedDetIds, 0, 0, 0, 0, 0, 0); // this adds rawDetIds to LayerSelectedDetIds
  }


  // remove any eventual zero elements - there should be none, but just in case
  for(std::vector<uint32_t>::iterator idets = LayerSelectedDetIds.begin(); idets != LayerSelectedDetIds.end(); idets++){
    if(*idets == 0) LayerSelectedDetIds.erase(idets);
  }
  

  //Histos for each detector, layer and module
  for (std::vector<uint32_t>::const_iterator detid_iter=LayerSelectedDetIds.begin();detid_iter!=LayerSelectedDetIds.end();detid_iter++){  //loop on all the active detid
    uint32_t detid = *detid_iter;

    if (detid < 1){
      edm::LogError("SiStripMonitorTrack")<< "[" <<__PRETTY_FUNCTION__ << "] invalid detid " << detid<< std::endl;
      continue;
    }
    if (DetectedLayers.find(folder_organizer.GetSubDetAndLayer(detid)) == DetectedLayers.end()){
      DetectedLayers[folder_organizer.GetSubDetAndLayer(detid)]=true;
    }    

    if(layerswitchncluson) {
      name="NumberOfClusters";
      if(MEMap.find(name)==MEMap.end()) {
	MEMap[name]=bookME1D("TH1nClusters", name.Data()); 
	name+="Trend";
	MEMap[name]=bookMETrend("TH1nClusters", name.Data());
      }
    }

    // book Layer plots      
    //    for (int j=0;j<2;j++){ 
    std::string flagtempo = "";
    folder_organizer.setLayerFolder(*detid_iter,folder_organizer.GetSubDetAndLayer(*detid_iter).second); 
    bookTrendMEs("layer",folder_organizer.GetSubDetAndLayer(*detid_iter).second,*detid_iter,flagtempo);
    //    }
  
  }//end loop on detector
  

}
  


void SiStripMonitorCluster::bookTrendMEs(TString name,int32_t layer,uint32_t id,std::string flag)//Histograms and Trends at LAYER LEVEL
{
  char rest[1024];
  int subdetid = ((id>>25)&0x7);
  if(       subdetid==3 ){
    // ---------------------------  TIB  --------------------------- //
    TIBDetId tib1 = TIBDetId(id);
    sprintf(rest,"TIB__layer__%d",tib1.layer());
  }else if( subdetid==4){
    // ---------------------------  TID  --------------------------- //
    TIDDetId tid1 = TIDDetId(id);
    sprintf(rest,"TID__side__%d__wheel__%d",tid1.side(),tid1.wheel());
  }else if( subdetid==5){
    // ---------------------------  TOB  --------------------------- //
    TOBDetId tob1 = TOBDetId(id);
    sprintf(rest,"TOB__layer__%d",tob1.layer());
  }else if( subdetid==6){
    // ---------------------------  TEC  --------------------------- //
    TECDetId tec1 = TECDetId(id);
    sprintf(rest,"TEC__side__%d__wheel__%d",tec1.side(),tec1.wheel());
  }else{
    // ---------------------------  ???  --------------------------- //
    edm::LogError("SiStripTkDQM|WrongInput")<<"no such subdetector type :"<<subdetid<<" no folder set!"<<std::endl;
    return;
  }

  SiStripHistoId hidmanager;
  std::string hid = hidmanager.createHistoLayer("",name.Data(),rest,flag);
  std::map<TString, ModMEs>::iterator iModME  = ModMEsMap.find(TString(hid));
  if(iModME==ModMEsMap.end()){
    ModMEs theModMEs; 

    //Cluster Width
    if(layerswitchcluswidthon) {
      theModMEs.LayerClusterWidth=bookME1D("TH1ClusterWidth", hidmanager.createHistoLayer("Summary_cWidth",name.Data(),rest,flag).c_str()); 
      dqmStore_->tag(theModMEs.LayerClusterWidth,layer); 
      theModMEs.LayerClusterWidthTrend=bookMETrend("TH1ClusterWidth", hidmanager.createHistoLayer("Trend_cWidth",name.Data(),rest,flag).c_str()); 
      dqmStore_->tag(theModMEs.LayerClusterWidthTrend,layer); 
    }

    //Cluster Noise
    if(layerswitchclusnoiseon) {
      theModMEs.LayerClusterNoise=bookME1D("TH1ClusterNoise", hidmanager.createHistoLayer("Summary_cNoise",name.Data(),rest,flag).c_str()); 
      dqmStore_->tag(theModMEs.LayerClusterNoise,layer); 
      theModMEs.LayerClusterNoiseTrend=bookMETrend("TH1ClusterNoise", hidmanager.createHistoLayer("Trend_cNoise",name.Data(),rest,flag).c_str()); 
      dqmStore_->tag(theModMEs.LayerClusterNoiseTrend,layer); 
    }

    //Cluster Charge
    if(layerswitchcluschargeon) {
      theModMEs.LayerClusterCharge=bookME1D("TH1ClusterCharge", hidmanager.createHistoLayer("Summary_cCharge",name.Data(),rest,flag).c_str());
      dqmStore_->tag(theModMEs.LayerClusterCharge,layer);
      theModMEs.LayerClusterChargeTrend=bookMETrend("TH1ClusterCharge", hidmanager.createHistoLayer("Trend_cCharge",name.Data(),rest,flag).c_str());
      dqmStore_->tag(theModMEs.LayerClusterChargeTrend,layer); 
    }

    //Cluster StoN
    if(layerswitchclusstonon) {
      theModMEs.LayerClusterStoN=bookME1D("TH1ClusterStoN", hidmanager.createHistoLayer("Summary_cStoN",name.Data(),rest,flag).c_str());
      dqmStore_->tag(theModMEs.LayerClusterStoN,layer); 
      theModMEs.LayerClusterStoNTrend=bookMETrend("TH1ClusterStoN", hidmanager.createHistoLayer("Trend_cStoN",name.Data(),rest,flag).c_str());
      dqmStore_->tag(theModMEs.LayerClusterStoNTrend,layer); 
    }

    //Cluster Position
    if(layerswitchclusposon) {
      theModMEs.LayerClusterPos=bookME1D("TH1ClusterPos", hidmanager.createHistoLayer("Summary_cPos",name.Data(),rest,flag).c_str());  
      dqmStore_->tag(theModMEs.LayerClusterPos,layer); 
    }


    //bookeeping
    ModMEsMap[hid]=theModMEs;
  }

}




void SiStripMonitorCluster::AllClusters( const edm::EventSetup& es)
{

  edm::DetSetVector<SiStripCluster>::const_iterator DSViter=cluster_detsetvektor->begin();
  for ( ; DSViter!=cluster_detsetvektor->end();++DSViter){

    //     uint32_t detid=DSViter->id();
    uint32_t detid=DSViter->id;

    if (find(ModulesToBeExcluded_.begin(),ModulesToBeExcluded_.end(),detid)!=ModulesToBeExcluded_.end()) continue;
    //Loop on Clusters
    LogDebug("SiStripMonitorTrack") << "on detid "<< detid << " N Cluster= " << DSViter->size();
    edm::DetSet<SiStripCluster>::const_iterator ClusIter = DSViter->begin();

    for(; ClusIter!=DSViter->end(); ClusIter++) {

      SiStripClusterInfo* SiStripClusterInfo_= new SiStripClusterInfo(detid,*ClusIter,es);
      //       LogDebug("SiStripMonitorTrack") << "ClusIter " << &*ClusIter << "\t " 
      // 				      << std::find(vPSiStripCluster.begin(),vPSiStripCluster.end(),&*ClusIter)-vPSiStripCluster.begin();
      //      if (std::find(vPSiStripCluster.begin(),vPSiStripCluster.end(),&*ClusIter) == vPSiStripCluster.end()){
      //       if ( clusterInfos(SiStripClusterInfo_,detid,"OffTrack",LV) ) countOff++;
      //std::cout << "created new sistripclusterinfo " << std::endl;

      bool clusterinforesult  = clusterInfos(SiStripClusterInfo_,detid);

      if ( clusterinforesult ) count++;

      delete SiStripClusterInfo_; 
    }//end of loop over clusters
    

  }//end of loop on detsetvectors


}




//------------------------------------------------------------------------
// bool SiStripMonitorCluster::clusterInfos(SiStripClusterInfo* cluster, const uint32_t& detid,std::string flag, const LocalVector LV)
// {
bool SiStripMonitorCluster::clusterInfos(SiStripClusterInfo* cluster, const uint32_t& detid)
{
  LogTrace("SiStripMonitorTrack") << "\n["<<__PRETTY_FUNCTION__<<"]" << std::endl;
  //folder_organizer.setDetectorFolder(0);
  if (cluster==0) return false;

  // if one imposes a cut on the clusters, apply it
  const  edm::ParameterSet ps = conf_.getParameter<edm::ParameterSet>("ClusterConditions");
  if( ps.getParameter<bool>("On") &&
      (cluster->getCharge()/cluster->getNoise() < ps.getParameter<double>("minStoN") ||
       cluster->getCharge()/cluster->getNoise() > ps.getParameter<double>("maxStoN") ||
       cluster->getWidth() < ps.getParameter<double>("minWidth") ||
       cluster->getWidth() > ps.getParameter<double>("maxWidth")                    )) return false;
  // start of the analysis
  
  int SubDet_enum = StripSubdetector(detid).subdetId()-3;

  NClus[SubDet_enum]++;
  std::stringstream ss;
  //  const_cast<SiStripClusterInfo*>(cluster)->print(ss);
  LogTrace("SiStripMonitorTrack") << "\n["<<__PRETTY_FUNCTION__<<"]\n" << ss.str() << std::endl;

  std::string name;
  char rest[1024];
  int subdetid = ((detid>>25)&0x7);
  if(       subdetid==3 ){
    // ---------------------------  TIB  --------------------------- //
    TIBDetId tib1 = TIBDetId(detid);
    sprintf(rest,"TIB__layer__%d",tib1.layer());
  }else if( subdetid==4){
    // ---------------------------  TID  --------------------------- //
    TIDDetId tid1 = TIDDetId(detid);
    sprintf(rest,"TID__side__%d__wheel__%d",tid1.side(),tid1.wheel());
  }else if( subdetid==5){
    // ---------------------------  TOB  --------------------------- //
    TOBDetId tob1 = TOBDetId(detid);
    sprintf(rest,"TOB__layer__%d",tob1.layer());
  }else if( subdetid==6){
    // ---------------------------  TEC  --------------------------- //
    TECDetId tec1 = TECDetId(detid);
    sprintf(rest,"TEC__side__%d__wheel__%d",tec1.side(),tec1.wheel());
  }else{
    // ---------------------------  ???  --------------------------- //
    edm::LogError("SiStripTkDQM|WrongInput")<<"no such subdetector type :"<<subdetid<<" no folder set!"<<std::endl;
    return 0;
  }

  SiStripHistoId hidmanager1;
   
  //Filling Layer Plots
  std::string flag = "";
  name= hidmanager1.createHistoLayer("","layer",rest,flag);


  fillTrendMEs(cluster,name);
  
  return true;

}//end of method clusterInfos







// void SiStripMonitorCluster::fillTrendMEs(SiStripClusterInfo* cluster,std::string name,float cos, std::string flag)
// { 
void SiStripMonitorCluster::fillTrendMEs(SiStripClusterInfo* cluster,std::string name)
{ 


  //   std::map<TString, ModMEs>::iterator iModME  = ModMEsMap.find(TString(name));
  std::map<TString, ModMEs>::iterator iModME  = ModMEsMap.find(name);

  if(iModME!=ModMEsMap.end()){


    if(layerswitchclusstonon) {
      fillME(iModME->second.LayerClusterStoN  ,cluster->getCharge()/cluster->getNoise());
      fillTrend(iModME->second.LayerClusterStoNTrend,cluster->getCharge()/cluster->getNoise());
    }
    
    if(layerswitchcluschargeon) {
      fillME(iModME->second.LayerClusterCharge,cluster->getCharge());
      fillTrend(iModME->second.LayerClusterChargeTrend,cluster->getCharge());
    }

    if(layerswitchclusnoiseon) {
      fillME(iModME->second.LayerClusterNoise ,cluster->getNoise());
      fillTrend(iModME->second.LayerClusterNoiseTrend,cluster->getNoise());
    }

    if(layerswitchcluswidthon) {
      fillME(iModME->second.LayerClusterWidth ,cluster->getWidth());
      fillTrend(iModME->second.LayerClusterWidthTrend,cluster->getWidth());
    }

    if(layerswitchclusposon) {
      fillME(iModME->second.LayerClusterPos   ,cluster->getPosition());
    }
    
  }

}//end of method fillTrendMEs





//--------------------------------------------------------------------------------
void SiStripMonitorCluster::fillTrend(MonitorElement* me ,float value)
{
  if(!me) return;
  //check the origin and check options
  int option = conf_.getParameter<edm::ParameterSet>("Trending").getParameter<int32_t>("UpdateMode");
  if(firstEvent==-1) firstEvent = eventNb;
  int CurrentStep = atoi(me->getAxisTitle(1).c_str()+8);
  int firstEventUsed = firstEvent;
  int presentOverflow = (int)me->getBinEntries(me->getNbinsX()+1);
  if(option==2) firstEventUsed += CurrentStep * int(me->getBinEntries(me->getNbinsX()+1));
  else if(option==3) firstEventUsed += CurrentStep * int(me->getBinEntries(me->getNbinsX()+1)) * me->getNbinsX();
  //fill
  me->Fill((eventNb-firstEventUsed)/CurrentStep,value);

  if(eventNb-firstEvent<1) return;
  // check if we reached the end
  if(presentOverflow == me->getBinEntries(me->getNbinsX()+1)) return;
  switch(option) {
  case 1:
    {
      // mode 1: rebin and change X scale
      int NbinsX = me->getNbinsX();
      float entries = 0.;
      float content = 0.;
      float error = 0.;
      int bin = 1;
      int totEntries = int(me->getEntries());
      for(;bin<=NbinsX/2;++bin) {
	content = (me->getBinContent(2*bin-1) + me->getBinContent(2*bin))/2.; 
	error   = pow((me->getBinError(2*bin-1)*me->getBinError(2*bin-1)) + (me->getBinError(2*bin)*me->getBinError(2*bin)),0.5)/2.; 
	entries = me->getBinEntries(2*bin-1) + me->getBinEntries(2*bin);
	me->setBinContent(bin,content*entries);
	me->setBinError(bin,error);
	me->setBinEntries(bin,entries);
      }
      for(;bin<=NbinsX+1;++bin) {
	me->setBinContent(bin,0);
	me->setBinError(bin,0);
	me->setBinEntries(bin,0); 
      }
      me->setEntries(totEntries);
      char buffer[256];
      sprintf(buffer,"EventId/%d",CurrentStep*2);
      me->setAxisTitle(std::string(buffer),1);
      break;
    }
  case 2:
    {
      // mode 2: slide
      int bin=1;
      int NbinsX = me->getNbinsX();
      for(;bin<=NbinsX;++bin) {
	me->setBinContent(bin,me->getBinContent(bin+1)*me->getBinEntries(bin+1));
	me->setBinError(bin,me->getBinError(bin+1));
	me->setBinEntries(bin,me->getBinEntries(bin+1));
      }
      break;
    }
  case 3:
    {
      // mode 3: reset
      int NbinsX = me->getNbinsX();
      for(int bin=0;bin<=NbinsX;++bin) {
	me->setBinContent(bin,0);
	me->setBinError(bin,0);
	me->setBinEntries(bin,0); 
      }
      break;
    }
  }
}










MonitorElement* SiStripMonitorCluster::bookMETrend(const char* ParameterSetLabel, const char* HistoName)
{
  Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
  edm::ParameterSet ParametersTrend =  conf_.getParameter<edm::ParameterSet>("Trending");
  MonitorElement* me = dqmStore_->bookProfile(HistoName,HistoName,
					      ParametersTrend.getParameter<int32_t>("Nbins"),
					      // 					      0,
					      ParametersTrend.getParameter<double>("xmin"),
					      ParametersTrend.getParameter<double>("xmax"),
					      // 					      ParametersTrend.getParameter<int32_t>("Nbins"),
					      100, //that parameter should not be there !?
					      ParametersTrend.getParameter<double>("ymin"),
					      ParametersTrend.getParameter<double>("ymax"),
					      "" );
  if(!me) return me;
  char buffer[256];
  sprintf(buffer,"EventId/%d",ParametersTrend.getParameter<int32_t>("Steps"));
  me->setAxisTitle(std::string(buffer),1);
  return me;
}

//------------------------------------------------------------------------------------------






MonitorElement* SiStripMonitorCluster::bookME1D(const char* ParameterSetLabel, const char* HistoName)
{
  Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
  return dqmStore_->book1D(HistoName,HistoName,
			   Parameters.getParameter<int32_t>("Nbinx"),
			   Parameters.getParameter<double>("xmin"),
			   Parameters.getParameter<double>("xmax")
			   );
}
