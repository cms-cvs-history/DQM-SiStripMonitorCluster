
// 2006.05.31 dkcira - plot some histograms written out from the TkDQM in a root file
{
  // some style parameters
  gStyle->SetCanvasColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetHistFillColor(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetMarkerStyle(23);
  gStyle->SetMarkerSize(2.5);
  gStyle->SetMarkerColor(1);

  // read in file with histograms
//  TString histo_file = "test_digi_cluster.root";
  TString histo_file = "~/scratch0/test_digi_cluster.root";
  TFile dqmf(histo_file);
//  dqmf.ls();


  //
  dqmf.cd("DQMData/SiStrip/MechanicalView/TOB/layer_3/forward_rods/rod_12/module_436440068");
  ADCsCoolestStrip__det__436440068->SetXTitle(ADCsCoolestStrip__det__436440068->GetTitle()); 
  ADCsHottestStrip__det__436440068->SetXTitle(ADCsHottestStrip__det__436440068->GetTitle()); 
  DigisPerDetector__det__436440068->SetXTitle(DigisPerDetector__det__436440068->GetTitle()); 
  ClusterCharge__det__436440068->SetXTitle(ClusterCharge__det__436440068->GetTitle()); 
  ClusterPosition__det__436440068->SetXTitle(ClusterPosition__det__436440068->GetTitle()); 
  ClusterWidth__det__436440068->SetXTitle(ClusterWidth__det__436440068->GetTitle()); 
  ClustersPerDetector__det__436440068->SetXTitle(ClustersPerDetector__det__436440068->GetTitle()); 
  NrOfClusterizedStrips__det__436440068->SetXTitle(NrOfClusterizedStrips__det__436440068->GetTitle()); 

  ADCsCoolestStrip__det__436440068->SetTitle("TOB/layer_3/forward_rods/rod_12");
  ADCsHottestStrip__det__436440068->SetTitle("TOB/layer_3/forward_rods/rod_12");
  DigisPerDetector__det__436440068->SetTitle("TOB/layer_3/forward_rods/rod_12");
  ClusterCharge__det__436440068->SetTitle("TOB/layer_3/forward_rods/rod_12");
  ClusterPosition__det__436440068->SetTitle("TOB/layer_3/forward_rods/rod_12");
  ClusterWidth__det__436440068->SetTitle("TOB/layer_3/forward_rods/rod_12");
  ClustersPerDetector__det__436440068->SetTitle("TOB/layer_3/forward_rods/rod_12");
  NrOfClusterizedStrips__det__436440068->SetTitle("TOB/layer_3/forward_rods/rod_12");


  char *s = new char[1];
  TCanvas *c1 = new TCanvas("c1","c1",500,0,700,700);
  ADCsCoolestStrip__det__436440068->Draw();
//  cout<<"hit return to continue"<<endl; gets(s);  // "return" is enough, root waits until you hit a key

  TCanvas *c2 = new TCanvas("c2","c2",500,0,700,700);
  ADCsHottestStrip__det__436440068->Draw();
//  cout<<"hit return to continue"<<endl; gets(s);

  TCanvas *c3 = new TCanvas("c3","c3",500,0,700,700);
  DigisPerDetector__det__436440068->Draw();
//  cout<<"hit return to continue"<<endl; gets(s);

  TCanvas *c4 = new TCanvas("c4","c4",500,0,700,700);
  ClusterCharge__det__436440068->Draw();
//  cout<<"hit return to continue"<<endl; gets(s);

  TCanvas *c5 = new TCanvas("c5","c5",500,0,700,700);
  ClusterPosition__det__436440068->Draw();
//  cout<<"hit return to continue"<<endl; gets(s);

  TCanvas *c6 = new TCanvas("c6","c6",500,0,700,700);
  ClusterWidth__det__436440068->Draw();
//  cout<<"hit return to continue"<<endl; gets(s);

  TCanvas *c7 = new TCanvas("c7","c7",500,0,700,700);
  ClustersPerDetector__det__436440068->Draw();
//  cout<<"hit return to continue"<<endl; gets(s);

  TCanvas *c8 = new TCanvas("c8","c8",500,0,700,700);
  NrOfClusterizedStrips__det__436440068->Draw();
//  cout<<"hit return to continue"<<endl; gets(s);

}
