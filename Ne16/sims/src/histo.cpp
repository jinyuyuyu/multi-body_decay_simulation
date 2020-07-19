#include "histo.h"

histo::histo()
{
  opf = new TFile("sim.root","recreate");
  opf->cd();

  hist_pcbeam = new TH1I("pcbeam","pcbeam",500,8950,9100);
  hist_Ebeam = new TH1I("Ebeam","Ebeam",400,2000,2140);
  hist_beamspot = new TH2I("beamspot","beamspot",100,-9,9,100,-9,9);
  hist_xTar = new TH1I("xTar","xTar",100,-9,9);
  hist_yTar = new TH1I("yTar","yTar",100,-9,9);

  hist_knockout_vel = new TH1I("knockout_vel","knockout_vel",200,11.,15.);
  hist_knockout_theta = new TH1I("knockout_theta","knockout_theta",200,-1,5);
  hist_knockout_phi = new TH1I("knockout_phi","knockpout_phi",200,-40,400);
  hist_mscat_vel = new TH1I("mscat_vel","mscat_vel",200,11.,15.);
  hist_mscat_theta = new TH1I("mscat_theta","mscat_theta",200,-1,5);
  hist_mscat_phi = new TH1I("mscat_phi","mscat_phi",200,-40,400);
  hist_mscat_pc = new TH1I("mscat_pc","mscat_pc",500,8000,9400);

  hist_ErelReal = new TH1I("ErelReal","ErelReal",1000,1.4,1.44);
  hist_ErelRecon = new TH1I("ErelRecon","ErelRecon",375,0,3);

  hist_p_range = new TH1I("p_range","p_range [mm]",1000,0,100);

  //in lab frame
  hist_p1Energy = new TH1I("p1_Energy","p1_energy",1000,70,130);
  hist_p1Theta = new TH1I("p1_Theta","p1_Theta",200,0,0.2);
  hist_p1Phi = new TH1I("p1_Phi","p1_Phi",200,0,8);
  hist_p2Energy = new TH1I("p2_Energy","p2_energy",1000,70,130);
  hist_p2Theta = new TH1I("p2_Theta","p2_Theta",200,0,0.2);
  hist_p2Phi = new TH1I("p2_Phi","p2_Phi",200,0,8);
  hist_coreEnergy = new TH1I("core_Energy","core_energy",1000,1200,2200);
  hist_coreTheta = new TH1I("core_Theta","core_Theta",200,0,0.2);
  hist_corePhi = new TH1I("core_Phi","core_Phi",200,0,8);
  //in center-of-mass frame
  //hist_p1Energy = new TH1I("p1_Energy","p1_energy",200,0,1);
  //hist_p1Theta = new TH1I("p1_Theta","p1_Theta",200,0,3.14);
  //hist_p1Phi = new TH1I("p1_Phi","p1_Phi",200,0,8);
  //hist_p2Energy = new TH1I("p2_Energy","p2_energy",200,0,1);
  //hist_p2Theta = new TH1I("p2_Theta","p2_Theta",200,0,3.14);
  //hist_p2Phi = new TH1I("p2_Phi","p2_Phi",200,0,8);
  //hist_coreEnergy = new TH1I("core_Energy","core_energy",200,0,1);
  //hist_coreTheta = new TH1I("core_Theta","core_Theta",200,0,3.14);
  //hist_corePhi = new TH1I("core_Phi","core_Phi",200,0,8);
}

void histo::write()
{
  opf->Write();
  cout << "histograms were written to sim.root" << endl;
  opf->Close();
}
