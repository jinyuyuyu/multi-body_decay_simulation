#include <iostream>
#include <TMath.h>
#include <TRandom3.h>
#include "constants.h"
#include "histo.h"
#include "loss.h"
#include "multScat.h"
#include "plf.h"
#include "frag.h"
#include "decay.h"
#include "range.h"
using namespace std;

int main()
{
  //variable parameters =========================
  double Qbreakup = 1.42;
  double P_acceptance = 0.01;
  double beamSize = 10.; //beam spot size [mm]
  double thickness = 1000.; //1mm 9Be [micron]
  double csiRes = 0.1275;

  //target and detector =========================
  double dist = 352;// [mm]
  double thickness_Al = 5920; //5.92mm Al [micron]
  double thickness_Si = 1014; //1014um Si [micron]
  double thickness_CsI = 50; //50mm CsI [mm]

  //beam ========================================
  double Abeam = 20.;
  double Zbeam = 12.;
  double mbeam = mass_20Mg;
  double EPAbeam = 103.428; //[MeV/A] Brho=2.5079
  string lossbeam = "Mg20inBe";
  double Ebeam0 = Abeam*EPAbeam;
  double pc0 = TMath::Sqrt(TMath::Power(Ebeam0+mbeam*amu,2)-TMath::Power(mbeam*amu,2));
  cout << "pc0=" << pc0 << endl;

  //decay properties ============================
  const int nFrags = 3;
  double mplf = mass_16Ne;
  double Zfrag[nFrags] = {1., 1., 8.};
  double mfrag[nFrags] = {mass_1H, mass_1H, mass_14O};
  string lossfrag[nFrags] = {"H1inBe","H1inBe","O14inBe"};

  //beam and proton loss file ===================
  CLoss loss_MginBe(lossbeam,mbeam);
  CLoss loss_HinAl("H1inAl",mass_1H);
  CLoss loss_HinSi("H1inSi",mass_1H);
  CRange range_HinCsI("Hydrogen");

  //other parameters ============================
  bool einstein = 1; //switch for newtonian (0) or relativistic (1)
  const int Nevents = 100000; //loop counts

  //decay relavant class ========================
  CPlf plf(Zbeam,mbeam,mplf,lossbeam,thickness,einstein);
  CFrag *frag[nFrags];
  for(int i=0;i<nFrags;i++)
  {
    frag[i] = new CFrag(Zfrag[i],mfrag[i],lossfrag[i],csiRes,thickness,dist,einstein);
  }
  CDecay decay(nFrags,frag,einstein);

  histo *Histo = new histo();
  TRandom3 ran;
  double xTar, yTar; //beam incident position on target

  int nDetProton[nFrags] = {0};
  int nDetCore = 0;
  int nDet = 0; //detection of all protons and the core
  int nDetClean = 0; //clean detection efficiency

  //=================================================================
  //loop begins =====================================================
  for(int ievent = 0;ievent<Nevents;ievent++)
  {
    if(ievent%3000==0) cout << "\r" << ievent << "/" << Nevents << flush;

    //beam momentum acceptance
    double delta_pc = (0.5-ran.Rndm())*P_acceptance*pc0; //jinyu
    double pc = delta_pc+pc0;
    double Ebeam = TMath::Sqrt(TMath::Power(pc,2)+TMath::Power(mbeam*amu,2))-mbeam*amu;
    Histo->hist_pcbeam->Fill(pc);
    Histo->hist_Ebeam->Fill(Ebeam);

    //beam spot size
    //xTar = (1.-2.*ran.Rndm())/2.*beamSize;
    //yTar = (1.-2.*ran.Rndm())/2.*beamSize;
    ran.Rannor(xTar,yTar);
    xTar *= beamSize / 6.; //plus minus 3sigma, beamSize=6*sigma
    yTar *= beamSize / 6.; //plus minus 3sigma, beamSize=6*sigma
    Histo->hist_beamspot->Fill(xTar,yTar);
    Histo->hist_xTar->Fill(xTar);
    Histo->hist_yTar->Fill(yTar);

    //reaction location in z-axis
    double dthick = ran.Rndm()*thickness;
    double dthick_beam = thickness-dthick;
    plf.getBeam(Ebeam,dthick_beam);
    plf.multiScat(Ebeam,dthick_beam/thickness); //multiple scattering
    //cout << "plf.beam->v[0]=" << plf.beam->v[0] << endl;
    Histo->hist_mscat_vel->Fill(plf.beam->velocity);
    Histo->hist_mscat_theta->Fill(plf.beam->theta*180./TMath::Pi());
    Histo->hist_mscat_phi->Fill(plf.beam->phi*180./TMath::Pi());
    Histo->hist_mscat_pc->Fill(plf.beam->pcTot);
    //plf.getPlf(v_beam); //neutron knockout
    plf.getPlf1(); //neutron knockout
    Histo->hist_knockout_vel->Fill(plf.frame->velocity);
    Histo->hist_knockout_theta->Fill(plf.frame->theta*180./TMath::Pi());
    Histo->hist_knockout_phi->Fill(plf.frame->phi*180./TMath::Pi());

    decay.ModeMicroCanonical(Qbreakup);
    //double ErelReal = decay.getErelReal();
    //cout << "ErelReal=" << ErelReal << endl;
    //Histo->hist_ErelReal->Fill(ErelReal);
    //remember to getJacobiPrimary

    for(int i=0;i<nFrags;i++)
    {
      frag[i]->AddVelocity(plf.frame->v);
    }
    double ErelReal = decay.getErelReal();
    //cout << "ErelReal=" << ErelReal << endl;
    Histo->hist_ErelReal->Fill(ErelReal);

    //fragment interaction in target ========
    bool stop = 0; //if the fragment is stopped in target
    for(int i=0;i<nFrags;i++)
    {
      if(frag[i]->targetInteraction(dthick,thickness))
      {
        stop = 1;
        break;
      }
    }
    if(stop) continue;

    Histo->hist_p1Energy->Fill(frag[0]->real->energy);
    Histo->hist_p1Theta->Fill(frag[0]->real->theta);
    Histo->hist_p1Phi->Fill(frag[0]->real->phi);
    Histo->hist_p2Energy->Fill(frag[1]->real->energy);
    Histo->hist_p2Theta->Fill(frag[1]->real->theta);
    Histo->hist_p2Phi->Fill(frag[1]->real->phi);
    Histo->hist_coreEnergy->Fill(frag[2]->real->energy);
    Histo->hist_coreTheta->Fill(frag[2]->real->theta);
    Histo->hist_corePhi->Fill(frag[2]->real->phi);

    ////detection =====================================================
    stop = 0;
    double p_range[nFrags-1];
    double p_straggle[nFrags-1];
    for(int i=0;i<nFrags-1;i++)
    {
      double thick_Al = thickness_Al/cos(frag[i]->real->theta);
      double E11 = loss_HinAl.getEout(frag[i]->real->energy,thick_Al);
      double thick_Si = thickness_Si/cos(frag[i]->real->theta);
      double E1 = loss_HinSi.getEout(E11,thick_Si);
      p_range[i] = range_HinCsI.getRange(E1);
      if(p_range[i]==0) {stop = 1; break;}
      if(p_range[i]*cos(frag[i]->real->theta)>thickness_CsI) {stop = 1; break;}
      Histo->hist_p_range->Fill(p_range[i]*cos(frag[i]->real->theta));
      p_straggle[i] = range_HinCsI.getLateralStraggle();
    }
    if(stop) continue;

    //proton hit on S4-CsI
    //cout << "before: nDetProton[0]=" << nDetProton[0] << " nDetProton[1]=" << nDetProton[1]
    //  << " nDetProton[2]=" << nDetProton[2] << endl;
    int nHitProton = 0;
    for(int i=0;i<nFrags-1;i++)
      nHitProton += frag[i]->protonHit(xTar,yTar,p_range[i],p_straggle[i]);
    nDetProton[nHitProton]++;
    //cout << "nHitProton=" << nHitProton << endl;
    //cout << "after: nDetProton[0]=" << nDetProton[0] << " nDetProton[1]=" << nDetProton[1]
    //  << " nDetProton[2]=" << nDetProton[2] << endl;
    //cout << endl;

    //core hit on fiber and S800
    int nHitCore = frag[nFrags-1]->coreHit(xTar,yTar);
    if(nHitCore==1) nDetCore++;

    if(nHitProton!=(nFrags-1) || nHitCore!=1) continue;
    nDet++;

    //protons can not hit on a same pie/ring/CsI
    stop = 0;
    for(int i=0;i<nFrags-2;i++)
    {
      for(int j=i+1;j<nFrags-1;j++)
      {
        if(frag[i]->ipie==frag[j]->ipie) {stop = 1; break;}
        if(frag[i]->iring==frag[j]->iring) {stop = 1; break;}
        if(frag[i]->icsi==frag[j]->icsi) {stop = 1; break;}
      }
      if(stop) break;
    }
    if(stop) continue;
    //if(frag[0]->ipie==frag[1]->ipie) continue;
    //if(frag[0]->iring==frag[1]->iring) continue;
    //if(frag[0]->icsi==frag[1]->icsi) continue;

    nDetClean++;

    //correct for energy loss in target
    for(int i=0;i<nFrags;i++)
      frag[i]->Egain(thickness/2.);

    //reconstruct invariant-mass
    double ErelRecon = decay.getErelRecon();
    Histo->hist_ErelRecon->Fill(ErelRecon);
  }
  //loop ends =======================================================
  //=================================================================

  cout << endl;
  cout << endl;

  cout << "simulated " << Nevents << " events" << endl;
  for(int i=0;i<nFrags;i++)
  {
    cout << "detect " << i << " proton: " << nDetProton[i] << endl;
  }
  cout << "detect core: " << nDetCore << endl;
  cout << "detection efficiency of " << nFrags-1 << " protons+core: " << (double)nDet/(double)Nevents << endl;
  cout << "clean efficiency: " << (double)nDetClean/(double)Nevents << endl;

  cout << endl;
  cout << endl;

  Histo->write();
}
