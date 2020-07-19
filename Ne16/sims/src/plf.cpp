#include "plf.h"

CPlf::CPlf(double Z,double mbeam,double mplf,string filename,double thick,bool einstein0)
{
  einstein = einstein0;
  beam = new CFrame(mbeam,einstein); //incident beam
  frame = new CFrame(mplf,einstein); //projectile-like nuclei after knockout reaction
  loss = new CLoss(filename,mbeam);
  momDist = new CMomDist();
  multScat = new CMultScat(Z,thick);
  factor = 1.3;
  mass = mbeam*amu;
}

CPlf::~CPlf()
{
  delete beam;
  delete frame;
  delete loss;
  delete momDist;
  delete multScat;
}

void CPlf::getPlf(double v_beam)
{
  double ptrx = momDist->getTransMom()*factor;
  double ptry = momDist->getTransMom()*factor;
  double pz = momDist->getLongMom()*factor;
  ptrx = 497.134;
  ptry = -99.5223;
  pz = -115.889;
  cout << "ptrx=" << ptrx << " ptry=" << ptry << " pz=" << pz << endl;

  double ptr = TMath::Sqrt(TMath::Power(ptrx,2)+TMath::Power(ptry,2));
  double ptot = TMath::Sqrt(TMath::Power(ptr,2)+TMath::Power(pz,2));
  double Etot = TMath::Sqrt(TMath::Power(ptot,2)+TMath::Power(mass,2));

  double pz_new = (pz+Etot*v_beam/c)/sqrt(1-TMath::Power(v_beam/c,2));
  double pmag = TMath::Sqrt(TMath::Power(pz_new,2)+TMath::Power(ptr,2));
  double phi = atan2(ptry,ptrx);
  double theta = atan2(ptr,pz_new);
  double vv = c*pmag/TMath::Sqrt(TMath::Power(pmag,2)+TMath::Power(mass,2));

  frame->v[0] = vv*sin(theta)*cos(phi);
  frame->v[1] = vv*sin(theta)*sin(phi);
  frame->v[2] = vv*cos(theta);
  frame->velocity = vv;
  frame->theta = theta;
  frame->phi = phi;
  cout << "pmag=" << pmag << endl;
  cout << "pz_new=" << pz_new << " ptr=" << ptr << endl;
  cout << "vv=" << vv << " theta=" << theta << " phi=" << phi << endl;
  frame->getEnergy();
  cout << "pmag=" << frame->pcTot << endl;
  cout << "pz_new=" << frame->pcTot*cos(frame->theta) << " ptr=" << frame->pcTot*sin(frame->theta) << endl;
  cout << "vv=" << frame->velocity << " theta=" << frame->theta << " phi=" << frame->phi << endl;
}

void CPlf::getPlf1()
{
  double ptrx = momDist->getTransMom()*factor;
  double ptry = momDist->getTransMom()*factor;
  double pz = momDist->getLongMom()*factor;
  //frame->v[0] = ptrx/mass*c;
  //frame->v[1] = ptry/mass*c;
  //frame->v[2] = pz/mass*c;
  frame->v[0] = ptrx/frame->mass*c;
  frame->v[1] = ptry/frame->mass*c;
  frame->v[2] = pz/frame->mass*c;
  frame->transformVelocity(beam->v);
  //cout << "pmag=" << frame->pcTot << endl;
  //cout << "pz_new=" << frame->pcTot*cos(frame->theta) << " ptr=" << frame->pcTot*sin(frame->theta) << endl;
  //cout << "vv=" << frame->velocity << " theta=" << frame->theta << " phi=" << frame->phi << endl;
}

void CPlf::multiScat(double Ebeam,double fractionalThick)
{
  if(fractionalThick == 0.) return;
  multScat->scatter(Ebeam,fractionalThick,beam->theta,beam->phi);
  double thetaNew = multScat->getThetaNew();
  double phiNew = multScat->getPhiNew();

  beam->theta = thetaNew;
  beam->phi = phiNew;
  //cout << "beam->theta=" << beam->theta << " beam->phi=" << beam->phi << " beam->energy=" << beam->energy << endl;
  beam->getVelocity();
}

void CPlf::getBeam(double Ebeam,double thick)
{
  beam->energy = Ebeam;
  beam->theta = 0;
  beam->phi = 0;
  
  if(thick>0)
    beam->energy = loss->getEout(Ebeam,thick/cos(beam->theta));
  beam->getVelocity();
}
