#include "multScat.h"

//angular straggling caused by multiple Coulomb scattering
//formula from R. Anne et. al., NIMB 34, 295 (1988)
//current settings are for 9Be target
//change Ztar, Atar, dtar for other targets

CMultScat::CMultScat(double Zpro0, double totalThick0)
{
  //constants;
  a0 = 0.529e-8;
  Ztar = 4.;
  Atar = 9.;
  dtar = 0.1848; //density of 9Be [mg/cm/cm/um]

  Zpro = Zpro0;
  N = dtar/Atar/1000*6.02e23; //number of target atoms in unit volume
  totalThick = totalThick0;

  a = 0.885*a0/TMath::Sqrt(TMath::Power(Zpro,2./3.)+TMath::Power(Ztar,2./3.));
  totalTau = TMath::Pi()*a*a*N*totalThick;
  factor = 16.26/Zpro/Ztar/TMath::Sqrt(TMath::Power(Zpro,2./3.)+TMath::Power(Ztar,2./3.))*1000;
  cout << "Zpro=" << Zpro << "  totalThick=" << totalThick << endl;
  cout << "a=" << a << "  totalTau=" << totalTau << "  factor=" << factor << endl;
}

void CMultScat::scatter(double energy, double fractionalThick, double theta, double phi)
{
  if(fractionalThick==0||totalThick==0)
  {
    thetaNew = theta;
    phiNew = phi;
    return;
  }

  double tau = fractionalThick*totalTau;
  double alphaBar = TMath::Power(tau,0.55);
  double alpha = alphaBar/energy/factor; //in [rad]
  //cout << "tau=" << tau << "  alphaBar=" << alphaBar << "  alpha=" << alpha << endl;
  double sigma = alpha*2./2.355; //alpha is the half width at half maximum

  double dTheta = ran.Gaus(0.,sigma);
  double dPhi = 2*TMath::Pi()*ran.Rndm();

  double x = sin(dTheta)*cos(dPhi);
  double y = sin(dTheta)*sin(dPhi);
  double z = cos(dTheta);

  //rotate in z-x plane by theta
  double xx = x*cos(theta) + z*sin(theta);
  double yy = y;
  double zz = z*cos(theta) - x*sin(theta);

  //rotate in x-y plane by phi
  double xxx = xx*cos(phi) - yy*sin(phi);
  double yyy = yy*cos(phi) + xx*sin(phi);
  double zzz = zz;

  thetaNew = acos(zzz);
  phiNew = atan2(yyy,xxx);
  if (phiNew < 0.) phiNew += 2.*pi;
}

double CMultScat::getThetaNew()
{
  return thetaNew;
}

double CMultScat::getPhiNew()
{
  return phiNew;
}
