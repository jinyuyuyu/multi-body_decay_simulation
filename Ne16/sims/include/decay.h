#ifndef DECAY_
#define DECAY_

#include <valarray>
#include <TRandom3.h>
#include <TGenPhaseSpace.h> 
#include <TLorentzVector.h> 
#include "frag.h"
#include "fragment.h"
#include "moscow.h"
#include "constants.h"

struct prop
{
  double Erel;
  double Dvelocity;
  double plfTheta;
  double plfPhi;
  double plfVel;
  double thetaEmission;
};


/**
 *!\brief selects the veloity vectors of the secondary fragments
 */

class CDecay
{
 public:
  CFrame* frag6Be;
  double CosTheta_T;
  double x_T;
  double x_Y[2];
  double CosTheta_Y[2];

  void getJacobi(CFrame**,bool);
  void getJacobiPrimary();
  void getJacobiSecondary();

  bool einstein;
  TRandom3 ran;
  CFrame *real[5]; 
  CFrame *recon[5];
  CFrame *plfRecon;
  CFrame *partCM[5];
  CFrag *frag[5];  //!< information about the decay fragments

  CDecay(int,CFrag**,bool einstein0);
  ~CDecay();
  double getErelReal();
  double getErelRecon();
  double getErel(CFrame**);
  double getErelNewton(CFrame**);
  double getErelRel(CFrame**);
  double getErelRel2(CFrame**);

  void ModeComplex();
  void ModeMoscow();
  void ModeMicroCanonical(double);
  void ModeMicroCanonicalBe();
  void micro(int,CFrame**,double);
  double micro2(int,CFrame**,double);
  void getEk6Be(CFrame**);
  double getEk3body(CFrame*,CFrame*,CFrame*);
  double getEk3bodyNewton(CFrame*,CFrame*,CFrame*);
  double getEk3bodyRel(CFrame*,CFrame*,CFrame*);
  double getEk4bodyRel(CFrame*,CFrame*,CFrame*,CFrame*);
  double getEk3();
  void getEk6BeRecon();
  void getEk6BeReal();
  prop Recon;

  int nFrags;
  double  sumA;
  double ErelRecon; //!<reconstructed relative kinetic energy


  double aaRel;
  double ppRel[6];
  double aapRelMin;
  double aapRelMax;
  double pEnergyMin;
  double pEnergyMax;
  double ppThetaRel;
  static double const EkTot8B;
  static double const gamma8B;
  double Ek6Be[6];
  double aRatio;

  int Nsolution; //!< number of p-p pairs with correct 6Be energy
  int Isolution; //!< solution #

  moscow * Moscow;
  
  TGenPhaseSpace genPhaseSpace;
};

#endif
