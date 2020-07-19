#ifndef FRAME_
#define FRAME_
#include <cmath>
#include "constants.h"
#include <iostream> 
using namespace std;

/**
 *!\brief Relativistic and Non-Relativistic stuff
 *
 * give velocity, energies, transforms to new reference frame
 * using either non-relativistic or Relativistics equations
 */

class CFrame
{
 public:
  CFrame(double,bool);

  bool einstein; //use newtonian or relativistic
  double mass;  //!< rest mass of fragment in MeV
  double energy;  //!< fragments energy in MeV
  double totEnergy; //!< rest mass plus kinetic energy
  double velocity; //!< fragments velocity in cm/ns 
  double v[3];  //!< velocity vector of fragment
  double pcTot;  //!< total momentum*c of fragment
  double pc[3]; //!< momentum*c vector of fragment MeV
  double theta; //!< theta nagle of fragment in radians
  double phi;  //!< phi angle
 
  double getVelocity();
  double getVelocityNewton();
  double getVelocityRel();
  double getEnergy();
  double getEnergyNewton();
  double getEnergyRel();
  void transformVelocity(double*);
  void transformVelocityNewton(double*);
  void transformVelocityRel(double*);
  void getVelocityFromMom();
  void getVelocityFromMom0();
  void getVelocityFromMomNewton();
  void getVelocityFromMomRel();
  void getAngle();

};


#endif
