#include "frag.h"

/**
 * Constructor
 \param Z0 is the atomic number of fragment
 \param mass0 is mass of fragment in amu
 \param filename specifies the name of the file with energy loss information
 \param CsI_res0 fractional energy resolution of this fragment in CsI
 \param thickness is target thickness in micron 
 \param dist0 is distance from target to detector in mm
 */

CFrag::CFrag(double Z0,double mass0, string filename, double CsI_res0,
    double thickness, double dist0, bool einstein0):fragment(Z0,mass0)
{
  if (filename != "") loss = new CLoss(filename,mass);
  CsI_res = CsI_res0;
  multScat = new CMultScat(Z,thickness);

  einstein = einstein0;
  real = new CFrame(mass,einstein);
  recon = new CFrame(mass,einstein);

  //S4-CsI setting
  dist = dist0; //352mm, distance from target to S4
  r_inner = 7.5; //mm
  r_outer = 62.5; //mm
  Nrings = 128;
  Npies = 128;
  r_CsI = 30; //radius of CsI division between inner and outer [mm]

  //fiber-S800 setting
  dist2Si = 76; //distance from Si to fiber
  r_S4Hole = 5.; //radius of S4 hole [mm]
  r_FiberPlate = 8.35; //radius of hole of fiber plate [mm]
  xmax = 8.;
  ymax = 8.;
  fiber_res = 0.25; //Kyle used 0.95 [mm]

  Plane = new plane_det(565.,20.,20.,0.95);//defining the "fibers"
  //Plane = new plane_det(dist+dist2Si,8.,8.,0.25);//defining the "fibers"
}
//*********************************************************
/**
 *Destructor
 */
CFrag::~CFrag()
{
  delete loss;
  delete multScat;
  delete real;
  delete recon;
  delete Plane;
}

//******************************************************************
/**
 * Add a velocity vector to the fragments velocity vector.
 * Used to transform between reference frames
 */
void CFrag::AddVelocity(double *Vplf)
{
  real->transformVelocity(Vplf);
}

//****************************************************************
/** 
 * returns the energy after the fragment has exited the target
 \param thick is the thickness of target that the particle has to traverse (mg/cm2)
 */
double CFrag::Eloss(double thick)
{
  real->energy = loss->getEout(real->energy,thick);
  return real->energy;
}

//*******************************************************************
/**
 * corrects energy of a detected particle for the energy loss
 * in the target.
 \param thick is the thickness of target material though which the particle passed (mg/cm2)
 */
double CFrag::Egain(double thick)
{
  if (thick > 0.)
    recon->energy = loss->getEin(recon->energy,thick/cos(recon->theta));

  recon->getVelocity();

  return recon->energy;
}

//***********************************************
//include multiple scattering
/**
 * Monte Carlo choice of small angle scattering due to passage through the target
 \param fractionalThick is the fractional thick of the target through which the particle passed
 */
void CFrag::multiScat(double fractionalThick)
{
  if(fractionalThick == 0.) return;
  multScat->scatter(real->energy,fractionalThick,real->theta,real->phi);
  double thetaNew = multScat->getThetaNew();
  double phiNew = multScat->getPhiNew();

  real->theta = thetaNew;
  real->phi = phiNew;
  real->getVelocity();
}

//*********************
/**
 * accounts for multiscattering and energy loss in the target
 \param dthick is thickness of target though the particle must pass (mg/cm2)
 \param thickness is total target thickness (mg/cm2)
 */
bool CFrag::targetInteraction(double dthick, double thickness)
{
  bool stopped = 0;
  if(dthick == 0.||thickness==0) return stopped;
  double thick = dthick/cos(real->theta);
  Eloss(thick);
  if(real->energy<=0.)
  {
    stopped = 1;
    return stopped;
  }
  multiScat(thick/thickness);
  return stopped;
}

int CFrag::protonHit(double xTarget,double yTarget,double p_range, double p_straggle)
{
  double x = dist*tan(real->theta)*cos(real->phi) + xTarget;
  double y = dist*tan(real->theta)*sin(real->phi) + yTarget;
  double r = sqrt(pow(x,2) + pow(y,2));

  //proton did not hit on S4
  if (r < r_inner) return 0;
  if (r > r_outer) return 0;

  bool inner = false;
  if (r<r_CsI) inner = true;

  //ring
  iring = (r-r_inner)/(r_outer-r_inner)*Nrings;
  //pie
  double phi = atan2(y,x);
  if (phi < 0.) phi += 2.*pi;
  ipie = phi/2./pi*Npies;
  //icsi
  if(inner) icsi = phi/2./pi*4.+16;
  else icsi = phi/2./pi*16;

  //hit position at the back/stop point of CsI
  x = (dist+p_range*cos(real->theta))*tan(real->theta)*cos(real->phi) 
    + xTarget + ran.Gaus(0.,p_straggle/sqrt(2.));
  y = (dist+p_range*cos(real->theta))*tan(real->theta)*sin(real->phi) 
    + yTarget + ran.Gaus(0.,p_straggle/sqrt(2.));
  r = sqrt(pow(x,2) + pow(y,2));

  if (r < r_inner) return 0.;
  if (r > r_outer) return 0.;

  bool innerBack = false;
  if(r<r_CsI) innerBack = true;
  //proton can have inner-outer energy sharing
  if(inner != innerBack) return 0.;

  int icsiBack;
  if(innerBack) icsiBack = phi/2./pi*4.+16;
  else icsiBack = phi/2./pi*16;
  //proton can not have energy sharing on neighboring CsIs
  if(icsi != icsiBack) return 0.;

  r = ((double)iring + ran.Rndm())*(r_outer-r_inner)/Nrings+r_inner;
  recon->theta = atan(r/dist);
  recon->phi = ((double)ipie + ran.Rndm())*2.*pi/Npies;
  recon->energy = real->energy + sqrt(real->energy)*CsI_res*ran.Gaus(0.,1.);

  recon->getVelocity();

  return 1;
}

// For heavy fragment in planar detector
int CFrag::hit4(double xtarget,double ytarget)
{
  bool is_hit = Plane->hit(real->theta,real->phi,xtarget,ytarget) ;

  if (is_hit) 
  {
    recon->theta = Plane->thetaHit;
    recon->phi = Plane->phiHit;

    recon->energy = real->energy;
    // recon->energy = real->energy + real->energy*0.02*(2*ran.Rndm() - 1);
    //recon->energy = real->energy + sqrt(real->energy)*CsI_res*ran.Gaus(0.,1.);
    //recon->energy = real->energy + real->energy*CsI_res*ran.Gaus(0.,1.);

    recon->getVelocity();
  }

  return is_hit;
}

int CFrag::coreHit(double xTar,double yTar)
{
  double r = dist*tan(real->theta);
  double x = r*cos(real->phi)+xTar;
  double y = r*sin(real->phi)+yTar;
  if(x*x+y*y>r_S4Hole*r_S4Hole) return 0; //exceed the 10mm hole of S4 detector

  r = (dist+dist2Si)*tan(real->theta);
  x = r*cos(real->phi)+xTar;
  y = r*sin(real->phi)+yTar;
  if(x*x+y*y>r_FiberPlate*r_FiberPlate) return 0; //exceed the 16.71mm hole of fiber plate
  if(fabs(x)>xmax || fabs(y)>ymax) return 0; //exceed the fiber width

  double det_x = x+fiber_res*ran.Gaus(0,1);
  double det_y = y+fiber_res*ran.Gaus(0,1);
  double det_r = sqrt(pow(det_x,2)+pow(det_y,2));
  recon->theta = atan2(det_r,dist+dist2Si);
  recon->phi = atan2(det_y,det_x);
  
  recon->energy = real->energy;
  recon->getVelocity();

  return 1; 
}
