#include "plane.h"

CPlane::CPlane(double dist0,double dist2Si0,double xmax0,double ymax0,double res0)
{
  dist = dist0;
  dist2Si = dist2Si0;
  xmax = xmax0;
  ymax = ymax0;
  res = res0;
}

int CPlane::hit(double theta,double phi,double xTar,double yTar)
{
  double r = dist*tan(theta);
  double x = r*cos(phi)+xTar;
  double y = r*sin(phi)+yTar;
  if(x*x+y*y>5.*5.) return 0; //exceed the 10mm hole of S4 detector

  r = (dist+dist2Si)*tan(theta);
  x = r*cos(phi)+xTar;
  y = r*sin(phi)+yTar;
  if(x*x+y*y>8.35*8.35) return 0; //exceed the 16.71mm hole of fiber plate
  if(fabs(x)>xmax || fabs(y)>ymax) return 0; //exceed the fiber width

  det_x = x+res*ran.Gaus(0,1);
  det_y = y+res*ran.Gaus(0,1);
  double det_r = sqrt(pow(det_x,2)+pow(det_y,2));
  det_theta = atan2(det_r,dist+dist2Si);
  det_phi = atan2(det_y,det_x);

  return 1;
}
