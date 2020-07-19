#ifndef PLANE_
#define PLANE_

#include <TRandom3.h>

class CPlane
{
  public:
    CPlane(double,double,double,double,double);
    int hit(double,double,double,double);

    TRandom3 ran;

    //factors of the planar detector
    double dist; //distance from target to Si
    double dist2Si; //distance from Si to target [mm]
    double xmax; //half width of fiber on x direction [mm]
    double ymax; //half width of fiber on y direction [mm]
    double res; //position resolution of fiber [mm]

    //detection output
    double det_x;
    double det_y;
    double det_theta;
    double det_phi;
};

#endif
