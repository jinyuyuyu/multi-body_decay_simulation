#ifndef FRAG_
#define FRAG_

#include <string>
#include <TRandom3.h>
#include "loss.h"
#include "multScat.h"
#include "frame.h"
#include "fragment.h"
#include "plane_det.h"
#include "constants.h"

using namespace std;

/**
 *!\brief information about a single fragment and its interection with the detector
 *
 */
class CFrag : public fragment
{
  public:

    CFrag(double,double,string,double,double,double,bool);
    ~CFrag();

    int hit4(double xtarget, double ytarget);
    int protonHit(double xTarget, double yTarget, double p_range, double p_straggle);
    int coreHit(double,double);

    void AddVelocity(double*);
    double Eloss(double); // energy loss in target
    double Egain(double); // corrects for energy loss in target
    void multiScat(double);
    bool targetInteraction(double,double);

    TRandom3 ran;
    CLoss *loss;
    CMultScat *multScat;
    plane_det * Plane;

    bool einstein;
    CFrame *real;      //<!real particles energy, velocity, etc
    CFrame *recon;    //<!reconstructed properties

    //S4-CsI
    double dist;
    double r_inner; //S4 [mm]
    double r_outer; //S4 [mm]
    double Nrings; //S4
    double Npies; //S4
    double r_CsI; //CsI [mm]
    double CsI_res;  
    int ipie, iring, icsi;
    //fiber-S800
    double dist2Si;
    double r_S4Hole;
    double r_FiberPlate;
    double xmax;
    double ymax;
    double fiber_res;
};

#endif
