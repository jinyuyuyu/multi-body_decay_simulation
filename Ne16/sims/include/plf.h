#ifndef PLF_
#define PLF_

#include <iostream>
#include <TMath.h>
#include "frame.h"
#include "loss.h"
#include "momDist.h"
#include "multScat.h"
#include "constants.h"

using namespace std;

class CPlf
{
  public:
    CPlf(double,double,double,string,double,bool);
    ~CPlf();
    void getPlf(double);
    void getPlf1();
    void multiScat(double,double);
    void getBeam(double,double);

    CFrame *beam;
    CFrame *frame;
    CLoss *loss;
    CMomDist *momDist;
    CMultScat *multScat;

    bool einstein;
    double factor;
    double mass;
};

#endif
