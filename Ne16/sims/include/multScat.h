#ifndef MULTSCAT_
#define MULTSCAT_

#include <iostream>
#include <TMath.h>
#include <TRandom3.h>
#include "constants.h"

using namespace std;

class CMultScat
{
  public:
    TRandom3 ran;

    double a0;
    double N;
    double Ztar;
    double Atar;
    double dtar;
    
    double Zpro;
    double totalThick;
    double a;
    double totalTau;
    double factor;

    double thetaNew;
    double phiNew;

    CMultScat(double, double);
    void scatter(double, double, double, double);
    double getThetaNew();
    double getPhiNew();
};

#endif
