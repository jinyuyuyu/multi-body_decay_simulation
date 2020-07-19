#ifndef LOSS_
#define LOSS_

#include <iostream>
#include <fstream>
#include <string>
#include <TGraph.h>

using namespace std;

class CLoss
{
  public:
    CLoss(string,double);
    ~CLoss();

    double getEout(double,double);
    double getEin(double,double);

    TGraph *grRE;
    TGraph *grER;

    double ionmass;
    double *x;
    double *y;
};

#endif
