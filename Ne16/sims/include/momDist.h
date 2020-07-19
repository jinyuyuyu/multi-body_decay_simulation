#ifndef MOMDIST_
#define MOMDIST_

#include <iostream>
#include <fstream>
#include <TRandom3.h>

using namespace std;

class CMomDist
{
  public:

    TRandom3 ran;
    int n;
    double *xtr;
    double *ytr;
    double *xz;
    double *yz;

    CMomDist();
    ~CMomDist();
    double getTransMom();
    double getLongMom();
};

#endif
