#ifndef RANGE_
#define RANGE_

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class CRange
{
  public:
    CRange(string filename);
    ~CRange();

    double getRange(double);
    double getLateralStraggle();

    double *x;
    double *y;
    double *z;
    int N;
    double straggle;
};

#endif
