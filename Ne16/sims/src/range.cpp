#include "range.h"

CRange::CRange(string filename)
{
  string ipfName = "lossfile/" + filename + ".range";
  cout << "Load range file: " << ipfName << endl;

  ifstream read(ipfName);
  if(!read.is_open())
  {
    cerr << "Could not open input file: " << ipfName << endl;
  }
  string dum;
  getline(read,dum);
  getline(read,dum);

  read >> N;
  x = new double[N];
  y = new double[N];
  z = new double[N];

  for(int i=0;i<N;i++)
  {
    read >> x[i] >> y[i] >> z[i];
  }
  read.close();
}

CRange::~CRange()
{
  delete []x;
  delete []y;
  delete []z;
}

double CRange::getRange(double E)
{
  if(E<x[0])
  {
    cout << "energy below tabulation in range" << endl;
    return 0;
  }

  if(E>x[N-1])
  {
    cout << "energy above tabulation in range" << endl;
    return 0;
  }

  int i=1;
  for(;;)
  {
    if(E<=x[i]) break;
    i++;
  }

  straggle = (z[i]-z[i-1])/(x[i]-x[i-1])*(E-x[i-1])+z[i-1];

  return (y[i]-y[i-1])/(x[i]-x[i-1])*(E-x[i-1])+y[i-1];
}

double CRange::getLateralStraggle()
{
  return straggle;
}
