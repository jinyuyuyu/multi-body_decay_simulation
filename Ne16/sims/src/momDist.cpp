//get momentum distribution of neutron knockout reaction
//calculated by the program: MomDis
//the input file is the differential cross section

#include "momDist.h"

CMomDist::CMomDist()
{
  n = 221;
  xtr = new double[n];
  ytr = new double[n];
  xz = new double[n];
  yz = new double[n];

  int i=0;
  ifstream read("lossfile/zn02_allM.out");
  if(!read.is_open())
  {
    cerr << "Could not open input file for z MomDist" << endl;
  }
  while(read.good())
  {
    read >> xz[i] >> yz[i];
    if(read.good())
    {
      if(i>0) yz[i] += yz[i-1];
      i++;
    }
  }
  if(i!=n) cerr << "in momDist i=" << i << " n=" << n << endl;
  for(int j=0;j<n;j++) yz[j] /= yz[n-1];
  read.close();
  read.clear();

  i=0;
  read.open("lossfile/trn02_allM.out");
  if(!read.is_open())
  {
    cerr << "Could not open input file for z MomDist" << endl;
  }
  while(read.good())
  {
    read >> xtr[i] >> ytr[i];
    if(read.good())
    {
      if(i>0) ytr[i] += ytr[i-1];
      i++;
    }
  }
  if(i!=n) cerr << "in momDist i=" << i << " n=" << n << endl;
  for(int j=0;j<n;j++) ytr[j] /= ytr[n-1];
  read.close();
}

CMomDist::~CMomDist()
{
  delete []xtr;
  delete []ytr;
  delete []xz;
  delete []yz;
}

double CMomDist::getTransMom()
{
  double probtr = ran.Rndm();
  int i=0;
  for(;;)
  {
    if(ytr[i]>probtr) break;
    i++;
    if(i==n) break;
  }

  double transMom;
  if(i==0) transMom = xtr[i];
  else transMom = xtr[i-1]+(xtr[i]-xtr[i-1])*ran.Rndm();

  return transMom;
}

double CMomDist::getLongMom()
{
  double probz = ran.Rndm();
  int i=0;
  for(;;)
  {
    if(yz[i]>probz) break;
    i++;
    if(i==n) break;
  }

  double longMom;
  if(i==0) longMom = xz[i];
  else longMom = xz[i-1]+(xz[i]-xz[i-1])*ran.Rndm();

  return longMom;
}
