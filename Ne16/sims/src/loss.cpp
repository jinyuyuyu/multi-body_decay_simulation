#include "loss.h"

//using range-energy curve obtained from LISE
//range is in [micron], energy is in [MeV/u]

CLoss::CLoss(string filename,double ionmass0)
{
  ionmass = ionmass0;

  cout << "Load lossfile: " << filename << endl;
  string ipfName1 = "lossfile/" + filename + "1.txt";
  string ipfName2 = "lossfile/" + filename + "2.txt";
  string ipfName3 = "lossfile/" + filename + "3.txt";

  int nColumns = 10;
  double ip[nColumns];
  int nPoints_RE = 0;
  x = new double[2000];
  y = new double[2000];
  for(int i=0;i<2000;i++)
  {
    x[i] = 0;
    y[i] = 0;
  }

  //read in ipfName1
  ifstream read(ipfName1);
  if(!read.is_open())
  {
    cerr << "Could not open input file: " << ipfName1 << endl;
  }
  read.ignore(1000,'\n').ignore(1000,'\n');
  while(read.good())
  {
    for(int i=0;i<nColumns;i++)
    {
      read >> ip[i];
    }
    x[nPoints_RE] = ip[0];
    y[nPoints_RE] = ip[2];
    if(read.good()) nPoints_RE++;
  }
  if(read.eof())
  {
    cout << "nPoints_RE=" << nPoints_RE << endl;
  }
  read.close();
  read.clear();

  //read in ipfName2
  double x_demarcation = x[nPoints_RE-1];
  read.open(ipfName2);
  if(!read.is_open())
  {
    cerr << "Could not open input file: " << ipfName2 << endl;
  }
  read.ignore(1000,'\n');
  while(read.good())
  {
    for(int i=0;i<nColumns;i++)
    {
      read >> ip[i];
    }
    if(ip[0]>x_demarcation)
    {
      x[nPoints_RE] = ip[0];
      y[nPoints_RE] = ip[2];
      if(read.good()) nPoints_RE++;
    }
  }
  if(read.eof())
  {
    cout << "nPoints_RE=" << nPoints_RE << endl;
  }
  read.close();
  read.clear();

  //read in ipfName3
  x_demarcation = x[nPoints_RE-1];
  read.open(ipfName3);
  if(!read.is_open())
  {
    cerr << "Could not open input file: " << ipfName3 << endl;
  }
  read.ignore(1000,'\n');
  while(read.good())
  {
    for(int i=0;i<nColumns;i++)
    {
      read >> ip[i];
    }
    if(ip[0]>x_demarcation)
    {
      x[nPoints_RE] = ip[0];
      y[nPoints_RE] = ip[2];
      if(read.good()) nPoints_RE++;
    }
  }
  if(read.eof())
  {
    cout << "nPoints_RE=" << nPoints_RE << endl;
  }
  read.close();

  grRE = new TGraph(nPoints_RE,x,y);
  grER = new TGraph(nPoints_RE,y,x);
  string name_grRE = "grRE_" + filename;
  string name_grER = "grER_" + filename;
  grRE->SetNameTitle(name_grRE.c_str(),name_grRE.c_str());
  grER->SetNameTitle(name_grER.c_str(),name_grER.c_str());
}

CLoss::~CLoss()
{
  delete grRE;
  delete grER;
  delete []x;
  delete []y;
}

double CLoss::getEout(double energy,double thick)
{
  if(thick==0) return energy;

  double Eout = 0.;
  double Ein = energy / ionmass;

  //range data is not valid for energy <=0.05 MeV/u
  if(Ein<=0.05)
  {
    Eout = 0.;
    return Eout;
  }

  double range = grRE->Eval(Ein);
  if(range <= thick)
  {
    Eout = 0.;
    return Eout;
  }

  Eout = grER->Eval(range-thick);
  if(Eout<=0.05) Eout = 0.;

  return Eout*ionmass;
}

double CLoss::getEin(double energy,double thick)
{
  if(thick==0) return energy;

  double Ein = 0.;
  double Eout = energy/ionmass;
  double range = 0.;

  if(Eout<=0.05)
  {
    range = 0.;
  }
  else
  {
    range = grRE->Eval(Eout);
  }

  Ein = grER->Eval(range+thick);

  return Ein*ionmass;
}
