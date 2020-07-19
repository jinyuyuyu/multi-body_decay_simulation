#include "decay.h"

double const CDecay::EkTot8B = 2.842; // energy released in 8B IAS decay to 6LI IAS
double const CDecay::gamma8B = 0.00075;

/**
 * Constructor
 \param part10 is a pointer to particle 1
 \param part20 is a pointer to particle 2
 \param part30 is a pointer to particle 3

*/
CDecay::CDecay(int nFrags0,CFrag **frag0,bool einstein0)
{

  ran.SetSeed(1);
  einstein = einstein0;
  nFrags = nFrags0;
  frag6Be = new CFrame (6,einstein);

  for(int i=0;i<nFrags;i++)
  {
    frag[i] = frag0[i];
    real[i] = frag0[i]->real;
    recon[i] = frag0[i]->recon;
  }

  sumA = 0.;
  for (int i=0;i<nFrags;i++) sumA += real[i]->mass;
  sumA /= amu;
  cout << "nFrags=" << nFrags << " sumA=" << sumA << endl;

  plfRecon = new CFrame(sumA,einstein);

  for (int i=0;i<nFrags;i++) 
  {
    partCM[i] = new CFrame(real[i]->mass/amu,einstein);
  }

  //string name("6be_p2_new.dat");
  //Moscow = new moscow(name);

  TLorentzVector W(0,0,0,mass_19Mg*amu-0.03303);
  double masses[3] = {real[0]->mass,real[1]->mass,real[2]->mass};
  genPhaseSpace.SetDecay(W,3,masses);
}

//*********************************************************
/**
 * Destructor
 */
CDecay::~CDecay()
{
  for (int i=0;i<3;i++)
  {
    delete partCM[i];
  }
  delete plfRecon;
  //delete Moscow;
}
//**********************************************************
/**
 * returns the reconstructed kinetical energy of the fragmens in their
 * center-of-mass frame using the real fragment velocities.
 */
double CDecay::getErelReal()
{
  return getErel(real);
}
//********************************************************
//*********************************************************
/**
 * returns the reconstructed kinetical energy of the fragmens in their
 * center-of-mass frame using the reconstructed or detected
 *  fragment velocities. 
 */
double CDecay::getErelRecon()
{
  return getErel(recon);
}
//*********************************************************
double CDecay::getErel(CFrame** part)
{
  if (einstein) return getErelRel(part);
  else
    return getErelNewton(part);
}
//**********************************************************
/** 
 * find the relative kinetic energy of the fragments in their
 * center-of-mass frame. Non-relativistic version
 \param part is a pointer to the fragments velocity vectors (real or reconstructed)
 */

double CDecay::getErelNewton(CFrame** part)
{

  for (int i=0;i<3;i++) 
  {
    plfRecon->v[i] = 0.;
    for (int j=0;j<5;j++)
      plfRecon->v[i] += part[j]->v[i]*part[j]->mass;
    plfRecon->v[i] /= sumA;  
  }

  plfRecon->getEnergy();

  ErelRecon = 0.;

  for (int j=0;j<5;j++)
  {
    partCM[j]->velocity = 0.;
    for (int i=0;i<3;i++)
    {
      partCM[j]->v[i] = part[j]->v[i] - plfRecon->v[i];
      partCM[j]->velocity += pow(partCM[j]->v[i],2);
    }
    partCM[j]->energy = real[j]->mass/2.*partCM[j]->velocity/pow(.9784,2);
    partCM[j]->velocity = sqrt(partCM[j]->velocity);
    ErelRecon += partCM[j]->energy;
  }

  return ErelRecon;

}

//**********************************************************
/** 
 * find the relative kinetic energy of the fragments in their
 * center-of-mass frame. Relativistic version
 \param part is a pointer to the fragments velocity vectors (real or reconstructed)
 */
double CDecay::getErelRel(CFrame **frame)
{
  for (int i=0;i<nFrags;i++) 
  {
    plfRecon->pc[i] = 0.;
    for (int j=0;j<3;j++) plfRecon->pc[i] += frame[j]->pc[i];
    //for (int j=0;j<3;j++) cout << frame[j]->pc[i] << "\t";
    //cout << plfRecon->pc[i] << endl;
  }
  plfRecon->totEnergy = 0;
  for (int i=0;i<3;i++) plfRecon->totEnergy += frame[i]->totEnergy;
  plfRecon->getVelocityFromMom();
  //cout << "plfRecon->totEnergy=" << plfRecon->totEnergy << endl;

  for (int j=0;j<nFrags;j++)
    for (int i=0;i<3;i++) partCM[j]->v[i] = frame[j]->v[i];

  ErelRecon = 0.;
  double vcm[3] = {-plfRecon->v[0],-plfRecon->v[1],-plfRecon->v[2]};
  for (int j=0;j<nFrags;j++)
  {
    partCM[j]->transformVelocity(vcm);
    ErelRecon += partCM[j]->energy;
    //cout << partCM[j]->energy << "\t";
  }
  //cout << "plfRecon->pcTot=" << plfRecon->pcTot << endl;
  //cout << "plfRecon->mass=" << plfRecon->mass << endl;
  //cout << "plfRecon->totEnergy=" << plfRecon->totEnergy << endl;
  //cout << "plfRecon->v" << " " << plfRecon->v[0] << " " << plfRecon->v[1] << " " << plfRecon->v[2] << endl;
  //cout << endl;

  //cout << "ErelRecon=" << ErelRecon << endl;
  return ErelRecon;

  //p-p relative
  int ii = 0;
  for (int i=0;i<1;i++)
  {
    for (int j=i+1;j<2;j++)
    {
      for (int k=0;k<3;k++) 
      {
        plfRecon->pc[k] = frame[i]->pc[k]+frame[j]->pc[k];
      }
      plfRecon->totEnergy = frame[i]->totEnergy
        + frame[j]->totEnergy;
      plfRecon->getVelocityFromMom();
      for (int k=0;k<3;k++) 
      {
        partCM[i]->v[k] = frame[i]->v[k];   
        partCM[j]->v[k] = frame[j]->v[k];   
      }
      ppRel[ii] = 0.;

      partCM[i]->transformVelocity(plfRecon->v);
      ppRel[ii] += partCM[i]->getEnergy();

      partCM[j]->transformVelocity(plfRecon->v);
      ppRel[ii] += partCM[j]->getEnergy();
      ii++;
    }
  }


  for (int i=0;i<3;i++) 
  {
    plfRecon->pc[i] = 0.;
    for (int j=0;j<4;j++)plfRecon->pc[i] += frame[j]->pc[i];
  }
  plfRecon->totEnergy = 0;
  for (int j=0;j<2;j++) plfRecon->totEnergy += frame[j]->totEnergy;
  plfRecon->getVelocityFromMom();


  for (int j=0;j<2;j++)
    for (int i=0;i<3;i++) partCM[j]->v[i] = frame[j]->v[i];

  double ErelReconA = 0.;
  for (int j=0;j<2;j++)
  {
    partCM[j]->transformVelocity(plfRecon->v);
    ErelReconA += partCM[j]->getEnergy();
  }

  aRatio = ErelReconA/ErelRecon;
}

double CDecay::getErelRel2(CFrame **frame)
{
  double Pm = 0;
  for (int i=0;i<3;i++) 
  {
    plfRecon->pc[i] = 0.;
    for (int j=0;j<nFrags;j++) plfRecon->pc[i] += frame[j]->pc[i];
    Pm += TMath::Power(plfRecon->pc[i],2);
    //for (int j=0;j<nFrags;j++) cout << frame[j]->pc[i] << "\t";
    //cout << plfRecon->pc[i] << endl;
  }
  plfRecon->totEnergy = 0;
  for (int i=0;i<nFrags;i++) plfRecon->totEnergy += frame[i]->totEnergy;
  //cout << "plfRecon->totEnergy=" << plfRecon->totEnergy << endl;

  double m = TMath::Power(plfRecon->totEnergy,2) - Pm;
  m = TMath::Sqrt(m);
  ErelRecon = m - sumA*amu;

  //cout << "ErelRecon=" << ErelRecon << endl;
  return ErelRecon;
}
//*************************************************************
/**
 * The momentum of the five fragments are chosen randomly from the 
 * allowable phase space.
 */
void CDecay::ModeMicroCanonical(double EkTot0)
{
  double EkTot = EkTot0;
  //for(;;)
  //{
  //  EkTot = ran.BreitWigner(EkTot8B,gamma8B);
  //  if( fabs(EkTot- EkTot8B) <  2.*gamma8B && EkTot > 0.) break;
  //}
  //micro2(nFrags,real,EkTot);
  micro(nFrags,real,EkTot);
}

//*****************************
void CDecay::micro(int N,CFrame **frame,double Ektot)
{
  valarray <double> vcm(3);
  for (int i=0;i<N;i++)
  {
    frame[i]->v[0] = ran.Gaus(0.,1.)/frame[i]->mass;
    frame[i]->v[1] = ran.Gaus(0.,1.)/frame[i]->mass;
    frame[i]->v[2] = ran.Gaus(0.,1.)/frame[i]->mass;

    for (int j=0;j<3;j++) vcm[j] += frame[i]->v[j]*frame[i]->mass;
  }

  vcm = vcm/sumA/amu;

  double testTotal= 0.;
  for (int i=0;i<N;i++)
  {
    frame[i]->velocity = 0.;
    for (int j=0;j<3;j++)
    {
      frame[i]->v[j] -= vcm[j];
      frame[i]->velocity += pow(frame[i]->v[j],2);
    }
    //frame[i]->energy = frame[i]->mass/2.*frame[i]->velocity/pow(.9784,2);
    frame[i]->energy = frame[i]->mass/2.*frame[i]->velocity/c/c;
    frame[i]->velocity = sqrt(frame[i]->velocity);
    testTotal += frame[i]->energy;
  }
  double ratio = sqrt(Ektot/testTotal);
  for (int i=0;i<N;i++)
  {
    frame[i]->velocity *= ratio;
    for (int j=0;j<3;j++) frame[i]->v[j] *= ratio;
    frame[i]->getEnergy();
  }

  //for(int i=0;i<N;i++)
  //{
  //  cout << "N=" << i << "  mass/amu=" << frame[i]->mass/amu << "  " << frame[i]->v[0] << "  " <<  frame[i]->v[1] << "  " << frame[i]->v[2] << endl;
  //  cout << "N=" << i << "  mass/amu=" << frame[i]->mass/amu << "  " << frame[i]->pc[0] << "  " <<  frame[i]->pc[1] << "  " << frame[i]->pc[2] << endl;
  //  cout << "N=" << i << "  mass/amu=" << frame[i]->mass/amu << "  " << frame[i]->totEnergy << endl;
  //}
  //cout << endl;
}

double CDecay::micro2(int N,CFrame **frame,double Ektot) {
  double weight = genPhaseSpace.Generate();
  TLorentzVector *pFrame[N];
  for(int np=0;np<N;np++) {
    pFrame[np] = genPhaseSpace.GetDecay(np);
    frame[np]->pc[0] = pFrame[np]->Px();
    frame[np]->pc[1] = pFrame[np]->Py();
    frame[np]->pc[2] = pFrame[np]->Pz();
    //cout << "np=" << np << " pc[0]=" << frame[np]->pc[0] << " pc[1]=" << frame[np]->pc[1] << " pc[2]=" << frame[np]->pc[2] << endl;
    frame[np]->getVelocityFromMom0();
  }
  

  return weight;
}

//*****************************
void CDecay::ModeComplex()
{
  double EkTot = 5.3;
  double Ek16Ne = 1.401;
  double Ek1 = (EkTot - Ek16Ne)/2.;
  //double Ek2 = EkTot - Ek16Ne - Ek1;


  double M1 = 1.;
  double M2 = 17.;
  double mu = M1*M2/(M1+M2);
  double vrel = sqrt(2.*Ek1/mu/amu)*c;
  double v1 = M2/(M1+M2)*vrel;
  double v2 = vrel - v1;

  double theta = acos(1.-2.*ran.Rndm());
  double phi = 2.*pi*ran.Rndm();


  real[0]->v[0] = v1*sin(theta)*cos(phi);
  real[0]->v[1] = v1*sin(theta)*sin(phi);
  real[0]->v[2] = v1*cos(theta);


  double recoil[3];

  for (int k=0;k<3;k++) recoil[k] = v2/v1*real[0]->v[k];


  //next decay
  M1 = 1.;
  M2 = 16.;
  mu = M1*M2/(M1+M2);
  vrel = sqrt(2.*Ek1/mu/amu)*c;
  v1 = M2/(M1+M2)*vrel;
  v2 = vrel - v1;

  theta = acos(1.-2.*ran.Rndm());
  phi = 2.*pi*ran.Rndm();

  real[1]->v[0] = v1*sin(theta)*cos(phi);
  real[1]->v[1] = v1*sin(theta)*sin(phi);
  real[1]->v[2] = v1*cos(theta);


  double recoil2[3];

  for (int k=0;k<3;k++) recoil2[k] = v2/v1*real[1]->v[k];


  for (int k=0;k<3;k++)
  {
    real[1]->v[k] += recoil[k];
    recoil2[k] += recoil[k];
  }

  //2p decay
  micro(3,&real[2],Ek16Ne);

  for (int k=0;k<3;k++)
  {
    real[2]->v[k] += recoil2[k];
    real[3]->v[k] += recoil2[k];
    real[4]->v[k] += recoil2[k];
  }

}


//**********************************************************
/**
 * Simulates the decay of 8C as two sequential prompt two-proton decays
 * passing through the 6Be ground state. The prompt decays are each 
 * treated as sampling available phase space.
 */
void CDecay::ModeMoscow()
{

  //avergae total kinetic energy for 6Be decay
  double const EkTot6Be = 1.3711;
  // width 
  double const gamma6Be = .092;

  //choose total kinetic 
  double EkTot;
  if (gamma6Be == 0.)EkTot= EkTot6Be;
  else
  {
    for (;;)
    {
      EkTot = ran.BreitWigner(EkTot6Be,gamma6Be);
      if( fabs(EkTot- EkTot6Be) <  2.*gamma6Be) break;
    }
  }


  CFrame ** Frag;
  Frag = new CFrame* [3];
  Frag[0] = new CFrame(1.,einstein);
  Frag[1] = new CFrame(1.,einstein);
  Frag[2] = new CFrame(4.,einstein);


  Moscow->getEvent(EkTot,Frag[2],Frag[1],Frag[0]);

  //look at correlations in Jacobi coordinates


  for (int i=0;i<3;i++)
  {
    real[0]->v[i] = Frag[0]->v[i];
    real[1]->v[i] = Frag[1]->v[i];
    real[2]->v[i] = Frag[2]->v[i];
  } 
  getJacobiPrimary();
}



//********************************************************************
double CDecay::getEk3body(CFrame* frag1, CFrame* frag2, CFrame* frag3)
{
  if (einstein) return getEk3bodyRel(frag1,frag2,frag3);
  else return getEk3bodyNewton(frag1,frag2,frag3);
}

//*********************************************************************
/**
 * returns the relative kinetic energy of three fragments in their
 * center-of-mass frame. Non-relativistic version
 \param frag1 is a pointer to fragment 1 
 \param frag2 is a pointer to fragment 2 
 \param frag3 is a pointer to fragment 3 
 */
double CDecay::getEk3bodyNewton(CFrame* frag1, CFrame* frag2, CFrame* frag3)
{
  double totMass = frag1->mass + frag2->mass + frag3->mass;

  CFrame frag1CM(frag1->mass,einstein);
  CFrame frag2CM(frag2->mass,einstein);
  CFrame frag3CM(frag3->mass,einstein);

  double vcm[3] = {0.};
  double ek = 0;
  for (int i=0;i<3;i++)
  {
    vcm[i] += frag1->v[i]*frag1->mass;
    vcm[i] += frag2->v[i]*frag2->mass;
    vcm[i] += frag3->v[i]*frag3->mass;
    vcm[i] /= totMass;

    frag1CM.v[i] = frag1->v[i] - vcm[i];
    frag2CM.v[i] = frag2->v[i] - vcm[i];
    frag3CM.v[i] = frag3->v[i] - vcm[i];

  }
  ek += frag1CM.getEnergy();
  ek += frag2CM.getEnergy();
  ek += frag3CM.getEnergy();

  return ek;
}
//*********************************************************
/**
 * returns the relative kinetic energy of three fragments in their
 * center-of-mass frame. Relativistic version
 \param frag1 is a pointer to fragment 1 
 \param frag2 is a pointer to fragment 2 
 \param frag3 is a pointer to fragment 3 
 */
double CDecay::getEk3bodyRel(CFrame* frag1, CFrame* frag2, CFrame* frag3)
{
  double totMass = frag1->mass + frag2->mass + frag3->mass;
  CFrame frame(totMass,einstein);

  for (int i=0;i<3;i++) 
  {
    frame.pc[i] = frag1->pc[i] + frag2->pc[i] + frag3->pc[i];
  }
  frame.totEnergy = frag1->totEnergy + frag2->totEnergy + frag3->totEnergy;
  frame.getVelocityFromMom();

  for (int i=0;i<3;i++)
  {
    frag6Be->pc[i] = frame.pc[i];
    frag6Be->v[i] = frame.pc[i]/frame.pcTot*frame.velocity;
  }
  frag6Be->totEnergy = frame.totEnergy;

  CFrame p1(frag1->mass,einstein);
  CFrame p2(frag2->mass,einstein);
  CFrame p3(frag3->mass,einstein);
  for (int i=0;i<3;i++)
  {
    p1.v[i] = frag1->v[i];
    p2.v[i] = frag2->v[i];
    p3.v[i] = frag3->v[i];

  }

  p1.transformVelocity(frame.v);
  p2.transformVelocity(frame.v);
  p3.transformVelocity(frame.v);


  ErelRecon = p1.getEnergy() + p2.getEnergy() + p3.getEnergy();

  return ErelRecon;

}


//*********************************************************
/**
 * returns the relative kinetic energy of three fragments in their
 * center-of-mass frame. Relativistic version
 \param frag1 is a pointer to fragment 1 
 \param frag2 is a pointer to fragment 2 
 \param frag3 is a pointer to fragment 3 
 */
double CDecay::getEk4bodyRel(CFrame* frag1, CFrame* frag2, CFrame* frag3, CFrame* frag4)
{
  double totMass = frag1->mass + frag2->mass + frag3->mass + frag4->mass;
  CFrame frame(totMass,einstein);

  for (int i=0;i<3;i++) 
  {
    frame.pc[i] = frag1->pc[i] + frag2->pc[i] + frag3->pc[i] + frag4->pc[i];
  }
  frame.totEnergy = frag1->totEnergy + frag2->totEnergy + frag3->totEnergy
    + frag4->totEnergy;
  frame.getVelocityFromMom();

  for (int i=0;i<3;i++)
  {
    frag6Be->pc[i] = frame.pc[i];
    frag6Be->v[i] = frame.pc[i]/frame.pcTot*frame.velocity;
  }
  frag6Be->totEnergy = frame.totEnergy;

  CFrame p1(frag1->mass,einstein);
  CFrame p2(frag2->mass,einstein);
  CFrame p3(frag3->mass,einstein);
  CFrame p4(frag4->mass,einstein);
  for (int i=0;i<3;i++)
  {
    p1.v[i] = frag1->v[i];
    p2.v[i] = frag2->v[i];
    p3.v[i] = frag3->v[i];
    p4.v[i] = frag4->v[i];

  }

  p1.transformVelocity(frame.v);
  p2.transformVelocity(frame.v);
  p3.transformVelocity(frame.v);
  p4.transformVelocity(frame.v);


  ErelRecon = p1.getEnergy() + p2.getEnergy() + p3.getEnergy() + p4.getEnergy();

  return ErelRecon;

}





//*******************
double CDecay::getEk3()
{
  return getEk3bodyRel(recon[0],recon[1],recon[2]);
}
//************************************************************


void CDecay::getJacobi(CFrame**part,bool com)
{
  double totMass = part[0]->mass + part[1]->mass + part[2]->mass;

  CFrame p1(part[0]->mass,einstein);
  CFrame p2(part[1]->mass,einstein);
  CFrame p3(part[2]->mass,einstein);
  double Etot;
  for (int i=0;i<3;i++)
  {
    p1.v[i] = part[0]->v[i];
    p2.v[i] = part[1]->v[i];
    p3.v[i] = part[2]->v[i];

  }

  if (com)
  {
    CFrame frame(totMass,einstein);


    for (int i=0;i<3;i++) 
    {
      frame.pc[i] = part[0]->pc[i] + part[1]->pc[i] + part[2]->pc[i];
    }
    frame.totEnergy = part[0]->totEnergy + part[1]->totEnergy + 
      part[2]->totEnergy;
    frame.getVelocityFromMom();



    p1.transformVelocity(frame.v);
    p2.transformVelocity(frame.v);
    p3.transformVelocity(frame.v);

    Etot = p1.getEnergy() + p2.getEnergy() + p3.getEnergy();
  }
  else
  {
    Etot = 0.;
    for (int i=0;i<3;i++) Etot += 0.5*p1.mass*pow(p1.v[i]/.9784,2)
      + 0.5*p2.mass*pow(p2.v[i]/.9784,2) + 0.5*p3.mass*pow(p3.v[i]/.9784,2);
  }

  //p-p relative velocity can use newton
  double vrel[3];
  double vrelative = 0.;
  double dot = 0.;
  double v3 = 0.;
  for (int i=0;i<3;i++)
  {
    vrel[i] = p1.v[i] - p2.v[i];
    vrelative += pow(vrel[i],2);
    dot += vrel[i]*p3.v[i];
    v3 += pow(p3.v[i],2);
  }
  vrelative = sqrt(vrelative);
  v3 = sqrt(v3);

  double Epp = 0.5*0.5*pow(vrelative/.9784,2);
  x_T = Epp/Etot;

  CosTheta_T = dot/v3/vrelative;


  //6Be-p1
  vrelative = 0.;
  dot = 0.;
  double v2 = 0.;
  for (int i=0;i<3;i++)
  {
    vrel[i] = p1.v[i] - p3.v[i];
    vrelative += pow(vrel[i],2);
    dot += vrel[i]*p2.v[i];
    v2 += pow(p2.v[i],2);
  }
  vrelative = sqrt(vrelative);
  v2 = sqrt(v2);

  double Epa = 0.5*28./29.*pow(vrelative/.9784,2);
  x_Y[0] = Epa/Etot;

  CosTheta_Y[0] = dot/v2/vrelative;


  //6Be-p2
  vrelative = 0.;
  dot = 0.;
  double v1 = 0.;
  for (int i=0;i<3;i++)
  {
    vrel[i] = p2.v[i] - p3.v[i];
    vrelative += pow(vrel[i],2);
    dot += vrel[i]*p1.v[i];
    v1 += pow(p1.v[i],2);
  }
  vrelative = sqrt(vrelative);
  v1 = sqrt(v1);

  Epa = 0.5*28./29.*pow(vrelative/.9784,2);
  x_Y[1] = Epa/Etot;

  CosTheta_Y[1] = dot/v1/vrelative;
}


void CDecay::getJacobiSecondary()
{
  getJacobi(recon,(bool)1);
}

void CDecay::getJacobiPrimary()
{
  getJacobi(real,(bool)0);
} 



