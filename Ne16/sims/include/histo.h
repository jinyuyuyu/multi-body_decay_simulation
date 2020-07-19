#ifndef HISTO_
#define HISTO_

#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

using namespace std;

class histo
{
  protected:
    TFile *opf;

  public:
    histo();
    ~histo();
    void write();

    TH1I *hist_pcbeam;
    TH1I *hist_Ebeam;
    TH2I *hist_beamspot;
    TH1I *hist_xTar;
    TH1I *hist_yTar;

    TH1I *hist_knockout_vel;
    TH1I *hist_knockout_theta;
    TH1I *hist_knockout_phi;

    TH1I *hist_mscat_vel;
    TH1I *hist_mscat_theta;
    TH1I *hist_mscat_phi;
    TH1I *hist_mscat_pc;

    TH1I *hist_ErelReal;
    TH1I *hist_ErelRecon;

    TH1I *hist_p_range;

    TH1I *hist_p1Energy;
    TH1I *hist_p1Theta;
    TH1I *hist_p1Phi;
    TH1I *hist_p2Energy;
    TH1I *hist_p2Theta;
    TH1I *hist_p2Phi;
    TH1I *hist_coreEnergy;
    TH1I *hist_coreTheta;
    TH1I *hist_corePhi;
};

#endif
