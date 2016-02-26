////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// ==================================================================         //
// Implemented by Hongtao Yang <Hongtao.Yang@cern.ch> on Feb. 6, 2016         //
// with theory inputs provided by Jan Stark <stark@in2p3.fr>                  //
// ==================================================================         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// class declaration include file below retrieved from workspace code storage
// class declaration include file below retrieved from workspace code storage
#include "HggGravitonLineShapePdf.h"

ClassImp(HggGravitonLineShapePdf);

//_____________________________________________________________________________
HggGravitonLineShapePdf:: HggGravitonLineShapePdf() {
}

//_____________________________________________________________________________
HggGravitonLineShapePdf::HggGravitonLineShapePdf(const char *name, const char *title,
					   RooAbsReal& _x, RooAbsReal& _mG, RooAbsReal& _GkM, bool isTopMassInfi) :
  RooAbsPdf(name, title),
  x("x", "Dependent", this, _x),
  mG("mG", "Mass", this, _mG),
  GkM("GkM", "GkM", this, _GkM)
{
}


//_____________________________________________________________________________
HggGravitonLineShapePdf::HggGravitonLineShapePdf(const HggGravitonLineShapePdf& other, const char* name) :
  RooAbsPdf(other, name),
  x("x", this, other.x),
  mG("mG", this, other.mG),
  GkM("GkM", this, other.GkM)
{
}


//_____________________________________________________________________________
Double_t HggGravitonLineShapePdf::evaluate() const {
  Double_t wRes=1.44*mG*GkM*GkM; // Width of the resonance
  // pre-defined variables
  Double_t m2Res=mG*mG;
  Double_t GamMRat=wRes/mG;
  Double_t sHat=x*x;
  // Breit-Wigner shape
  Double_t weightBW=1 / ((sHat - m2Res)*(sHat - m2Res) + (sHat * GamMRat)*(sHat * GamMRat));
  Double_t p0=0, p1=0, p2=0, p3=0, p4=0, p5=0, p6=0;
  if (x >= 50 && x < 2500){
    p0 = 1.06963e+01;
    p1 = -1.86534e-01;
    p2 = 9.00278e-04;
    p3 = 1.01576e-06;
    p4 = -1.29332e-09;
    p5 = 4.64752e-13;
    p6 = -5.68297e-17;
  } else if (x >= 2500 && x < 6000){
    p0 = -1.77843e+03;
    p1 = 3.70114e+00;
    p2 = -1.06644e-03;
    p3 = 3.77245e-08;
    p4 = 2.49922e-11;
    p5 = -3.96985e-15;
    p6 = 1.88504e-19;
  } else if (x >= 6000 && x < 8000){
    p0 = 1.43663e+03;
    p1 = 3.40467e-01;
    p2 = -7.44453e-05;
    p3 = -1.08578e-08;
    p4 = 6.74486e-13;
    p5 = 2.82952e-16;
    p6 = -2.20149e-20;
  } else if (x >= 8000 && x <= 10500){
    p0 = 1.21245e+03;
    p1 = -1.08389e-01;
    p2 = -1.02834e-05;
    p3 = 4.11376e-12;
    p4 = 8.55312e-14;
    p5 = 6.98307e-18;
    p6 = -6.52683e-22;
  }
  Double_t weightPL = p0+(p1*x)+(p2*(x*x))+(p3*(x*x*x))+(p4*(x*x*x*x))+(p5*(x*x*x*x*x))+(p6*(x*x*x*x*x*x));

  Double_t w=weightBW*weightPL;            
  return w;
}


