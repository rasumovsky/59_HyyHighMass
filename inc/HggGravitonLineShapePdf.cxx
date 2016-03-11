// class declaration include file below retrieved from workspace code storage
// class declaration include file below retrieved from workspace code storage
#include "HggGravitonLineShapePdf.h"

ClassImp(HggGravitonLineShapePdf);

//_____________________________________________________________________________
HggGravitonLineShapePdf:: HggGravitonLineShapePdf() {
}

//_____________________________________________________________________________
HggGravitonLineShapePdf::HggGravitonLineShapePdf(const char *name, const char *title,
					   RooAbsReal& _x, RooAbsReal& _mG, RooAbsReal& _GkM, int _cme) :
  RooAbsPdf(name, title),
  x("x", "Dependent", this, _x),
  mG("mG", "Mass", this, _mG),
  GkM("GkM", "GkM", this, _GkM),
  cme(_cme)
{
}


//_____________________________________________________________________________
HggGravitonLineShapePdf::HggGravitonLineShapePdf(const HggGravitonLineShapePdf& other, const char* name) :
  RooAbsPdf(other, name),
  x("x", this, other.x),
  mG("mG", this, other.mG),
  GkM("GkM", this, other.GkM),
  cme(other.cme)
{
}

//_____________________________________________________________________________
Double_t HggGravitonLineShapePdf::evaluate() const {
  // pre-defined variables
  Double_t m2Res=mG*mG;
  Double_t sHat=x*x;
  
  if(cme==13){
    Double_t wRes=1.44*mG*GkM*GkM; // Width of the resonance
    Double_t GamMRat=wRes/mG;
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
  else if(cme==8){
    Double_t wRes = 1.43813*mG*GkM*GkM;
    Double_t GamMRat=wRes/mG;
    Double_t weightBW = 1 / ((sHat - m2Res)*(sHat - m2Res) + (sHat * GamMRat)*(sHat * GamMRat));
    // Parameters derived from Bernstein polynomial, in the range from 200 GeV to 3.5 TeV.
    // p0 fixed to 0
    // p1        = 0.0607621    +/-  0.015913  (limited)
    // p2        = 0.327415     +/-  0.080904  (limited)
    // p3        = 0.274664     +/-  0.0737216 (limited)
    // p4        = 0.134802     +/-  0.0385843 (limited)
    // p5        = 0.0690092    +/-  0.0198695 (limited)
    // p6        = 0.039084     +/-  0.00980382        (limited)
    
    double p[7]={0, 0.0607621, 0.327415, 0.274664, 0.134802, 0.0690092, 0.039084};
    Double_t weightPL = 1;
    if(x<200) weightPL=bernstein(6, p, 200, 200, 3500);
    else if(x>3500) weightPL=bernstein(6, p, 3500, 200, 3500);
    else weightPL=bernstein(6, p, x, 200, 3500);
    Double_t w=weightBW*weightPL;          
    return w;
  }
  else{
    std::cerr<<"Uknown center of mass energy "<<cme<<std::endl;
    abort();
  }
}

//_____________________________________________________________________________
Double_t HggGravitonLineShapePdf::bernstein(int degree, double* p, double mgg, double xmin, double xmax) const 
{
  Double_t mggx = (mgg - xmin) / (xmax - xmin); // rescale to [0,1]

  if(degree == 0) {

    return p[0];

  } else if(degree == 1) {

    Double_t a0 = p[0]; // c0
    Double_t a1 = p[1] - a0; // c1 - c0
    return a1 * mggx + a0;

  } else if(degree == 2) {

    Double_t a0 = p[0]; // c0
    Double_t a1 = 2 * (p[1] - a0); // 2 * (c1 - c0)
    Double_t a2 = p[2] - a1 - a0; // c0 - 2 * c1 + c2
    return (a2 * mggx + a1) * mggx + a0;

  } else if(degree > 2) {

    Double_t t = mggx;
    Double_t s = 1 - mggx;

    Double_t result = p[0] * s;    
    for(Int_t i = 1; i < degree; i++) {
      result = (result + t * TMath::Binomial(degree, i) * (p[i])) * s;
      t *= mggx;
    }
    result += t * p[degree]; 

    return result;
  }

  // in case list of arguments passed is empty
  return TMath::SignalingNaN();
}
