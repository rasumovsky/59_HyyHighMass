////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// ==================================================================         //
// Implemented by Hongtao Yang <Hongtao.Yang@cern.ch> on Feb. 6, 2016         //
// With theory inputs provided by Yee Yap <yee.yap@lpnhe.in2p3.fr>            //
// and Lydia Roos <lroos@lpnhe.in2p3.fr>                                      //
// Everything should be in GeV to be consistent with the implementation       //
// ==================================================================         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// class declaration include file below retrieved from workspace code storage
// class declaration include file below retrieved from workspace code storage
#include "HggScalarLineShapePdf.h"

ClassImp(HggScalarLineShapePdf);

//_____________________________________________________________________________
HggScalarLineShapePdf:: HggScalarLineShapePdf() {
}

//_____________________________________________________________________________
HggScalarLineShapePdf::HggScalarLineShapePdf(const char *name, const char *title,
					     RooAbsReal& x, RooAbsReal& mean, RooAbsReal& width, RooAbsReal& alpha, int cme, bool performTruncation) :
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _mean("mean", "Mass", this, mean),
  _width("width", "Width", this, width),
  _alpha("alpha", "Position of truncation", this, alpha),
  _cme(cme),
  _performTruncation(performTruncation)
{
  _nu=1;
  _sigma=2;
  _epsilon=1e-5;
}


//_____________________________________________________________________________
HggScalarLineShapePdf::HggScalarLineShapePdf(const HggScalarLineShapePdf& other, const char* name) :
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _mean("mean", this, other._mean),
  _width("width", this, other._width),
  _alpha("alpha", this, other._alpha),
  _cme(other._cme),
  _performTruncation(other._performTruncation),
  _nu(other._nu),
  _sigma(other._sigma),
  _epsilon(other._epsilon)
{
}

//_____________________________________________________________________________
// Matrix element correction
Double_t HggScalarLineShapePdf::ME(Double_t mgg) const {
  Double_t topMass=173.34;
  Double_t tau_Q=4*topMass*topMass/(mgg*mgg);
  std::complex<Double_t> f;

  if(tau_Q>=1){
    f=asin(1./sqrt(tau_Q));
    f=f*f;
  }
  else{
    f=log((1.+sqrt(1.-tau_Q))/(1.-sqrt(1.-tau_Q)))-std::complex<Double_t>(0, TMath::Pi());
    f=-0.25*f*f;
  }
  
  std::complex<Double_t> Bgg_1=1.5*tau_Q*(1.+(1.-tau_Q)*f);
  Double_t Bgg=abs(Bgg_1)*abs(Bgg_1)*mgg*mgg;

  return Bgg;
}

//_____________________________________________________________________________
// Parton luminosity correction
Double_t HggScalarLineShapePdf::PL(Double_t mgg) const {
  Double_t mgg_norm=0;

  if(_cme==13){
    // Using formula provided by Yee: RooGenericPdf("lumipdf", "lumipdf", " (1-(mgg/13000))^28 * (mgg/13000)^-3.3 * (5e-9 - 3e-8*(mgg/13000) + 3e-7*(mgg/13000)^2)", *mgg);
    // Double_t Lgg=pow(1-mgg_norm, 28)*pow(mgg_norm, -3.3)*(5e-9 - 3e-8*mgg_norm + 3e-7*mgg_norm*mgg_norm);
  
    // Update on Feb 10, 2016
    // RooGenericPdf lumipdf = new RooGenericPdf("lumipdf", "lumipdf", " (1-(mgg/13000)^(1./3))^8.95 * (mgg/13000)^-4.1  (-2.95e-10 + 1.32e-7*(mgg/13000) -1.7e-7*(mgg/13000)^2) ", mgg);
    mgg_norm=mgg/13000.;	// Normalize mgg by sqrt(s)
    return pow(1-pow(mgg_norm,1./3.), 8.95)*pow(mgg_norm, -4.1)*(-2.95e-10 + 1.32e-7*mgg_norm - 1.7e-7*mgg_norm*mgg_norm);
  }
  else if(_cme==8){
    // 8 TeV parameterization
    // RooGenericPdf lumipdf8TeV = new RooGenericPdf("lumipdf8TeV", "lumipdf8TeV", " (1-(mgg/8000)^(1./3))^10.55 * (mgg/8000)^-2.76 * (2.04e-6 - 1.9e-6*(mgg/8000) + 1.25e-5*(mgg/8000)^2 - 1.52e-5*(mgg/8000)^3)", mgg);
    mgg_norm=mgg/8000.;	// Normalize mgg by sqrt(s)
    return pow((1-pow(mgg_norm, 1./3.)), 10.55)*pow(mgg_norm, -2.76)*(2.04e-6 - 1.9e-6*mgg_norm + 1.25e-5*mgg_norm*mgg_norm - 1.52e-5*mgg_norm*mgg_norm*mgg_norm);
  }
  else{
    std::cerr<<"\tERROR: unknown center-of-mass-energy type "<<_cme
	     <<". Choose between 13 (TeV) and 8 (TeV). Aborting..."
	     <<std::endl;
    exit(-1);
  }
}
  
//_____________________________________________________________________________
// Matrix element times parton luminosity
Double_t HggScalarLineShapePdf::MEPL(Double_t mgg) const {
  return ME(mgg)*PL(mgg);
}

//_____________________________________________________________________________
// Relativistic Breit-Wigner
Double_t HggScalarLineShapePdf::BW(Double_t mgg) const {
  return (mgg*mgg*_width/_mean)/((mgg*mgg-_mean*_mean)*(mgg*mgg-_mean*_mean)+(mgg*mgg*_width/_mean)*(mgg*mgg*_width/_mean));
}

//_____________________________________________________________________________
Double_t HggScalarLineShapePdf::evaluate() const {
  Double_t BggLgg=MEPL(_x);
  Double_t RBW=BW(_x);

  if(_performTruncation){
    Double_t cutoff=_mean-_alpha*_width;
    if(_nu>0){			// Extrapolate ME*Lgg with 2nd order polynomial
      if(_x<cutoff){
	Double_t deriv1st=dMEPLdmgg(cutoff);
	Double_t deriv2nd=d2MEPLdmgg2(cutoff);
	
	Double_t a=deriv2nd*0.5;
	Double_t b=deriv1st;
	Double_t c=MEPL(cutoff);
	Double_t root1=(-b+sqrt(b*b-4*a*c))/(2*a)+cutoff, root2=(-b-sqrt(b*b-4*a*c))/(2*a)+cutoff;
	// 2nd polynomial is dangerous. We need to make sure it does not go negative in the domain of x
	// Otherwise we would rather do nothing...
	if(!((root1>_x.min()&&root1<cutoff)||(root2>_x.min()&&root2<cutoff)))
	  BggLgg=deriv2nd*0.5*(_x-cutoff)*(_x-cutoff)+deriv1st*(_x-cutoff)+MEPL(cutoff);
	// Double_t c=MEPL(cutoff);
	// Double_t a=deriv1st/c;
	// Double_t b=0.5*(deriv2nd/c-a*a);
	
	// BggLgg=c*exp(a*(x-cutoff)+b*(x-cutoff)*(x-cutoff));
      }
    }
    else{
      BggLgg=MEPL(_x)*ROOT::Math::gaussian_cdf(_x, _sigma*_width, cutoff-2*_sigma*_width);
      // else BggLgg*=1/pow(1+exp(-(_x-(_mean-_alpha*_width))/(_sigma*_width)),_nu); // Sigmoid truncation
    }
  }

  return BggLgg*RBW*_x;
  
}


Double_t HggScalarLineShapePdf::dMEPLdmgg(Double_t mgg) const
{
  Double_t h=_epsilon*(_x.max()-_x.min());
  Double_t xx;
  xx = mgg+h;     Double_t f1 = MEPL(xx);
  xx = mgg-h;     Double_t f2 = MEPL(xx);
   
  xx = mgg+h/2;   Double_t g1 = MEPL(xx);
  xx = mgg-h/2;   Double_t g2 = MEPL(xx);

  //compute the central differences
  Double_t h2    = 1/(2.*h);
  Double_t d0    = f1 - f2;
  Double_t d2    = g1 - g2;
  Double_t deriv = h2*(8*d2 - d0)/3.;

  return deriv;
}

Double_t HggScalarLineShapePdf::d2MEPLdmgg2(Double_t mgg) const
{
  Double_t h=_epsilon*(_x.max()-_x.min());
  Double_t xx;
  xx = mgg+h;     Double_t f1 = MEPL(xx);
  xx = mgg;       Double_t f2 = MEPL(xx);
  xx = mgg-h;     Double_t f3 = MEPL(xx);

  xx = mgg+h/2;   Double_t g1 = MEPL(xx);
  xx = mgg-h/2;   Double_t g3 = MEPL(xx);

  //compute the central differences
  Double_t hh    = 1/(h*h);
  Double_t d0    = f3 - 2*f2 + f1;
  Double_t d2    = 4*g3 - 8*f2 +4*g1;

  Double_t deriv = hh*(4*d2 - d0)/3.;
  return deriv;
}

