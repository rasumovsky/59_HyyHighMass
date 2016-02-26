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

#ifndef ROOT_HggTruthLineShapePdf
#define ROOT_HggTruthLineShapePdf

#include <math.h>
#include "Math/ProbFuncMathCore.h"
#include "Riostream.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooFit.h"
#include "RooMath.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "Math/Functor.h"
#include "Math/RichardsonDerivator.h"
 
class RooRealVar;

class HggScalarLineShapePdf : public RooAbsPdf {
  
 public:
  
  HggScalarLineShapePdf();
  HggScalarLineShapePdf(const char *name, const char *title, RooAbsReal& x,
			RooAbsReal& mean, RooAbsReal& width, RooAbsReal& alpha,
			int cme=13, bool performTruncation=false);
  
  HggScalarLineShapePdf(const HggScalarLineShapePdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new HggScalarLineShapePdf(*this,newname); }
  inline virtual ~HggScalarLineShapePdf() { }

  void setTopMassInfinite(bool flag=true){_performTruncation=flag;}
  void setNu(Double_t nu){_nu=nu;}
  void setSigma(Double_t sigma){_sigma=sigma;}
  void setEpsilon(Double_t epsilon){_epsilon=epsilon;}
  Double_t getNu(){return _nu;}
  Double_t getSigma(){return _sigma;}
  Double_t getEpsilon(){return _epsilon;}
 protected:

  RooRealProxy _x;
  RooRealProxy _mean;
  RooRealProxy _width;
  RooRealProxy _alpha;		// Position of cut-off, quantified as mX-alpha*wX
  
  int _cme;			// Center-of-mass energy
  bool _performTruncation;	// Flag for damping away the low-mass bump
  Double_t _nu;			// Controling the type of damping
  Double_t _sigma;		// Sigma of the Gaussian CDF (in unit of wX)
  Double_t _epsilon;		// Small number for calculating derivative using Richardson algorithm
  
  Double_t ME(Double_t mgg) const;		// Matrix element part
  Double_t MEPL(Double_t mgg) const;		// Matrix element times parton luminosity
  Double_t dMEPLdmgg(Double_t mgg) const;	// First derivative of matrix element
  Double_t d2MEPLdmgg2(Double_t mgg) const;	// Second derivative of matrix element
  Double_t PL(Double_t mgg) const;		// Parton luminosity part
  Double_t BW(Double_t mgg) const;		// Breit-Wigner part
  Double_t evaluate() const;

 private:
  
  ClassDef(HggScalarLineShapePdf,1); // Crystal Ball lineshape PDF
    
};
  
#endif
