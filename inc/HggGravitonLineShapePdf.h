////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// ==================================================================         //
// Implemented by Hongtao Yang <Hongtao.Yang@cern.ch> on Feb. 6, 2016         //
// with theory inputs provided by Jan Stark <stark@in2p3.fr>                  //
// ==================================================================         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_HggGravitonLineShapePdf
#define ROOT_HggGravitonLineShapePdf

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

class RooRealVar;

class HggGravitonLineShapePdf : public RooAbsPdf {
  
 public:
  
  HggGravitonLineShapePdf();
  HggGravitonLineShapePdf(const char *name, const char *title, RooAbsReal& _x,
		       RooAbsReal& _mG, RooAbsReal& _GkM, bool isTopMassInfi=false);
  
  HggGravitonLineShapePdf(const HggGravitonLineShapePdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new HggGravitonLineShapePdf(*this,newname); }
  inline virtual ~HggGravitonLineShapePdf() { }

 protected:

  RooRealProxy x ;
  RooRealProxy mG ;
  RooRealProxy GkM ;

  Double_t evaluate() const;
  
 private:

  ClassDef(HggGravitonLineShapePdf,1); // Crystal Ball lineshape PDF
    
};
  
#endif
