// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SudakovFormFactor class.
//

#include "SudakovFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig++/Shower/ShowerHandler.h"

using namespace Herwig;

void SudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _splittingFn << _alpha << _pdfmax << _particles;
}

void SudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _splittingFn >> _alpha >> _pdfmax >> _particles;
}

AbstractClassDescription<SudakovFormFactor> SudakovFormFactor::initSudakovFormFactor;
// Definition of the static class description member.

void SudakovFormFactor::Init() {

  static ClassDocumentation<SudakovFormFactor> documentation
    ("The SudakovFormFactor class is the base class for the implementation of Sudakov"
     " form factors in Herwig++");

  static Reference<SudakovFormFactor,SplittingFunction>
    interfaceSplittingFunction("SplittingFunction",
			       "A reference to the SplittingFunction object",
			       &Herwig::SudakovFormFactor::_splittingFn,
			       false, false, true, false);

  static Reference<SudakovFormFactor,ShowerAlpha>
    interfaceAlpha("Alpha",
		   "A reference to the Alpha object",
		   &Herwig::SudakovFormFactor::_alpha,
		   false, false, true, false);

  static Parameter<SudakovFormFactor,double> interfacePDFmax
    ("PDFmax",
     "Maximum value of PDF weight. ",
     &SudakovFormFactor::_pdfmax, 35.0, 1.0, 4000.0,
     false, false, Interface::limited);

}

bool SudakovFormFactor::
PDFVeto(const Energy2 t, const double x,
	const tcPDPtr parton0, const tcPDPtr parton1,
	Ptr<BeamParticleData>::transient_const_pointer beam) const {
  tcPDFPtr pdf = tcPDFPtr();
  //using the pdf's associated with the ShowerHandler assures, that
  //modified pdf's are used for the secondary interactions via 
  //CascadeHandler::resetPDFs(...)
  if(ShowerHandler::currentHandler()->firstPDF().particle() == beam)
    pdf = ShowerHandler::currentHandler()->firstPDF().pdf();
  if(ShowerHandler::currentHandler()->secondPDF().particle() == beam)
    pdf = ShowerHandler::currentHandler()->secondPDF().pdf();

  assert(pdf);
  // remember: pdf's q is cut in pdf class.  should probably be done here! 
  // this would correspond to QSPAC in F-HERWIG.     

  double newpdf(0.0), oldpdf(0.0);
  //different treatment of MPI ISR by hand:
  if( ShowerHandler::currentHandler()->FirstInt() ){
    newpdf=pdf->xfx(beam,parton0,t,x/z());
    oldpdf=pdf->xfx(beam,parton1,t,x);
  }else{
    newpdf=pdf->xfx(beam,parton0,t,x/z()) - pdf->xfvx(beam,parton0,t,x/z());
    oldpdf=pdf->xfx(beam,parton1,t,x) - pdf->xfvx(beam,parton1,t,x);
  }
  
  if(newpdf<=0.) return true;
  if(oldpdf<=0.) return false;
  double ratio = newpdf/oldpdf;

  // ratio / PDFMax must be a probability <= 1.0
  if (ratio > _pdfmax) {
    generator()->log() << "PDFVeto warning: Ratio (" << ratio 
			    << ") > " << name() << ":PDFmax ("
		       <<_pdfmax <<") for " 
		       << parton0->PDGName() << " to " 
		       << parton1->PDGName() << "\n";
  }
  return ratio < UseRandom::rnd()*_pdfmax;
}

void SudakovFormFactor::addSplitting(const IdList & in) {
  bool add=true;
  for(unsigned int ix=0;ix<_particles.size();++ix) {
    if(_particles[ix].size()==in.size()) {
      bool match=true;
      for(unsigned int iy=0;iy<in.size();++iy) {
	if(_particles[ix][iy]!=in[iy]) {
	  match=false;
	  break;
	}
      }
      if(match) {
	add=false;
	break;
      }
    }
  }
  if(add) _particles.push_back(in);
}
