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

using namespace Herwig;

void SudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _splittingFn << _alpha << _pdfmax << _particles
     << _a << _b << _c << _kinCutoffScale;
}

void SudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _splittingFn >> _alpha >> _pdfmax >> _particles 
     >> _a >> _b >> _c >> _kinCutoffScale;
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

  static Parameter<SudakovFormFactor,double> interfaceaParameter
    ("aParameter",
     "The a parameter for the kinematic cut-off",
     &SudakovFormFactor::_a, 0.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<SudakovFormFactor,double> interfacebParameter
    ("bParameter",
     "The b parameter for the kinematic cut-off",
     &SudakovFormFactor::_b, 2.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<SudakovFormFactor,Energy> interfacecParameter
    ("cParameter",
     "The c parameter for the kinematic cut-off",
     &SudakovFormFactor::_c, GeV, 0.3*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<SudakovFormFactor,Energy>
    interfaceKinScale ("cutoffKinScale",
		       "kinematic cutoff scale for the parton shower phase"
		       " space (unit [GeV])",
		       &SudakovFormFactor::_kinCutoffScale, GeV, 
		       0.75*GeV, 0.001*GeV, 10.0*GeV,false,false,false);
}

bool SudakovFormFactor::
PDFVeto(const Energy2 t, const double x,
	const tcPDPtr parton0, const tcPDPtr parton1,
	Ptr<BeamParticleData>::transient_const_pointer beam) const {
  tcPDFPtr pdf=beam->pdf();
  assert(pdf);
  // remember: pdf's q is cut in pdf class.  should probably be done here! 
  // this would correspond to QSPAC in F-HERWIG. 
  double newpdf=pdf->xfx(beam,parton0,t,x/z());
  double oldpdf=pdf->xfx(beam,parton1,t,x);
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
