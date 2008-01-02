// -*- C++ -*-
//
// SudakovFormFactor.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig++/Shower/ShowerHandler.h"

using namespace Herwig;

void SudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _splittingFn << _alpha << _pdfmax << _particles
     << _a << _b << ounit(_c,GeV) << ounit(_kinCutoffScale,GeV)
     << _pdffactor;
}

void SudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _splittingFn >> _alpha >> _pdfmax >> _particles 
     >> _a >> _b >> iunit(_c,GeV) >> iunit(_kinCutoffScale,GeV)
     >> _pdffactor;
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
		       2.3*GeV, 0.001*GeV, 10.0*GeV,false,false,false);

  static Switch<SudakovFormFactor,unsigned int> interfacePDFFactor
    ("PDFFactor",
     "Include additional factors in the overestimate for the PDFs",
     &SudakovFormFactor::_pdffactor, 0, false, false);
  static SwitchOption interfacePDFFactorOff
    (interfacePDFFactor,
     "Off",
     "Don't include any factors",
     0);
  static SwitchOption interfacePDFFactorOverZ
    (interfacePDFFactor,
     "OverZ",
     "Include an additional factor of 1/z",
     1);
  static SwitchOption interfacePDFFactorOverOneMinusZ
    (interfacePDFFactor,
     "OverOneMinusZ",
     "Include an additional factor of 1/(1-z)",
     2);
  static SwitchOption interfacePDFFactorOverZOneMinusZ
    (interfacePDFFactor,
     "OverZOneMinusZ",
     "Include an additional factor of 1/z/(1-z)",
     3);
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
  //different treatment of MPI ISR is done via CascadeHandler::resetPDFs()
  newpdf=pdf->xfx(beam,parton0,t,x/z());
  oldpdf=pdf->xfx(beam,parton1,t,x);

  
  if(newpdf<=0.) return true;
  if(oldpdf<=0.) return false;
  double ratio = newpdf/oldpdf;

  double maxpdf(_pdfmax);
  switch (_pdffactor) {
  case 1:
    maxpdf /= z();
    break;
  case 2:
    maxpdf /= 1.-z();
    break;
  case 3:
    maxpdf /= (z()*(1.-z()));
    break;
  }

  // ratio / PDFMax must be a probability <= 1.0
  if (ratio > maxpdf) {
    generator()->log() << "PDFVeto warning: Ratio > " << name() 
		       << ":PDFmax (by a factor of"
		       << ratio/maxpdf <<") for " 
		       << parton0->PDGName() << " to " 
		       << parton1->PDGName() << "\n";
  }
  return ratio < UseRandom::rnd()*maxpdf;
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
