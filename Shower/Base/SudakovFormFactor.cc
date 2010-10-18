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
#include "ShowerKinematics.h"

using namespace Herwig;

void SudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << splittingFn_ << alpha_ << pdfmax_ << particles_ << pdffactor_
     << a_ << b_ << ounit(c_,GeV) << ounit(kinCutoffScale_,GeV) << cutOffOption_
     << ounit(vgcut_,GeV) << ounit(vqcut_,GeV) 
     << ounit(pTmin_,GeV) << ounit(pT2min_,GeV2);
}

void SudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> splittingFn_ >> alpha_ >> pdfmax_ >> particles_ >> pdffactor_
     >> a_ >> b_ >> iunit(c_,GeV) >> iunit(kinCutoffScale_,GeV) >> cutOffOption_
     >> iunit(vgcut_,GeV) >> iunit(vqcut_,GeV) 
     >> iunit(pTmin_,GeV) >> iunit(pT2min_,GeV2);
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
			       &Herwig::SudakovFormFactor::splittingFn_,
			       false, false, true, false);

  static Reference<SudakovFormFactor,ShowerAlpha>
    interfaceAlpha("Alpha",
		   "A reference to the Alpha object",
		   &Herwig::SudakovFormFactor::alpha_,
		   false, false, true, false);

  static Parameter<SudakovFormFactor,double> interfacePDFmax
    ("PDFmax",
     "Maximum value of PDF weight. ",
     &SudakovFormFactor::pdfmax_, 35.0, 1.0, 100000.0,
     false, false, Interface::limited);

  static Switch<SudakovFormFactor,unsigned int> interfacePDFFactor
    ("PDFFactor",
     "Include additional factors in the overestimate for the PDFs",
     &SudakovFormFactor::pdffactor_, 0, false, false);
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


  static Switch<SudakovFormFactor,unsigned int> interfaceCutOffOption
    ("CutOffOption",
     "The type of cut-off to use to end the shower",
     &SudakovFormFactor::cutOffOption_, 0, false, false);
  static SwitchOption interfaceCutOffOptionDefault
    (interfaceCutOffOption,
     "Default",
     "Use the standard Herwig++ cut-off on virtualities with the minimum"
     " virtuality depending on the mass of the branching particle",
     0);
  static SwitchOption interfaceCutOffOptionFORTRAN
    (interfaceCutOffOption,
     "FORTRAN",
     "Use a FORTRAN-like cut-off on virtualities",
     1);
  static SwitchOption interfaceCutOffOptionpT
    (interfaceCutOffOption,
     "pT",
     "Use a cut on the minimum allowed pT",
     2);

  static Parameter<SudakovFormFactor,double> interfaceaParameter
    ("aParameter",
     "The a parameter for the kinematic cut-off",
     &SudakovFormFactor::a_, 0.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<SudakovFormFactor,double> interfacebParameter
    ("bParameter",
     "The b parameter for the kinematic cut-off",
     &SudakovFormFactor::b_, 2.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<SudakovFormFactor,Energy> interfacecParameter
    ("cParameter",
     "The c parameter for the kinematic cut-off",
     &SudakovFormFactor::c_, GeV, 0.3*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<SudakovFormFactor,Energy>
    interfaceKinScale ("cutoffKinScale",
		       "kinematic cutoff scale for the parton shower phase"
		       " space (unit [GeV])",
		       &SudakovFormFactor::kinCutoffScale_, GeV, 
		       2.3*GeV, 0.001*GeV, 10.0*GeV,false,false,false);
  
  static Parameter<SudakovFormFactor,Energy> interfaceGluonVirtualityCut
    ("GluonVirtualityCut",
     "For the FORTRAN cut-off option the minimum virtuality of the gluon",
     &SudakovFormFactor::vgcut_, GeV, 0.85*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<SudakovFormFactor,Energy> interfaceQuarkVirtualityCut
    ("QuarkVirtualityCut",
     "For the FORTRAN cut-off option the minimum virtuality added to"
     " the mass for particles other than the gluon",
     &SudakovFormFactor::vqcut_, GeV, 0.85*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<SudakovFormFactor,Energy> interfacepTmin
    ("pTmin",
     "The minimum pT if using a cut-off on the pT",
     &SudakovFormFactor::pTmin_, GeV, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);
}

bool SudakovFormFactor::
PDFVeto(const Energy2 t, const double x,
	const tcPDPtr parton0, const tcPDPtr parton1,
	Ptr<BeamParticleData>::transient_const_pointer beam) const {
  assert(pdf_);

  Energy2 theScale = t;
  if (theScale < sqr(freeze_)) theScale = sqr(freeze_);

  double newpdf(0.0), oldpdf(0.0);
  //different treatment of MPI ISR is done via CascadeHandler::resetPDFs()
  newpdf=pdf_->xfx(beam,parton0,theScale,x/z());
  oldpdf=pdf_->xfx(beam,parton1,theScale,x);

  
  if(newpdf<=0.) return true;
  if(oldpdf<=0.) return false;
  double ratio = newpdf/oldpdf;

  double maxpdf(pdfmax_);
  switch (pdffactor_) {
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
		       << ":PDFmax (by a factor of "
		       << ratio/maxpdf <<") for " 
		       << parton0->PDGName() << " to " 
		       << parton1->PDGName() << "\n";
  }
  return ratio < UseRandom::rnd()*maxpdf;
}

void SudakovFormFactor::addSplitting(const IdList & in) {
  bool add=true;
  for(unsigned int ix=0;ix<particles_.size();++ix) {
    if(particles_[ix].size()==in.size()) {
      bool match=true;
      for(unsigned int iy=0;iy<in.size();++iy) {
	if(particles_[ix][iy]!=in[iy]) {
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
  if(add) particles_.push_back(in);
}

void SudakovFormFactor::removeSplitting(const IdList & in) {
  for(vector<IdList>::iterator it=particles_.begin();
      it!=particles_.end();++it) {
    if(it->size()==in.size()) {
      bool match=true;
      for(unsigned int iy=0;iy<in.size();++iy) {
	if((*it)[iy]!=in[iy]) {
	  match=false;
	  break;
	}
      }
      if(match) {
	vector<IdList>::iterator itemp=it;
	--itemp;
	particles_.erase(it);
	it = itemp;
      }
    }
  }
}

Energy2 SudakovFormFactor::guesst(Energy2 t1,unsigned int iopt,
					 const IdList &ids,
					 double enhance,bool ident) const {
  unsigned int pdfopt = iopt!=1 ? 0 : pdffactor_;
  double c =
    1./((splittingFn_->integOverP(zlimits_.second,ids,pdfopt) -
	 splittingFn_->integOverP(zlimits_.first ,ids,pdfopt))* 
	alpha_->overestimateValue()/Constants::twopi*enhance);
  assert(iopt<=2);
  if(iopt==1) {
    c/=pdfmax_;
    if(ident) c*=0.5;
  }
  else if(iopt==2) c*=-1.;
  if(splittingFn_->interactionOrder()==1) {
    return t1*pow(UseRandom::rnd(),c);
  }
  else {
    assert(false && "Units are dubious here.");
    int nm(splittingFn()->interactionOrder()-1);
    c/=Math::powi(alpha_->overestimateValue()/Constants::twopi,nm);
    return t1 /       pow (1. - nm*c*log(UseRandom::rnd()) 
			   * Math::powi(t1*UnitRemoval::InvE2,nm)
			   ,1./double(nm));
  }
}

double SudakovFormFactor::guessz (unsigned int iopt, const IdList &ids) const {
  unsigned int pdfopt = iopt!=1 ? 0 : pdffactor_;
  double lower = splittingFn_->integOverP(zlimits_.first,ids,pdfopt);
  return splittingFn_->invIntegOverP
    (lower + UseRandom::rnd()*(splittingFn_->integOverP(zlimits_.second,ids,pdfopt) - 
			       lower),ids,pdfopt);
}

void SudakovFormFactor::doinit() {
  Interfaced::doinit();
  pT2min_ = cutOffOption()==2 ? sqr(pTmin_) : ZERO; 
}
