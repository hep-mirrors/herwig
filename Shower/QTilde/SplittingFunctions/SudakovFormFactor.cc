// -*- C++ -*-
//
// SudakovFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig/Shower/QTilde/Kinematics/ShowerKinematics.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"
#include "Herwig/Shower/QTilde/Kinematics/FS_QTildeShowerKinematics1to2.h"
#include "Herwig/Shower/QTilde/Kinematics/IS_QTildeShowerKinematics1to2.h"
#include "Herwig/Shower/QTilde/Kinematics/Decay_QTildeShowerKinematics1to2.h"
#include "Herwig/Shower/QTilde/Kinematics/KinematicHelpers.h"
#include "SudakovCutOff.h"
#include "ThePEG/Utilities/EnumIO.h"

#include <array>
using std::array;

using namespace Herwig;

DescribeAbstractClass<SudakovFormFactor,Interfaced>
describeSudakovFormFactor ("Herwig::SudakovFormFactor","");

void SudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << alpha_ << pdfmax_ << particles_ << pdffactor_ << cutoff_
     << oenum(interactionType_) << oenum(colourStructure_) 
     << angularOrdered_ << scaleChoice_ << strictAO_ << colourFactor_;
}

void SudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> alpha_ >> pdfmax_ >> particles_ >> pdffactor_ >> cutoff_
     >> ienum(interactionType_) >> ienum(colourStructure_) 
     >> angularOrdered_ >> scaleChoice_ >> strictAO_ >> colourFactor_;
}

void SudakovFormFactor::Init() {

  static ClassDocumentation<SudakovFormFactor> documentation
    ("The SudakovFormFactor class is the base class for the implementation of Sudakov"
     " form factors in Herwig");

  static Reference<SudakovFormFactor,ShowerAlpha>
    interfaceAlpha("Alpha",
		   "A reference to the Alpha object",
		   &Herwig::SudakovFormFactor::alpha_,
		   false, false, true, false);

  static Reference<SudakovFormFactor,SudakovCutOff>
    interfaceCutoff("Cutoff",
		   "A reference to the SudakovCutOff object",
		   &Herwig::SudakovFormFactor::cutoff_,
		   false, false, true, false);

  static Parameter<SudakovFormFactor,double> interfacePDFmax
    ("PDFmax",
     "Maximum value of PDF weight. ",
     &SudakovFormFactor::pdfmax_, 35.0, 1.0, 1000000.0,
     false, false, Interface::limited);

  static Switch<SudakovFormFactor,unsigned int> interfacePDFFactor
    ("PDFFactor",
     "Include additional factors in the overestimate for the PDFs",
     &SudakovFormFactor::pdffactor_, 0, false, false);
  static SwitchOption interfacePDFFactorNo
    (interfacePDFFactor,
     "No",
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
  static SwitchOption interfacePDFFactorOverRootZ
    (interfacePDFFactor,
     "OverRootZ",
     "Include an additional factor of 1/sqrt(z)",
     4);
  static SwitchOption interfacePDFFactorRootZ
    (interfacePDFFactor,
     "RootZ",
     "Include an additional factor of sqrt(z)",
     5);

  static Switch<SudakovFormFactor,bool> interfaceAngularOrdered
    ("AngularOrdered",
     "Whether or not this interaction is angular ordered, "
     "normally only g->q qbar and gamma-> f fbar are the only ones which aren't.",
     &SudakovFormFactor::angularOrdered_, true, false, false);
  static SwitchOption interfaceAngularOrderedYes
    (interfaceAngularOrdered,
     "Yes",
     "Interaction is angular ordered",
     true);
  static SwitchOption interfaceAngularOrderedNo
    (interfaceAngularOrdered,
     "No",
     "Interaction isn't angular ordered",
     false);

  static Switch<SudakovFormFactor,unsigned int> interfaceScaleChoice
    ("ScaleChoice",
     "The scale choice to be used",
     &SudakovFormFactor::scaleChoice_, 2, false, false);
  static SwitchOption interfaceScaleChoicepT
    (interfaceScaleChoice,
     "pT",
     "pT of the branching",
     0);
  static SwitchOption interfaceScaleChoiceQ2
    (interfaceScaleChoice,
     "Q2",
     "Q2 of the branching",
     1);
  static SwitchOption interfaceScaleChoiceFromAngularOrdering
    (interfaceScaleChoice,
     "FromAngularOrdering",
     "If angular order use pT, otherwise Q2",
     2);

  static Switch<SudakovFormFactor,bool> interfaceStrictAO
    ("StrictAO",
     "Whether or not to apply strict angular-ordering,"
     " i.e. for QED even in QCD emission, and vice versa",
     &SudakovFormFactor::strictAO_, true, false, false);
  static SwitchOption interfaceStrictAOYes
    (interfaceStrictAO,
     "Yes",
     "Apply strict ordering",
     true);
  static SwitchOption interfaceStrictAONo
    (interfaceStrictAO,
     "No",
     "Don't apply strict ordering",
     false);

  static Switch<SudakovFormFactor,ColourStructure> interfaceColourStructure
    ("ColourStructure",
     "The colour structure for the splitting function",
     &SudakovFormFactor::colourStructure_, Undefined, false, false);
  static SwitchOption interfaceColourStructureTripletTripletOctet
    (interfaceColourStructure,
     "TripletTripletOctet",
     "3 -> 3 8",
     TripletTripletOctet);
  static SwitchOption interfaceColourStructureOctetOctetOctet
    (interfaceColourStructure,
     "OctetOctetOctet",
     "8 -> 8 8",
     OctetOctetOctet);
  static SwitchOption interfaceColourStructureOctetTripletTriplet
    (interfaceColourStructure,
     "OctetTripletTriplet",
     "8 -> 3 3bar",
     OctetTripletTriplet);
  static SwitchOption interfaceColourStructureTripletOctetTriplet
    (interfaceColourStructure,
     "TripletOctetTriplet",
     "3 -> 8 3",
     TripletOctetTriplet); 
  static SwitchOption interfaceColourStructureSextetSextetOctet
    (interfaceColourStructure,
     "SextetSextetOctet",
     "6 -> 6 8",
     SextetSextetOctet);
  static SwitchOption interfaceColourStructureChargedChargedNeutral
    (interfaceColourStructure,
     "ChargedChargedNeutral",
     "q -> q 0",
     ChargedChargedNeutral);
  static SwitchOption interfaceColourStructureNeutralChargedCharged
    (interfaceColourStructure,
     "NeutralChargedCharged",
     "0 -> q qbar",
     NeutralChargedCharged);
  static SwitchOption interfaceColourStructureChargedNeutralCharged
    (interfaceColourStructure,
     "ChargedNeutralCharged",
     "q -> 0 q",
     ChargedNeutralCharged);

  static Switch<SudakovFormFactor,ShowerInteraction> 
    interfaceInteractionType
    ("InteractionType",
     "Type of the interaction",
     &SudakovFormFactor::interactionType_, 
     ShowerInteraction::UNDEFINED, false, false);
  static SwitchOption interfaceInteractionTypeQCD
    (interfaceInteractionType,
     "QCD","QCD",ShowerInteraction::QCD);
  static SwitchOption interfaceInteractionTypeQED
    (interfaceInteractionType,
     "QED","QED",ShowerInteraction::QED);

}

bool SudakovFormFactor::alphaSVeto(Energy2 pt2) const {
  double ratio=alphaSVetoRatio(pt2,1.);
  return UseRandom::rnd() > ratio;
}

double SudakovFormFactor::alphaSVetoRatio(Energy2 pt2, double factor) const {
  factor *= ShowerHandler::currentHandler()->renormalizationScaleFactor();
  return alpha_->ratio(pt2, factor);
}


bool SudakovFormFactor::PDFVeto(const Energy2 t, const double x,
				const tcPDPtr parton0, const tcPDPtr parton1,
				Ptr<BeamParticleData>::transient_const_pointer beam) const {
  double ratio=PDFVetoRatio(t,x,parton0,parton1,beam,1.);
  return UseRandom::rnd() > ratio;
}

double SudakovFormFactor::PDFVetoRatio(const Energy2 t, const double x,
				       const tcPDPtr parton0, const tcPDPtr parton1,
				       Ptr<BeamParticleData>::transient_const_pointer beam,double factor) const {
  assert(pdf_);
  Energy2 theScale = t * sqr(ShowerHandler::currentHandler()->factorizationScaleFactor()*factor);
  if (theScale < sqr(freeze_)) theScale = sqr(freeze_);

  const double newpdf=pdf_->xfx(beam,parton0,theScale,x/z());
  if(newpdf<=0.) return 0.;

  const double oldpdf=pdf_->xfx(beam,parton1,theScale,x);
  if(oldpdf<=0.) return 1.;
  
  const double ratio = newpdf/oldpdf;
  double maxpdf = pdfmax_;

  switch (pdffactor_) {
  case 0: break;
  case 1: maxpdf /= z(); break;
  case 2: maxpdf /= 1.-z(); break;
  case 3: maxpdf /= (z()*(1.-z())); break;
  case 4: maxpdf /= sqrt(z()); break;
  case 5: maxpdf *= sqrt(z()); break;
  default :
    throw Exception() << "SudakovFormFactor::PDFVetoRatio invalid PDFfactor = "
		      << pdffactor_ << Exception::runerror;
    
  }

  if (ratio > maxpdf) {
    generator()->log() << "PDFVeto warning: Ratio > " << name()
                       << ":PDFmax (by a factor of "
                       << ratio/maxpdf <<") for "
                       << parton0->PDGName() << " to "
                       << parton1->PDGName() << "\n";
  }
  return ratio/maxpdf ;
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

void SudakovFormFactor::guesstz(Energy2 t1,unsigned int iopt,
				  const IdList &ids,
				  double enhance,bool ident,
				  double detune, 
				  Energy2 &t_main, double &z_main) {
  unsigned int pdfopt = iopt!=1 ? 0 : pdffactor_;
  double lower = integOverP(zlimits_.first ,ids,pdfopt);
  double upper = integOverP(zlimits_.second,ids,pdfopt);
  double c = 1./((upper - lower) * colourFactor()
           * alpha_->overestimateValue()/Constants::twopi*enhance*detune);
  double r = UseRandom::rnd();
  assert(iopt<=2);
  if(iopt==1) {
    c/=pdfmax_;
    //symmetry of FS gluon splitting
    if(ident) c*=0.5;
  }
  else if(iopt==2) c*=-1.;
  // guessing t
  if(iopt!=2 || c*log(r)<log(Constants::MaxEnergy2/t1)) {
    t_main = t1*pow(r,c);
  }
  else
    t_main = Constants::MaxEnergy2;
  // guessing z
  z_main = invIntegOverP(lower + UseRandom::rnd()
         *(upper - lower),ids,pdfopt);
}

bool SudakovFormFactor::guessTimeLike(Energy2 &t,Energy2 tmin,double enhance,
				  double detune) {
  Energy2 told = t;
  // calculate limits on z and if lower>upper return
  if(!computeTimeLikeLimits(t)) return false;
  // guess values of t and z
  guesstz(told,0,ids_,enhance,ids_[1]==ids_[2],detune,t,z_);
  // actual values for z-limits
  if(!computeTimeLikeLimits(t)) return false;
  if(t<tmin) {
    t=-1.0*GeV2;
    return false;
  }
  else
    return true; 
} 

bool SudakovFormFactor::guessSpaceLike(Energy2 &t, Energy2 tmin, const double x,
				   double enhance,
				   double detune) {
  Energy2 told = t;
  // calculate limits on z if lower>upper return
  if(!computeSpaceLikeLimits(t,x)) return false;
  // guess values of t and z
  guesstz(told,1,ids_,enhance,ids_[1]==ids_[2],detune,t,z_);
  // actual values for z-limits
  if(!computeSpaceLikeLimits(t,x)) return false;
  if(t<tmin) {
    t=-1.0*GeV2;
    return false;
  }
  else
    return true; 
} 

bool SudakovFormFactor::PSVeto(const Energy2 t) {
  // still inside PS, return true if outside
  // check vs overestimated limits
  if (z() < zlimits_.first || z() > zlimits_.second) return true;

  Energy2 m02 = (ids_[0]->id()!=ParticleID::g && ids_[0]->id()!=ParticleID::gamma) ?
  	masssquared_[0] : Energy2();
  
  Energy2 pt2 = QTildeKinematics::pT2_FSR(t,z(),m02,masssquared_[1],masssquared_[2],
					  masssquared_[1],masssquared_[2]);

  // if pt2<0 veto
  if (pt2<cutoff_->pT2min()) return true;
  // otherwise calculate pt and return
  pT_ = sqrt(pt2);
  return false;
}

ShoKinPtr SudakovFormFactor::generateNextTimeBranching(const Energy startingScale,
						   const IdList &ids,
						   const RhoDMatrix & rho,
						   double enhance,
						   double detuning) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  q_ = ZERO;
  z_ = 0.;
  phi_ = 0.; 
  // perform initialization
  Energy2 tmax(sqr(startingScale)),tmin;
  initialize(ids,tmin);
  // check max > min
  if(tmax<=tmin) return ShoKinPtr();
  // calculate next value of t using veto algorithm
  Energy2 t(tmax);
  // no shower variations to calculate
  if(ShowerHandler::currentHandler()->showerVariations().empty()){
    // Without variations do the usual Veto algorithm
    // No need for more if-statements in this loop.
    do {
      if(!guessTimeLike(t,tmin,enhance,detuning)) break;
    }
    while(PSVeto(t) ||
        SplittingFnVeto(z()*(1.-z())*t,ids,true,rho,detuning) || 
        alphaSVeto(pTScale() ? sqr(z()*(1.-z()))*t : z()*(1.-z())*t));
  }
  else {
    bool alphaRew(true),PSRew(true),SplitRew(true);
    do {
      if(!guessTimeLike(t,tmin,enhance,detuning)) break;
      PSRew=PSVeto(t);
      if (PSRew) continue;
      SplitRew=SplittingFnVeto(z()*(1.-z())*t,ids,true,rho,detuning);
      alphaRew=alphaSVeto(pTScale() ? sqr(z()*(1.-z()))*t : z()*(1.-z())*t);
      double factor=alphaSVetoRatio(pTScale() ? sqr(z()*(1.-z()))*t : z()*(1.-z())*t,1.)*
                    SplittingFnVetoRatio(z()*(1.-z())*t,ids,true,rho,detuning);

      tShowerHandlerPtr ch = ShowerHandler::currentHandler();

      if( !(SplitRew || alphaRew) ) {
        //Emission
        q_ = t > ZERO ? Energy(sqrt(t)) : -1.*MeV;
        if (q_ <= ZERO) break;
      }

        for ( map<string,ShowerVariation>::const_iterator var =
	          ch->showerVariations().begin();
	          var != ch->showerVariations().end(); ++var ) {
          if ( ( ch->firstInteraction() && var->second.firstInteraction ) ||
	           ( !ch->firstInteraction() && var->second.secondaryInteractions ) ) {

                double newfactor = alphaSVetoRatio(pTScale() ?
                                        sqr(z()*(1.-z()))*t :
                                        z()*(1.-z())*t,var->second.renormalizationScaleFactor)
		  * SplittingFnVetoRatio(z()*(1.-z())*t,ids,true,rho,detuning);

                double varied;
                if ( SplitRew || alphaRew ) {
                  // No Emission
                  varied = (1. - newfactor) / (1. - factor);
                } else {
                  // Emission
                  varied = newfactor / factor;
                }

                map<string,double>::iterator wi = ch->currentWeights().find(var->first);
	        if ( wi != ch->currentWeights().end() )
	          wi->second *= varied;
	        else {
                  assert(false);
                  //ch->currentWeights()[var->first] = varied;
                }
	  }
        }
      
    }
    while(PSRew || SplitRew || alphaRew);
  }
  q_ = t > ZERO ? Energy(sqrt(t)) : -1.*MeV;
  if(q_ < ZERO) return ShoKinPtr();
  
  // return the ShowerKinematics object
  return new_ptr(FS_QTildeShowerKinematics1to2(q_,z(),phi(),pT(),this)); 
}

ShoKinPtr SudakovFormFactor::
generateNextSpaceBranching(const Energy startingQ,
			   const IdList &ids,
			   double x,
			   const RhoDMatrix & rho,
			   double enhance,
			   Ptr<BeamParticleData>::transient_const_pointer beam,
			   double detuning) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  q_ = ZERO;
  z_ = 0.;
  phi_ = 0.;
  // perform the initialization
  Energy2 tmax(sqr(startingQ)),tmin;
  initialize(ids,tmin);
  // check max > min
  if(tmax<=tmin) return ShoKinPtr();
  // calculate next value of t using veto algorithm
  Energy2 t(tmax),pt2(ZERO);
  // no shower variations
  if(ShowerHandler::currentHandler()->showerVariations().empty()) {
    // Without variations do the usual Veto algorithm
    // No need for more if-statements in this loop.
    do {
      if(!guessSpaceLike(t,tmin,x,enhance,detuning)) break;
      pt2 = QTildeKinematics::pT2_ISR(t,z(),masssquared_[2]);
    }
    while(pt2 < cutoff_->pT2min()||
        z() > zlimits_.second||
	  SplittingFnVeto((1.-z())*t/z(),ids,false,rho,detuning)||
        alphaSVeto(pTScale() ? sqr(1.-z())*t : (1.-z())*t)||
        PDFVeto(t,x,ids[0],ids[1],beam));
  }
  // shower variations
  else {
    bool alphaRew(true),PDFRew(true),ptRew(true),zRew(true),SplitRew(true);
    do {
      if(!guessSpaceLike(t,tmin,x,enhance,detuning)) break;
      pt2 = QTildeKinematics::pT2_ISR(t,z(),masssquared_[2]);
      ptRew=pt2 < cutoff_->pT2min();
      zRew=z() > zlimits_.second;
      if (ptRew||zRew) continue;
      SplitRew=SplittingFnVeto((1.-z())*t/z(),ids,false,rho,detuning);
      alphaRew=alphaSVeto(pTScale() ? sqr(1.-z())*t : (1.-z())*t);
      PDFRew=PDFVeto(t,x,ids[0],ids[1],beam);
      double factor=PDFVetoRatio(t,x,ids[0],ids[1],beam,1.)*
                    alphaSVetoRatio(pTScale() ? sqr(1.-z())*t : (1.-z())*t,1.)*
	SplittingFnVetoRatio((1.-z())*t/z(),ids,false,rho,detuning);

      tShowerHandlerPtr ch = ShowerHandler::currentHandler();

      if( !(PDFRew || SplitRew || alphaRew) ) {
        //Emission
        q_ = t > ZERO ? Energy(sqrt(t)) : -1.*MeV;
        if (q_ <= ZERO) break;
      }

        for ( map<string,ShowerVariation>::const_iterator var =
	          ch->showerVariations().begin();
	          var != ch->showerVariations().end(); ++var ) {
          if ( ( ch->firstInteraction() && var->second.firstInteraction ) ||
	           ( !ch->firstInteraction() && var->second.secondaryInteractions ) ) {



            double newfactor = PDFVetoRatio(t,x,ids[0],ids[1],beam,var->second.factorizationScaleFactor)*
                           alphaSVetoRatio(pTScale() ?
                           sqr(1.-z())*t : (1.-z())*t,var->second.renormalizationScaleFactor)
	      *SplittingFnVetoRatio((1.-z())*t/z(),ids,false,rho,detuning);

            double varied;
            if( PDFRew || SplitRew || alphaRew) {
                // No Emission
                varied = (1. - newfactor) / (1. - factor);
            } else {
                // Emission
                varied = newfactor / factor;
            }


            map<string,double>::iterator wi = ch->currentWeights().find(var->first);
            if ( wi != ch->currentWeights().end() )
	      wi->second *= varied;
	    else {
	      assert(false);
	      //ch->currentWeights()[var->first] = varied;
            }
	  }
        }
      
    }
    while( PDFRew || SplitRew || alphaRew);
  }
  if(t > ZERO && zlimits_.first < zlimits_.second)  q_ = sqrt(t);
  else return ShoKinPtr();
  
  pT_ = sqrt(pt2);
  // create the ShowerKinematics and return it
  return new_ptr(IS_QTildeShowerKinematics1to2(q_,z(),phi(),pT(),this)); 
}

void SudakovFormFactor::initialize(const IdList & ids, Energy2 & tmin) {
  ids_=ids;
  tmin = 4.*cutoff_->pT2min();
  masses_ = cutoff_->virtualMasses(ids);
  masssquared_.clear();
  for(unsigned int ix=0;ix<masses_.size();++ix) {
    masssquared_.push_back(sqr(masses_[ix]));
    if(ix>0) tmin=max(masssquared_[ix],tmin);
  }
}

ShoKinPtr SudakovFormFactor::generateNextDecayBranching(const Energy startingScale,
						    const Energy stoppingScale,
						    const Energy minmass,
						    const IdList &ids,
						    const RhoDMatrix & rho, 
						    double enhance,
						    double detuning) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to this method.
  q_ = Constants::MaxEnergy;
  z_ = 0.;
  phi_ = 0.;
  // perform initialisation
  Energy2 tmax(sqr(stoppingScale)),tmin;
  initialize(ids,tmin);
  tmin=sqr(startingScale);
  // check some branching possible
  if(tmax<=tmin) return ShoKinPtr();
  // perform the evolution
  Energy2 t(tmin),pt2(-MeV2);
  do {
    if(!guessDecay(t,tmax,minmass,enhance,detuning)) break;
    pt2 = QTildeKinematics::pT2_Decay(t,z(),masssquared_[0],masssquared_[2]);
  }
  while(SplittingFnVeto((1.-z())*t/z(),ids,true,rho,detuning)|| 
	alphaSVeto(pTScale() ? sqr(1.-z())*t : (1.-z())*t ) ||
	pt2<cutoff_->pT2min() ||
	t*(1.-z())>masssquared_[0]-sqr(minmass));
  if(t > ZERO) {
    q_ = sqrt(t);
    pT_ = sqrt(pt2);
  }
  else return ShoKinPtr();
  phi_ = 0.;
  // create the ShowerKinematics object
  return new_ptr(Decay_QTildeShowerKinematics1to2(q_,z(),phi(),pT(),this)); 
}

bool SudakovFormFactor::guessDecay(Energy2 &t,Energy2 tmax, Energy minmass,
				   double enhance, double detune) {
  minmass = max(minmass,GeV);
  // previous scale
  Energy2 told = t;
  // overestimated limits on z
  if(tmax<masssquared_[0]) {
    t=-1.0*GeV2;
    return false;
  }
  Energy2 tm2 = tmax-masssquared_[0];
  Energy tm  = sqrt(tm2); 
  zlimits_ = make_pair(sqr(minmass/masses_[0]),
		       1.-sqrt(masssquared_[2]+cutoff_->pT2min()+
			       0.25*sqr(masssquared_[2])/tm2)/tm
		       +0.5*masssquared_[2]/tm2);
  if(zlimits_.second<zlimits_.first) {
    t=-1.0*GeV2;
    return false;
  }
  // guess values of t and z
  guesstz(told,2,ids_,enhance,ids_[1]==ids_[2],detune,t,z_);
  // actual values for z-limits
  if(t<masssquared_[0])  {
    t=-1.0*GeV2;
    return false;
  }
  tm2 = t-masssquared_[0];
  tm  = sqrt(tm2); 
  zlimits_ = make_pair(sqr(minmass/masses_[0]),
		   1.-sqrt(masssquared_[2]+cutoff_->pT2min()+
			   0.25*sqr(masssquared_[2])/tm2)/tm
		   +0.5*masssquared_[2]/tm2);
  if(t>tmax||zlimits_.second<zlimits_.first) {
    t=-1.0*GeV2;
    return false;
  }
  else
    return true; 
} 

bool SudakovFormFactor::computeTimeLikeLimits(Energy2 & t) {
  if (t < 1e-20 * GeV2) {
    t=-1.*GeV2;
    return false;
  }
  const Energy2 pT2min = cutoff_->pT2min();
  // special case for gluon radiating
  if(ids_[0]->id()==ParticleID::g||ids_[0]->id()==ParticleID::gamma) {
    // no emission possible
    if(t<16.*(masssquared_[1]+pT2min)) {
      t=-1.*GeV2;
      return false;
    }
    // overestimate of the limits
    zlimits_.first  = 0.5*(1.-sqrt(1.-4.*sqrt((masssquared_[1]+pT2min)/t)));
    zlimits_.second = 1.-zlimits_.first;
  }
  // special case for radiated particle is gluon 
  else if(ids_[2]->id()==ParticleID::g||ids_[2]->id()==ParticleID::gamma) {
    zlimits_.first  =    sqrt((masssquared_[1]+pT2min)/t);
    zlimits_.second = 1.-sqrt((masssquared_[2]+pT2min)/t);
  }
  else if(ids_[1]->id()==ParticleID::g||ids_[1]->id()==ParticleID::gamma) {
    zlimits_.second  =    sqrt((masssquared_[2]+pT2min)/t);
    zlimits_.first   = 1.-sqrt((masssquared_[1]+pT2min)/t);
  }
  else {
    zlimits_.first  =    (masssquared_[1]+pT2min)/t;
    zlimits_.second = 1.-(masssquared_[2]+pT2min)/t; 
  }
  if(zlimits_.first>=zlimits_.second) {
    t=-1.*GeV2;
    return false;
  }
  return true;
}

bool SudakovFormFactor::computeSpaceLikeLimits(Energy2 & t, double x) {
  if (t < 1e-20 * GeV2) {
    t=-1.*GeV2;
    return false;
  }
  // compute the limits
  zlimits_.first = x;
  double yy = 1.+0.5*masssquared_[2]/t;
  zlimits_.second = yy - sqrt(sqr(yy)-1.+cutoff_->pT2min()/t); 
  // return false if lower>upper
  if(zlimits_.second<zlimits_.first) {
    t=-1.*GeV2;
    return false;
  }
  else
    return true;
}

namespace {

tShowerParticlePtr findCorrelationPartner(ShowerParticle & particle,
					  bool forward,
					  ShowerInteraction inter) {
  tPPtr child = &particle;
  tShowerParticlePtr mother;
  if(forward) {
    mother = !particle.parents().empty() ? 
      dynamic_ptr_cast<tShowerParticlePtr>(particle.parents()[0]) : tShowerParticlePtr();
  }
  else {
    mother = particle.children().size()==2 ?
      dynamic_ptr_cast<tShowerParticlePtr>(&particle) : tShowerParticlePtr();
  }
  tShowerParticlePtr partner;
  while(mother) {
    tPPtr otherChild;
    if(forward) {
      for (unsigned int ix=0;ix<mother->children().size();++ix) {
	if(mother->children()[ix]!=child) {
	  otherChild = mother->children()[ix];
	  break;
	}
      }
    }
    else {
      otherChild = mother->children()[1];
    }
    tShowerParticlePtr other = dynamic_ptr_cast<tShowerParticlePtr>(otherChild);
    if((inter==ShowerInteraction::QCD && otherChild->dataPtr()->coloured()) ||
       (inter==ShowerInteraction::QED && otherChild->dataPtr()->charged())) {
      partner = other;
      break;
    }
    if(forward && !other->isFinalState()) {
      partner = dynamic_ptr_cast<tShowerParticlePtr>(mother);
      break;
    }
    child = mother;
    if(forward) {
      mother = ! mother->parents().empty() ?
	dynamic_ptr_cast<tShowerParticlePtr>(mother->parents()[0]) : tShowerParticlePtr();
    }
    else {
      if(mother->children()[0]->children().size()!=2)
	break;
      tShowerParticlePtr mtemp = 
	dynamic_ptr_cast<tShowerParticlePtr>(mother->children()[0]);
      if(!mtemp)
	break;
      else
	mother=mtemp;
    }
  }
  if(!partner) {
    if(forward) {
      partner = dynamic_ptr_cast<tShowerParticlePtr>( child)->partner();
    }
    else {
      if(mother) {
	tShowerParticlePtr parent;
	if(!mother->children().empty()) {
	  parent = dynamic_ptr_cast<tShowerParticlePtr>(mother->children()[0]);
	}
	if(!parent) {
	  parent = dynamic_ptr_cast<tShowerParticlePtr>(mother);
	}
	partner = parent->partner();
      }
      else {
	partner = dynamic_ptr_cast<tShowerParticlePtr>(&particle)->partner();
      }
    }
  }
  return partner;
}

pair<double,double> softPhiMin(double phi0, double phi1, double A, double B, double C, double D) {
  double c01 = cos(phi0 - phi1);
  double s01 = sin(phi0 - phi1);
  double s012(sqr(s01)), c012(sqr(c01));
  double A2(A*A), B2(B*B), C2(C*C), D2(D*D);
  if(abs(B/A)<1e-10 && abs(D/C)<1e-10) return make_pair(phi0,phi0+Constants::pi);
  double root = sqr(B2)*C2*D2*sqr(s012) + 2.*A*B2*B*C2*C*D*c01*s012 + 2.*A*B2*B*C*D2*D*c01*s012
		     + 4.*A2*B2*C2*D2*c012 - A2*B2*C2*D2*s012 - A2*B2*sqr(D2)*s012 - sqr(B2)*sqr(C2)*s012 
		     - sqr(B2)*C2*D2*s012 - 4.*A2*A*B*C*D2*D*c01 - 4.*A*B2*B*C2*C*D*c01 + sqr(A2)*sqr(D2)
		     + 2.*A2*B2*C2*D2 + sqr(B2)*sqr(C2);
  if(root<0.) return make_pair(phi0,phi0+Constants::pi);
  root = sqrt(root);
  double denom  = (-2.*A*B*C*D*c01 + A2*D2 + B2*C2);
  double denom2 = (-B*C*c01 + A*D);

  double num    = B2*C*D*s012;
  double y1 = B*s01*(-C*(num + root) + D*denom) / denom2;
  double y2 = B*s01*(-C*(num - root) + D*denom) / denom2;
  double x1 = -(num + root );
  double x2 = -(num - root );
  if(denom<0.) {
    y1*=-1.;
    y2*=-1.;
    x1*=-1.;
    x2*=-1.;
  }
  return make_pair(atan2(y1,x1) + phi0,atan2(y2,x2) + phi0);
}

}

double SudakovFormFactor::generatePhiForward(ShowerParticle & particle,
					 const IdList & ids,
					 ShoKinPtr kinematics,
					 const RhoDMatrix & rho) {
  // no correlations, return flat phi
  if(! dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->correlations())
    return Constants::twopi*UseRandom::rnd();
  // get the kinematic variables
  double  z = kinematics->z();
  Energy2 t = z*(1.-z)*sqr(kinematics->scale());
  Energy pT = kinematics->pT();
  // if soft correlations
  Energy2 pipj,pik;
  bool canBeSoft[2] = {ids[1]->id()==ParticleID::g || ids[1]->id()==ParticleID::gamma,
		       ids[2]->id()==ParticleID::g || ids[2]->id()==ParticleID::gamma };
  array<Energy2,3> pjk;
  array<Energy,3> Ek;
  Energy Ei,Ej;
  Energy2 m12(ZERO),m22(ZERO);
  InvEnergy2 aziMax(ZERO);
  bool softAllowed = dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()&&
    (canBeSoft[0] || canBeSoft[1]);
  if(softAllowed) {
    // find the partner for the soft correlations
    tShowerParticlePtr partner=findCorrelationPartner(particle,true,interactionType());
    // remember we want the softer gluon 
    bool swapOrder = !canBeSoft[1] || (canBeSoft[0] && canBeSoft[1] && z < 0.5);
    double zFact = !swapOrder ? (1.-z) : z;
    // compute the transforms to the shower reference frame
    // first the boost
    Lorentz5Momentum pVect = particle.showerBasis()->pVector();
    Lorentz5Momentum nVect = particle.showerBasis()->nVector();
    Boost beta_bb;
    if(particle.showerBasis()->frame()==ShowerBasis::BackToBack) {
      beta_bb = -(pVect + nVect).boostVector();
    }
    else if(particle.showerBasis()->frame()==ShowerBasis::Rest) {
      beta_bb = -pVect.boostVector();
    }
    else
      assert(false);
    pVect.boost(beta_bb);
    nVect.boost(beta_bb);
    Axis axis;
    if(particle.showerBasis()->frame()==ShowerBasis::BackToBack) {
      axis = pVect.vect().unit();
    }
    else if(particle.showerBasis()->frame()==ShowerBasis::Rest) {
      axis = nVect.vect().unit();
    }
    else
      assert(false);
    // and then the rotation
    LorentzRotation rot;
    if(axis.perp2()>0.) {
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      rot.rotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    }
    else if(axis.z()<0.) {
      rot.rotate(Constants::pi,Axis(1.,0.,0.));
    }
    rot.invert();
    pVect *= rot;
    nVect *= rot;
    // shower parameters
    Energy2 pn = pVect*nVect, m2 = pVect.m2();
    double alpha0 = particle.showerParameters().alpha;
    double  beta0 = 0.5/alpha0/pn*
      (sqr(particle.dataPtr()->mass())-sqr(alpha0)*m2+sqr(particle.showerParameters().pt));
    Lorentz5Momentum qperp0(particle.showerParameters().ptx,
			    particle.showerParameters().pty,ZERO,ZERO);
    assert(partner);
    Lorentz5Momentum pj = partner->momentum();
    pj.boost(beta_bb);
    pj *= rot;
    // compute the two phi independent dot products
    pik = 0.5*zFact*(sqr(alpha0)*m2 - sqr(particle.showerParameters().pt) + 2.*alpha0*beta0*pn )
      +0.5*sqr(pT)/zFact;
    Energy2 dot1 = pj*pVect;
    Energy2 dot2 = pj*nVect;
    Energy2 dot3 = pj*qperp0;
    pipj = alpha0*dot1+beta0*dot2+dot3;
    // compute the constants for the phi dependent dot product
    pjk[0] = zFact*(alpha0*dot1+dot3-0.5*dot2/pn*(alpha0*m2-sqr(particle.showerParameters().pt)/alpha0))
      +0.5*sqr(pT)*dot2/pn/zFact/alpha0;
    pjk[1] = (pj.x() - dot2/alpha0/pn*qperp0.x())*pT;
    pjk[2] = (pj.y() - dot2/alpha0/pn*qperp0.y())*pT;
    m12 = sqr(particle.dataPtr()->mass());
    m22 = sqr(partner->dataPtr()->mass());
    if(swapOrder) {
      pjk[1] *= -1.;
      pjk[2] *= -1.;
    }
    Ek[0] = zFact*(alpha0*pVect.t()-0.5*nVect.t()/pn*(alpha0*m2-sqr(particle.showerParameters().pt)/alpha0))
      +0.5*sqr(pT)*nVect.t()/pn/zFact/alpha0;
    Ek[1] = -nVect.t()/alpha0/pn*qperp0.x()*pT;
    Ek[2] = -nVect.t()/alpha0/pn*qperp0.y()*pT;
    if(swapOrder) {
      Ek[1] *= -1.;
      Ek[2] *= -1.;
    }
    Energy mag2=sqrt(sqr(Ek[1])+sqr(Ek[2]));
    Ei = alpha0*pVect.t()+beta0*nVect.t();
    Ej = pj.t();
    double phi0 = atan2(-pjk[2],-pjk[1]);
    if(phi0<0.) phi0 += Constants::twopi;
    double phi1 = atan2(-Ek[2],-Ek[1]);
    if(phi1<0.) phi1 += Constants::twopi;
    double xi_min = pik/Ei/(Ek[0]+mag2), xi_max = pik/Ei/(Ek[0]-mag2), xi_ij = pipj/Ei/Ej;
    if(xi_min>xi_max) swap(xi_min,xi_max);
    if(xi_min>xi_ij) softAllowed = false;
    Energy2 mag = sqrt(sqr(pjk[1])+sqr(pjk[2]));
    if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==1) {
      aziMax = -m12/sqr(pik) -m22/sqr(pjk[0]+mag) +2.*pipj/pik/(pjk[0]-mag);
    }
    else if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==2) {
      double A = (pipj*Ek[0]- Ej*pik)/Ej/sqr(Ej);
      double B = -sqrt(sqr(pipj)*(sqr(Ek[1])+sqr(Ek[2])))/Ej/sqr(Ej);
      double C = pjk[0]/sqr(Ej);
      double D = -sqrt(sqr(pjk[1])+sqr(pjk[2]))/sqr(Ej);
      pair<double,double> minima = softPhiMin(phi0,phi1,A,B,C,D);
      aziMax = 0.5/pik/(Ek[0]-mag2)*(Ei-m12*(Ek[0]-mag2)/pik  + max(Ej*(A+B*cos(minima.first -phi1))/(C+D*cos(minima.first -phi0)),
								    Ej*(A+B*cos(minima.second-phi1))/(C+D*cos(minima.second-phi0))));
    }
    else
      assert(false);
  }
  // if spin correlations
  vector<pair<int,Complex> > wgts;     
  if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->spinCorrelations()) {
    // calculate the weights
    wgts = generatePhiForward(z,t,ids,rho);
  }
  else {
    wgts = {{ {0, 1.} }};
  }
  // generate the azimuthal angle
  double phi,wgt;
  static const Complex ii(0.,1.);
  unsigned int ntry(0);
  double phiMax(0.),wgtMax(0.);
  do {
    phi = Constants::twopi*UseRandom::rnd();
    // first the spin correlations bit (gives 1 if correlations off)
    Complex spinWgt = 0.;
    for(unsigned int ix=0;ix<wgts.size();++ix) {
      if(wgts[ix].first==0)
  	spinWgt += wgts[ix].second;
      else
  	spinWgt += exp(double(wgts[ix].first)*ii*phi)*wgts[ix].second;
    }
    wgt = spinWgt.real();
    if(wgt-1.>1e-10) {
      generator()->log() << "Forward spin weight problem " << wgt << " " << wgt-1. 
			 << " " << ids[0]->id() << " " << ids[1]->id() << " " << ids[2]->id() << " " << " " << phi << "\n";
      generator()->log() << "Weights \n";
      for(unsigned int ix=0;ix<wgts.size();++ix)
	generator()->log() << wgts[ix].first << " " << wgts[ix].second << "\n";
    }
    // soft correlations bit
    double aziWgt = 1.;
    if(softAllowed) {
      Energy2 dot = pjk[0]+pjk[1]*cos(phi)+pjk[2]*sin(phi);
      Energy  Eg  = Ek[0]+Ek[1]*cos(phi)+Ek[2]*sin(phi);
      if(pipj*Eg>pik*Ej) {
	if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==1) {
	  aziWgt = (-m12/sqr(pik) -m22/sqr(dot) +2.*pipj/pik/dot)/aziMax;
	}
	else if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==2) {
	  aziWgt = max(ZERO,0.5/pik/Eg*(Ei-m12*Eg/pik  + (pipj*Eg - Ej*pik)/dot)/aziMax);
	}
	if(aziWgt-1.>1e-10||aziWgt<-1e-10) {
	  generator()->log() << "Forward soft weight problem " << aziWgt << " " << aziWgt-1. 
			     << " " << ids[0]->id() << " " << ids[1]->id() << " " << ids[2]->id() << " " << " " << phi << "\n";
	}
      }
      else {
	aziWgt = 0.;
      }
    }
    wgt *= aziWgt;
    if(wgt>wgtMax) {
      phiMax = phi;
      wgtMax = wgt;
    }
    ++ntry;
  }
  while(wgt<UseRandom::rnd()&&ntry<10000);
  if(ntry==10000) {
    generator()->log() << "Too many tries to generate phi in forward evolution\n";
    phi = phiMax;
  }
  // return the azimuthal angle
  return phi;
}

double SudakovFormFactor::generatePhiBackward(ShowerParticle & particle,
					  const IdList & ids,
					  ShoKinPtr kinematics,
					  const RhoDMatrix & rho) {
  // no correlations, return flat phi
  if(! dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->correlations())
    return Constants::twopi*UseRandom::rnd();
  // get the kinematic variables
  double z = kinematics->z();
  Energy2 t = (1.-z)*sqr(kinematics->scale())/z;
  Energy pT = kinematics->pT();
  // if soft correlations 
  bool softAllowed = dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations() &&
    (ids[2]->id()==ParticleID::g || ids[2]->id()==ParticleID::gamma);
  Energy2 pipj,pik,m12(ZERO),m22(ZERO);
  array<Energy2,3> pjk;
  Energy Ei,Ej,Ek;
  InvEnergy2 aziMax(ZERO);
  if(softAllowed) {
    // find the partner for the soft correlations
    tShowerParticlePtr partner=findCorrelationPartner(particle,false,interactionType());
    double zFact = (1.-z);
    // compute the transforms to the shower reference frame
    // first the boost
    Lorentz5Momentum pVect = particle.showerBasis()->pVector();
    Lorentz5Momentum nVect = particle.showerBasis()->nVector();
    assert(particle.showerBasis()->frame()==ShowerBasis::BackToBack);
    Boost beta_bb = -(pVect + nVect).boostVector();
    pVect.boost(beta_bb);
    nVect.boost(beta_bb);
    Axis axis = pVect.vect().unit();
    // and then the rotation
    LorentzRotation rot;
    if(axis.perp2()>0.) {
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      rot.rotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    }
    else if(axis.z()<0.) {
      rot.rotate(Constants::pi,Axis(1.,0.,0.));
    }
    rot.invert();
    pVect *= rot;
    nVect *= rot;
    // shower parameters
    Energy2 pn = pVect*nVect;
    Energy2 m2 = pVect.m2();
    double alpha0 = particle.x();
    double  beta0 = -0.5/alpha0/pn*sqr(alpha0)*m2;
    Lorentz5Momentum pj = partner->momentum();
    pj.boost(beta_bb);
    pj *= rot;
    double beta2 = 0.5*(1.-zFact)*(sqr(alpha0*zFact/(1.-zFact))*m2+sqr(pT))/alpha0/zFact/pn; 
    // compute the two phi independent dot products
    Energy2 dot1 = pj*pVect;
    Energy2 dot2 = pj*nVect;
    pipj = alpha0*dot1+beta0*dot2;
    pik  = alpha0*(alpha0*zFact/(1.-zFact)*m2+pn*(beta2+zFact/(1.-zFact)*beta0));
    // compute the constants for the phi dependent dot product
    pjk[0] = alpha0*zFact/(1.-zFact)*dot1+beta2*dot2;
    pjk[1] = pj.x()*pT;
    pjk[2] = pj.y()*pT;
    m12 = ZERO;
    m22 = sqr(partner->dataPtr()->mass());
    Energy2 mag = sqrt(sqr(pjk[1])+sqr(pjk[2]));
    if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==1) {
      aziMax = -m12/sqr(pik) -m22/sqr(pjk[0]+mag) +2.*pipj/pik/(pjk[0]-mag);
    }
    else if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==2) {
      Ek = alpha0*zFact/(1.-zFact)*pVect.t()+beta2*nVect.t();
      Ei = alpha0*pVect.t()+beta0*nVect.t();
      Ej = pj.t();
      if(pipj*Ek> Ej*pik) {
	aziMax = 0.5/pik/Ek*(Ei-m12*Ek/pik  + (pipj*Ek- Ej*pik)/(pjk[0]-mag));
      }
      else {
	aziMax = 0.5/pik/Ek*(Ei-m12*Ek/pik);
      }
    }
    else {
      assert(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==0);
    }
  }
  // if spin correlations
  vector<pair<int,Complex> > wgts;
  if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->spinCorrelations()) {
    // get the weights
    wgts = generatePhiBackward(z,t,ids,rho);
  }
  else {
    wgts = {{ {0, 1.} }};
  }
  // generate the azimuthal angle
  double phi,wgt;
  static const Complex ii(0.,1.);
  unsigned int ntry(0);
  double phiMax(0.),wgtMax(0.);
  do {
    phi = Constants::twopi*UseRandom::rnd();
    Complex spinWgt = 0.;
    for(unsigned int ix=0;ix<wgts.size();++ix) {
      if(wgts[ix].first==0)
  	spinWgt += wgts[ix].second;
      else
  	spinWgt += exp(double(wgts[ix].first)*ii*phi)*wgts[ix].second;
    }
    wgt = spinWgt.real();
    if(wgt-1.>1e-10) {
      generator()->log() << "Backward weight problem " << wgt << " " << wgt-1. 
			 << " " << ids[0]->id() << " " << ids[1]->id() << " " << ids[2]->id() << " " << " " << z << " " << phi << "\n";
      generator()->log() << "Weights \n";
      for(unsigned int ix=0;ix<wgts.size();++ix)
  	generator()->log() << wgts[ix].first << " " << wgts[ix].second << "\n";
    }
    // soft correlations bit
    double aziWgt = 1.;
    if(softAllowed) {
      Energy2 dot = pjk[0]+pjk[1]*cos(phi)+pjk[2]*sin(phi);
      if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==1) {
	aziWgt = (-m12/sqr(pik) -m22/sqr(dot) +2.*pipj/pik/dot)/aziMax;
      }
      else if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==2) {
	aziWgt = max(ZERO,0.5/pik/Ek*(Ei-m12*Ek/pik  + pipj*Ek/dot - Ej*pik/dot)/aziMax);
      }
      if(aziWgt-1.>1e-10||aziWgt<-1e-10) {
 	generator()->log() << "Backward soft weight problem " << aziWgt << " " << aziWgt-1. 
 			   << " " << ids[0]->id() << " " << ids[1]->id() << " " << ids[2]->id() << " " << " " << phi << "\n";
      }
    }
    wgt *= aziWgt;
    if(wgt>wgtMax) {
      phiMax = phi;
      wgtMax = wgt;
    }
    ++ntry;
  }
  while(wgt<UseRandom::rnd()&&ntry<10000);
  if(ntry==10000) {
    generator()->log() << "Too many tries to generate phi in backward evolution\n";
    phi = phiMax;
  }
  // return the azimuthal angle
  return phi;
}

double SudakovFormFactor::generatePhiDecay(ShowerParticle & particle,
				       const IdList & ids,
				       ShoKinPtr kinematics,
				       const RhoDMatrix &) {
  // only soft correlations in this case
  // no correlations, return flat phi
  if( !(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations() &&
	(ids[2]->id()==ParticleID::g || ids[2]->id()==ParticleID::gamma )))
    return Constants::twopi*UseRandom::rnd();
  // get the kinematic variables
  double  z = kinematics->z();
  Energy pT = kinematics->pT();
  // if soft correlations
  // find the partner for the soft correlations
  tShowerParticlePtr partner = findCorrelationPartner(particle,true,interactionType());
  double zFact(1.-z);
  // compute the transforms to the shower reference frame
  // first the boost
  Lorentz5Momentum pVect = particle.showerBasis()->pVector();
  Lorentz5Momentum nVect = particle.showerBasis()->nVector();
  assert(particle.showerBasis()->frame()==ShowerBasis::Rest);
  Boost beta_bb = -pVect.boostVector();
  pVect.boost(beta_bb);
  nVect.boost(beta_bb);
  Axis axis = nVect.vect().unit();
  // and then the rotation
  LorentzRotation rot;
  if(axis.perp2()>0.) {
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    rot.rotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  }
  else if(axis.z()<0.) {
    rot.rotate(Constants::pi,Axis(1.,0.,0.));
  }
  rot.invert();
  pVect *= rot;
  nVect *= rot;
  // shower parameters
  Energy2 pn = pVect*nVect;
  Energy2 m2 = pVect.m2();
  double alpha0 = particle.showerParameters().alpha;
  double  beta0 = 0.5/alpha0/pn*
    (sqr(particle.dataPtr()->mass())-sqr(alpha0)*m2+sqr(particle.showerParameters().pt));
  Lorentz5Momentum qperp0(particle.showerParameters().ptx,
			  particle.showerParameters().pty,ZERO,ZERO);
  Lorentz5Momentum pj = partner->momentum();
  pj.boost(beta_bb);
  pj *= rot;
  // compute the two phi independent dot products
  Energy2 pik = 0.5*zFact*(sqr(alpha0)*m2 - sqr(particle.showerParameters().pt) + 2.*alpha0*beta0*pn )
    +0.5*sqr(pT)/zFact;
  Energy2 dot1 = pj*pVect;
  Energy2 dot2 = pj*nVect;
  Energy2 dot3 = pj*qperp0;
  Energy2 pipj = alpha0*dot1+beta0*dot2+dot3;
  // compute the constants for the phi dependent dot product
  array<Energy2,3> pjk;
  pjk[0] = zFact*(alpha0*dot1+dot3-0.5*dot2/pn*(alpha0*m2-sqr(particle.showerParameters().pt)/alpha0))
    +0.5*sqr(pT)*dot2/pn/zFact/alpha0;
  pjk[1] = (pj.x() - dot2/alpha0/pn*qperp0.x())*pT;
  pjk[2] = (pj.y() - dot2/alpha0/pn*qperp0.y())*pT;
  Energy2 m12 = sqr(particle.dataPtr()->mass());
  Energy2 m22 = sqr(partner->dataPtr()->mass());
  Energy2 mag = sqrt(sqr(pjk[1])+sqr(pjk[2]));
  InvEnergy2 aziMax;
  array<Energy,3> Ek;
  Energy Ei,Ej;
  if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==1) {
    aziMax = -m12/sqr(pik) -m22/sqr(pjk[0]+mag) +2.*pipj/pik/(pjk[0]-mag);
  }
  else if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==2) {
    Ek[0] = zFact*(alpha0*pVect.t()+-0.5*nVect.t()/pn*(alpha0*m2-sqr(particle.showerParameters().pt)/alpha0))
      +0.5*sqr(pT)*nVect.t()/pn/zFact/alpha0;
    Ek[1] = -nVect.t()/alpha0/pn*qperp0.x()*pT;
    Ek[2] = -nVect.t()/alpha0/pn*qperp0.y()*pT;
    Energy mag2=sqrt(sqr(Ek[1])+sqr(Ek[2]));
    Ei = alpha0*pVect.t()+beta0*nVect.t();
    Ej = pj.t();
    aziMax = 0.5/pik/(Ek[0]-mag2)*(Ei-m12*(Ek[0]-mag2)/pik  + pipj*(Ek[0]+mag2)/(pjk[0]-mag) - Ej*pik/(pjk[0]-mag) );
  }
  else
    assert(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==0);
  // generate the azimuthal angle
  double phi,wgt(0.);
  unsigned int ntry(0);
  double phiMax(0.),wgtMax(0.);
  do {
    phi = Constants::twopi*UseRandom::rnd();
    Energy2 dot = pjk[0]+pjk[1]*cos(phi)+pjk[2]*sin(phi);
    if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==1) {
      wgt = (-m12/sqr(pik) -m22/sqr(dot) +2.*pipj/pik/dot)/aziMax;
    }
    else if(dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->softCorrelations()==2) {
      if(qperp0.m2()==ZERO) {
	wgt = 1.;
      }
      else {
	Energy Eg = Ek[0]+Ek[1]*cos(phi)+Ek[2]*sin(phi);
	wgt = max(ZERO,0.5/pik/Eg*(Ei-m12*Eg/pik  + (pipj*Eg - Ej*pik)/dot)/aziMax);
      }	
    }
    if(wgt-1.>1e-10||wgt<-1e-10) {
      generator()->log() << "Decay soft weight problem " << wgt << " " << wgt-1. 
			 << " " << ids[0]->id() << " " << ids[1]->id() << " " << ids[2]->id() << " " << " " << phi << "\n";
    }
    if(wgt>wgtMax) {
      phiMax = phi;
      wgtMax = wgt;
    }
    ++ntry;
  }
  while(wgt<UseRandom::rnd()&&ntry<10000);
  if(ntry==10000) {
    phi = phiMax;
    generator()->log() << "Too many tries to generate phi\n";
  }
  // return the azimuthal angle
  return phi;
}


Energy SudakovFormFactor::calculateScale(double zin, Energy pt, IdList ids,
				     unsigned int iopt) {
  Energy2 tmin;
  initialize(ids,tmin);
  // final-state branching
  if(iopt==0) {
    Energy2 scale=(sqr(pt)+masssquared_[1]*(1.-zin)+masssquared_[2]*zin);
    if(ids[0]->id()!=ParticleID::g) scale -= zin*(1.-zin)*masssquared_[0];
    scale /= sqr(zin*(1-zin));
    return scale<=ZERO ? sqrt(tmin) : sqrt(scale);
  }
  else if(iopt==1) {
    Energy2 scale=(sqr(pt)+zin*masssquared_[2])/sqr(1.-zin);
    return scale<=ZERO ? sqrt(tmin) : sqrt(scale);
  }
  else if(iopt==2) {
    Energy2 scale = (sqr(pt)+zin*masssquared_[2])/sqr(1.-zin)+masssquared_[0];
    return scale<=ZERO ? sqrt(tmin) : sqrt(scale);
  }
  else {
    throw Exception() << "Unknown option in SudakovFormFactor::calculateScale() "
		      << "iopt = " << iopt << Exception::runerror;
  }
}

void SudakovFormFactor::colourConnection(tShowerParticlePtr parent,
                                         tShowerParticlePtr first,
                                         tShowerParticlePtr second,
					 ShowerPartnerType partnerType, 
                                         const bool back) const {
  if(colourStructure_==TripletTripletOctet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
                                      parent->antiColourLine());
      // ensure input consistency
      assert((  cparent.first && !cparent.second && 
		partnerType==ShowerPartnerType::QCDColourLine) || 
             ( !cparent.first &&  cparent.second && 
		partnerType==ShowerPartnerType::QCDAntiColourLine));
      // q -> q g
      if(cparent.first) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(second);
        newline->addColoured     ( first);
        newline->addAntiColoured (second);
      }
      // qbar -> qbar g
      else {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.second->addAntiColoured(second);
        newline->addColoured(second);
        newline->addAntiColoured(first);
      }
      // Set progenitor
      first->progenitor(parent->progenitor());
      second->progenitor(parent->progenitor());
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
                                     first->antiColourLine());
      // ensure input consistency
      assert((  cfirst.first && !cfirst.second && 
		partnerType==ShowerPartnerType::QCDColourLine) || 
             ( !cfirst.first &&  cfirst.second && 
		partnerType==ShowerPartnerType::QCDAntiColourLine));
      // q -> q g
      if(cfirst.first) {
        ColinePtr newline=new_ptr(ColourLine());
        cfirst.first->addAntiColoured(second);
        newline->addColoured(second);
        newline->addColoured(parent);
      }
      // qbar -> qbar g
      else {
        ColinePtr newline=new_ptr(ColourLine());
        cfirst.second->addColoured(second);
        newline->addAntiColoured(second);
        newline->addAntiColoured(parent);
      }
      // Set progenitor
      parent->progenitor(first->progenitor());
      second->progenitor(first->progenitor());
    }
  }
  else if(colourStructure_==OctetOctetOctet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
                                      parent->antiColourLine());
      // ensure input consistency
      assert(cparent.first&&cparent.second);
      // ensure first gluon is hardest
      if( first->id()==second->id() && parent->showerKinematics()->z()<0.5 ) 
	swap(first,second);
      // colour line radiates
      if(partnerType==ShowerPartnerType::QCDColourLine) {
	// The colour line is radiating
	ColinePtr newline=new_ptr(ColourLine());
	cparent.first->addColoured(second);
	cparent.second->addAntiColoured(first);
	newline->addColoured(first);
	newline->addAntiColoured(second);
      }
      // anti colour line radiates
      else if(partnerType==ShowerPartnerType::QCDAntiColourLine) {
	ColinePtr newline=new_ptr(ColourLine());
	cparent.first->addColoured(first);
	cparent.second->addAntiColoured(second);
	newline->addColoured(second);
	newline->addAntiColoured(first);
      }
      else
	assert(false);
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
                                     first->antiColourLine());
      // ensure input consistency
      assert(cfirst.first&&cfirst.second);
      // The colour line is radiating
      if(partnerType==ShowerPartnerType::QCDColourLine) {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.first->addAntiColoured(second);
	cfirst.second->addAntiColoured(parent);
	newline->addColoured(parent);
	newline->addColoured(second);
      }
      // anti colour line radiates
      else if(partnerType==ShowerPartnerType::QCDAntiColourLine) {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.first->addColoured(parent);
	cfirst.second->addColoured(second);
	newline->addAntiColoured(second);
	newline->addAntiColoured(parent);
      }
      else
	assert(false);
    }    
  }
  else if(colourStructure_ == OctetTripletTriplet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
                                      parent->antiColourLine());
      // ensure input consistency
      assert(cparent.first&&cparent.second);
      cparent.first ->addColoured    ( first);
      cparent.second->addAntiColoured(second);
      // Set progenitor
      first->progenitor(parent->progenitor());
      second->progenitor(parent->progenitor());
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
                                     first->antiColourLine());
      // ensure input consistency
      assert(( cfirst.first && !cfirst.second) ||
             (!cfirst.first &&  cfirst.second));
      // g -> q qbar
      if(cfirst.first) {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.first->addColoured(parent);
	newline->addAntiColoured(second);
	newline->addAntiColoured(parent);	
      }
      // g -> qbar q
      else {
        ColinePtr newline=new_ptr(ColourLine());
        cfirst.second->addAntiColoured(parent);
        newline->addColoured(second);
        newline->addColoured(parent);
      }
      // Set progenitor
      parent->progenitor(first->progenitor());
      second->progenitor(first->progenitor());
    }
  }
  else if(colourStructure_ == TripletOctetTriplet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
                                      parent->antiColourLine());
      // ensure input consistency
      assert(( cparent.first && !cparent.second) || 
             (!cparent.first &&  cparent.second));
      // q -> g q
      if(cparent.first) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(first);
        newline->addColoured    (second);
        newline->addAntiColoured( first);
      }
      // qbar -> g qbar
      else {
	ColinePtr newline=new_ptr(ColourLine());
	cparent.second->addAntiColoured(first);
	newline->addColoured    ( first);
	newline->addAntiColoured(second);	
      }
      // Set progenitor
      first->progenitor(parent->progenitor());
      second->progenitor(parent->progenitor());
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
                                     first->antiColourLine());
      // ensure input consistency
      assert(cfirst.first&&cfirst.second);
      // q -> g q
      if(parent->id()>0) {
        cfirst.first ->addColoured(parent);
        cfirst.second->addColoured(second);
      }
      else {
        cfirst.first ->addAntiColoured(second);
        cfirst.second->addAntiColoured(parent);
      }
      // Set progenitor
      parent->progenitor(first->progenitor());
      second->progenitor(first->progenitor());
    }
  }
  else if(colourStructure_==SextetSextetOctet) {
    //make sure we're not doing backward evolution
    assert(!back);

    //make sure something sensible
    assert(parent->colourLine() || parent->antiColourLine());
   
    //get the colour lines or anti-colour lines
    bool isAntiColour=true;
    ColinePair cparent;
    if(parent->colourLine()) {
      cparent = ColinePair(const_ptr_cast<tColinePtr>(parent->colourInfo()->colourLines()[0]), 
			   const_ptr_cast<tColinePtr>(parent->colourInfo()->colourLines()[1]));
      isAntiColour=false;
    }
    else {
      cparent = ColinePair(const_ptr_cast<tColinePtr>(parent->colourInfo()->antiColourLines()[0]), 
			   const_ptr_cast<tColinePtr>(parent->colourInfo()->antiColourLines()[1]));
    }
    
    //check for sensible input
    //    assert(cparent.first && cparent.second);

    // sextet has 2 colour lines
    if(!isAntiColour) {
      //pick at random which of the colour topolgies to take
      double topology = UseRandom::rnd();
      if(topology < 0.25) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(second);
        cparent.second->addColoured(first);
        newline->addColoured(first);
        newline->addAntiColoured(second);
      }
      else if(topology >=0.25 && topology < 0.5) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(first);
        cparent.second->addColoured(second);
        newline->addColoured(first);
        newline->addAntiColoured(second); 
      }
      else if(topology >= 0.5 && topology < 0.75) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(second);
        cparent.second->addColoured(first); 
        newline->addColoured(first); 
        newline->addAntiColoured(second); 
      }
      else {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(first);
        cparent.second->addColoured(second);
        newline->addColoured(first);
        newline->addAntiColoured(second);
      }
    }
    // sextet has 2 anti-colour lines
    else {
      double topology = UseRandom::rnd();
      if(topology < 0.25){
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addAntiColoured(second);
        cparent.second->addAntiColoured(first);
        newline->addAntiColoured(first);
        newline->addColoured(second);
      }
      else if(topology >=0.25 && topology < 0.5) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addAntiColoured(first);
        cparent.second->addAntiColoured(second);
        newline->addAntiColoured(first);
        newline->addColoured(second); 
      }
      else if(topology >= 0.5 && topology < 0.75) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addAntiColoured(second);
        cparent.second->addAntiColoured(first);
        newline->addAntiColoured(first);
        newline->addColoured(second); 
      }
      else {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addAntiColoured(first);
        cparent.second->addAntiColoured(second);
        newline->addAntiColoured(first);
        newline->addColoured(second);
      }
    }   
  }
  else if(colourStructure_ == ChargedChargedNeutral) {
    if(!parent->data().coloured()) return;
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
				      parent->antiColourLine());
      // q -> q g
      if(cparent.first) {
	cparent.first->addColoured(first);
      }
      // qbar -> qbar g
      if(cparent.second) {
	cparent.second->addAntiColoured(first);
      }
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
				     first->antiColourLine());
      // q -> q g
      if(cfirst.first) {
	cfirst.first->addColoured(parent);
      }
      // qbar -> qbar g
      if(cfirst.second) {
	cfirst.second->addAntiColoured(parent);
      }
    }
  }
  else if(colourStructure_ == ChargedNeutralCharged) {
    if(!parent->data().coloured()) return;
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
				      parent->antiColourLine());
      // q -> q g
      if(cparent.first) {
	cparent.first->addColoured(second);
      }
      // qbar -> qbar g
      if(cparent.second) {
	cparent.second->addAntiColoured(second);
      }
    }
    else {
      if (second->dataPtr()->iColour()==PDT::Colour3 ) {
	ColinePtr newline=new_ptr(ColourLine());
	newline->addColoured(second);
	newline->addColoured(parent);
      }
      else if (second->dataPtr()->iColour()==PDT::Colour3bar ) {
	ColinePtr newline=new_ptr(ColourLine());
	newline->addAntiColoured(second);
	newline->addAntiColoured(parent);
      }
    }
  }
  else if(colourStructure_ == NeutralChargedCharged ) {
    if(!back) {
      if(first->dataPtr()->coloured()) {
	ColinePtr newline=new_ptr(ColourLine());
	if(first->dataPtr()->iColour()==PDT::Colour3) {
	  newline->addColoured    (first );
	  newline->addAntiColoured(second);
	}
	else if (first->dataPtr()->iColour()==PDT::Colour3bar) {
	  newline->addColoured    (second);
	  newline->addAntiColoured(first );
	}
	else if(parent->dataPtr()->coloured()||
		first ->dataPtr()->coloured()||
		second->dataPtr()->coloured())
	  assert(false);
      }
    }
    else {   
      ColinePair cfirst = ColinePair(first->colourLine(), 
				     first->antiColourLine());
      // gamma -> q qbar
      if(cfirst.first) {
	cfirst.first->addAntiColoured(second);
      }
      // gamma -> qbar q
      else if(cfirst.second) {
	cfirst.second->addColoured(second);
      }
      else if(parent->dataPtr()->coloured()||
	      first ->dataPtr()->coloured()||
	      second->dataPtr()->coloured())
	assert(false);
    }
  }
  else {
    assert(false);
  }
}

namespace {

  bool hasColour(tPPtr p) {
    PDT::Colour colour = p->dataPtr()->iColour();
    return colour==PDT::Colour3 || colour==PDT::Colour8 || colour == PDT::Colour6;
  }

  bool hasAntiColour(tPPtr p) {
    PDT::Colour colour = p->dataPtr()->iColour();
    return colour==PDT::Colour3bar || colour==PDT::Colour8 || colour == PDT::Colour6bar;
  }
  
}

void SudakovFormFactor::evaluateFinalStateScales(ShowerPartnerType partnerType,
						 Energy scale, double z,
						 tShowerParticlePtr parent,
						 tShowerParticlePtr emitter,
						 tShowerParticlePtr emitted) {
  // identify emitter and emitted
  double zEmitter = z, zEmitted = 1.-z;
  bool bosonSplitting(false);
  // special for g -> gg, particle highest z is emitter
  if(emitter->id() == emitted->id() && emitter->id() == parent->id() &&
     zEmitted > zEmitter) {
    swap(zEmitted,zEmitter);
    swap( emitted, emitter);
  }
  // otherwise if particle ID same
  else if(emitted->id()==parent->id()) {
    swap(zEmitted,zEmitter);
    swap( emitted, emitter);
  }
  // no real emitter/emitted
  else if(emitter->id()!=parent->id()) {
    bosonSplitting = true;
  }
  // may need to add angularOrder flag here
  // now the various scales
  // QED
  if(partnerType==ShowerPartnerType::QED) {
    assert(colourStructure()==ChargedChargedNeutral ||
	   colourStructure()==ChargedNeutralCharged ||
	   colourStructure()==NeutralChargedCharged );
    // normal case
    if(!bosonSplitting) {
      assert(colourStructure()==ChargedChargedNeutral ||
	     colourStructure()==ChargedNeutralCharged );
      // set the scales
      // emitter
      emitter->scales().QED         = zEmitter*scale;
      emitter->scales().QED_noAO    =          scale;
      if(strictAO_)
	emitter->scales().QCD_c       = min(zEmitter*scale,parent->scales().QCD_c      );
      else
	emitter->scales().QCD_c       = min(         scale,parent->scales().QCD_c      );
      emitter->scales().QCD_c_noAO  = min(scale,parent->scales().QCD_c_noAO );
      if(strictAO_)
	emitter->scales().QCD_ac      = min(zEmitter*scale,parent->scales().QCD_ac     );
      else
	emitter->scales().QCD_ac      = min(         scale,parent->scales().QCD_ac     );
      emitter->scales().QCD_ac_noAO = min(scale,parent->scales().QCD_ac_noAO);
      // emitted 
      emitted->scales().QED         = zEmitted*scale;
      emitted->scales().QED_noAO    =          scale;
      emitted->scales().QCD_c       = ZERO;
      emitted->scales().QCD_c_noAO  = ZERO;
      emitted->scales().QCD_ac      = ZERO;
      emitted->scales().QCD_ac_noAO = ZERO;
    }
    // gamma -> f fbar
    else {
      assert(colourStructure()==NeutralChargedCharged );
      // emitter
      emitter->scales().QED         = zEmitter*scale;
      emitter->scales().QED_noAO    =          scale;
      if(hasColour(emitter)) {
	emitter->scales().QCD_c       = zEmitter*scale;
	emitter->scales().QCD_c_noAO  =          scale;
      }
      if(hasAntiColour(emitter)) {
	emitter->scales().QCD_ac      = zEmitter*scale;
	emitter->scales().QCD_ac_noAO =          scale;
      }
      // emitted 
      emitted->scales().QED         = zEmitted*scale;
      emitted->scales().QED_noAO    =          scale;
      if(hasColour(emitted)) {
	emitted->scales().QCD_c       = zEmitted*scale;
	emitted->scales().QCD_c_noAO  =          scale;
      }
      if(hasAntiColour(emitted)) {
	emitted->scales().QCD_ac      = zEmitted*scale;
	emitted->scales().QCD_ac_noAO =          scale;
      }
    }
  }
  // QCD
  else {
   // normal case eg q -> q g and g -> g g
    if(!bosonSplitting) {
      if(strictAO_)
	emitter->scales().QED         = min(zEmitter*scale,parent->scales().QED     );
      else
	emitter->scales().QED         = min(         scale,parent->scales().QED     );
      emitter->scales().QED_noAO    = min(scale,parent->scales().QED_noAO);
      if(partnerType==ShowerPartnerType::QCDColourLine) {
	emitter->scales().QCD_c       = zEmitter*scale;
	emitter->scales().QCD_c_noAO  =          scale;
	emitter->scales().QCD_ac      = min(zEmitter*scale,parent->scales().QCD_ac     );
	emitter->scales().QCD_ac_noAO = min(         scale,parent->scales().QCD_ac_noAO);
      }
      else {
	emitter->scales().QCD_c       = min(zEmitter*scale,parent->scales().QCD_c      );
	emitter->scales().QCD_c_noAO  = min(         scale,parent->scales().QCD_c_noAO );
	emitter->scales().QCD_ac      = zEmitter*scale;
	emitter->scales().QCD_ac_noAO =          scale;
      }
      // emitted 
      emitted->scales().QED         = ZERO;
      emitted->scales().QED_noAO    = ZERO;
      emitted->scales().QCD_c       = zEmitted*scale;
      emitted->scales().QCD_c_noAO  =          scale;
      emitted->scales().QCD_ac      = zEmitted*scale;
      emitted->scales().QCD_ac_noAO =          scale;
    }
    // g -> q qbar
    else {
      // emitter
      if(emitter->dataPtr()->charged()) {
	emitter->scales().QED         = zEmitter*scale;
	emitter->scales().QED_noAO    =          scale;
      }
      emitter->scales().QCD_c       = zEmitter*scale;
      emitter->scales().QCD_c_noAO  =          scale;
      emitter->scales().QCD_ac      = zEmitter*scale;
      emitter->scales().QCD_ac_noAO =          scale;
      // emitted 
      if(emitted->dataPtr()->charged()) {
	emitted->scales().QED         = zEmitted*scale;
	emitted->scales().QED_noAO    =          scale;
      }
      emitted->scales().QCD_c       = zEmitted*scale;
      emitted->scales().QCD_c_noAO  =          scale;
      emitted->scales().QCD_ac      = zEmitted*scale;
      emitted->scales().QCD_ac_noAO =          scale;
    }
  }
}

void SudakovFormFactor::evaluateInitialStateScales(ShowerPartnerType partnerType,
						   Energy scale, double z,
						   tShowerParticlePtr parent,
						   tShowerParticlePtr spacelike,
						   tShowerParticlePtr timelike) {
  // scale for time-like child
  Energy AOScale = (1.-z)*scale;
  // QED
  if(partnerType==ShowerPartnerType::QED) {
    if(parent->id()==spacelike->id()) {
      // parent
      parent   ->scales().QED         =   scale;
      parent   ->scales().QED_noAO    =   scale;
      parent   ->scales().QCD_c       = min(scale,spacelike->scales().QCD_c      );
      parent   ->scales().QCD_c_noAO  = min(scale,spacelike->scales().QCD_c_noAO );
      parent   ->scales().QCD_ac      = min(scale,spacelike->scales().QCD_ac     );
      parent   ->scales().QCD_ac_noAO = min(scale,spacelike->scales().QCD_ac_noAO);
      // timelike
      timelike->scales().QED         = AOScale;
      timelike->scales().QED_noAO    =   scale;
      timelike->scales().QCD_c       =    ZERO;
      timelike->scales().QCD_c_noAO  =    ZERO;
      timelike->scales().QCD_ac      =    ZERO;
      timelike->scales().QCD_ac_noAO =    ZERO;
    }
    else if(parent->id()==timelike->id()) {
      parent   ->scales().QED         =   scale;
      parent   ->scales().QED_noAO    =   scale;
      if(hasColour(parent)) {
	parent   ->scales().QCD_c       = scale;
	parent   ->scales().QCD_c_noAO  = scale;
      }
      if(hasAntiColour(parent)) {
	parent   ->scales().QCD_ac      = scale;
	parent   ->scales().QCD_ac_noAO = scale;
      }
      // timelike 
      timelike->scales().QED         = AOScale;
      timelike->scales().QED_noAO    =   scale;
      if(hasColour(timelike)) {
	timelike->scales().QCD_c       = AOScale;
	timelike->scales().QCD_c_noAO  =   scale;
      }
      if(hasAntiColour(timelike)) {
	timelike->scales().QCD_ac      = AOScale;
	timelike->scales().QCD_ac_noAO =   scale;
      }
    }
    else {
      parent   ->scales().QED         = scale;
      parent   ->scales().QED_noAO    = scale;
      parent   ->scales().QCD_c       = ZERO ;
      parent   ->scales().QCD_c_noAO  = ZERO ;
      parent   ->scales().QCD_ac      = ZERO ;
      parent   ->scales().QCD_ac_noAO = ZERO ;
      // timelike 
      timelike->scales().QED         = AOScale;
      timelike->scales().QED_noAO    =   scale;
      if(hasColour(timelike)) {
	timelike->scales().QCD_c       = min(AOScale,spacelike->scales().QCD_ac     );
	timelike->scales().QCD_c_noAO  = min(  scale,spacelike->scales().QCD_ac_noAO);
      }
      if(hasAntiColour(timelike)) {
	timelike->scales().QCD_ac      = min(AOScale,spacelike->scales().QCD_c      );
	timelike->scales().QCD_ac_noAO = min(  scale,spacelike->scales().QCD_c_noAO );
      }
    }
  }
  // QCD
  else {
    // timelike 
    if(timelike->dataPtr()->charged()) {
      timelike->scales().QED         = AOScale;
      timelike->scales().QED_noAO    =   scale;
    }
    if(hasColour(timelike)) {
      timelike->scales().QCD_c       = AOScale;
      timelike->scales().QCD_c_noAO  =   scale;
    }
    if(hasAntiColour(timelike)) {
      timelike->scales().QCD_ac      = AOScale;
      timelike->scales().QCD_ac_noAO =   scale;
    }
    if(parent->id()==spacelike->id()) {
      parent   ->scales().QED         = min(scale,spacelike->scales().QED        );
      parent   ->scales().QED_noAO    = min(scale,spacelike->scales().QED_noAO   );
      parent   ->scales().QCD_c       = min(scale,spacelike->scales().QCD_c      );
      parent   ->scales().QCD_c_noAO  = min(scale,spacelike->scales().QCD_c_noAO );
      parent   ->scales().QCD_ac      = min(scale,spacelike->scales().QCD_ac     );
      parent   ->scales().QCD_ac_noAO = min(scale,spacelike->scales().QCD_ac_noAO);
    }
    else {
      if(parent->dataPtr()->charged()) {
	parent   ->scales().QED         = scale;
	parent   ->scales().QED_noAO    = scale;
      }
      if(hasColour(parent)) {
	parent   ->scales().QCD_c      = scale;
	parent   ->scales().QCD_c_noAO  = scale;
      }
      if(hasAntiColour(parent)) {
	parent   ->scales().QCD_ac      = scale;
	parent   ->scales().QCD_ac_noAO = scale;
      }
    }
  }
}

void SudakovFormFactor::evaluateDecayScales(ShowerPartnerType partnerType,
					    Energy scale, double z,
					    tShowerParticlePtr parent,
					    tShowerParticlePtr spacelike,
					    tShowerParticlePtr timelike) {
  assert(parent->id()==spacelike->id());
  // angular-ordered scale for 2nd child
  Energy AOScale = (1.-z)*scale;
  // QED
  if(partnerType==ShowerPartnerType::QED) {
    // timelike
    timelike->scales().QED         = AOScale;
    timelike->scales().QED_noAO    =   scale;
    timelike->scales().QCD_c       =    ZERO;
    timelike->scales().QCD_c_noAO  =    ZERO;
    timelike->scales().QCD_ac      =    ZERO;
    timelike->scales().QCD_ac_noAO =    ZERO;
    // spacelike
    spacelike->scales().QED         =   scale;
    spacelike->scales().QED_noAO    =   scale;
  }
  // QCD
  else {
    // timelike 
    timelike->scales().QED         = ZERO;
    timelike->scales().QED_noAO    = ZERO;
    timelike->scales().QCD_c       = AOScale;
    timelike->scales().QCD_c_noAO  =   scale;
    timelike->scales().QCD_ac      = AOScale;
    timelike->scales().QCD_ac_noAO =   scale;
    // spacelike
    spacelike->scales().QED         = max(scale,parent->scales().QED        );
    spacelike->scales().QED_noAO    = max(scale,parent->scales().QED_noAO   );
  }
  spacelike->scales().QCD_c       = max(scale,parent->scales().QCD_c      );
  spacelike->scales().QCD_c_noAO  = max(scale,parent->scales().QCD_c_noAO );
  spacelike->scales().QCD_ac      = max(scale,parent->scales().QCD_ac     );
  spacelike->scales().QCD_ac_noAO = max(scale,parent->scales().QCD_ac_noAO);
}

void SudakovFormFactor::doinit() {
  Interfaced::doinit();
  assert(interactionType_!=ShowerInteraction::UNDEFINED);
  assert((colourStructure_>0&&interactionType_==ShowerInteraction::QCD) ||
   	 (colourStructure_<0&&interactionType_==ShowerInteraction::QED) );
  // compute the colour factors if need
  if(colourStructure_==TripletTripletOctet) {
    colourFactor_ = 4./3.;
  }
  else if(colourStructure_==OctetOctetOctet) {
    colourFactor_ = 3.;
  }
  else if(colourStructure_==OctetTripletTriplet) {
    colourFactor_ = 0.5;
  }
  else if(colourStructure_==TripletOctetTriplet) {
    colourFactor_ = 4./3.;
  }
  else if(colourStructure_==SextetSextetOctet) {
    colourFactor_ = 10./3.;
  }
  else if(colourStructure_<0) {
    colourFactor_ = 1.;
  }
  else {
    assert(false);
  }
}

double SudakovFormFactor::colourFactor() const {
  if(colourStructure_>0)
    return colourFactor_;
  else if(colourStructure_<0) {
    if(colourStructure_==ChargedChargedNeutral ||
       colourStructure_==ChargedNeutralCharged) {
      return sqr(double(ids_[0]->iCharge())/3.);
    }
    else if(colourStructure_==NeutralChargedCharged) {
      double fact = sqr(double(ids_[1]->iCharge())/3.);
      if(ids_[1]->coloured())
	fact *= abs(double(ids_[1]->iColour()));
      return fact;
    }
    else {
      assert(false);
      return 0.;
    }
  }
  else {
    assert(false);
    return 0.;
  }
}

bool SudakovFormFactor::checkColours(const IdList & ids) const {
  if(colourStructure_==TripletTripletOctet) {
    if(ids[0]!=ids[1]) return false;
    if((ids[0]->iColour()==PDT::Colour3||ids[0]->iColour()==PDT::Colour3bar) &&
       ids[2]->iColour()==PDT::Colour8) return true;
    return false;
  }
  else if(colourStructure_==OctetOctetOctet) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(ids[ix]->iColour()!=PDT::Colour8) return false;
    }
    return true;
  }
  else if(colourStructure_==OctetTripletTriplet) {
    if(ids[0]->iColour()!=PDT::Colour8) return false;
    if(ids[1]->iColour()==PDT::Colour3&&ids[2]->iColour()==PDT::Colour3bar)
      return true;
    if(ids[1]->iColour()==PDT::Colour3bar&&ids[2]->iColour()==PDT::Colour3)
      return true;
    return false;
  }
  else if(colourStructure_==TripletOctetTriplet) {
    if(ids[0]!=ids[2]) return false;
    if((ids[0]->iColour()==PDT::Colour3||ids[0]->iColour()==PDT::Colour3bar) &&
       ids[1]->iColour()==PDT::Colour8) return true;
    return false;
  }
  else if(colourStructure_==SextetSextetOctet) {
    if(ids[0]!=ids[1]) return false;
    if((ids[0]->iColour()==PDT::Colour6 || ids[0]->iColour()==PDT::Colour6bar) &&
       ids[2]->iColour()==PDT::Colour8) return true;
    return false;
  }
  else if(colourStructure_==ChargedChargedNeutral) {
    if(ids[0]!=ids[1]) return false;
    if(ids[2]->iCharge()!=0) return false;
    if(ids[0]->iCharge()==ids[1]->iCharge()) return true;
    return false;
  }
  else if(colourStructure_==ChargedNeutralCharged) {
    if(ids[0]!=ids[2]) return false;
    if(ids[1]->iCharge()!=0) return false;
    if(ids[0]->iCharge()==ids[2]->iCharge()) return true;
    return false;
  }
  else if(colourStructure_==NeutralChargedCharged) {
    if(ids[1]->id()!=-ids[2]->id()) return false;
    if(ids[0]->iCharge()!=0) return false;
    if(ids[1]->iCharge()==-ids[2]->iCharge()) return true;
    return false;
  }
  else {
    assert(false);
  }
  return false;
}
