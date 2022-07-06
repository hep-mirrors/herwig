// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SudakovFormFactor class.
//

#include "SudakovFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"

using namespace Herwig;


void SudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << alpha_ << pdfMax_ << pdfFactor_ << angularOrdered_
     << oenum(interactionType_) << oenum(colourStructure_) << particles_ ;
}

void SudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> alpha_ >> pdfMax_ >> pdfFactor_ >> angularOrdered_ 
     >> ienum(interactionType_) >> ienum(colourStructure_) >> particles_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<SudakovFormFactor,Interfaced>
describeHerwigSudakovFormFactor("Herwig::SudakovFormFactor", "HwShower.so");

void SudakovFormFactor::Init() {

  static ClassDocumentation<SudakovFormFactor> documentation
    ("There is no documentation for the SudakovFormFactor class");

  static Reference<SudakovFormFactor,ShowerAlpha>
    interfaceAlpha("Alpha",
		   "A reference to the Alpha object",
		   &Herwig::SudakovFormFactor::alpha_,
		   false, false, true, false);

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
  static SwitchOption interfaceInteractionTypeEW
    (interfaceInteractionType,
     "EW","EW",ShowerInteraction::EW);

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
  static SwitchOption interfaceColourStructureTripletTripletSinglet
    (interfaceColourStructure,
     "TripletTripletSinglet",
     "3 -> 3 1",
     TripletTripletSinglet);
  static SwitchOption interfaceColourStructureOctetOctetSinglet
    (interfaceColourStructure,
     "OctetOctetSinglet",
     "8 -> 8 1",
     OctetOctetSinglet);
  static SwitchOption interfaceColourStructureEpsilon
    (interfaceColourStructure,
     "Epsilon",
     "3 -> 3 3",
     Epsilon);
  static SwitchOption interfaceColourStructureOctetSinglet
    (interfaceColourStructure,
     "OctetSinglet",
     "8 -> 8 1",
     OctetSinglet);
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
  static SwitchOption interfaceColourStructureEW
    (interfaceColourStructure,
     "EW",
     "q -> q W/Z, q -> q h0, V -> V' V'', V -> V H",
     EW);

  static Parameter<SudakovFormFactor,double> interfacePdfMax
    ("PDFmax",
     "Maximum value of PDF weight. ",
     &SudakovFormFactor::pdfMax_, 35.0, 1.0, 1000000.0,
     false, false, Interface::limited);

  static Switch<SudakovFormFactor,unsigned int> interfacePDFFactor
    ("PDFFactor",
     "Include additional factors in the overestimate for the PDFs",
     &SudakovFormFactor::pdfFactor_, 0, false, false);
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
  else if(colourStructure_==TripletTripletSinglet) {
    if(ids[0]->iColour()!=ids[1]->iColour()) return false;
    if((ids[0]->iColour()==PDT::Colour3||ids[0]->iColour()==PDT::Colour3bar) &&
       ids[2]->iColour()==PDT::Colour0) return true;
    return false;
  }
  else if(colourStructure_==OctetOctetSinglet) {
    if(ids[0]!=ids[1]) return false;
    if(ids[0]->iColour()==PDT::Colour8 && ids[2]->iColour()==PDT::Colour0) return true;
    return false;
  }
  else if(colourStructure_==Epsilon) {
    if(ids[0]->iColour()!=PDT::Colour3&&ids[0]->iColour()!=PDT::Colour3bar) return false;
    if(ids[0]->iColour()!=-ids[1]->iColour()) return false;
    if(ids[0]->iColour()!=-ids[2]->iColour()) return false;
    return true;
  }
  else if(colourStructure_==OctetSinglet) {
    if(ids[0]->iColour()==PDT::Colour8 && ids[1]->iColour()==PDT::Colour0) return true;
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
  else if(colourStructure_==EW) {
    return true;
  }
  else {
    assert(false);
  }
  return false;
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

bool SudakovFormFactor::PDFVeto(const Energy2 t, const double x, const double z,
				const tcPDPtr parton0, const tcPDPtr parton1,
				Ptr<BeamParticleData>::transient_const_pointer beam) const {
  double ratio=PDFVetoRatio(t,x,z,parton0,parton1,beam,1.);
  return UseRandom::rnd() > ratio;
}

double SudakovFormFactor::PDFVetoRatio(const Energy2 t, const double x, const double z,
					   const tcPDPtr parton0, const tcPDPtr parton1,
					   Ptr<BeamParticleData>::transient_const_pointer beam,double factor) const {
  assert(pdf_);
  Energy2 theScale = t * sqr(ShowerHandler::currentHandler()->factorizationScaleFactor()*factor);
  if (theScale < sqr(freeze_)) theScale = sqr(freeze_);

  const double newpdf=pdf_->xfx(beam,parton0,theScale,x/z);
  if(newpdf<=0.) return 0.;

  const double oldpdf=pdf_->xfx(beam,parton1,theScale,x);
  if(oldpdf<=0.) return 1.;
  
  const double ratio = newpdf/oldpdf;
  double maxpdf = pdfMax_;

  switch (pdfFactor_) {
  case 0: break;
  case 1: maxpdf /= z; break;
  case 2: maxpdf /= 1.-z; break;
  case 3: maxpdf /= (z*(1.-z)); break;
  case 4: maxpdf /= sqrt(z); break;
  case 5: maxpdf *= sqrt(z); break;
  default :
    throw Exception() << "SudakovFormFactor::PDFVetoRatio invalid PDFfactor = "
		      << pdfFactor_ << Exception::runerror;
    
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

bool SudakovFormFactor::alphaSVeto(Energy2 pt2) const {
  double ratio=alphaSVetoRatio(pt2,1.);
  return UseRandom::rnd() > ratio;
}

double SudakovFormFactor::alphaSVetoRatio(Energy2 pt2, double factor) const {
  if(ShowerHandler::currentHandlerIsSet())
    factor *= ShowerHandler::currentHandler()->renormalizationScaleFactor();
  return alpha_->ratio(pt2, factor);
}
