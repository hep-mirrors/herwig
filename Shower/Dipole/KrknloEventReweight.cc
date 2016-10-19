#include "KrknloEventReweight.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;
using namespace ThePEG;


IBPtr KrknloEventReweight::clone() const {
  return new_ptr(*this);
}

IBPtr KrknloEventReweight::fullclone() const {
  return new_ptr(*this);
}

void KrknloEventReweight::persistentOutput(PersistentOStream & os) const {
  os << alphaSScale_R_ << alphaSScale_V_ << scaleFactor_ << mcPDF_ << mode_;
}

void KrknloEventReweight::persistentInput(PersistentIStream & is, int) {
  is >> alphaSScale_R_ >> alphaSScale_V_ >> scaleFactor_ >> mcPDF_ >> mode_;
}


bool KrknloEventReweight::secondaryInteractions() const {
  return false;
}

bool KrknloEventReweight::firstInteraction() const {
  return true;

}

double KrknloEventReweight::weightCascade(const PPair& in, const PList& , const PList& , Ptr<AlphaSBase>::tptr as) const { 

  // This weight will be applied only once per event even if there was no radiation
  // we use it to is to implemnt VS correction 
  //===============================================================================
  // Warning: it won't be apply if we cut the emission after 1 emission!

  using Constants::pi;
  const double CF=4.0/3.0;
  
  double w = 1.0;

  const Lorentz5Momentum & qa = in.first->momentum();
  const Lorentz5Momentum & qb = in.second->momentum();  
  
  Energy2 q2;
  if ( alphaSScale_V_ == 0 ) { 
  	q2 = (qa+qb)*(qa+qb); /* better way? */ 
  }
  else { 
    const PID theType = ( mode_ == 0 ) ? ParticleID::Z0 : ParticleID::h0;
    q2 = sqr(getParticleData(theType)->mass());
  }
  
  // AlphaS
  const double alpha_s = as->value(q2) / (2*pi);
  w *= alpha_s;

  if ( mode_ == 0 ) { //Z boson
    if ( mcPDF_ == 1 ) { // fullPDF
      w *= CF*(4.0/3.0*sqr(pi) + 1.0/2.0);
    }
    else if ( mcPDF_ == 0 ){ //DY PDF
      w *= CF*(4.0/3.0*sqr(pi) - 5.0/2.0);
    }
  }
  else if ( mode_ == 1 ){ //Higgs
    const double CA = 3.0;
    const double Tf = 0.5 * 5.0; // Tf=NfTR
    w *= CA * ( 4.0/3.0*sqr(pi) + 473.0/36.0 - 59.0/18.0*Tf/CA );
  }
  else { 
    w=0.; 
  }
  
  return 1.+w; // VS weight
}

double KrknloEventReweight::weight(const PPair& in, const PList& out, const PList& , Ptr<AlphaSBase>::tptr as) const {

  double w = 1.0;
 
  using Constants::pi;
    
  // Weight first emission
  if ( out.size() ==  1 ) {
    Lorentz5Momentum qa = in.first->momentum();
    Lorentz5Momentum qb = in.second->momentum();  
    const Lorentz5Momentum & q = out.front()->momentum();
    
    if(out.front()->parents()[0]->id() == in.second->id() ) 
    	swap(qa, qb); // qA should be emitter

    const double alpha = (qb*q)/(qa*qb);
    const double beta  = (qa*q)/(qa*qb);
    const double x = 1. - alpha - beta;
  
    if ( mode_ == 0 ) { // Z Boson
      //=============== real weight for qq ===================================
      if ( out.front()->id() == 21 ) {
        w *=  1. - 2.*alpha*beta/(1. + sqr(x));
      }
      //=============== real weight for qg =================================== 
      else if( abs(out.front()->id()) <= 6 ) { 
	      w *= 1. + beta*(beta+2.*x)/(sqr(x) + sqr(1. - x));
      }
      //=============== TEST we should have no other emissions weighted=======
      else {
        throw Exception() << "Error: unsupported particle with ID "
               << out.front()->id()
               << " in KrknloEventReweight::weight()"
               << Exception::runerror;
      }
    }
    else if ( mode_ == 1 ) { // Higgs 
	    //=============== real weight for qq =================================== 
      if (out.front()->id() == 21) {
	       w *= (1.0 + pow(x,4.0) + pow(alpha,4.0) + pow(beta,4.0))/(1.0 + pow(x,4) + pow(1-x,4) );
      }
      //=============== real weight for qg ===================================                                    
      else if( abs(out.front()->id()) <= 6 ) {
        w *= (1.0 + sqr(alpha))/(1.0 + sqr(1.0-x));
      }
      // TEST we should have no other emissions weighted                              
      else {
        throw Exception() << "Error: unsupported particle with ID "
               << out.front()->id()
               << " in KrknloEventReweight::weight()"
               << Exception::runerror;
      }

      }
    else {
      throw Exception() << "Error: unsupported run mode " 
             << mode_
             << "in KrknloEventReweight::weight()"
             << Exception::runerror; // This should not be possible via interface
    }
    
    // AlphaS reweight
    const Energy2 M2 
      = sqr(getParticleData((mode_ == 0) ? 
            ParticleID::Z0 : ParticleID::h0)->mass());
    const Energy2 pt2 = alpha*beta*2*qa*qb; //scale used in the PS evolution                                       
    
    switch( alphaSScale_R_ ) {
        case 0  :
          w *= 1.0;
          break;
        case 1  :
          w *= as->value(M2) / as->value(pt2);
          break;
        case 2:
          if (pt2 > M2){ 
            w *= as->value(M2) / as->value(pt2); 
          }
          break;
    }

  }

  return w;
}

DescribeClass<KrknloEventReweight, DipoleEventReweight>
  describeKrknloEventReweight("Herwig::KrknloEventReweight", 
                              "HwKrknloEventReweight.so");



void KrknloEventReweight::Init() {
  
  static ClassDocumentation<KrknloEventReweight> documentation
    ("There is no documentation for the KrknloEventReweight class");

    static Switch<KrknloEventReweight,unsigned int> interfaceAlphaS_V
      ("AlphaS_V",
       "The argument used to evaluate alphaS for virtual correction",
       &KrknloEventReweight::alphaSScale_V_, 1, false, false);
       
    static SwitchOption interfaceAlphaS_V_Q2
     (interfaceAlphaS_V, "Q2", "as(Q2)", 0);
    static SwitchOption interfaceAlphaS_V_M2
      (interfaceAlphaS_V, "M2", "as(M2)", 1);


    static Switch<KrknloEventReweight,unsigned int> interfaceAlphaS_R
     ("AlphaS_R",
      "The argument used to evaluate alphaS for real emission",
      &KrknloEventReweight::alphaSScale_R_, 0, false, false);

    static SwitchOption interfaceAlphaS_R_Q2
      (interfaceAlphaS_R, "Q2", "as(Q2)", 0);
    static SwitchOption interfaceAlphaS_R_M2
     (interfaceAlphaS_R, "M2", "as(M2)", 1);
    static SwitchOption interfaceAlphaS_R_M2_Freeze
      (interfaceAlphaS_R, "Freeze", "as(Q2 < M2 ? Q2 : MZ)", 2);

    static Switch<KrknloEventReweight, unsigned int> interfaceMode
      ("Mode",
       "Mode of the method: Higgs or Z boson for now",
       &KrknloEventReweight::mode_, 0, false, false);

    static SwitchOption interfaceModeZ
      (interfaceMode,
       "Z", "Zboson", 0);
    static SwitchOption interfaceModeH
      (interfaceMode, "H", "Higgs", 1);

    static Parameter<KrknloEventReweight, double> interfaceAlphaSScaleFactor
     ("ScaleFactor",
      "The scale factor used to locally rescale the argument of the running coupling [OPTION CURRENTLY DISABLED]",
        &KrknloEventReweight::scaleFactor_, 1.0, 0.0, 0,
        false, false, Interface::lowerlim);

    static Switch<KrknloEventReweight,unsigned int> interfacePDF
      ("PDF",
       "The full MC PDF or just quark MC PDF",
       &KrknloEventReweight::mcPDF_, 0, false, false);

    static SwitchOption interfacePDFcomplete
      (interfacePDF,
       "MC",
       "Full MC Scheme PDFs",
       1);
    static SwitchOption interfacePDFquarkonly
      (interfacePDF,
       "MCDY",
       "MC-DY Scheme PDFs",
       0);

}
