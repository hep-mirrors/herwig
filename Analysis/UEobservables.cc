// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEobservables class.
//

#include "UEobservables.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "Herwig++/Interfaces/KtJetInterface.h"
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

UEobservables::UEobservables() {}

UEobservables::UEobservables(const UEobservables & x)
  : AnalysisHandler(x), theShowerHandler(x.theShowerHandler) {}

UEobservables::~UEobservables() {}

void UEobservables::dofinish() {
  AnalysisHandler::dofinish();
  theRootData.finish();
}

void UEobservables::doinitrun() {
  AnalysisHandler::doinitrun();
  string gen(generator()->filename());
  gen = gen.substr(2, gen.size());
  //string pt(gen.substr(gen.rfind("-")+1, gen.size()));
  theRootData.init((gen+".root").c_str(),gen.c_str());
}

void UEobservables::analyze(tEventPtr event, long , int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  /** get the final-state particles */
  tPVector particles=event->getFinalState();
  tPVector perp_particles;
  int nch(0), nch_t(0);
  double dphi(0), eta(0);
  Lorentz5Momentum p;

  for (tPVector::const_iterator pit = particles.begin(); pit != particles.end(); ++pit){
      /** Select only the charged particles  */
      if( ChargedSelector::Check(**pit) ){
	  p = (**pit).momentum();
	  /** write a TClonesArray of TLorentzVector's to the file */
//	  cout << "debug: " << p.px() << " " << p.pz() << " " << p.e() << endl;
	  theRootData.fillarray("all", nch, p.x()/GeV, p.y()/GeV, p.z()/GeV, p.e()/GeV );
	  /**
	     Select the particles according the noted selection cuts and do the
	     analysis (towards away and transverse region) only on this set of particles.
	  */
//	  cerr << p.perp();
	  eta = (**pit).eta();
	  if ( fabs( eta ) < 1 && p.perp() > 0.5*GeV ){
	      perp_particles.push_back( *pit );
	  }

	  nch++;
      }


  }
  
  /*
    get the "transverse" region, defined by the azimuthal angle difference
    to the reconstructed jet with larges pt and write all the particles that 
    are in there to the ROOT file in the same manner as above.
  */
  if(perp_particles.size()){

      KtJetInterface jet;
      KtJet::KtEvent ev(jet.convert(perp_particles), 4, 2, 2, 0.7);//pt recom scheme
      /**
	 cerr << "Number of final state jets in event " << ieve << " is " << ev.getNJets() << endl;
      */

      /** Get the leading jet (largest pt) */
      vector<KtJet::KtLorentzVector> jets = ev.getJetsPt();
      vector<KtJet::KtLorentzVector>::const_iterator itr = jets.begin();
  
      /** and write its pt, eta, phi to the file */
      theRootData.fill("pt1", (*itr).perp()/1000.);//is in MeV convert to GeV
      theRootData.fill("eta1", (*itr).eta());
      theRootData.fill("phi1", (*itr).phi());


	  
      for (tPVector::const_iterator pit = perp_particles.begin(); pit != perp_particles.end(); ++pit){

	  p = (**pit).momentum();

	  dphi = fabs( p.phi() - (*itr).phi() )*180.0/Constants::pi;
	  if (dphi > 180) dphi = fabs( dphi - 360.0 );
      
	  if( dphi > 60 && dphi < 120 ){
	      theRootData.fillarray("perp", nch_t, p.x()/GeV, p.y()/GeV, p.z()/GeV, p.e()/GeV );
	      nch_t++;
	  }
      }

/**
      cerr << endl << "event " << ieve << " enthielt " << nch << " geladene Teilchen und "
	   << nch_t << " davon im transversalen Bereich, pt:" << (*itr).perp()/GeV << endl << endl;
*/

      theRootData.write();

  }else{
//      cerr << "no particle survived the selection cuts for the UEAnalysis!" << endl;
  }
}

LorentzRotation UEobservables::transform(tEventPtr ) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void UEobservables::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void UEobservables::analyze(tPPtr) {}

void UEobservables::persistentOutput(PersistentOStream & os) const {
  os << theShowerHandler;
}

void UEobservables::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> theShowerHandler;
}

ClassDescription<UEobservables> UEobservables::initUEobservables;
// Definition of the static class description member.

void UEobservables::Init() {

  static ClassDocumentation<UEobservables> documentation
    ("There is no documentation for the UEobservables class");

  
  static Reference<UEobservables,ShowerHandler> interfaceShowerHandler
    ("ShowerHandler",
     "A reference to the ShowerHandler",
     &UEobservables::theShowerHandler, true, false, true, false, false);

}

