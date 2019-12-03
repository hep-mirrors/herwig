// -*- C++ -*-
//
// ConstituentReshuffler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ConstituentReshuffler class.
//

#include <config.h>
#include "ConstituentReshuffler.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include <limits>

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "DipolePartonSplitter.h"

#include "Herwig/Utilities/GSLBisection.h"

#include "Herwig/Shower/Dipole/DipoleShowerHandler.h"

#include "Herwig/Shower/ShowerHandler.h"

using namespace Herwig;

ConstituentReshuffler::ConstituentReshuffler() 
  : HandlerBase() {}

ConstituentReshuffler::~ConstituentReshuffler() {}

IBPtr ConstituentReshuffler::clone() const {
  return new_ptr(*this);
}

IBPtr ConstituentReshuffler::fullclone() const {
  return new_ptr(*this);
}

double ConstituentReshuffler::ReshuffleEquation::aUnit() {
  return 1.;
}

double ConstituentReshuffler::ReshuffleEquation::vUnit() {
  return 1.;
}

double ConstituentReshuffler::DecayReshuffleEquation::aUnit() {
  return 1.;
}

double ConstituentReshuffler::DecayReshuffleEquation::vUnit() {
  return 1.;
}


double ConstituentReshuffler::ReshuffleEquation::operator() (double xi) const {

  double r = - w/GeV;

  for (PList::iterator p = p_begin; p != p_end; ++p) {
    r += sqrt(sqr((**p).dataPtr()->constituentMass()) +
	      xi*xi*(sqr((**p).momentum().t())-sqr((**p).dataPtr()->mass()))) / GeV;
  }
  return r;  

}

double ConstituentReshuffler::DecayReshuffleEquation::operator() (double xi) const {
  double r = - w/GeV;

  for (PList::iterator pIt = p_begin; pIt != p_end; ++pIt) {
    r += sqrt(sqr((**pIt).dataPtr()->constituentMass()) +
 	      xi*xi*(sqr((**pIt).momentum().t())-sqr((**pIt).dataPtr()->mass()))) / GeV;
  }

  for (PList::iterator rIt = r_begin; rIt != r_end; ++rIt) {
    r +=  sqrt(sqr((**rIt).momentum().m()) +
	       xi*xi*(sqr((**rIt).momentum().t())-sqr((**rIt).momentum().m()))) / GeV;
  }

  return r;  

}

void ConstituentReshuffler::reshuffle(PList& out,
				      PPair& in,
				      PList& intermediates,
				      const bool decay,
				      PList& decayPartons,
				      PList& decayRecoilers) {

  assert(ShowerHandler::currentHandler()->retConstituentMasses());

  if ( !decay ) {
  
    if (out.size() == 0)
      return;

    if (out.size() == 1) {

      PPtr recoiler;
      PPtr parton = out.front();

      if (DipolePartonSplitter::colourConnected(parton,in.first) &&
	  DipolePartonSplitter::colourConnected(parton,in.second)) {
	if (UseRandom::rnd() < .5)
	  recoiler = in.first;
	else
	  recoiler = in.second;
      } else if (DipolePartonSplitter::colourConnected(parton,in.first)) {
	recoiler = in.first;
      } else if (DipolePartonSplitter::colourConnected(parton,in.second)) {
	recoiler = in.second;
      } else assert(false);

      assert(abs(recoiler->momentum().vect().perp2()/GeV2) < 1e-6);

      double sign = recoiler->momentum().z() < 0.*GeV ? -1. : 1.;

      Energy2 qperp2 = parton->momentum().perp2();

      if (qperp2/GeV2 < Constants::epsilon) {
	// no emission off a 2 -> singlet process which
	// needed a single forced splitting: should never happen (?)
	assert(false);
	throw Veto();
      }

      Energy2 m2 = sqr(parton->dataPtr()->constituentMass());
      
      Energy abs_q = parton->momentum().vect().mag();
      Energy qz = parton->momentum().z();
      Energy abs_pz = recoiler->momentum().t();
      assert(abs_pz > 0.*GeV);

      Energy xi_pz = sign*(2.*qperp2*abs_pz + m2*(abs_q + sign*qz))/(2.*qperp2);
      Energy x_qz = (2.*qperp2*qz + m2*(qz+sign*abs_q))/(2.*qperp2);

      Lorentz5Momentum recoiler_momentum 
	(0.*GeV,0.*GeV,xi_pz,xi_pz < 0.*GeV ? - xi_pz : xi_pz);

      recoiler_momentum.rescaleMass();

      Lorentz5Momentum parton_momentum 
	(parton->momentum().x(),parton->momentum().y(),x_qz,sqrt(m2+qperp2+x_qz*x_qz));

      parton_momentum.rescaleMass();

      PPtr n_parton = new_ptr(Particle(parton->dataPtr()));
      n_parton->set5Momentum(parton_momentum);

      DipolePartonSplitter::change(parton,n_parton,false);

      out.pop_front();
      intermediates.push_back(parton);
      out.push_back(n_parton);

      PPtr n_recoiler = new_ptr(Particle(recoiler->dataPtr()));
      n_recoiler->set5Momentum(recoiler_momentum);

      DipolePartonSplitter::change(recoiler,n_recoiler,true);

      intermediates.push_back(recoiler);

      if (recoiler == in.first) {
	in.first = n_recoiler;
      }

      if (recoiler == in.second) {
	in.second = n_recoiler;
      }

      return;

    }

  }

  Energy zero (0.*GeV);
  Lorentz5Momentum Q (zero,zero,zero,zero);
    
  for (PList::iterator p = out.begin();
       p != out.end(); ++p) {
    Q += (**p).momentum();
  }

  Boost beta = Q.findBoostToCM();

  list<Lorentz5Momentum> mbackup;

  bool need_boost = (beta.mag2() > Constants::epsilon);

  if (need_boost) {

    for (PList::iterator p = out.begin();
	 p != out.end(); ++p) {
      Lorentz5Momentum mom = (**p).momentum();
      mbackup.push_back(mom);
      (**p).set5Momentum(mom.boost(beta));
    }

  }

  double xi;

  // Only partons
  if ( decayRecoilers.size()==0 ) {
    ReshuffleEquation solve (Q.m(),out.begin(),out.end());

    GSLBisection solver(1e-10,1e-8,10000);

    try {
      xi = solver.value(solve,0.0,1.1);
    } catch (GSLBisection::GSLerror) {
      throw DipoleShowerHandler::RedoShower();
    } catch (GSLBisection::IntervalError) {
      throw DipoleShowerHandler::RedoShower();
    }
  }

  // Partons and decaying recoilers
  else {
    DecayReshuffleEquation solve (Q.m(),decayPartons.begin(),decayPartons.end(),decayRecoilers.begin(),decayRecoilers.end());

    GSLBisection solver(1e-10,1e-8,10000);

    try {
      xi = solver.value(solve,0.0,1.1);
    } catch (GSLBisection::GSLerror) {
      throw DipoleShowerHandler::RedoShower();
    } catch (GSLBisection::IntervalError) {
      throw DipoleShowerHandler::RedoShower();
    }
  }


  PList reshuffled;

  list<Lorentz5Momentum>::const_iterator backup_it;
  if (need_boost)
    backup_it = mbackup.begin();
    

  // Reshuffling of non-decaying partons only
  if ( decayRecoilers.size()==0 ) {

    for (PList::iterator p = out.begin();
	 p != out.end(); ++p) {

      PPtr rp = new_ptr(Particle((**p).dataPtr()));

      DipolePartonSplitter::change(*p,rp,false);

      Lorentz5Momentum rm;

      rm = Lorentz5Momentum (xi*(**p).momentum().x(),
			     xi*(**p).momentum().y(),
			     xi*(**p).momentum().z(),
			     sqrt(sqr((**p).dataPtr()->constituentMass()) +
				  xi*xi*(sqr((**p).momentum().t())-sqr((**p).dataPtr()->mass()))));
     
      rm.rescaleMass();

      if (need_boost) {
	(**p).set5Momentum(*backup_it);
	++backup_it;
	rm.boost(-beta);
      }

      rp->set5Momentum(rm);

      intermediates.push_back(*p);
      reshuffled.push_back(rp);
      

    }
  }

  // For the case of a decay process with non-partonic recoilers
  else {
    assert ( decay );

    for (PList::iterator p = out.begin();
	 p != out.end(); ++p) {

      // Flag to update spinInfo
      bool updateSpin = false;

      PPtr rp = new_ptr(Particle((**p).dataPtr()));

      DipolePartonSplitter::change(*p,rp,false);

      Lorentz5Momentum rm;

      // If the particle is a parton and not a recoiler
      if ( find( decayRecoilers.begin(), decayRecoilers.end(), *p ) == decayRecoilers.end() ) {
	rm = Lorentz5Momentum (xi*(**p).momentum().x(),
			       xi*(**p).momentum().y(),
			       xi*(**p).momentum().z(),
			       sqrt(sqr((**p).dataPtr()->constituentMass()) +
				    xi*xi*(sqr((**p).momentum().t())-sqr((**p).dataPtr()->mass()))));
      }


      // Otherwise the parton is a recoiler 
      // and its invariant mass must be preserved
      else {
        if ( (*p)-> spinInfo() )
          updateSpin = true;
	rm = Lorentz5Momentum (xi*(**p).momentum().x(),
			       xi*(**p).momentum().y(),
			       xi*(**p).momentum().z(),
			       sqrt(sqr((**p).momentum().m()) +
				    xi*xi*(sqr((**p).momentum().t())-sqr((**p).momentum().m()))));
      }
 
      rm.rescaleMass();

      if (need_boost) {
	(**p).set5Momentum(*backup_it);
	++backup_it;
	rm.boost(-beta);
      }

      rp->set5Momentum(rm);

      // Update SpinInfo if required
      if ( updateSpin )
        updateSpinInfo(*p, rp);
      
      intermediates.push_back(*p);
      reshuffled.push_back(rp);
      

    }
  }
 
  out.clear();
  out.splice(out.end(),reshuffled);

}


void ConstituentReshuffler::hardProcDecayReshuffle(PList& decaying,
						   PList& eventOutgoing,
						   PList& eventHard,
						   PPair& eventIncoming,
						   PList& eventIntermediates) {
  
  // Note, when this function is called, the particle pointers
  // in theDecays/decaying are those prior to the showering.
  // Here we find the newest pointers in the outgoing.
  // The update of the PPtrs in theDecays is done in DipoleShowerHandler::constituentReshuffle()
  // as this needs to be done if ConstituentReshuffling is switched off.

  
  //Make sure the shower should return constituent masses:
  assert(ShowerHandler::currentHandler()->retConstituentMasses());

  // Find the outgoing decaying particles
  PList recoilers;
  for ( PList::iterator decIt = decaying.begin(); decIt != decaying.end(); ++decIt) {
    
    // First find the particles in the intermediates
    PList::iterator pos = find(eventIntermediates.begin(),eventIntermediates.end(), *decIt);

    // Colourless particle or coloured particle that did not radiate.
    if(pos==eventIntermediates.end()) {
      
      // Check that this is not a particle from a subsequent decay.
      // e.g. the W from a top decay from an LHE file.
      if ( find( eventHard.begin(), eventHard.end(), *decIt ) == eventHard.end() &&
	   find( eventOutgoing.begin(), eventOutgoing.end(), *decIt ) == eventOutgoing.end() )
      continue;
    
      else
	recoilers.push_back( *decIt );

    }
    
    // Coloured decaying particle that radiated
    else {
      PPtr unstable = *pos;
      while(!unstable->children().empty()) {
	unstable = unstable->children()[0];
      }
      assert( find( eventOutgoing.begin(),eventOutgoing.end(), unstable ) != eventOutgoing.end() );

      recoilers.push_back( unstable );
    }
  }

  // Make a list of partons
  PList partons;
  for ( PList::iterator outPos = eventOutgoing.begin(); outPos != eventOutgoing.end(); ++outPos ) {
    if ( find (recoilers.begin(), recoilers.end(), *outPos ) == recoilers.end() ) {
      partons.push_back( *outPos );
    }
  }


  // If no outgoing partons, do nothing
  if ( partons.size() == 0 ){
    return;
  }
  
  // Otherwise reshuffling needs to be done.

  // If there is only one parton, attempt to reshuffle with 
  // the incoming to be consistent with the reshuffle for a
  // hard process with no decays.
  else if ( partons.size() == 1 &&
            ( DipolePartonSplitter::colourConnected(partons.front(),eventIncoming.first) || 
				     DipolePartonSplitter::colourConnected(partons.front(),eventIncoming.second) ) ) {
      
    // Erase the parton from the event outgoing
    eventOutgoing.erase( find( eventOutgoing.begin(), eventOutgoing.end(), partons.front() ) );
    // Perform the reshuffle, this update the intermediates and the incoming
    reshuffle(partons, eventIncoming, eventIntermediates);
    // Update the outgoing
    eventOutgoing.push_back(partons.front());
    return;
  }
    
  // If reshuffling amongst the incoming is not possible
  // or if we have multiple outgoing partons.
  else {

    // Create a complete list of the outgoing from the process
    PList out;
    // Make an empty list for storing the new intermediates
    PList intermediates;
    // Empty incoming particles pair
    PPair in;

    // A single parton which cannot be reshuffled 
    // with the incoming.
    if ( partons.size() == 1 ) {

      // Populate the out for the reshuffling
      out.insert(out.end(),partons.begin(),partons.end());
      out.insert(out.end(),recoilers.begin(),recoilers.end());
      assert( out.size() > 1 );
    
      // Perform the reshuffle with the temporary particle lists
      reshuffle(out, in, intermediates, true, partons, recoilers);    
    }
  
    // If there is more than one parton, reshuffle only
    // amongst the partons
    else {
      assert(partons.size() > 1);

      // Populate the out for the reshuffling
      out.insert(out.end(),partons.begin(),partons.end());
      assert( out.size() > 1 );
    
      // Perform the reshuffle with the temporary particle lists
      reshuffle(out, in, intermediates, true);
    }
  
    // Update the dipole event record
    updateEvent(intermediates, eventIntermediates, out, eventOutgoing, eventHard );
    return;
  }
}
    

void ConstituentReshuffler::decayReshuffle(PerturbativeProcessPtr& decayProc,
					   PList& eventOutgoing,
					   PList& eventHard,
					   PList& eventIntermediates ) {
  
  // Separate particles into those to be assigned constituent masses
  // i.e. non-decaying coloured partons
  // and those which must only absorb recoil
  // i.e. non-coloured and decaying particles
  PList partons;
  PList recoilers;

  //Make sure the shower should return constituent masses:
  assert(ShowerHandler::currentHandler()->retConstituentMasses());
  

  // Populate the particle lists from the outgoing of the decay process
  for( unsigned int ix = 0; ix<decayProc->outgoing().size(); ++ix) {

    // Identify recoilers
    if ( !decayProc->outgoing()[ix].first->coloured() || 
	 ShowerHandler::currentHandler()->decaysInShower(decayProc->outgoing()[ix].first->id() ) )
      recoilers.push_back(decayProc->outgoing()[ix].first);

    else 
      partons.push_back(decayProc->outgoing()[ix].first);
  }

  // If there are no outgoing partons, then no reshuffling
  // needs to be done
  if ( partons.size() == 0 )
    return;

  // Reshuffling needs to be done:
  else {

    // Create a complete list of the outgoing from the process
    PList out;
    // Make an empty list for storing the new intermediates
    PList intermediates;
    // Empty incoming particles pair
    PPair in;


    // SW - 15/06/2018, 31/01/2019 - Always include 'recoilers' in
    // reshuffling, regardless of the number of partons to be put on their
    // constituent mass shell. This is because reshuffling between 2 partons
    // frequently leads to a redoShower exception. This treatment is
    // consistent with the AO shower

      // Populate the out for the reshuffling
      out.insert(out.end(),partons.begin(),partons.end());
      out.insert(out.end(),recoilers.begin(),recoilers.end());
      assert( out.size() > 1 );
    
      // Perform the reshuffle with the temporary particle lists
      reshuffle(out, in, intermediates, true, partons, recoilers);    

    // Update the dipole event record and the decay process
    updateEvent(intermediates, eventIntermediates, out, eventOutgoing, eventHard, decayProc );
    return;
  }
}


void ConstituentReshuffler::updateEvent( PList& intermediates,
					 PList& eventIntermediates,
#ifndef NDEBUG
					 PList& out,
#else
					 PList&,
#endif
					 PList& eventOutgoing,
					 PList& eventHard,
					 PerturbativeProcessPtr decayProc ) {

  // Loop over the new intermediates following the reshuffling
  for (PList::iterator p = intermediates.begin();
       p != intermediates.end(); ++p) {

    // Update the event record intermediates
    eventIntermediates.push_back(*p);

    // Identify the reshuffled particle
    assert( (*p)->children().size()==1 );
    PPtr reshuffled = (*p)->children()[0];
    assert( find(out.begin(), out.end(), reshuffled) != out.end() );

    // Update the event record outgoing
    PList::iterator posOut = find(eventOutgoing.begin(), eventOutgoing.end(), *p);

    if ( posOut != eventOutgoing.end() ) {
      eventOutgoing.erase(posOut);
      eventOutgoing.push_back(reshuffled);
    }
      
    else {
      PList::iterator posHard = find(eventHard.begin(), eventHard.end(), *p);
      assert( posHard != eventHard.end() );
      eventHard.erase(posHard);
      eventHard.push_back(reshuffled);
    }
    
    // Replace the particle in the the decay process outgoing
    if ( decayProc ) {
      vector<pair<PPtr,PerturbativeProcessPtr> >::iterator decayOutIt = decayProc->outgoing().end();
      for ( decayOutIt = decayProc->outgoing().begin();
	    decayOutIt!= decayProc->outgoing().end(); ++decayOutIt ) {
	
	if ( decayOutIt->first == *p ){
	  break;
	}
      }
      assert( decayOutIt != decayProc->outgoing().end() );
      decayOutIt->first = reshuffled;
    }
  }  
}
 
void ConstituentReshuffler::updateSpinInfo( PPtr& oldPart,
                                            PPtr& newPart ) {

  const Lorentz5Momentum& oldMom = oldPart->momentum();
  const Lorentz5Momentum& newMom = newPart->momentum();

  // Rotation from old momentum to +ve z-axis
  LorentzRotation oldToZAxis;
  Axis axisOld(oldMom.vect().unit());
  if( axisOld.perp2() > 1e-12 ) {
    double sinth(sqrt(1.-sqr(axisOld.z())));
    oldToZAxis.rotate( -acos(axisOld.z()),Axis(-axisOld.y()/sinth,axisOld.x()/sinth,0.));
  }

  // Rotation from new momentum to +ve z-axis
  LorentzRotation newToZAxis;
  Axis axisNew(newMom.vect().unit());
  if( axisNew.perp2() > 1e-12 ) {
    double sinth(sqrt(1.-sqr(axisNew.z())));
    newToZAxis.rotate( -acos(axisNew.z()),Axis(-axisNew.y()/sinth,axisNew.x()/sinth,0.));
  }

  // Boost from old momentum to new momentum along z-axis
  Lorentz5Momentum momOldRotated = oldToZAxis*Lorentz5Momentum(oldMom);
  Lorentz5Momentum momNewRotated = newToZAxis*Lorentz5Momentum(newMom);
  
  Energy2 a = sqr(momOldRotated.z()) + sqr(momNewRotated.t());
  Energy2 b = 2.*momOldRotated.t()*momOldRotated.z();
  Energy2 c = sqr(momOldRotated.t()) - sqr(momNewRotated.t());
  double beta;
  
  // The rotated momentum should always lie along the +ve z-axis
  if ( momOldRotated.z() > ZERO )
    beta = (-b + sqrt(sqr(b)-4.*a*c)) / 2. / a;
  else
    beta = (-b - sqrt(sqr(b)-4.*a*c)) / 2. / a;
  
  LorentzRotation boostOldToNew(0., 0., beta);

  // Total transform
  LorentzRotation transform = (newToZAxis.inverse())*boostOldToNew*oldToZAxis;

  // Assign the same spin info to the old and new particles
  newPart->spinInfo(oldPart->spinInfo());
  newPart->spinInfo()->transform(oldMom, transform);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ConstituentReshuffler::persistentOutput(PersistentOStream &) const {
}

void ConstituentReshuffler::persistentInput(PersistentIStream &, int) {
}

ClassDescription<ConstituentReshuffler> ConstituentReshuffler::initConstituentReshuffler;
// Definition of the static class description member.

void ConstituentReshuffler::Init() {

  static ClassDocumentation<ConstituentReshuffler> documentation
    ("The ConstituentReshuffler class implements reshuffling "
     "of partons on their nominal mass shell to their constituent "
     "mass shells.");

}

