// -*- C++ -*-
//
// PartonSplitter.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonSplitter class.
//

#include "PartonSplitter.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/EventRecord/Step.h>
#include "ThePEG/Interface/Parameter.h"
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Utilities/DescribeClass.h>
#include "ClusterHadronizationHandler.h"
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/PDT/PDT.h>
#include "CheckId.h"


using namespace Herwig;

IBPtr PartonSplitter::clone() const {
  return new_ptr(*this);
}

IBPtr PartonSplitter::fullclone() const {
  return new_ptr(*this);
}

void PartonSplitter::persistentOutput(PersistentOStream & os) const {
  os << _quarkSelector << ounit(_gluonDistance,femtometer)
     << _splitGluon << _splitPwtUquark << _splitPwtDquark << _splitPwtSquark
     << _enhanceSProb << ounit(_m0,GeV) << _massMeasure;
}

void PartonSplitter::persistentInput(PersistentIStream & is, int) {
  is >> _quarkSelector >> iunit(_gluonDistance,femtometer)
     >> _splitGluon >> _splitPwtUquark >> _splitPwtDquark >> _splitPwtSquark
     >>_enhanceSProb >> iunit(_m0,GeV) >> _massMeasure;
}

DescribeClass<PartonSplitter,Interfaced>
describePartonSplitter("Herwig::PartonSplitter","");


void PartonSplitter::Init() {

  static ClassDocumentation<PartonSplitter> documentation
    ("This class is reponsible of the nonperturbative splitting of partons");

  static Switch<PartonSplitter,int> interfaceSplit
    ("Split",
     "Option for different splitting options",
     &PartonSplitter::_splitGluon, 0, false, false);
  static SwitchOption interfaceSplitDefault
    (interfaceSplit,
     "ud",
     "Normal cluster splitting where only u and d  quarks are drawn is used.",
     0);
  static SwitchOption interfaceSplitAll
    (interfaceSplit,
     "uds",
     "Alternative cluster splitting where all light quark pairs (u, d, s) can be drawn.",
     1);

  static Parameter<PartonSplitter,double> interfaceSplitPwtUquark
    ("SplitPwtUquark",
     "Weight for splitting in U quarks",
     &PartonSplitter::_splitPwtUquark, 1, 0.0, 1.0,
     false, false, Interface::limited);
  static Parameter<PartonSplitter,double> interfaceSplitPwtDquark
    ("SplitPwtDquark",
     "Weight for splitting in D quarks",
     &PartonSplitter::_splitPwtDquark, 1, 0.0, 1.0,
     false, false, Interface::limited);
  static Parameter<PartonSplitter,double> interfaceSplitPwtSquark
    ("SplitPwtSquark",
     "Weight for splitting in S quarks",
     &PartonSplitter::_splitPwtSquark, 0.5, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<PartonSplitter,int> interfaceEnhanceSProb
    ("EnhanceSProb",
     "Option to enhance the strangeness weight (MassSplit Switch needs to be on)",
     &PartonSplitter::_enhanceSProb,0,false,false);
  static SwitchOption interfaceEnhanceSProbNo
    (interfaceEnhanceSProb,
     "No",
     "No scaling for strangeness",
     0);
  static SwitchOption interfaceEnhanceSProbScale
    (interfaceEnhanceSProb,
     "Scaled",
     "Power-law scaling for strangeness",
     1);
  static SwitchOption interfaceEnhanceSProbExp
    (interfaceEnhanceSProb,
     "Exponential",
     "Exponential suppression for strangeness",
     2);

   static Switch<PartonSplitter,int> interfaceMassMeasure
     ("MassMeasure",
      "Option to use different mass measures",
      &PartonSplitter::_massMeasure,0,false,false);
   static SwitchOption interfaceMassMeasureMass
     (interfaceMassMeasure,
      "Mass",
      "Mass Measure",
      0);
   static SwitchOption interfaceMassMeasureLambda
     (interfaceMassMeasure,
      "Lambda",
      "Lambda Measure",
      1);

  static Parameter<PartonSplitter,Energy> interfaceMassScale
   ("MassScale",
    "Mass scale for g->qqb strangeness enhancement",
    &PartonSplitter::_m0, GeV, 20.*GeV, 1.*GeV, 1000.*GeV,
    false, false, Interface::limited);
}

void PartonSplitter::split(PVector & tagged) {
  // set the gluon c tau once and for all
  static bool first = true;

  if(first) {
    _gluonDistance = hbarc*getParticleData(ParticleID::g)->constituentMass()/
      ClusterHadronizationHandler::currentHandler()->minVirtuality2();
    first = false;
  }
  // Copy incoming for the (possible sorting and) splitting
  PVector particles = tagged;
  // Switch to enhance strangeness
  if (_enhanceSProb >= 1){
    colourSorted(particles);
  }

  PVector newtag;
  Energy2 Q02 = 0.99*sqr(getParticleData(ParticleID::g)->constituentMass());
  // Loop over all of the particles in the event.
  // Loop counter for colourSorted
  for(PVector::iterator pit = particles.begin(); pit!=particles.end(); ++pit) {
    // only considering gluons so add other particles to list of particles
    if( (**pit).data().id() != ParticleID::g ) {
      newtag.push_back(*pit);
      continue;
    }
    // should not have been called for massless or space-like gluons
    if((**pit).momentum().m2() <= 0.0*sqr(MeV) ) {
      throw Exception()
	<< "Spacelike or massless gluon m2= " << (**pit).momentum().m2()/GeV2
	<< "GeV2 in PartonSplitter::split()"
	<< Exception::eventerror;
    }
    // time like gluon gets split
    PPtr ptrQ = PPtr();
    PPtr ptrQbar = PPtr();
    if (_enhanceSProb == 0){
      splitTimeLikeGluon(*pit,ptrQ,ptrQbar);
    }
    else {
      size_t i = pit - particles.begin();
      massSplitTimeLikeGluon(*pit, ptrQ, ptrQbar, i);
    }

    ptrQ->scale(Q02);
    ptrQbar->scale(Q02);

    (*pit)->colourLine()->addColoured(ptrQ);
    (*pit)->addChild(ptrQ);
    newtag.push_back(ptrQ);

    (*pit)->antiColourLine()->addAntiColoured(ptrQbar);
    (*pit)->addChild(ptrQbar);
    newtag.push_back(ptrQbar);

    // set the life length of gluon
    Length distance = UseRandom::rndExp(_gluonDistance);
    (**pit).setLifeLength((distance/(**pit).mass())*(**pit).momentum());
    // assume quarks same position as gluon
    ptrQ   ->setVertex((**pit).decayVertex());
    ptrQ   ->setLifeLength(Lorentz5Distance());
    ptrQbar->setVertex((**pit).decayVertex());
    ptrQbar->setLifeLength(Lorentz5Distance());
  }
  swap(tagged,newtag);
}

void PartonSplitter::splitTimeLikeGluon(tcPPtr ptrGluon,
					PPtr & ptrQ,
					PPtr & ptrQbar){
  // select the quark flavour
  tPDPtr quark;
  long idNew=0;
  switch(_splitGluon){
    case 0:
      quark = _quarkSelector.select(UseRandom::rnd());
      break;
    case 1:
      if ( ptrGluon->momentum().m() <
	   2.0 *getParticle(ThePEG::ParticleID::s)->data().constituentMass() ) {
	throw Exception() << "Impossible Kinematics in PartonSplitter::splitTimeLikeGluon()"
			  << Exception::runerror;
      }
      // Only allow light quarks u,d,s with the probabilities
      double prob_d = _splitPwtDquark;
      double prob_u = _splitPwtUquark;
      double prob_s = _splitPwtSquark;

      int choice = UseRandom::rnd3(prob_u, prob_d, prob_s);
      switch(choice) {
        case 0: idNew = ThePEG::ParticleID::u; break;
        case 1: idNew = ThePEG::ParticleID::d; break;
        case 2: idNew = ThePEG::ParticleID::s; break;
      }
      ptrQ = getParticle(idNew);
      ptrQbar = getParticle(-idNew);
  break;
  }
  // Solve the kinematics of the two body decay  G --> Q + Qbar
  Lorentz5Momentum momentumQ;
  Lorentz5Momentum momentumQbar;
  double cosThetaStar = UseRandom::rnd( -1.0 , 1.0 );
  using Constants::pi;
  double phiStar = UseRandom::rnd( -pi , pi );

  Energy constituentQmass;
  if(_splitGluon==0) {
    constituentQmass = quark->constituentMass();
  }else{
    constituentQmass = ptrQ->data().constituentMass();
  }

 if (ptrGluon->momentum().m() < 2.0*constituentQmass) {
    throw Exception() << "Impossible Kinematics in PartonSplitter::splitTimeLikeGluon()"
		      << Exception::eventerror;
  }
  Kinematics::twoBodyDecay(ptrGluon->momentum(), constituentQmass,
			   constituentQmass, cosThetaStar, phiStar, momentumQ,
			   momentumQbar );
  // Create quark and anti-quark particles of the chosen flavour
  // and set they 5-momentum (the mass is the constituent one).
  if(_splitGluon==0) {
    ptrQ    = new_ptr(Particle(quark      ));
    ptrQbar = new_ptr(Particle(quark->CC()));
  }

  ptrQ    ->set5Momentum( momentumQ    );
  ptrQbar ->set5Momentum( momentumQbar );
}

void PartonSplitter::doinit() {
  Interfaced::doinit();
  // calculate the probabilties for the gluon to branch into each quark type
  // based on the available phase-space, as in fortran.
  Energy mg=getParticleData(ParticleID::g)->constituentMass();
  for( int ix=1; ix<6; ++ix ) {
    PDPtr quark = getParticleData(ix);
    Energy pcm = Kinematics::pstarTwoBodyDecay(mg,quark->constituentMass(),
					       quark->constituentMass());
    if(pcm>ZERO) _quarkSelector.insert(pcm/GeV,quark);
  }
  if(_quarkSelector.empty())
    throw InitException() << "At least one quark must have constituent mass less "
			  << "then the constituent mass of the gluon in "
			  << "PartonSplitter::doinit()" << Exception::runerror;
}

// Method to colour sort the event and calculate the masses of the
// pre-clusters
// Convention is to have the triplet of the colour singlet first,
// then all gluons, then the antitriplet (and repeat for all the
// colour-singlets in the event)
void PartonSplitter::colourSorted(PVector& tagged) {
  // Set up the output
  PVector sorted;
  // Reset the storage variables for doing the mass-based strangeness
  // enhancement
  _colSingletSize.resize(0);
  _colSingletm2.resize(0);
  // Variable to exit out of gluon loops
  bool gluonLoop = false;
  // Partons left to consider
  PVector left = tagged;
  // Loop over all the tagged particles
  while (int(left.size()) > 0){
    // Pick the first particle available
    PPtr p = left[0];
    // Temporary holding pen for momenta
    Lorentz5Momentum tempMom(ZERO, ZERO, ZERO, ZERO);
    // Don't do anything if the particle has no colour
    // Simply add it back into the list of particles
    if ( !p->coloured() ) {
      sorted.push_back(p);
      // nparts is the index of the particle after adding it to the sorted list
      int nparts = 1;
      // Add on the last entry of colSingletSize if the vector is not empty
      // This is essentially the index the particle will have once it
      // Get placed into the new colour sorted event
      nparts += (_colSingletSize.empty()) ? 0 : _colSingletSize.back();
      tempMom += p->momentum();
      Energy2 singletm2 = tempMom.m2();
      // Store the number of particles and mass of the colour-singlet
      _colSingletSize.push_back(nparts);
      _colSingletm2.push_back(singletm2);
      // Remove the particle from the list of particles left
      left.erase(remove(left.begin(),left.end(), p), left.end());
    }
    // Temporary holding pen for partons
    PVector temp;
    // Variable to sum end-point masses i.e. triplets and anti-triplets
    Energy endPointMass = ZERO;
    // If the particle in question is a gluon, search up it's antiColourLine
    // until we get to the triplet.
    // Note there are situations where we have a gluon loop
    if ( p->hasColour() && p->hasAntiColour() ){
      // Search up its anticolour line until we get to the start particle
      tPPtr partner = p->antiColourLine()->endParticle();
      // Dummy index used to loop over the anticolour line trace
      tPPtr dummy = partner;
      // Store the partner particle
      temp.push_back(partner);
      tempMom += partner->momentum();
      left.erase(remove(left.begin(),left.end(), partner), left.end());
      // While loop continues while we still reach a particle with with
      // anti-colour, i.e. a gluon
      while ( dummy->hasAntiColour() ){
        dummy = dummy->antiColourLine()->endParticle();
        // Check that we haven't already included it via colour indices
        // If we have, it is a gluon loop.
        if ( find(left.begin(), left.end(), dummy) == left.end() ) {
          gluonLoop = true;
          break;
        }
        // Store the dummy partons in reverse
        temp.push_back(dummy);
        tempMom += dummy->momentum();
        // Remove counted ones from the remaining list
        left.erase(remove(left.begin(),left.end(), dummy), left.end());
      }
      // Number of particles in this colour singlets so far
      int nparts = int(temp.size());
      // Insert the new particles in the reverse order
      sorted.insert(sorted.end(), temp.rbegin(), temp.rend());
      endPointMass += ((temp.back())->mass());
      // If it is a gluon loop we've already looped over the colour-singlet
      // in its entirety, so we need to end early
      if (gluonLoop){
        // Store the index of the final entry
        nparts += (_colSingletSize.empty()) ? 0 : _colSingletSize.back();
        // Insert the new particles in the correct order
        // i.e. triplet first, then the gluons we have seen so far
        Energy2 singletm2 = tempMom.m2();
        _colSingletSize.push_back(nparts);
        _colSingletm2.push_back(singletm2);
        continue;
      }
      // If it is not a gluon loop, we now need to trace the colourLine
      // down, until we reach the triplet.
      // Works similarly to the first half
      // Reset the temp PVector
      temp.resize(0);
      // Push back the particle in question
      temp.push_back(p);
      // tempMom hasn't been reset, add the particle we started with
      tempMom += p->momentum();
      left.erase(remove(left.begin(),left.end(), p), left.end());
      // Search down its colour line until we get to the end particle
      dummy = p->colourLine()->startParticle();
      temp.push_back(dummy);
      tempMom += dummy->momentum();
      left.erase(remove(left.begin(),left.end(), dummy), left.end());
      while ( dummy->hasColour() ){
        dummy = dummy->colourLine()->startParticle();
        temp.push_back(dummy);
        tempMom += dummy->momentum();
        left.erase(remove(left.begin(),left.end(), dummy), left.end());
      }
      endPointMass += ((temp.back())->mass());
      // Update size of colour singlet
      nparts += int(temp.size());
      nparts += (_colSingletSize.empty()) ? 0 : _colSingletSize.back();
      // Insert the new particles in the correct order
      Energy2 singletm2 = tempMom.m2();
      sorted.insert(sorted.end(), temp.begin(), temp.end());
      endPointMass += ((sorted.back())->mass());
      _colSingletSize.push_back(nparts);
      // Chooses whether to use invariant mass of singlet
      // or to use the lambda measure i.e. m^2 - (\sum m_i)^2
      Energy2 m2  = (_massMeasure == 0) ? singletm2 : singletm2 - sqr(endPointMass);
      _colSingletm2.push_back(m2);
    }
    // Else if it's a quark
    else if ( p->hasColour() ) {
      // Search up its colour line until we get to the start particle
      // Works the same way as the second half of the gluon handling
      tPPtr partner = p->colourLine()->startParticle();
      tPPtr dummy = partner;
      temp.push_back(p);
      endPointMass += ((temp.back())->mass());
      temp.push_back(partner);
      tempMom += p->momentum();
      tempMom += partner->momentum();
      left.erase(remove(left.begin(),left.end(), p), left.end());
      left.erase(remove(left.begin(),left.end(), partner), left.end());
      while ( dummy->hasColour() ){
        dummy = dummy->colourLine()->startParticle();
        temp.push_back(dummy);
        tempMom += dummy->momentum();
        left.erase(remove(left.begin(),left.end(), dummy), left.end());
      }
      // Number of particles in this colour singlets
      int nparts = int(temp.size());
      nparts += (_colSingletSize.empty()) ? 0 : _colSingletSize.back();
      // Insert the new particles in the correct order
      Energy2 singletm2 = tempMom.m2();
      sorted.insert(sorted.end(), temp.begin(), temp.end());
      endPointMass += ((sorted.back())->mass());
      _colSingletSize.push_back(nparts);
      Energy2 m2  = (_massMeasure == 0) ? singletm2 : singletm2 - sqr(endPointMass);
      _colSingletm2.push_back(m2);
    }
    // Else it's an antiquark
    else if ( p->hasAntiColour() ) {
      // Search along anti-colour line, storing particles, and reversing the order
      // at the end
      // Works in the same way as the first half of the gluon handling
      tPPtr partner = p->antiColourLine()->endParticle();
      tPPtr dummy = partner;
      temp.push_back(p);
      endPointMass += ((temp.back())->mass());
      temp.push_back(partner);
      tempMom += p->momentum();
      tempMom += partner->momentum();
      left.erase(remove(left.begin(),left.end(), p), left.end());
      left.erase(remove(left.begin(),left.end(), partner), left.end());
      while ( dummy->hasAntiColour() ){
        dummy = dummy->antiColourLine()->endParticle();
        temp.push_back(dummy);
        tempMom += dummy->momentum();
        left.erase(remove(left.begin(),left.end(), dummy), left.end());
      }
      // Number of particles in this colour singlets
      int nparts = int(temp.size());
      nparts += (_colSingletSize.empty()) ? 0 : _colSingletSize.back();
      // Insert the particles in the reverse order
      Energy2 singletm2 = tempMom.m2();
      sorted.insert(sorted.end(), temp.rbegin(), temp.rend());
      endPointMass += ((temp.back())->mass());
      _colSingletSize.push_back(nparts);
      Energy2 m2  = (_massMeasure == 0) ? singletm2 : singletm2 - sqr(endPointMass);
      _colSingletm2.push_back(m2);
    }
  }

  // Check that the algorithm hasn't missed any particles.
  assert( sorted.size() == tagged.size() );

  swap(sorted,tagged);

}

void PartonSplitter::massSplitTimeLikeGluon(tcPPtr ptrGluon,
					PPtr & ptrQ,
					PPtr & ptrQbar, size_t i){
  // select the quark flavour
  tPDPtr quark;
  long idNew=0;

  if ( ptrGluon->momentum().m() <
 2.0 *getParticle(ThePEG::ParticleID::s)->data().constituentMass() ) {
throw Exception() << "Impossible Kinematics in PartonSplitter::massSplitTimeLikeGluon()"
    << Exception::runerror;
  }
  // Only allow light quarks u,d,s with the probabilities
  double prob_d = _splitPwtDquark;
  double prob_u = _splitPwtUquark;
  double prob_s = enhanceStrange(i);
  int choice = UseRandom::rnd3(prob_u, prob_d, prob_s);
  switch(choice) {
    case 0: idNew = ThePEG::ParticleID::u; break;
    case 1: idNew = ThePEG::ParticleID::d; break;
    case 2: idNew = ThePEG::ParticleID::s; break;
  }
  ptrQ = getParticle(idNew);
  ptrQbar = getParticle(-idNew);

  // Solve the kinematics of the two body decay  G --> Q + Qbar
  Lorentz5Momentum momentumQ;
  Lorentz5Momentum momentumQbar;
  double cosThetaStar = UseRandom::rnd( -1.0 , 1.0 );
  using Constants::pi;
  double phiStar = UseRandom::rnd( -pi , pi );

  Energy constituentQmass;
  constituentQmass = ptrQ->data().constituentMass();

 if (ptrGluon->momentum().m() < 2.0*constituentQmass) {
    throw Exception() << "Impossible Kinematics in PartonSplitter::massSplitTimeLikeGluon()"
          << Exception::eventerror;
  }
  Kinematics::twoBodyDecay(ptrGluon->momentum(), constituentQmass,
         constituentQmass, cosThetaStar, phiStar, momentumQ,
         momentumQbar );
  // Create quark and anti-quark particles of the chosen flavour
  // and set they 5-momentum (the mass is the constituent one).
  ptrQ    ->set5Momentum( momentumQ    );
  ptrQbar ->set5Momentum( momentumQbar );

}

double PartonSplitter::enhanceStrange(size_t i){

  // Get the m2 of the relevant colour singlet
  // First we need to get the index of the colour-singlet the indexed i
  // parton has
  auto const it = lower_bound(_colSingletSize.begin(), _colSingletSize.end(), i);

  // Get the index of the colourSinglet mass
  int indx = distance(_colSingletSize.begin(), it);

  Energy2 mass2 = _colSingletm2[indx];
  Energy2 m2 = _m0*_m0;

  // Scaling strangeness enhancement
  if (_enhanceSProb == 1){
    // If the mass is significantly smaller than the characteristic mass,
    // just set the prob to 0
    double scale = double(m2/mass2);
    return (_maxScale < scale || scale < 0.) ? 0. : pow(_splitPwtSquark,scale);
  }
  // Exponential strangeness enhancement
  else if (_enhanceSProb == 2){
    double scale = double(m2/mass2);
    return (_maxScale < scale || scale < 0.) ? 0. : exp(-scale);
  }
  else
    return _splitPwtSquark;
}
