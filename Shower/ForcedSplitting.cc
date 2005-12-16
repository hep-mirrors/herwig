#include "BackwardEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ShowerParticle.h"
#include "ShowerKinematics.h"
#include "IS_QtildaShowerKinematics1to2.h"

using namespace Herwig;

tEGPtr ForcedSplitting::generator() const { return _backEvolve->generator(); }

Ptr<SplittingGenerator>::pointer ForcedSplitting::splittingGenerator() { 
  return _backEvolve->splittingGenerator(); 
}

void ForcedSplitting::setParticleList(ShowerParticleVector &p) {
  _particles = &p;
}


void ForcedSplitting::persistentOutput(PersistentOStream & os) const {
  os << _backEvolve;
}


void ForcedSplitting::persistentInput(PersistentIStream & is, int) {
  is >> _backEvolve;
}


ClassDescription<ForcedSplitting> ForcedSplitting::initForcedSplitting;
// Definition of the static class description member.

void ForcedSplitting::Init() {
  static ClassDocumentation<ForcedSplitting> documentation
    ("This class is responsible for correctly tying the parton shower to "
     "the remaining flavours in the hadron and producing the correct remnant");

  static Reference<ForcedSplitting,BackwardEvolver> 
    interfaceBackEvolve("BackwardEvolver", 
		      "A reference to the BackwardEvolver object", 
		      &Herwig::ForcedSplitting::_backEvolve,
		      false, false, true, false);
}

/****
 * Now we need to force the final (up to two) splittings so that we are left
 * with only the valence quarks that make up the incoming hadron. If we have
 * terminated the shower on a sea quark, then we need two splittings, one
 * to a gluon and one to a valence quark. If we are on a gluon we just go
 * to a valence quark. We must also choose qtilda and z for each splitting so
 * that the kinematics can be reconstructed properly. This is no longer 
 * sampled according to the splitting functions, as they no longer have space
 * in the virtuality (since the shower has terminated). Instead we use a new
 * distribution.
 * NOTE: temporarily chosen linearly in z and logarithmically in qtilda, this
 * may be changed later.
 ****/
int ForcedSplitting::split(tShowerParticlePtr &part, ShowerParticleVector &p,
			   tEHPtr ch)
{
  setParticleList(p);
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  Forced splittings:" << endl;
  }

  long hadronId = part->parents()[0]->id();
  Energy oldQ;
  long quarks[3];
  int maxIdx = 3;
  int idx = 0;
  int hasEmitted;
  long lg = ParticleID::g;
  if(abs(hadronId) > 99) { // We have a hadron
    quarks[0] = hadronId % 10;
    quarks[1] = (hadronId/10)%10;
    quarks[2] = (hadronId/100)%10;
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log()<< "  Constituents are " << quarks[0] << ", " 
			<< quarks[1] << ", " << quarks[2] << endl;
    }
    // NOTE TODO: Must make sure that the sign is correct for the meson quarks
    if(quarks[2] == 0) maxIdx = 2; // we have a meson
    oldQ = part->evolutionScales()[ShowerIndex::QCD];
 
    // Look first at sea quarks, these must go to a gluon, we then handle
    // the gluons in the next step
    if(part->id() != quarks[0] && part->id() != quarks[1] && 
       part->id() != quarks[2] && part->id() != ParticleID::g) { 
      hasEmitted = forceSplit(part, part->id(), lg, -part->id(), oldQ);
      if(hasEmitted == -1) return hasEmitted;
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "Created gluon splitting, gluon has scale " 
			   << oldQ/GeV << endl;
      }      
    }
    // We now handle the gluons, either it is where the shower terminated or
    // it has been created by splitting a sea quark
    if(part->id() == ParticleID::g) { // gluon
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "  Particle is a gluon\n";
      }      
      // Create new particles, splitting is q->g q
      // First choose which q
      idx = UseRandom::irnd(maxIdx);
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "  Chosen to split into " << quarks[idx] << endl;
      }
      hasEmitted = forceSplit(part,part->id(),quarks[idx],quarks[idx], oldQ);
      if(hasEmitted == -1) return hasEmitted;
    } else {
      // Otherwise figure out which particle we have ended on so we ignore it
      // in the remnant
      for(int i = 0; i<3; i++) if(part->id() == quarks[i]) idx = i;
      _particles->insert(_particles->end(), part);
    }
    makeRemnant(part,maxIdx,quarks,idx,ch);
  }
  return 1;
}

int ForcedSplitting::forceSplit(tShowerParticlePtr &part, long p, long np, 
				long child, Energy &oldQ) {
  int hasEmitted;
  ShowerIndex::InteractionType inter = ShowerIndex::QCD;
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  Particle is a sea quark\n" << flush; 
  }
  // determine some bounds for qtilde
  Energy minQ = splittingGenerator()->showerVariables()->kinScale();
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  delta = " << minQ/GeV 
		       << " and oldQ = " << oldQ/GeV << endl;
  }
  // Create Shower Kinematics
  part->setSplittingFn(splittingGenerator()->getSplittingFunction(p,np));
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  got SplitFun " << part->splitFun() << endl;
  }
  ShoKinPtr kin = forcedSplitting(*part,oldQ,minQ);
  if (kin) {
    part->setShowerKinematics(kin);
    hasEmitted = 1;
  } else {
    hasEmitted = -1; 
    return hasEmitted;
  }
  
  // Create new particles, splitting is g->q qbar
  ShowerParticlePtr newParent = new_ptr(ShowerParticle(getParticleData(np)));
  ShowerParticlePtr childs = new_ptr(ShowerParticle(getParticleData(child)));
  newParent->setFinalState(false);
  childs->setFinalState(true);
  
  // make sure, otherChild is included in TL shower.
  _particles->insert(_particles->end(), childs);
  _particles->insert(_particles->end(), newParent);
  
  _backEvolve->createBranching(part,newParent,childs, kin->qtilde(),inter);
  
  // Store the old data so we can do the gluon splitting
  oldQ = kin->qtilde();
  part = newParent;
  
  return hasEmitted;
}

void ForcedSplitting::makeRemnant(tShowerParticlePtr &part, int maxIdx, 
				  long quarks[3], int idx, tEHPtr ch) {
  tPPtr hadron;
  if(part->parents().size() == 1) hadron = part->parents()[0];
  else if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log()  << "  Created final splitting " 
			<< "no remnant present!" << endl;
  }
  
  // First decide what the remnant is
  long remId;
  int sign, spin;
  if(maxIdx == 2) { // Meson hadronic state
    remId = quarks[(idx+1)%2];
  } else { // Baryonic hadron
    // Get the other 2 elements of the array
    // use modulus to simplify things. idx is the array entry for the
    // parton which eventually leads to the hard process, the other two
    // elements of the array constitute the remnant.
    long id1 = quarks[(idx+1)%3];
    long id2 = quarks[(idx+2)%3];
    
    if (abs(id1) > abs(id2)) swap(id1, id2);
    sign = (id1 < 0) ? -1 : 1; // Needed for the spin 0/1 part
    remId = id2*1000+id1*100;
    
    // Now decide if we have spin 0 diquark or spin 1 diquark
    if(id1 == id2 || UseRandom::rndbool()) spin = 3; // spin 1
    else spin = 1; // otherwise spin 0
    remId += sign*spin;
  }
   
  // Create the remnant and set its momentum, also reset all of the decay 
  // products from the hadron
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  will create remnant with id = " 
		       << remId << ", " << flush 
		       << getParticleData(remId) << flush << endl;
  } 
  ShowerParticlePtr remnant = new_ptr(ShowerParticle(getParticleData(remId)));
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  newRemnant-> = "
		       << remnant << ", hadron-> = " 
		       << hadron << ", p_hadron = " 
		       << hadron->momentum() << endl << flush
		       << "  remId = " << remId << ", " << flush << endl
		       << "  part = " << part << ", part->x() = " 
		       << part->x() << flush << endl; 
  }
  remnant->setMomentum((1-part->x())*hadron->momentum());

  for(int i = hadron->children().size()-1; i!= -1; i--) {
    ShowerParticlePtr child = 
	dynamic_ptr_cast<ShowerParticlePtr>(hadron->children()[i]);
    if(child) {
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower) {
	generator()->log() << "  Hadron = " << hadron 
			   << "(id=" << hadron->id() << ") " 
			   << " has child " << child
			   << ", i = " << i << ", id = " << child->id() 
			   << ", " << (part == child ? " =  part":"") 
			   << endl << flush;
      }
      //if(child->getThePEGBase()) 
      //hadron->abandonChild(child->getThePEGBase());
      hadron->abandonChild(child);
      if(part != child) ch->currentStep()->removeParticle(child);
    } else {
      hadron->abandonChild(hadron->children()[i]);
      ch->currentStep()->removeParticle(hadron->children()[i]);
    }
  }
  // Add the remnant to the step, this will be changed again if the
  // shower is vetoed. Set the colour connections as well
  hadron->addChild(remnant);
  hadron->abandonChild(part);
  hadron->addChild(part);
  if (part->id() < 0) part->colourNeighbour(remnant);      
  else part->antiColourNeighbour(remnant);      
}

ShoKinPtr ForcedSplitting::forcedSplitting(const ShowerParticle &particle,
					   Energy lastQ, Energy minQ) {
  // Now generate the new z and qtilde
  Energy newQ;
  double newZ,z0,z1;
  Energy kinCutoff;
  ShowerVarsPtr vars = splittingGenerator()->showerVariables();;
  if ( particle.id() == ParticleID::g ) {
    kinCutoff = (vars->kinScale() - 0.3*(particle.children()[0])->mass())/2.3;
  } else {
    kinCutoff = (vars->kinScale() - 0.3*particle.data().mass())/2.3;
  }
  // Generate z with the same distributions as for regular splittings
  tSplittingFnPtr sf = particle.splitFun();
  // Bounds on z
  z0 = particle.x();
  double yy = 1.+sqr(kinCutoff/lastQ)/2.;
  z1 = yy - sqrt(sqr(yy)-1.); 
  //  z1 = 1.;
  bool cant;
  do {
    cant = false;
    double randQ = UseRandom::rnd();
    double randZ = UseRandom::rnd();
    if(!sf) {
      if(HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Shower) {
	generator()->log() << "  The particle has no splitting function!"
			   << " Will use flat in z distribution." << endl;
      }
      newZ = z0 + (z1-z0)*randZ;
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "  Particle has no SplitFun()!, " 
			   << z0 << " < " << newZ << " < " << z1 << endl;
      }
    } else {
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "  Trying " << z0 << " < z < " << z1 << " " 
			   << "sf = " << sf << endl;
      }
      newZ = sf->invIntegOverP(sf->integOverP(z0) 
			       + randZ*(sf->integOverP(z1) - 
					sf->integOverP(z0)));
    }
    // For the qtilde lets just start with a simple distribution
    // weighted towards the lower value: dP/dQ = 1/Q -> Q(R) =
    // Q0^(1-R) Qmax^R
    //    newQ = pow(minQ,1-randQ)*pow(lastQ,randQ);
    newQ = pow(kinCutoff, 1-randQ)*pow(lastQ,randQ);
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "  newQ = " << newQ/GeV 
			 << ", newZ = " << newZ << endl;
    }
    // find out whether a next splitting would be possible
    
    if (particle.id() != ParticleID::g) { 
      //      Energy q0g = (vars->kinScale()-0.0015*GeV)/2.3;
      Energy q0g = (vars->kinScale())/2.3;
      double zm;
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	yy = 1.+sqr(q0g/lastQ)/2.;
	zm = yy - sqrt(sqr(yy)-1.); 
	generator()->log() 
	  << "  Checking: x = " << particle.x()
	  << " < zm^2 = " << zm*zm
	  << (particle.x() < zm*zm ? " y":" NO! Can never split twice!") 
	  << endl;
      }
      yy = 1.+sqr(q0g/newQ)/2.; // maximum z if next Q is max as as well 
      zm = yy - sqrt(sqr(yy)-1.); 
      yy = 1.+sqr(q0g/lastQ)/2.; 
      double zm2 = yy - sqrt(sqr(yy)-1.); 
      double xp = particle.x()/newZ;
      //if (xp > zm*zm) {
      //      if (xp > zm) {
      if (xp > zm) {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  yy = 1.+sqr(q0g/lastQ)/2.;
	  zm = yy - sqrt(sqr(yy)-1.); 
	  generator()->log() 
	    << "  Forced: Can't split again! xp = " 
	    << xp << " > zm = " << zm << endl
	    << "  Even with lastQ: zm = " << zm << endl;
	}
	cant = true;
      }
    }
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << (sqr((1.-newZ)*newQ) > newZ*sqr(kinCutoff) ? 
			     "  pt ok":"  pt not ok") 
			 << " for a further branching." << endl;
    }
    // check kinematics...
  } while(sqr((1.-newZ)*newQ) < newZ*sqr(kinCutoff));
  //  } while(sqr((1.-newZ)*newQ) < newZ*sqr(kinCutoff) || cant);
  
  if (cant) return ShoKinPtr();

  Lorentz5Momentum p, n, pthis, ppartner, pcm;
  if(particle.isFromHardSubprocess()) {
//     pthis = particle.momentum();
//     cout << "parent is " << particle.parents()[0]->id() << endl;
//     ppartner = particle.partners()[ShowerIndex::QCD]->momentum();
//     pcm = pthis; 
//     pcm.boost((pthis + ppartner).findBoostToCM());	  
//     p = Lorentz5Momentum(0.0, pcm.vect());
//     n = Lorentz5Momentum(0.0, -pcm.vect()); 
//     p.boost( -(pthis + ppartner).findBoostToCM() );
//     n.boost( -(pthis + ppartner).findBoostToCM() );
    pcm = particle.parents()[0]->momentum();
    p = Lorentz5Momentum(0.0, pcm.vect());
    n = Lorentz5Momentum(0.0, -pcm.vect()); 
  } else {
    p = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
      ->showerKinematics()->getBasis()[0];
    n = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
      ->showerKinematics()->getBasis()[1];
  } 
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  create ShowerKinematics with " 
		       << endl 
		       << "  p = " << p << endl 
		       << "  n = " << n << endl;
  }
  
  Ptr<IS_QtildaShowerKinematics1to2>::pointer showerKin = 
    new_ptr(IS_QtildaShowerKinematics1to2(p, n));

  // Phi is uniform
  showerKin->qtilde(newQ);
  showerKin->setResScale(vars->cutoffQScale(ShowerIndex::QCD));
  showerKin->setKinScale(vars->kinScale()); 
  showerKin->z(newZ);
  showerKin->phi(2.*pi*UseRandom::rnd());

  return showerKin;
}
