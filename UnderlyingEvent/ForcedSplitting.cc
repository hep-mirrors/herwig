#include "ForcedSplitting.h"
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Handlers/DecayHandler.h>

using namespace std;
using namespace ThePEG;
using namespace Herwig;


// This is the routine that is called to start the algorithm. 
void ForcedSplitting::handle(EventHandler &ch, const tPVector &tagged,
			     const Hint &hint) throw(Veto,Stop,Exception) {
  // Find beam particles
  PPair beam = ch.currentCollision()->incoming();
  PVector::const_iterator it;
  tPPtr rem1,rem2;
  tShowerParticlePtr part1,part2;
  for(it = beam.first->children().begin(); it != beam.first->children().end();
      it++) {
    if((*it)->children().size()==0) rem1 = *it;
    else part1 = dynamic_ptr_cast<ShowerParticlePtr>(*it);
  }
  for(it = beam.second->children().begin(); it !=beam.second->children().end();
      it++) {
    if((*it)->children().size()==0) rem2 = *it;
    else part2 = dynamic_ptr_cast<ShowerParticlePtr>(*it);
  }
  
  tStepPtr step = ch.newStep();
  if(rem1) split(rem1,part1,step);
  if(rem2) split(rem2,part2,step);
}

void ForcedSplitting::persistentOutput(PersistentOStream & os) const {
  os << kinCutoff;
}


void ForcedSplitting::persistentInput(PersistentIStream & is, int) {
  is >> kinCutoff;
}


ClassDescription<ForcedSplitting> ForcedSplitting::initForcedSplitting;
// Definition of the static class description member.

void ForcedSplitting::Init() {
  static ClassDocumentation<ForcedSplitting> documentation
    ("This class is responsible for correctly tying the parton shower to "
     "the remaining flavours in the hadron and producing the correct remnant");

  static Parameter<ForcedSplitting,double>
    interfacekinCut("KinCutoff", "Parameter kinCutoff used to constrain qtilde"                    " distribution",
		    &Herwig::ForcedSplitting::kinCutoff,1.,0.,5.,false,false,
		    false);
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
void ForcedSplitting::split(const tPPtr rem, tShowerParticlePtr part, 
			    const tStepPtr step) {
  long hadronId = rem->parents()[0]->id();
  Energy oldQ;
  double oldx;
  long quarks[3];
  int maxIdx = 3;
  int idx = 0;
  long lg = ParticleID::g;
  long currentPart = part->id();
  Lorentz5Momentum usedMomentum = 0;
  Lorentz5Momentum lastp = rem->momentum();
  PPtr lastColour = part;
  PPtr newPart;
  ColinePtr x1, x2;
  x1 = x2 = ColinePtr();
  if(abs(hadronId) > 99) { // We have a hadron
    quarks[0] = hadronId % 10;
    quarks[1] = (hadronId/10)%10;
    quarks[2] = (hadronId/100)%10;
    // NOTE TODO: Must make sure that the sign is correct for the meson quarks
    if(quarks[2] == 0) maxIdx = 2; // we have a meson

    //    oldQ = part->evolutionScales()[ShowerIndex::QCD];
    Energy Q2, mb2, mc2;
    double lambda;
    mb2 = sqr(rem->momentum());
    Q2 = sqr(rem->parents()[0]->momentum());
    mc2 = sqr(part->momentum());
    if(sqr(Q2+mb2-mc2)>4.*mb2*Q2) lambda = sqrt(sqr(Q2+mb2-mc2) - 4.*mb2*Q2);
    else                         lambda = 0.;
    // A qtilde and x value for the remnant
    oldQ = sqrt((Q2+mb2-mc2+lambda)/2.);
    oldx = 1.-part->x();

    // Look first at sea quarks, these must go to a gluon, we then handle
    // the gluons in the next step
    if(currentPart != quarks[0] && currentPart != quarks[1] && 
       currentPart != quarks[2] && currentPart != ParticleID::g) { 
      // Create the new parton with its momentum and parent/child relationship
      // set
      newPart = forceSplit(rem, -currentPart, oldQ, oldx, lastp, 
			   usedMomentum, step);
      // Set the proper colour connections
      x1 = new_ptr(ColourLine());
      if(lastColour->colourLine()) x1->addAntiColoured(newPart);
      else if(lastColour->antiColourLine()) x1->addColoured(newPart);
      currentPart = lg;
    }
    // We now handle the gluons, either it is where the shower terminated or
    // it has been created by splitting a sea quark
    if(currentPart == ParticleID::g) { // gluon
      // Create new particles, splitting is q->g q
      // First choose which q
      idx = UseRandom::irnd(maxIdx);
      Lorentz5Momentum s = lastp-usedMomentum;
      // Generate the next parton, with s momentum remaining in the remnant.
      newPart = forceSplit(rem, quarks[idx], oldQ, oldx, s, usedMomentum,step);
      // Several colour connection cases...
      if(x1) {
	bool npC = newPart->hasColour();
	if(npC) 
	  {
	    x2 = lastColour->antiColourLine();
	    if(!x1->antiColoured().empty())   x1->addColoured(newPart);
	    else if(!x1->coloured().empty())  x2->addColoured(newPart);
	  } 
	else 
	  {
	    x2 = lastColour->colourLine();
	    if(!x1->coloured().empty()) x1->addAntiColoured(newPart);
	    else if(!x1->antiColoured().empty()) x2->addAntiColoured(newPart);
	  }
      } 
      else 
	{
	  x1 = lastColour->colourLine();
	  x2 = lastColour->antiColourLine();
	  if(newPart->hasColour()) x2->addColoured(newPart);
	  else if(newPart->hasAntiColour()) x1->addAntiColoured(newPart);
	  //cerr << "testng new part " << *newPart << endl;
	}
      currentPart = quarks[idx];
    } 
    // Lastly, do the final split into the (di)quark and a parton
    newPart = finalSplit(rem,maxIdx,quarks,currentPart,usedMomentum, step);
    // Set colour connections, case 1, no other forced splittings
    if(!x1 || !x2) 
      {
	if(rem->colourLine()) rem->colourLine()->addColoured(newPart);
	else if(rem->antiColourLine()) 
	  rem->antiColourLine()->addAntiColoured(newPart);
      } 
    else
      {
	if(getParticleData(currentPart)->coloured()) x1->addAntiColoured(newPart);
	else x2->addColoured(newPart);
      }
  }
}

// This creates the parton to split and sets it momentum and parent/child
// relationships
PPtr ForcedSplitting::forceSplit(const tPPtr rem, long child, Energy &oldQ, 
				 double &oldx, Lorentz5Momentum &pf, 
				 Lorentz5Momentum &p,const tStepPtr step) {
  Lorentz5Momentum beam = rem->parents()[0]->momentum();
  PPtr parton = new_ptr(Particle(getParticleData(child)));
  double m2 = sqr(getParticleData(child)->mass());
  Lorentz5Momentum partonp = emit(beam,oldQ,oldx,m2,pf);
  p += partonp;
  parton->set5Momentum(partonp);
  step->addDecayProduct(rem,parton);
  return parton;
}

// This forces the final output of the remnant ((di)quark) and sets the
// momentum and parent/child relationships
PPtr ForcedSplitting::finalSplit(const tPPtr rem, int maxIdx, 
				 long quarks[3], int idx, 
				 Lorentz5Momentum usedMomentum, 
				 const tStepPtr step) {
  tPPtr hadron;
  if(rem->parents().size() == 1) hadron = rem->parents()[0];
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
  PPtr remnant = new_ptr(Particle(getParticleData(remId)));
  remnant->setMomentum(rem->momentum()-usedMomentum);

  // Add the remnant to the step, this will be changed again if the
  // shower is vetoed. Set the colour connections as well
  step->addDecayProduct(rem,remnant);
  return remnant;
}

// This defines the momentum for an emitted parton, currently no pt is
// given to the produced partons, z is generated uniformly.
Lorentz5Momentum ForcedSplitting::emit(const Lorentz5Momentum &par,
				       Energy &lastQ, double &lastx, 
				       double emittedm2,
				       Lorentz5Momentum &pf) {
  // Now generate the new z and qtilde
  Energy q;
  double z,z0,z1;
  kinCutoff = 0.75*GeV;

  //cout << "Emit, kincutoff = " << kinCutoff << ", lastQ = " << lastQ 
    //     << ", emittedm2 = " << emittedm2 << endl;
  //cout << "beam momentum = " << par << endl;
  //cout << "pf = " << pf << endl;
  // Bounds on z
  z0 = 0.;//lastx;
  double yy = 1.+sqr(kinCutoff/lastQ)/2.;
  z1 = yy - sqrt(sqr(yy)-1.); 
  //cout << "Bounds on z " << z0 << ":" << z1 << endl;
  do {
    double randQ = UseRandom::rnd();
    double randZ = UseRandom::rnd();
    z = z0 + (z1-z0)*randZ;
    
    // For the qtilde lets just start with a simple distribution
    // weighted towards the lower value: dP/dQ = 1/Q -> Q(R) =
    // Q0^(1-R) Qmax^R
    q = pow(kinCutoff, 1-randQ)*pow(lastQ,randQ);
    q = kinCutoff/sqrt(z)/(1.-z);
    // check kinematics...
  } while(sqr(z*(1.-z)*q) <= z*sqr(kinCutoff));

  //cout << "qtilde = " << q << " and z = " << z << endl;
  double phi = 2.*pi*UseRandom::rnd();
  //cout << "qtilde = " << q << " and z = " << z << endl;
  Lorentz5Momentum p, n, pcm, pfboost;
  //cout << "qtilde = " << q << " and z = " << z << endl;
  LorentzMomentum newp, qthat;
  pcm = par;

  // Now boost pcm and pf to z only frame
  pfboost = pf;

  p = Lorentz5Momentum(0.0, pcm.vect());
  n = Lorentz5Momentum(0.0, -pcm.vect());
  qthat = LorentzMomentum(cos(phi), sin(phi), 0.0, 0.0);

  double p_dot_n = p.dot(n);
  //cout << "p_dot_n = " << p_dot_n << endl;
  //cout << "p_dot_qt = " << p.dot(qthat) << endl
  //    << "n_dot_qt = " << n.dot(qthat) << endl;

  double lastalpha = pfboost.dot(n)/p_dot_n;
  double lastqt2 = pfboost.dot(qthat);
  lastqt2 = sqr(lastqt2);
  //cout << "lasta = " << lastalpha << ", lastqt2 = " << lastqt2 << endl;

  // Ignore the mass of the quark, if it is here as it must be a valance quark,
  // and assume it is u or d
  // Compute first transverse components using q,z
  double pt2 = sqr(z*(1.-z)*q) - z*sqr(kinCutoff);
  double qt02 = lastqt2;
  double qt2 = z*z*qt02+pt2;
  double q02 = pfboost.dot(pfboost);
  //cout << "q02 = " << q02 << ", pt2 = " << pt2 << ",qt2 = " << qt2 << endl;

  // Now compute alpha and beta
  double alpha = z*lastalpha;
  double q2 = z*q02 - z*emittedm2/(1.-z) - pt2/(1.-z);
  double beta = (qt2 + q2)/(2.*alpha*p_dot_n);

  //cout << "alpha = " << alpha << " beta = " << beta << ", q2 = " << q2 << endl;
  // Store results for iterative calls
  lastx = z*lastx;
  lastalpha = alpha;
  // Compute momentum
  newp = alpha*p + beta*n + qt2*qthat;

  //cout << "newp = " << newp << endl;

  // Now boost newp to the original frame

  // Return the momentum of the emitted parton
  LorentzMomentum k = pf-newp;
  pf = newp;

  //cout << "k = " << k << ", " << k.dot(k) << endl;
  return k;
}

