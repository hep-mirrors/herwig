// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Remnant class.
//

#include "Remnant.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/VSelector.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Remnant.tcc"
#endif


using namespace Herwig;
using namespace ThePEG;

Remnant::~Remnant() {}

ClassDescription<Remnant> Remnant::initRemnant;
// Definition of the static class description member.

void Remnant::Init() {

  static ClassDocumentation<Remnant> documentation
    ("There is no documentation for the Remnant class");

}

Remnant::Remnant(tcEventPDPtr x) :Particle(x)
{}

Remnant::Remnant(PartonBinInstance & pb,const LorentzMomentum & p) : 
Particle(CurrentGenerator::current().getParticleData(ExtraParticleID::Remnant))
{ 
  // set the momentum of the remnant
  set5Momentum(p);
  // id of the particle
  _parent=pb.particleData();
  int pid(_parent->id());
  // get the valence flavours
  _sign = pid < 0 ? -1 : 1;
  // beam particle is a baryon
  if(BaryonMatcher::Check(*(pb.particleData())))
    {
      // get the valence flavours
      _valence.resize(3);
      _valence[0]=(abs(pid)/1000)%10;
      _valence[1]=(abs(pid)/100)%10;
      _valence[2]=(abs(pid)/10)%10;
    }
  // beam particle is a meson
  else if(MesonMatcher::Check(*(pb.particleData())))
    {throw Exception() << "Meson requested in Remant::Remnant() but not implemented "
		       << Exception::runerror;}
  // unknown type of beam particle
  else
    {throw Exception() << " requested in Remant::Remnant() but not implemented "
			<< Exception::runerror;}
  // work out the flavours of the remnants
  obtainConstituents(pb.partonData()->id());
}

PPtr Remnant::clone() const {
  return dynamic_ptr_cast<PPtr>(ptr_new<RemnantPtr>(*this));
}

PPtr Remnant::fullclone() const {
  return clone();
}

void Remnant::obtainConstituents(int extracted)
{
  // set the code of the extracted particle
  _extracted=extracted;
  // construct the remnant
  _constituents.resize(0);
  // copy of the valence partons to construct the remnant
  vector<int> vtemp(_valence);
  // see if the parton is one of the valence ones
  vector<int>::iterator v=find(vtemp.begin(),vtemp.end(),_sign*_extracted);
  // if it is
  bool isvalence(false);
  if(v!=vtemp.end())
    {
      vtemp.erase(v);
      isvalence=true;
    }
  // if valence then the remnant is a diquark
  if(isvalence)
    {
      // this is the spin 1 diquark
      int idqr = 1000*max(vtemp[0],vtemp[1])+100*min(vtemp[0],vtemp[1])+3;
      // if flavours the same could be spin-0 (makes no difference in Hw++)
      if(vtemp[0]!=vtemp[1] && UseRandom::rnd() < 0.25) idqr-=2;
      _constituents.push_back(CurrentGenerator::current().getParticleData(_sign*idqr)->produceParticle());
    }
  // otherwise all constituents
  else
    {
      // obtain the possible quarks and diquarks for the valence bit
      VSelector< pair< int, int > > valenceselector;
      int iq1,iq2,iq3;
      for(iq1 = 0; iq1 < 3; iq1++)
	{
	  iq2 = (iq1+1)%3;
	  iq3 = (iq2+1)%3;
	  // This is the id of the diquark (spin 1) that accompanies iq1
	  int idq = 1000*max(vtemp[iq2], vtemp[iq3]) +
	    100*min(vtemp[iq2], vtemp[iq3]) + 3;
	  valenceselector.insert(3.0, make_pair(vtemp[iq1], idq));
	  if(vtemp[iq2] == vtemp[iq3]) continue;
	  // If they are different, we have spin 0 combination too
	  valenceselector.insert(1.0, make_pair(vtemp[iq1], idq-2));
	}
      // select a quark-diquark pair and add to remnant
      pair<int,int> rr = valenceselector.select(UseRandom::current());
      _constituents.push_back(CurrentGenerator::current().getParticleData(rr.first *_sign)->produceParticle());
      _constituents.push_back(CurrentGenerator::current().getParticleData(rr.second*_sign)->produceParticle());
      // if we extracted a sea quark/antiquark then we to add the antiparticle
      // as well
      if(_extracted!=ParticleID::g)
	{_constituents.push_back(CurrentGenerator::current().getParticleData(-_extracted)->produceParticle());}
    }
}

void Remnant::regenerate(tPPtr extracted,Lorentz5Momentum ptotal)
{
  // change the momentum
  set5Momentum(ptotal);
  // change the constituents
  obtainConstituents(extracted->id());
  // remake the colour connections
  // remove old colour connection
  if(this->colourLine())
    this->colourLine()->removeColoured(this,false);
  if(this->antiColourLine())
    this->antiColourLine()->removeColoured(this,true);
  // make the new colour lines
  if(extracted->colourLine())
    extracted->colourLine()->addColoured(this,true);
  if(extracted->antiColourLine())
    extracted->antiColourLine()->addColoured(this,false);
}

void Remnant::createRemnant(tStepPtr pstep)
{
  cout << "testing in create remnant\n";
  for(unsigned int ix=0;ix<_valence.size();++ix)
    {cout << "testing valence " << ix << " " << _valence[ix] << '\n';}
  for(unsigned int ix=0;ix<_constituents.size();++ix)
    {cout << "testing consitituents " << _constituents[ix]->PDGName() << '\n';}
  // if only one constituent just add it
  if(_constituents.size()==1)
    {
      _constituents[0]->set5Momentum(momentum());
      // set up the colours
      if(this->colourLine()    )
	this->colourLine()->addColoured(_constituents[0],true);
      if(this->antiColourLine())
	this->antiColourLine()->addColoured(_constituents[0],true);
      this->addChild(_constituents[0]);
      pstep->addDecayProduct(_constituents[0]);
      return;
    }
  // if two constituents
  else if(_constituents.size()==2)
    {
      cerr << "Remnant::createRemnant() testing two constituents " << _extracted << '\n';
      exit(1);
    }
  else if(_constituents.size()==3)
    {
      cerr << "Remnant::createRemnant() testing three constituents " << _extracted << '\n';
      // first we need a forced splitting of the 

      exit(1);
    }
  else
    {
      cerr << "Remnant::createRemnant() testing #constituents != 1,2 or 3 " << _extracted << '\n';
      exit(1);
    }
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


//     oldQ = part->evolutionScales()[ShowerIndex::QCD];

//     // Look first at sea quarks, these must go to a gluon, we then handle
//     // the gluons in the next step
//     if(part->id() != quarks[0] && part->id() != quarks[1] && 
//        part->id() != quarks[2] && part->id() != ParticleID::g) { 


//       // determine some bounds for qtilde
//       minQ = _splittingGenerator->showerVariables()->kinScale();
//       // Create Shower Kinematics
//       part->setSplittingFn(_splittingGenerator->
// 			   getSplittingFunction(part->id(),
// 						ParticleID::g));
//       ShoKinPtr kinematics = forcedSplitting(*part,oldQ,minQ);
//       if (kinematics) {
// 	part->setShowerKinematics(kinematics);
// 	hasEmitted = 1;
//       } else {
// 	hasEmitted = -1; 
// 	return hasEmitted;
//       }

//       // Create new particles, splitting is g->q qbar
//       ShowerParticlePtr newParent = new_ptr(
// 		   ShowerParticle(getParticleData(ParticleID::g)));
//       ShowerParticlePtr otherChild = new_ptr(
// 		   ShowerParticle(getParticleData(-part->id())));
//       newParent->setFinalState(false);
//       otherChild->setFinalState(true);
//       //newParent->setFromHardSubprocess(true);

//       // make sure, otherChild is included in TL shower.
//       allShowerParticles.insert(allShowerParticles.end(), otherChild);
//       allShowerParticles.insert(allShowerParticles.end(), newParent);

//       createBranching(part,newParent,otherChild,kinematics->qtilde(),inter);

//       // Store the old data so we can do the gluon splitting
//       oldQ = kinematics->qtilde();
//       part = newParent;

//       // Put into list so it will be final showered
//       //allShowerParticles.push_back(otherChild);
//       //allShowerParticles.push_back(newParent);
//     }






















//     // We now handle the gluons, either it is where the shower terminated or
//     // it has been created by splitting a sea quark
//     int idx = 0;
//         // determine some bounds for qtilde
//       minQ = _splittingGenerator->showerVariables()->kinScale();

//       // Create new particles, splitting is q->g q
//       // First choose which q
//       idx = UseRandom::irnd(maxIdx);

//       part->setSplittingFn(_splittingGenerator->
// 			   getSplittingFunction(ParticleID::g, quarks[idx]));


//       // Create Shower Kinematics
//       ShoKinPtr kinematics = forcedSplitting(*part,oldQ,minQ);
//       if (kinematics) {
// 	part->setShowerKinematics(kinematics);
// 	hasEmitted = 1;
//       } else {
// 	hasEmitted = -1; 
// 	return hasEmitted;
//       }

//       ShowerParticlePtr newParent = new_ptr(
// 		   ShowerParticle(getParticleData(quarks[idx])));
//       ShowerParticlePtr otherChild = new_ptr(
// 		   ShowerParticle(getParticleData(quarks[idx])));
//       newParent->setFinalState(false);
//       otherChild->setFinalState(true);
//       //newParent->setFromHardSubprocess(true);

//       // Set the colour and parent/child relationships
//       createBranching(part,newParent,otherChild,kinematics->qtilde(),inter);

//       // Add these so that they will be treated properly later
//       allShowerParticles.push_back(otherChild);
//       allShowerParticles.push_back(newParent);
//       //newParent->setFromHardSubprocess(true);
//       part = newParent;
//     } else {
//       // Otherwise figure out which particle we have ended on so we ignore it
//       // in the remnant
//       for(int i = 0; i<3; i++) if(part->id() == quarks[i]) idx = i;
//       allShowerParticles.push_back(part);
//       //part->setFromHardSubprocess(true);
//     }
//     // set remnant. 
//     tPPtr hadron;
// //     cout << part << '\n' 
// // 	 << part->parents().size()  << '\n';
//     if(part->parents().size() == 1) hadron = part->parents()[0];

//     // First decide what the remnant is
//     long remId;
//     int sign, spin;
//     if(maxIdx == 2) { // Meson hadronic state
//       remId = quarks[(idx+1)%2];
//     } else { // Baryonic hadron
//       // Get the other 2 elements of the array
//       long id1 = quarks[(idx+1)%2];
//       long id2 = quarks[(idx+2)%2];
//       // hack... ask Phil about that assignment!
//       if (abs(id1) > abs(id2)) swap(id1, id2);
//       sign = (id1 < 0) ? -1 : 1; // Needed for the spin 0/1 part
//       remId = id2*1000+id1*100;
//       // Now decide if we have spin 0 diquark or spin 1 diquark
//       if(id1 == id2 || UseRandom::rndbool()) spin = 3; // spin 1
//       else spin = 1; // otherwise spin 0
//       remId += sign*spin;

//       // Create the remnant and set its momentum, also reset all of the decay 
//       // products from the hadron
//       //      PPtr newRemnant = new_ptr(Particle(getParticleData(remId)));
//       ShowerParticlePtr newRemnant = new_ptr(ShowerParticle(getParticleData(remId)));

//       newRemnant->setMomentum((1-part->x())*hadron->momentum());

//       for(int i = hadron->children().size()-1; i!= -1; i--) {
// 	PPtr child = hadron->children()[i];		
// 	hadron->abandonChild(child);
// 	if(part != child) ch->currentStep()->removeParticle(child);
//       }
      
//       // Add the remnant to the step, this will be changed again if the
//       // shower is vetoed. Set the colour connections as well
//       hadron->addChild(newRemnant);
//       hadron->addChild(part);
//       if (part->id() < 0) part->colourNeighbour(newRemnant);      
//       else part->antiColourNeighbour(newRemnant);      

//       for(int i = hadron->children().size()-1; i!= -1; i--) {
// 	PPtr child = hadron->children()[i];
//       }
//     }
//   }
//   // Do we veto the whole shower after the final state showering or do we
//   // seperately veto the initial state shower and final state shower?

//   // do timelike evolution of new particles here? 
//   // I think not before 1st reconstruction! 
// //   while(!particlesYetToShower.empty()) {
// //     tShowerParticlePtr part = particlesYetToShower.back();
// //     particlesYetToShower.pop_back();
// //     hasEmitted = hasEmitted || 
// //       _forwardEvolver->timeLikeShower(ch, showerVars, part, allShowerParticles);
// //   } 
  

// ShoKinPtr BackwardEvolver::forcedSplitting(const ShowerParticle &particle,
// 					   Energy lastQ, Energy minQ) {
//   // Now generate the new z and qtilde
//   Energy newQ;
//   double newZ,z0,z1;
//   Energy kinCutoff;
//   ShowerVarsPtr vars = _splittingGenerator->showerVariables();;
//   if ( particle.id() == ParticleID::g ) {
//     kinCutoff = (vars->kinScale() - 0.3*(particle.children()[0])->mass())/2.3;
//   } else {
//     kinCutoff = (vars->kinScale() - 0.3*particle.data().mass())/2.3;
//   }
//   // Generate z with the same distributions as for regular splittings
//   tSplittingFnPtr sf = particle.splitFun();
//   // Bounds on z
//   z0 = particle.x();
//   double yy = 1.+sqr(kinCutoff/lastQ)/2.;
//   z1 = yy - sqrt(sqr(yy)-1.); 
//   //  z1 = 1.;
//   bool cant;
//   do {
//     cant = false;
//     double randQ = UseRandom::rnd();
//     double randZ = UseRandom::rnd();
//     if(!sf) {
//       newZ = z0 + (z1-z0)*randZ;
//     } else {
//       newZ = sf->invIntegOverP(sf->integOverP(z0) 
// 			       + randZ*(sf->integOverP(z1) - 
// 					sf->integOverP(z0)));
//     }
//     // For the qtilde lets just start with a simple distribution
//     // weighted towards the lower value: dP/dQ = 1/Q -> Q(R) =
//     // Q0^(1-R) Qmax^R
//     //    newQ = pow(minQ,1-randQ)*pow(lastQ,randQ);
//     newQ = pow(kinCutoff, 1-randQ)*pow(lastQ,randQ);


//     // find out whether a next splitting would be possible
    
//     if (particle.id() != ParticleID::g) { 
//       //      Energy q0g = (vars->kinScale()-0.0015*GeV)/2.3;
//       Energy q0g = (vars->kinScale())/2.3;
//       double zm;

//       yy = 1.+sqr(q0g/newQ)/2.; // maximum z if next Q is max as as well 
//       zm = yy - sqrt(sqr(yy)-1.); 
//       yy = 1.+sqr(q0g/lastQ)/2.; 
//       double zm2 = yy - sqrt(sqr(yy)-1.); 
//       double xp = particle.x()/newZ;
//       if (xp > zm) cant = true;
//     }
//     // check kinematics...
//   } while(sqr((1.-newZ)*newQ) < newZ*sqr(kinCutoff));
//   //  } while(sqr((1.-newZ)*newQ) < newZ*sqr(kinCutoff) || cant);
  
//   if (cant) return ShoKinPtr();

//   Lorentz5Momentum p, n, pthis, ppartner, pcm;
//   if(particle.isFromHardSubprocess()) {
// //     pthis = particle.momentum();
// //     cout << "parent is " << particle.parents()[0]->id() << '\n';
// //     ppartner = particle.partners()[ShowerIndex::QCD]->momentum();
// //     pcm = pthis; 
// //     pcm.boost((pthis + ppartner).findBoostToCM());	  
// //     p = Lorentz5Momentum(0.0, pcm.vect());
// //     n = Lorentz5Momentum(0.0, -pcm.vect()); 
// //     p.boost( -(pthis + ppartner).findBoostToCM() );
// //     n.boost( -(pthis + ppartner).findBoostToCM() );
//     pcm = particle.parents()[0]->momentum();
//     p = Lorentz5Momentum(0.0, pcm.vect());
//     n = Lorentz5Momentum(0.0, -pcm.vect()); 
//   } else {
//     p = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
//       ->showerKinematics()->getBasis()[0];
//     n = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
//       ->showerKinematics()->getBasis()[1];
//   }
  
//   Ptr<IS_QtildaShowerKinematics1to2>::pointer showerKin = 
//     new_ptr(IS_QtildaShowerKinematics1to2(p, n));

//   // Phi is uniform
//   showerKin->qtilde(newQ);
//   showerKin->setResScale(vars->cutoffQScale(ShowerIndex::QCD));
//   showerKin->setKinScale(vars->kinScale()); 
//   showerKin->z(newZ);
//   showerKin->phi(2.*pi*UseRandom::rnd());

//   return showerKin;
// }
