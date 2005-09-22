// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonRemnants class.
//

#include "BaryonRemnants.h"
#include <ThePEG/CLHEPWrap/LorentzVector.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/PDT/StandardMatchers.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/PDF/PDFBase.h>
#include <ThePEG/Utilities/UtilityBase.h>
#include <ThePEG/Utilities/Direction.h>
#include <ThePEG/PDT/EnumParticles.h>

using namespace ThePEG;
using namespace Herwig;

BaryonRemnants::~BaryonRemnants() {}

// This handles all partons (udscbtg)
bool BaryonRemnants::canHandle(tcPDPtr particle, 
			       const cPDVector &partons) const {
  if(!BaryonMatcher::Check(*particle)) return false;
  for(int i = 0, N = partons.size(); i < N; ++i)
    if(!StandardQCDPartonMatcher::Check(*partons[i])) return false;
  return true;
}

void BaryonRemnants::createRemnants(PartonBinInstance &pb) const {
  const LorentzMomentum pr = pb.particle()->momentum() - 
                             pb.parton()->momentum();
  const PVector &rem = pb.remnants();
  Energy2 s = Utilities::sumMomentum(rem).m2();
  Utilities::setMomentumFromCMS(rem.begin(), rem.end(), s, pr);
  for(int i = 0, N=rem.size(); i<N; i++)
    rem[i]->scale(pb.particle()->momentum().mass2());
}

Lorentz5Momentum BaryonRemnants::generate(PartonBinInstance & pb, 
					  const double * r,
					  Energy2 scale, 
					  const LorentzMomentum &parent) 
  const {
  LorentzMomentum p(0.0, 0.0, parent.rho(), parent.e());
  double x = pb.xi();
  //double eps = pb.eps();
  pb.remnantWeight(1.0);
  PVector rem;
  
  typedef Ptr<BaryonRemInfo>::pointer BRemIPtr;
  BRemIPtr ip;
  if(!(ip=dynamic_ptr_cast<BRemIPtr>(pb.remnantInfo()))) {
    // This is the first time for this parton bin, fill in info
    ip = new_ptr(BaryonRemInfo());
    pb.remnantInfo(ip);
    BaryonRemInfo &info = *ip;
    info.iq = pb.partonData()->id();
    int pid = pb.particleData()->id();
    info.sign = pid < 0 ? -1 : 1;
   
    // Get valence flavours
    info.flav = info.valenceFlav = vector<int>(3);
    info.flav[0] = info.valenceFlav[0] = (pid = abs(pid)/10)%10;
    info.flav[1] = info.valenceFlav[1] = (pid/=10)%10;
    info.flav[2] = info.valenceFlav[2] = (pid/=10)%10;
    vector<int>::iterator v = find(info.valenceFlav.begin(), 
				   info.valenceFlav.end(),
				   info.sign*info.iq);
    // If we found it, then extracted parton is a valence parton
    if(v!=info.valenceFlav.end()) {
      info.valenceFlav.erase(v);
      info.maybeValence = true;
    } else info.maybeValence = false;

    int iq1,iq2,iq3;
    for(iq1 = 0; iq1 < 3; iq1++) {
      iq2 = (iq1+1)%3;
      iq3 = (iq2+1)%3;
      // This is the id of the diquark (spin 1) that accompanies iq1
      int idq = 1000*max(info.flav[iq2], info.flav[iq3]) +
 	         100*min(info.flav[iq2], info.flav[iq3]) + 3;
      info.flavsel.insert(3.0, make_pair(info.flav[iq1], idq));
      if(info.flav[iq2] == info.flav[iq3]) continue;
      // If they are different, we have spin 0 combination too
      info.flavsel.insert(1.0, make_pair(info.flav[iq1], idq-2));
    }
  }
  BaryonRemInfo &info = *ip;
  
  if(info.maybeValence){
    int idqr = 1000*max(info.valenceFlav[0], info.valenceFlav[1]) +
                100*min(info.valenceFlav[0], info.valenceFlav[1]) + 3;
    // Also include the chance for spin 0 combination
    if(info.valenceFlav[0] != info.valenceFlav[1] && rnd() < 0.25) idqr-=2;
    rem.push_back(getParticleData(info.sign*idqr)->produceParticle());
  } else {
    pair<int,int> rr = info.flavsel.select(UseRandom::rnd());
    int iqr = rr.first * info.sign;
    int idqr = rr.second * info.sign;
    /*** 
     * We just want to randomly pick a diquark and put it in as the remnant.
     * This needs to be handled properly at the end of the initial state shower
     * as this method assumes that the parton drawn from the pdf must initial 
     * state shower and the remnant is associated with the final element of 
     * that shower, not the one involved in the hard subprocess.
     * For the Parton Extractor we must provide valid colour connections, 
     * this is very annoying!
     ***/
    if(info.iq == ParticleID::g) {
      rem.push_back(getParticleData(info.sign*idqr)->produceParticle());
      rem.push_back(getParticleData(info.sign*iqr)->produceParticle());
    } else { // Must be a sea quark
      if(info.iq < 0) 
	rem.push_back(getParticleData(info.sign*iqr)->produceParticle());
      else
	rem.push_back(getParticleData(info.sign*idqr)->produceParticle());
    }
  }

  pb.remnants(rem);
  p = parent*x;
  return p;
}

NoPIOClassDescription<BaryonRemnants> BaryonRemnants::initBaryonRemnants;

void BaryonRemnants::Init() {}

