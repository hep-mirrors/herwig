// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonRemnants class.
//

#include "BaryonRemnants.h"
#include "Herwig++/Hadronization/Remnant.h"
#include <ThePEG/Vectors/LorentzVector.h>
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

void BaryonRemnants::boostRemnants(PartonBinInstance &pb) const {
  const LorentzMomentum pr = pb.particle()->momentum() - 
                             pb.parton()->momentum();
  const PVector &rem = pb.remnants();
  Energy2 s = Utilities::sumMomentum(rem).m2();
  Utilities::setMomentumFromCMS(rem.begin(), rem.end(), s, pr);
  for(int i = 0, N=rem.size(); i<N; i++)
    rem[i]->scale(pb.particle()->momentum().mass2());
}

Lorentz5Momentum BaryonRemnants::generate(PartonBinInstance & pb, 
					  const double *,
					  Energy2, 
					  const LorentzMomentum &parent) 
  const {
  
  LorentzMomentum p(0*GeV, 0*GeV, parent.rho(), parent.e());
  // set the weight for the remnant
  pb.remnantWeight(1.0);
  // set an empty remnant info as this feature of ThePEG is useless
  if(!pb.remnantInfo()) {
    RemIPtr Rem=new_ptr(RemInfoBase());
    pb.remnantInfo(Rem);
  }
  // calculate the momentum of the extracted parton
  double x = pb.xi();
  p = LorentzMomentum(0*GeV, 0*GeV, parent.rho(), parent.e());
  p = parent*x;
  // create the remnant
  PPtr remnant=new_ptr(Remnant(pb,parent-p));
  // insert the remnant in the parton bin
  //  cout << '\n' << *remnant << '\n';

  PVector rem(1,remnant);
  pb.remnants(rem);
  // return the momentum of the particle

  return p;
}

NoPIOClassDescription<BaryonRemnants> BaryonRemnants::initBaryonRemnants;

void BaryonRemnants::Init() {}
