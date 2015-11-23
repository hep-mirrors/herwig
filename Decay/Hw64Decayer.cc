// -*- C++ -*-
//
// Hw64Decayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Hw64Decayer class.
//

#include "Hw64Decayer.h"
#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include "Herwig/PDT/GenericMassGenerator.h"

using namespace Herwig;

void Hw64Decayer::Init() {

   static ClassDocumentation<Hw64Decayer> documentation
     ("Class to decay all particles in HERWIG by the algorithms used in HERWIG 6.4",
      "Some decays used the Fortran HERWIG decay algorithm \\cite{Corcella:2000bw}.",
      "%\\cite{Corcella:2000bw}\n"
      "\\bibitem{Corcella:2000bw}\n"
      "  G.~Corcella {\\it et al.},\n"
      "  %``HERWIG 6.5: an event generator for Hadron Emission Reactions With\n"
      "  %Interfering Gluons (including supersymmetric processes),''\n"
      "  JHEP {\\bf 0101} (2001) 010\n"
      "  [arXiv:hep-ph/0011363].\n"
      "  %%CITATION = JHEPA,0101,010;%%\n"
      );

  static Switch<Hw64Decayer,int> interfaceMECode
    ("MECode",
     "The code for the ME type to use in the decay",
     &Hw64Decayer::MECode, 0, false, false);
  static SwitchOption interfaceMECodePhaseSpace
    (interfaceMECode,
     "PhaseSpace",
     "Use a phase-space distribution",
     0);
  static SwitchOption interfaceMECodeFreeVA
    (interfaceMECode,
     "FreeVA",
     "Free Massless (V-A)*(V-A) ME",
     100);
  static SwitchOption interfaceMECodeBoundVA
    (interfaceMECode,
     "BoundVA",
     "Bound Massless (V-A)*(V-A) ME",
     101);

  static Parameter<Hw64Decayer,unsigned int> interfaceMassTry
    ("MassTry",
     "The maximum number of attempts to generate the off-shell masses of the"
     " decay products.",
     &Hw64Decayer::_masstry, 50, 1, 1000,
     false, false, Interface::limited);

}

ClassDescription<Hw64Decayer> Hw64Decayer::initHw64Decayer;

bool Hw64Decayer::accept(tcPDPtr, const tPDVector & children) const  {
  return children.size() <= 3;
}

ParticleVector Hw64Decayer::decay(const Particle & p, 
				  const tPDVector & children) const {
  useMe();
  // storage for the decay products and number of decay products
  ParticleVector rval;
  unsigned int numProds(children.size());
  // check that it is possible to kinematically perform the decay
  Energy minmass(ZERO);
  vector<Energy> minmasses(numProds);
  vector<tcGenericMassGeneratorPtr> massgen(numProds,tcGenericMassGeneratorPtr());
  tcMassGenPtr mtemp;
  for(unsigned int ix=0;ix<numProds;++ix) {
    minmasses[ix]=children[ix]->massMin();
    minmass+=minmasses[ix];
    mtemp=children[ix]->massGenerator();
    if(mtemp) massgen[ix]=dynamic_ptr_cast<tcGenericMassGeneratorPtr>(mtemp);
  }
  // check not decaying a massless particle
  if(p.mass() < 0.000001*MeV) 
    throw Exception() << "HwDecayer called on a particle with no mass " 
		      << p.PDGName() << ", " << p.mass()/MeV 
		      << Exception::eventerror;
  // throw a veto if not kinematically possible
  if(minmass>p.mass()&&numProds!=1) throw Veto();
  // Create a vectors for momenta and masses
  vector<Lorentz5Momentum> products(numProds);
  vector<Energy> masses(numProds);
  // now generate the masses of the particles starting with a random one
  // to avoid bias
  unsigned int ntry(0);
  Energy outmass;
  if(numProds!=1) {
    do {
      unsigned int istart=UseRandom::irnd(numProds);
      outmass=ZERO;
      for(unsigned int ix=istart;ix<numProds;++ix) { 
	masses[ix] = massgen[ix] ?
	  massgen[ix]->mass(*(children[ix]),minmasses[ix],
			    p.mass()-minmass+minmasses[ix]) :
	  children[ix]->generateMass();
	outmass+=masses[ix];
	if(outmass>p.mass()) break;
      }
      for(unsigned int ix=0;ix<istart;++ix) {
	masses[ix] = massgen[ix] ?
	  massgen[ix]->mass(*(children[ix]),minmasses[ix],
			    p.mass()-minmass+minmasses[ix]) :
	  children[ix]->generateMass();
	outmass+=masses[ix];
	if(outmass>p.mass()) break;
      }
      ++ntry;
    }
    while(ntry<_masstry&&outmass>p.mass());
    if(outmass>p.mass()) throw Veto();
  }
  else masses[0]=p.mass();
  for(unsigned int ix=0;ix<numProds;++ix) products[ix].setMass(masses[ix]);
  // The K -> KL0 and KS0
  if(numProds == 1) {
    products[0]=p.momentum();
  }
  // 2 Body Decay
  else if(numProds == 2) {
    double CosAngle, AzmAngle; 
    Kinematics::generateAngles(CosAngle, AzmAngle);
    Kinematics::twoBodyDecay(p.momentum(), masses[0], masses[1],
			     CosAngle, AzmAngle, products[0], products[1]);
  }
  // Three Body Decay
  else if(numProds == 3) {
    // Free Massless (V-A)*(V-A) ME
    if(MECode == 100) {
      Kinematics::threeBodyDecay(p.momentum(), products[0], products[1], 
				 products[2], &VAWt);
    } 
    // Bound Massless (V-A)*(V-A) ME
    else if(MECode == 101) {
      Energy2 wtmx, dot1, dot2;
      Energy4 wtmx2;
      double xs;
      wtmx = ( (p.mass() - products[2].mass())*(p.mass() + products[2].mass())
	       + (products[1].mass() - products[0].mass())*
	       (products[1].mass() + products[0].mass()))/2.0;
      wtmx2 = sqr(wtmx);
      // Find sum of masses of constituent particles 
      int IPDG = abs(p.id());
      Energy m1, m2, m3;
      if(IPDG >= 1000)
	m1 = generator()->getParticleData((IPDG/1000)%10)->constituentMass();
      else
	m1 = ZERO;
      m2 = generator()->getParticleData((IPDG/100)%10)->constituentMass();
      m3 = generator()->getParticleData((IPDG/10)%10)->constituentMass();
      xs = 1.0 - Math::absmax<Energy>(m1, Math::absmax<Energy>(m2, m3))/(m1+m2+m3);
      // Do decay, repeat until meets condition
      do {
	Kinematics::threeBodyDecay(p.momentum(), products[1], products[2], 
				   products[0], &VAWt);
	dot1 = p.momentum().dot(products[1]);
	dot2 = p.momentum().dot(products[0]);
      } 
      while(dot1*(wtmx-dot1-xs*dot2) < UseRandom::rnd()*wtmx2);
    }
    // Three Body via phase space
    else {
      Kinematics::threeBodyDecay(p.momentum(), products[0], products[1], 
				 products[2]);
    }
  }  
  if(products[0] == Lorentz5Momentum()) return ParticleVector();
  cPDVector productParticles(numProds);
  for(unsigned int ix=0;ix<numProds;++ix) {
    productParticles[ix]=children[ix];
  }
  // set the momenta of the particles and return the answer
  setParticleMomentum(rval, productParticles, products);
  return rval;
}
   
void Hw64Decayer::persistentOutput(PersistentOStream &os) const  {
  os << MECode << _masstry;
}

void Hw64Decayer::persistentInput(PersistentIStream &is, int) {
  is >> MECode >> _masstry;
}

double Hw64Decayer::VAWt(Energy2 t0, Energy2 t1, Energy2 t2, InvEnergy4 t3) { 
  return (t1-t0)*(t0-t2)*t3;
}

void Hw64Decayer::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  output << "newdef " << name() << ":MECode "  << MECode << " \n";
  output << "newdef " << name() << ":MassTry " << _masstry << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}


