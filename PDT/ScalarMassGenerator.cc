// -*- C++ -*-
//
// ScalarMassGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMassGenerator class.
//

#include "ScalarMassGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "GenericWidthGenerator.h"

using namespace Herwig;

void ScalarMassGenerator::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,GeV) << ounit(_mass1,GeV) << ounit(_mass2,GeV)
     << ounit(_m2plus,GeV2) << ounit(_m2minus,GeV2);
}

void ScalarMassGenerator::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,GeV) >> iunit(_mass1,GeV) >> iunit(_mass2,GeV)
     >> iunit(_m2plus,GeV2) >> iunit(_m2minus,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarMassGenerator,GenericMassGenerator>
describeHerwigScalarMassGenerator("Herwig::ScalarMassGenerator", "Herwig.so");

void ScalarMassGenerator::Init() {

  static ClassDocumentation<ScalarMassGenerator> documentation
    ("The ScalarMassGenerator class inherits from the "
     "GenericMassGenerator class and includes finite width effects "
     "for the scalar f_0 and a_0 mesons.",
     "Finite width effects for the scalar $f_0$ and $a_0$ mesons follow \\cite{Flatte:1976xu}.",
     "%\\cite{Flatte:1976xu}\n"
     "\\bibitem{Flatte:1976xu}\n"
     "  S.~M.~Flatte,\n"
     "   ``Coupled - Channel Analysis Of The Pi Eta And K Anti-K Systems Near K Anti-K\n"
     "  Threshold,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 63}, 224 (1976).\n"
     "  %%CITATION = PHLTA,B63,224;%%\n"
     );

  static ParVector<ScalarMassGenerator,Energy> interfacecoupling
    ("Coupling",
     "The coupling",
     &ScalarMassGenerator::_coupling,
     GeV, 0, ZERO, ZERO, Constants::MaxEnergy, false, false, true);

  static ParVector<ScalarMassGenerator,Energy> interfacemass1
    ("Mass1",
     "The mass for first particle",
     &ScalarMassGenerator::_mass1,
     GeV, 0, ZERO, ZERO, Constants::MaxEnergy, false, false, true);

  static ParVector<ScalarMassGenerator,Energy> interfacemass2
    ("Mass2",
     "The mass for second particle",
     &ScalarMassGenerator::_mass2,
     GeV, 0, ZERO, ZERO, Constants::MaxEnergy, false, false, true);

}

void ScalarMassGenerator::dataBaseOutput(ofstream & output,bool header) {
  if(header) output << "update Mass_Generators set parameters=\"";
  GenericMassGenerator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_coupling.size();++ix)
    output << "insert " << fullName() << ":Coupling " 
	   << ix << " " << _coupling[ix]/GeV << "\n";
  for(unsigned int ix=0;ix<_mass1.size();++ix)
    output << "insert " << fullName() << ":Mass1 " 
	   << ix << " " << _mass1[ix]/GeV << "\n";
  for(unsigned int ix=0;ix<_mass2.size();++ix)
    output << "insert " << fullName() 
	   << ":Mass2 " << ix << " " << _mass2[ix]/GeV << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void ScalarMassGenerator::doinit() {
  if(_coupling.size()!=_mass1.size()||
     _coupling.size()!=_mass2.size()) 
    throw InitException() << "Parameter vectors have inconsistent sizes in "
			  << "ScalarMassGenerator::doinit()"
			  << Exception::runerror;
  // initialise the local variables
  for(unsigned int ix=0;ix<_mass1.size();++ix) {
    _m2plus .push_back(sqr(_mass1[ix]+_mass2[ix]));
    _m2minus.push_back(sqr(_mass1[ix]-_mass2[ix]));
  }
  // rest of the initialisation is handled in the base class
  GenericMassGenerator::doinit();
}

double ScalarMassGenerator::weight(Energy q,int shape) const {
  useMe();
  Energy2 q2    = sqr(q);
  Energy2 mass2 = sqr(nominalMass());
  Energy2 mwidth= nominalMass()*nominalWidth();
  Energy2 gamma[2]={0.*MeV2,ZERO};
  if(shape==1) {
    for(unsigned int ix=0;ix<_coupling.size();++ix) {
      double lambda = (mass2-_m2plus[ix])*(mass2-_m2minus[ix])/sqr(mass2);
      if(lambda>=0.) gamma[0]+=sqr(_coupling[ix])*sqrt( lambda)*q/nominalMass();
      else           gamma[1]+=sqr(_coupling[ix])*sqrt(-lambda)*q/nominalMass();
    }
  }
  else {
    for(unsigned int ix=0;ix<_coupling.size();++ix) {
      double lambda = (q2-_m2plus[ix])*(q2-_m2minus[ix])/sqr(q2);
      if(lambda>=0.) gamma[0]+=sqr(_coupling[ix])*sqrt( lambda);
      else           gamma[1]+=sqr(_coupling[ix])*sqrt(-lambda);
    }
  }
  Energy2 numer(ZERO);
  if(shape==2)      numer = nominalMass()*gamma[0]/q;
  else if(shape==3) numer = nominalMass()*nominalWidth();
  else              numer = gamma[0];
  complex<Energy2> denom = (shape==2) ?
    mass2 - q2 +mass2/q2*complex<Energy2>(gamma[1],-gamma[0]) :
    mass2 - q2 +complex<Energy2>(gamma[1],-gamma[0]); 
  // complex denominantor
  Energy4 den = real(denom*conj(denom));
  return numer/den*(sqr(mass2-q2)+sqr(mwidth))/Constants::pi/mwidth;
}

InvEnergy2 ScalarMassGenerator::BreitWignerWeight(Energy q, int shape) const {
  useMe();
  Energy2 q2    = sqr(q);
  Energy2 mass2 = sqr(nominalMass());
  Energy2 gamma[2]={0.*MeV2,ZERO};
  if(shape==1) {
    for(unsigned int ix=0;ix<_coupling.size();++ix) {
      double lambda = (mass2-_m2plus[ix])*(mass2-_m2minus[ix])/sqr(mass2);
      if(lambda>=0.) gamma[0]+=sqr(_coupling[ix])*sqrt( lambda)*q/nominalMass();
      else           gamma[1]+=sqr(_coupling[ix])*sqrt(-lambda)*q/nominalMass();
    }
  }
  else {
    for(unsigned int ix=0;ix<_coupling.size();++ix) {
      double lambda = (q2-_m2plus[ix])*(q2-_m2minus[ix])/sqr(q2);
      if(lambda>=0.) gamma[0]+=sqr(_coupling[ix])*sqrt( lambda);
      else           gamma[1]+=sqr(_coupling[ix])*sqrt(-lambda);
    }
  }
  Energy2 numer(ZERO);
  if(shape==2)      numer = nominalMass()*gamma[0]/q;
  else if(shape==3) numer = nominalMass()*nominalWidth();
  else              numer = gamma[0];
  complex<Energy2> denom = (shape==2) ?
    mass2 - q2 +mass2/q2*complex<Energy2>(gamma[1],-gamma[0]) :
    mass2 - q2 +complex<Energy2>(gamma[1],-gamma[0]); 
  // complex denominantor
  Energy4 den = real(denom*conj(denom));
  return numer/den/Constants::pi;
}
