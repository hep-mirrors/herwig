// -*- C++ -*-
//
// SMHiggsMassGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsMassGenerator class.
//

#include "SMHiggsMassGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void SMHiggsMassGenerator::persistentOutput(PersistentOStream & os) const {
  os << _shape << _hwidth;
}

void SMHiggsMassGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _shape >> _hwidth;
}

ClassDescription<SMHiggsMassGenerator> SMHiggsMassGenerator::initSMHiggsMassGenerator;
// Definition of the static class description member.

void SMHiggsMassGenerator::Init() {

  static ClassDocumentation<SMHiggsMassGenerator> documentation
    ("The SMHiggsMassGenerator class implements the mass distribution for the"
     " Higgs boson as in hep-ph/9505211.",
     "The Higgs mass was distributed as in \\cite{Seymour:1995qg}.",
     "%\\cite{Seymour:1995qg}\n"
     "\\bibitem{Seymour:1995qg}\n"
     "  M.~H.~Seymour,\n"
     "  %``The Higgs boson line shape and perturbative unitarity,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 354}, 409 (1995)\n"
     "  [arXiv:hep-ph/9505211].\n"
     "  %%CITATION = PHLTA,B354,409;%%\n"
     );

  static Switch<SMHiggsMassGenerator,unsigned int> interfaceHiggsShape
    ("HiggsShape",
     "The distribution for the Higgs mass",
     &SMHiggsMassGenerator::_shape, 1, false, false);
  static SwitchOption interfaceHiggsShapeNormal
    (interfaceHiggsShape,
     "Normal",
     "The standard running width distribution",
     0);
  static SwitchOption interfaceHiggsShapeImproved
    (interfaceHiggsShape,
     "Improved",
     "The improved shape of hep-ph/9505211",
     1);
}

bool SMHiggsMassGenerator::accept(const ParticleData & part) const {
  if(part.id()!=ParticleID::h0) return false;
  return GenericMassGenerator::accept(part);
}

void SMHiggsMassGenerator::doinit() {
  if(particle()->massGenerator()==this) { 
    if(particle()->widthGenerator()) {
      _hwidth=dynamic_ptr_cast<GenericWidthGeneratorPtr>(particle()->widthGenerator());
    }
    if(!_hwidth) throw InitException() 
		   << "Must be using the Herwig::GenericWidthGenerator in "
		   << "SMHiggsMassGenerator::doinit()" << Exception::runerror;
  }
  GenericMassGenerator::doinit();
}

void SMHiggsMassGenerator::dataBaseOutput(ofstream & output,bool header) {
  if(header) output << "update Mass_Generators set parameters=\"";
  output << "newdef " << fullName() << ":BreitWignerShape "   << _shape << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
