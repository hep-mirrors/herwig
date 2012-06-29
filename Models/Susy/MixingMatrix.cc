// -*- C++ -*-
//
// MixingMatrix.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MixingMatrix class.
//

#include "MixingMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void MixingMatrix::persistentOutput(PersistentOStream & os) const {
  os << _theMixingMatrix << _theIds << _theSize;
}

void MixingMatrix::persistentInput(PersistentIStream & is, int) {
  is >> _theMixingMatrix >> _theIds >> _theSize;
}

ClassDescription<MixingMatrix> MixingMatrix::initMixingMatrix;
// Definition of the static class description member.

void MixingMatrix::Init() {

  static ClassDocumentation<MixingMatrix> documentation
    ("The MixingMatrix class implements the storage of the SUSY mixing "
     "matrices.");

}

void MixingMatrix::adjustPhase(long id) { 
  unsigned int irow(0);
  while(irow < size().first && _theIds[irow] != id) 
    ++irow;
  for(unsigned int c = 0; c < _theSize.second; ++c)
    _theMixingMatrix[irow][c] *= Complex(0., 1.);
}

ostream & Herwig::operator<<(ostream & os,const MixingMatrix & mix) {
  for(unsigned int ix=0;ix<mix.size().first;++ix) {
    for(unsigned int iy=0;iy<mix.size().second;++iy) {
      os << mix(ix,iy) << "\t";
    }
    os << "\n";
  }
  return os;
}
