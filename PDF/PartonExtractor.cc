// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonExtractor class.
//

#include "PartonExtractor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

void PartonExtractor::persistentOutput(ThePEG::PersistentOStream & os) const {
}

void PartonExtractor::persistentInput(ThePEG::PersistentIStream & is, int) {
}

ThePEG::ClassDescription<PartonExtractor> PartonExtractor::initPartonExtractor;
// Definition of the static class description member.

void PartonExtractor::Init() {

  static ThePEG::ClassDocumentation<PartonExtractor> documentation
    ("The PartonExtractor class of Herwig++ uses all the functionality from ThePEG"
     " with the only change being we integrate over the off-shell mass of the "
     " photon in processes where a photon entering the hard process is emitted"
     "from an incoming lepton beam");

}

ThePEG::pair<int,int> PartonExtractor::nDims(const ThePEG::PBPair & pbins) {
  // if photon from a lepton generate scale
  bool genscale[2]={false,false};
  for(unsigned int ix=0;ix<2;++ix) {
    ThePEG::PBPtr bin = ix==0 ? pbins.first : pbins.second;
    int id = abs(bin->particle()->id()); 
    if( ( id==ThePEG::ParticleID::eminus || id==ThePEG::ParticleID::muminus ) &&
	bin->parton()->id()==ThePEG::ParticleID::gamma )
      genscale[ix]=true;
  }
  return ThePEG::make_pair(pbins.first ->nDim(genscale[0]),
			   pbins.second->nDim(genscale[1]));
}
