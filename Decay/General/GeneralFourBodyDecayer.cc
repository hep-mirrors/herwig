// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralFourBodyDecayer class.
//

#include "GeneralFourBodyDecayer.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

GeneralFourBodyDecayer::GeneralFourBodyDecayer() {}

void GeneralFourBodyDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoing << _diagrams << _reftag << _reftagcc;
}

void GeneralFourBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoing >> _diagrams >> _reftag >> _reftagcc;
}

DescribeAbstractClass<GeneralFourBodyDecayer,DecayIntegrator>
describeGeneralFourBodyDecayer("Herwig::GeneralFourBodyDecayer",
			       "Herwig.so");

void GeneralFourBodyDecayer::Init() {

  static ClassDocumentation<GeneralFourBodyDecayer> documentation
    ("There is no documentation for the GeneralFourBodyDecayer class");

}

ParticleVector GeneralFourBodyDecayer::decay(const Particle & parent,
					     const tPDVector & children) const {
  assert(false);
  return ParticleVector();
}

int GeneralFourBodyDecayer::modeNumber(bool & cc, tcPDPtr parent,
				       const tPDVector & children) const {
  assert(false);
  return 0;
}

bool GeneralFourBodyDecayer::setDecayInfo(PDPtr incoming,vector<PDPtr> outgoing,
					  const vector<PrototypeVertexPtr> & process) {
  // set the member variables from the info supplied
  // external particles
  _incoming        = incoming;
  _outgoing        = outgoing;
  assert( _outgoing.size() == 4 );
  // convert and store the diagrams
  for(unsigned int ix=0;ix<process.size();++ix)
    _diagrams.push_back(NBDiagram(process[ix]));
  // Construct reference tags for testing in modeNumber function
  OrderedParticles refmode(_outgoing.begin(), _outgoing.end());
  OrderedParticles::const_iterator dit = refmode.begin();
  _reftag = _incoming->name() + "->";
  for( ; dit != refmode.end(); ++dit) {
    if( dit != refmode.begin() )  _reftag += string(",");
    _reftag += (**dit).name();
  }
  //CC-mode
  refmode.clear();
  _reftagcc = _incoming->CC() ? _incoming->CC()->name() : _incoming->name();
  _reftagcc += "->";
  for( unsigned int i = 0;  i < 3; ++i ) {
    if( _outgoing[i]->CC() ) refmode.insert( _outgoing[i]->CC() );
    else refmode.insert( _outgoing[i] );
  }
  dit = refmode.begin();
  for( ; dit != refmode.end(); ++dit) {
    if( dit != refmode.begin() )  _reftag += string(",");
    _reftagcc += (**dit).name();
  }
  // set the color factors
  cerr << "testing in setDecayInfo\n";
//   return setColourFactors();
  exit(0);
}
