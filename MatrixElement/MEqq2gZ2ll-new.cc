// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2gZ2ll class.
//

#include "MEqq2gZ2ll-new.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEqq2gZ2ll.tcc"
#endif

#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
using namespace Herwig;

MEqq2gZ2ll::~MEqq2gZ2ll() {}

void MEqq2gZ2ll::getDiagrams() const {
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  for(int i = 1; i <= _maxflavour; ++i) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();
    int step = (_withNeutrinos ? 1 : 2);
    for(long leptonID = ParticleID::eminus; 
	leptonID <= ParticleID::nu_tau; leptonID+=step) {
      tcPDPtr lm = getParticleData(leptonID);
      tcPDPtr lp = lm->CC();
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0, 3, lm, 3, lp, -1)));
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, gamma, 3, lm, 3, lp, -2)));
    }
  }
}

Energy2 MEqq2gZ2ll::scale() const {
  return sHat();
}


double MEqq2gZ2ll::me2() const {
  // declare the variables we need
  unsigned int inhel1,inhel2,outhel1,outhel2;
  VectorWaveFunction inter;
  double me[3]={0.,0.,0.};
  complex<double> diag1,diag2;
  int iq,iqbar,ilp,ilm;
  // get the order right
  if(mePartonData()[0]->id()>0) {iq=0;iqbar=1;}
  else {iq=1;iqbar=0;}
  if(mePartonData()[2]->id()>0) {ilm=2;ilp=3;}
  else {ilm=3;ilp=2;}
  // setup momenta and particle data for the wavefunctions 
  SpinorWaveFunction    ein( meMomenta()[iq       ],mePartonData()[iq       ],incoming);
  SpinorBarWaveFunction pin( meMomenta()[iqbar    ],mePartonData()[iqbar    ],incoming);
  SpinorBarWaveFunction fout(meMomenta()[ilm      ],mePartonData()[ilm      ],outgoing);
  SpinorWaveFunction    aout(meMomenta()[ilp      ],mePartonData()[ilp      ],outgoing);
  // sum over helicities to get the matrix element
  for(inhel1=0;inhel1<2;++inhel1)
    {
      ein.reset(inhel1);
      for(inhel2=0;inhel2<2;++inhel2)
	{
	  pin.reset(inhel2);
	  for(outhel1=0;outhel1<2;++outhel1)
	    {
	      fout.reset(outhel1);
	      for(outhel2=0;outhel2<2;++outhel2)
		{
		  aout.reset(outhel2);
		  // first the Z exchange diagram
		  inter = _theFFZVertex->evaluate(sHat(),1,_Z0,ein,pin);
		  diag1 = _theFFZVertex->evaluate(sHat(),aout,fout,inter);
		  // first the photon exchange diagram
		  inter = _theFFPVertex->evaluate(sHat(),1,_gamma,ein,pin);
		  diag2 = _theFFPVertex->evaluate(sHat(),aout,fout,inter);
		  // add up squares of individual terms
		  me[1] += real(diag1*conj(diag1));
		  me[2] += real(diag2*conj(diag2));
		  // the full thing including interference
		  diag1 +=diag2;
		  me[0] += real(diag1*conj(diag1));
		}
	    }
	}
    }
  // results
  // factor 12 from 4 helicity and 3 colour
  for(int ix=0;ix<3;++ix){me[ix]/=12.0;}
  DVector save;
  save.push_back(me[1]);
  save.push_back(me[2]);
  meInfo(save);
  //  meInfo(save << lastCont << lastBW);
  return me[0];
}

unsigned int MEqq2gZ2ll::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2gZ2ll::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEqq2gZ2ll::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
  }
  return sel;
}

Selector<const ColourLines *>
MEqq2gZ2ll::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c("1 -2");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c);
  return sel;
}


void MEqq2gZ2ll::persistentOutput(PersistentOStream & os) const {
  os <<  _maxflavour << _theFFZVertex << _theFFPVertex 
     << _gamma << _Z0 << _withNeutrinos; 

}

void MEqq2gZ2ll::persistentInput(PersistentIStream & is, int) {
  is >> _maxflavour >> _theFFZVertex >> _theFFPVertex 
     >> _gamma >> _Z0 >> _withNeutrinos; 
}

ClassDescription<MEqq2gZ2ll> MEqq2gZ2ll::initMEqq2gZ2ll;
// Definition of the static class description member.

void MEqq2gZ2ll::Init() {

  static Parameter<MEqq2gZ2ll,int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEqq2gZ2ll::_maxflavour, 5, 0, 8, false, false, true);

  static Switch<MEqq2gZ2ll,bool> interfaceWithNeutrinos
    ("WithNeutrinos",
     "Should this ME produce neutrinos?",
     &MEqq2gZ2ll::_withNeutrinos, false, false, false);
  static SwitchOption interfaceWithNeutrinosNeutrinosOff
    (interfaceWithNeutrinos,
     "NeutrinosOff",
     "MEqq2gZ2ll does not produce neutrinos.",
     false);
  static SwitchOption interfaceWithNeutrinosNeutrinosOn
    (interfaceWithNeutrinos,
     "NeutrinosOn",
     "MEqq2gZ2ll also produces neutrinos.",
     true);

  static ClassDocumentation<MEqq2gZ2ll> documentation
    ("The \\classname{MEqq2gZ2ll} class implements the matrix element for"
     "q qbar to leptons via Z and photon exchange using helicity amplitude"
     "techniques");

}

