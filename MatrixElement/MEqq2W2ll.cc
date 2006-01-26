// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2W2ll class.
//

#include "MEqq2W2ll.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEqq2W2ll.tcc"
#endif
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

MEqq2W2ll::~MEqq2W2ll() {}

void MEqq2W2ll::getDiagrams() const {
  typedef std::vector<pair<ParticleID::ParticleCodes, 
    ParticleID::ParticleCodes> > Pairvector;

  Pairvector parentpair, childpair;
  parentpair.reserve(6);
  childpair.reserve(3);

  // don't even think of putting 'break' in here!
  switch(_maxflavour) {
  case 5:
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::ubar));
  case 4:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::cbar));
  case 3:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::ubar));
  case 2:
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::ubar));
  default:
    ;
  }
  
  childpair.push_back(make_pair(ParticleID::eminus,   ParticleID::nu_ebar));
  childpair.push_back(make_pair(ParticleID::muminus,  ParticleID::nu_mubar));
  childpair.push_back(make_pair(ParticleID::tauminus, ParticleID::nu_taubar));

  tcPDPtr Wminus = getParticleData(ParticleID::Wminus);
  tcPDPtr Wplus = getParticleData(ParticleID::Wplus);
  
  Pairvector::const_iterator parent = parentpair.begin();
  for (; parent != parentpair.end(); ++parent) {
    tcPDPtr qNeg1 = getParticleData(parent->first);
    tcPDPtr qNeg2 = getParticleData(parent->second);
    tcPDPtr qPos1 = qNeg1->CC();
    tcPDPtr qPos2 = qNeg2->CC();
    
    Pairvector::const_iterator child = childpair.begin();
    for (; child != childpair.end(); ++child) {
      tcPDPtr lNeg1 = getParticleData(child->first);
      tcPDPtr lNeg2 = getParticleData(child->second);
      tcPDPtr lPos1 = lNeg1->CC();
      tcPDPtr lPos2 = lNeg2->CC();
 
      add(new_ptr((Tree2toNDiagram(2), 
		   qNeg1, qNeg2, 
		   1, Wminus, 
		   3, lNeg1, 3, lNeg2, -1)));

      add(new_ptr((Tree2toNDiagram(2), 
		   qPos1, qPos2, 
		   1, Wplus, 
		   3, lPos1, 3, lPos2, -2)));
    }
  }
}

Energy2 MEqq2W2ll::scale() const {
  return sHat();
}

double MEqq2W2ll::me2() const {
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

  bool positive = mePartonData()[iq]->iCharge() 
    + mePartonData()[iqbar]->iCharge() == 3;

  bool negative = mePartonData()[iq]->iCharge() 
    + mePartonData()[iqbar]->iCharge() == -3;

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

		  if (positive) {
		    // the Wp exchange
		    inter = _theFFWVertex->evaluate(sHat(),1,_Wp,ein,pin);
		    diag1 = _theFFWVertex->evaluate(sHat(),aout,fout,inter);
		    // add up squares of individual terms
		    me[1] += real(diag1*conj(diag1));
		  } else if (negative) {
		    // the Wm exchange
		    inter = _theFFWVertex->evaluate(sHat(),1,_Wm,ein,pin);
		    diag2 = _theFFWVertex->evaluate(sHat(),aout,fout,inter);
		    // add up squares of individual terms
		    me[2] += real(diag2*conj(diag2));
		  } else {
		    cerr << "Error in MEqq2W2ll::me2() charge setup.\n";
		  }
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

unsigned int MEqq2W2ll::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2W2ll::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEqq2W2ll::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
  }
  return sel;
}

Selector<const ColourLines *>
MEqq2W2ll::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c;
  if (diag->id() == -1)
    c=ColourLines("1 -2");
  else 
    c=ColourLines("-1 2");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c);
  return sel;
}

void MEqq2W2ll::persistentOutput(PersistentOStream & os) const {
  os << _maxflavour << _theFFWVertex << _Wp << _Wm;
}

void MEqq2W2ll::persistentInput(PersistentIStream & is, int) {
  is >> _maxflavour >> _theFFWVertex >> _Wp >> _Wm;
  
}

ClassDescription<MEqq2W2ll> MEqq2W2ll::initMEqq2W2ll;
// Definition of the static class description member.

void MEqq2W2ll::Init() {

  static Parameter<MEqq2W2ll,int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEqq2W2ll::_maxflavour, 5, 0, 5, false, false, true);
  
  static ClassDocumentation<MEqq2W2ll> documentation
    ("The \\classname{MEqq2W2ll} class implements the matrix element for"
     "q qbar to leptons via W exchange using helicity amplitude"
     "techniques");

}

