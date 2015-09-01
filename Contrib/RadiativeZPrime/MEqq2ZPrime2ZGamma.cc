// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2ZPrime2ZGamma class.
//

#include "MEqq2ZPrime2ZGamma.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "RadiativeZPrimeModel.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "AnomalousVVVVertex.h"

using namespace RadiativeZPrime;

MEqq2ZPrime2ZGamma::MEqq2ZPrime2ZGamma()  : _maxflavour(5) {
  vector<unsigned int> mopt(2,1);
  mopt[0]=2;
  massOption(mopt);
  rescalingOption(2);
}

void MEqq2ZPrime2ZGamma::doinit() {
  HwMEBase::doinit();
  _zPrime = getParticleData(32); 
  tcSMPtr sm = generator()->standardModel();
  tcRadiativeZPrimeModelPtr model = 
    dynamic_ptr_cast<tcRadiativeZPrimeModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the RadiativeZPrimeModel in "
			       << "MEqq2ZPrime2ZGamma::doinit()" << Exception::abortnow;
  _theFFZPrimeVertex = model->vertexFFZPrime();
  _theGammaZPrimeZVertex = model->vertexGammaZPrimeZ();
}

void MEqq2ZPrime2ZGamma::getDiagrams() const {
  tcPDPtr Z0    = getParticleData(ParticleID::Z0);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  for(unsigned int i = 1; i <= _maxflavour; ++i) {
    tcPDPtr q  = getParticleData(long(i));
    tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, _zPrime, 3, Z0, 3, gamma, -1)));
  }
}

Energy2 MEqq2ZPrime2ZGamma::scale() const {
  return sHat();
}

double MEqq2ZPrime2ZGamma::me2() const {
  vector<SpinorWaveFunction>    fin;
  vector<SpinorBarWaveFunction> ain;
  vector<VectorWaveFunction> gammaout,Zout;
  SpinorWaveFunction       q(meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbar(meMomenta()[1],mePartonData()[1],incoming);
  VectorWaveFunction   Zwave(meMomenta()[2],mePartonData()[2],outgoing);
  VectorWaveFunction   Gwave(meMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix!=1) { 
      Gwave.reset(ix);
      gammaout.push_back(Gwave);
    }
    if(ix!=2) {
      q.reset(ix)   ; fin.push_back(q);
      qbar.reset(ix); ain.push_back(qbar);
    }
    Zwave.reset(ix);
    Zout.push_back(Zwave);
  }
  return qqME(fin,ain,Zout,gammaout,false);
}

unsigned int MEqq2ZPrime2ZGamma::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2ZPrime2ZGamma::orderInAlphaEW() const {
  return 4;
}

Selector<MEBase::DiagramIndex>
MEqq2ZPrime2ZGamma::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) sel.insert(1., i);
  return sel;
}

Selector<const ColourLines *>
MEqq2ZPrime2ZGamma::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines c1("1 -2");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c1);
  return sel;
}


void MEqq2ZPrime2ZGamma::persistentOutput(PersistentOStream & os) const {
  os << _theFFZPrimeVertex << _theGammaZPrimeZVertex << _zPrime << _maxflavour;
}

void MEqq2ZPrime2ZGamma::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZPrimeVertex >> _theGammaZPrimeZVertex >> _zPrime >> _maxflavour;
}

ClassDescription<MEqq2ZPrime2ZGamma> MEqq2ZPrime2ZGamma::initMEqq2ZPrime2ZGamma;
// Definition of the static class description member.

void MEqq2ZPrime2ZGamma::Init() {

  static ClassDocumentation<MEqq2ZPrime2ZGamma> documentation
    ("The MEqq2ZPrime2ZGamma class implements the matrix element for q qbar -> Z gamma"
     " via a resonant Z' in the RadiativeZPrime Model");

  static Parameter<MEqq2ZPrime2ZGamma,unsigned int> interfaceMaxFlavour
    ("MaxFlavour",
     "The heaviest incoming quark flavour this matrix element is allowed to handle",
     &MEqq2ZPrime2ZGamma::_maxflavour, 5, 1, 6,
     false, false, Interface::limited);

}

double MEqq2ZPrime2ZGamma::qqME(vector<SpinorWaveFunction>    & fin ,
				vector<SpinorBarWaveFunction> & ain ,
				vector<VectorWaveFunction> & Zout,
				vector<VectorWaveFunction> & gammaout,
				bool calc) const {
  // scale
  Energy2 mb2(scale());
  // matrix element to be stored
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1,PDT::Spin1);
  // declare the variables we need
  unsigned int ihel1,ihel2,ohel1,ohel2;
  VectorWaveFunction inter;
  double me(0.);
  Complex diag;
  // sum over helicities to get the matrix element
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      // intermediate for Z'
      inter=_theFFZPrimeVertex->evaluate(mb2,1,_zPrime,fin[ihel1],ain[ihel2]);
      for(ohel1=0;ohel1<3;++ohel1) {
	for(ohel2=0;ohel2<2;++ohel2) {
	  diag = _theGammaZPrimeZVertex->evaluate(mb2,gammaout[ohel2],inter,Zout[ohel1]);
	  me += norm(diag);
	  if(calc) menew(ihel1,ihel2,ohel1,2*ohel2) = diag;
	}
      }
    }
  }
  if(calc) _me.reset(menew);
  // spin and colour factor
  return me/12.;
}

void MEqq2ZPrime2ZGamma::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  if(hard[0]->id()<0) swap(order[0],order[1]);
  if(hard[2]->id()==ParticleID::gamma) swap(order[2],order[3]);
  vector<SpinorWaveFunction>    fin;
  vector<SpinorBarWaveFunction> ain;
  vector<VectorWaveFunction> gammaout,Zout;
  SpinorWaveFunction(   fin  ,hard[order[0]],incoming,false,true);
  SpinorBarWaveFunction(ain  ,hard[order[1]],incoming,false,true);
  VectorWaveFunction(Zout    ,hard[order[2]],outgoing,true,false,true);
  VectorWaveFunction(gammaout,hard[order[3]],outgoing,true,true ,true);
  gammaout[1]=gammaout[2];
  qqME(fin,ain,Zout,gammaout,true);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) 
    hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
}
