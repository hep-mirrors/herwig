// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEChargedCurrentDIS class.
//

#include "MEChargedCurrentDIS.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"

using namespace Herwig;

MEChargedCurrentDIS::MEChargedCurrentDIS() 
  : _maxflavour(5), _massopt(0) {
  vector<unsigned int> mopt(2,1);
  mopt[1] = _massopt;
  massOption(mopt);
}

void MEChargedCurrentDIS::doinit() {
  DISBase::doinit();
  _wp = getParticleData(ThePEG::ParticleID::Wplus );
  _wm = getParticleData(ThePEG::ParticleID::Wminus);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() 
    << "Must be the Herwig StandardModel class in "
    << "MEChargedCurrentDIS::doinit" << Exception::abortnow;
  // vertices
  _theFFWVertex = hwsm->vertexFFW();
}


void MEChargedCurrentDIS::getDiagrams() const {
  // possible quarks
  typedef std::vector<pair<long,long> > Pairvector;
  Pairvector quarkpair;
  quarkpair.reserve(6);
  // don't even think of putting 'break' in here!
  switch(_maxflavour) {
  case 6:
    quarkpair.push_back(make_pair(ParticleID::s, ParticleID::t));
    quarkpair.push_back(make_pair(ParticleID::d, ParticleID::t));
    quarkpair.push_back(make_pair(ParticleID::b, ParticleID::t));
    [[fallthrough]];
  case 5:
    quarkpair.push_back(make_pair(ParticleID::b, ParticleID::c));
    quarkpair.push_back(make_pair(ParticleID::b, ParticleID::u));
    [[fallthrough]];
  case 4:
    quarkpair.push_back(make_pair(ParticleID::s, ParticleID::c));
    quarkpair.push_back(make_pair(ParticleID::d, ParticleID::c));
    [[fallthrough]];
  case 3:
    quarkpair.push_back(make_pair(ParticleID::s, ParticleID::u));
    [[fallthrough]];
  case 2:
    quarkpair.push_back(make_pair(ParticleID::d, ParticleID::u));
    [[fallthrough]];
  default:
    ;
  }
  // create the diagrams
  for(int il1=11;il1<=14;++il1) {
    int il2 = il1%2==0 ? il1-1 : il1+1;
    for(unsigned int iz=0;iz<2;++iz) {
      tcPDPtr lepin  = iz==1 ? getParticleData(il1) : getParticleData(-il1);
      tcPDPtr lepout = iz==1 ? getParticleData(il2) : getParticleData(-il2);
      tcPDPtr inter  = lepin->iCharge()-lepout->iCharge()==3 ? _wp : _wm;
      for(unsigned int iq=0;iq<quarkpair.size();++iq) {
	tcPDPtr first  = getParticleData(quarkpair[iq].first );
	tcPDPtr second = getParticleData(quarkpair[iq].second);
	if(inter==_wp) {
	  add(new_ptr((Tree2toNDiagram(3), lepin, inter, first       , 
		       1, lepout, 2, second     , -1)));
	  add(new_ptr((Tree2toNDiagram(3), lepin, inter, second->CC(), 
		       1, lepout, 2, first->CC(), -2)));
	}
	else {
	  add(new_ptr((Tree2toNDiagram(3), lepin, inter, second     , 
	  	       1, lepout, 2, first       , -1)));
	  add(new_ptr((Tree2toNDiagram(3), lepin, inter, first->CC(), 
	  	       1, lepout, 2, second->CC(), -2)));
	}
      }
    }
  }
}

unsigned int MEChargedCurrentDIS::orderInAlphaS() const {
  return 0;
}

unsigned int MEChargedCurrentDIS::orderInAlphaEW() const {
  return 2;
}

Selector<const ColourLines *>
MEChargedCurrentDIS::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c1("3 5");
  static ColourLines c2("-3 -5");
  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 )
    sel.insert(1.0, &c1);
  else
    sel.insert(1.0, &c2);
  return sel;
}

void MEChargedCurrentDIS::persistentOutput(PersistentOStream & os) const {
  os << _theFFWVertex << _maxflavour << _wp << _wm << _massopt;
}

void MEChargedCurrentDIS::persistentInput(PersistentIStream & is, int) {
  is >> _theFFWVertex >> _maxflavour >> _wp >> _wm >> _massopt;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEChargedCurrentDIS,DISBase>
describeHerwigMEChargedCurrentDIS("Herwig::MEChargedCurrentDIS", "HwMEDIS.so");

void MEChargedCurrentDIS::Init() {

  static ClassDocumentation<MEChargedCurrentDIS> documentation
    ("The MEChargedCurrentDIS class implements the matrix elements "
     "for leading-order charged current deep inelastic scattering");

  static Parameter<MEChargedCurrentDIS,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEChargedCurrentDIS::_maxflavour, 5, 2, 6, false, false, true);

  static Switch<MEChargedCurrentDIS,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the mass of the outgoing quarks",
     &MEChargedCurrentDIS::_massopt, 0, false, false);
  static SwitchOption interfaceMassOptionMassless
    (interfaceMassOption,
     "Massless",
     "Treat the outgoing quarks as massless",
     0);
  static SwitchOption interfaceMassOptionMassive
    (interfaceMassOption,
     "Massive",
     "Treat the outgoing quarks as massive",
     1);

}

Selector<MEBase::DiagramIndex>
MEChargedCurrentDIS::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) sel.insert(1., i);
  return sel;
}

double MEChargedCurrentDIS::helicityME(vector<SpinorWaveFunction>    & f1,
				       vector<SpinorWaveFunction>    & f2,
				       vector<SpinorBarWaveFunction> & a1,
				       vector<SpinorBarWaveFunction> & a2,
				       bool lorder, bool qorder, bool calc) const {
  // scale
  Energy2 mb2(scale());
  // matrix element to be stored
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // pick a W boson
  tcPDPtr ipart = (mePartonData()[0]->iCharge()-mePartonData()[1]->iCharge())==3 ?
    _wp : _wm;
  // declare the variables we need
  VectorWaveFunction inter;
  double me(0.);
  Complex diag;
  // sum over helicities to get the matrix element
  unsigned int hel[4];
  unsigned int lhel1,lhel2,qhel1,qhel2;
  for(lhel1=0;lhel1<2;++lhel1) {
    for(lhel2=0;lhel2<2;++lhel2) {
      // intermediate W
      inter = _theFFWVertex->evaluate(mb2,3,ipart,f1[lhel1],a1[lhel2]);
      for(qhel1=0;qhel1<2;++qhel1) {
	for(qhel2=0;qhel2<2;++qhel2) {
	  hel[0] = lhel1;
	  hel[1] = qhel1;
	  hel[2] = lhel2;
	  hel[3] = qhel2;
	  if(!lorder) swap(hel[0],hel[2]);
	  if(!qorder) swap(hel[1],hel[3]);
	  diag = _theFFWVertex->evaluate(mb2,f2[qhel1],a2[qhel2],inter);
	  me += norm(diag);
	  menew(hel[0],hel[1],hel[2],hel[3]) = diag;
	}
      }
    }
  }
  // spin and colour factor
  me *= 0.25;
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
			 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    me = menew.average(rho[0],rho[1]);
  }
  if(calc) _me.reset(menew);
  return me;
}

double MEChargedCurrentDIS::me2() const {
  vector<SpinorWaveFunction>    f1,f2;
  vector<SpinorBarWaveFunction> a1,a2;
  bool lorder,qorder;
  SpinorWaveFunction    l1,q1;
  SpinorBarWaveFunction l2,q2;
  // lepton wave functions
  if(mePartonData()[0]->id()>0) {
    lorder=true;
    l1 = SpinorWaveFunction   (meMomenta()[0],mePartonData()[0],incoming);
    l2 = SpinorBarWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
  }
  else {
    lorder=false;
    l1 = SpinorWaveFunction   (meMomenta()[2],mePartonData()[2],outgoing);
    l2 = SpinorBarWaveFunction(meMomenta()[0],mePartonData()[0],incoming);
  }
  // quark wave functions
  if(mePartonData()[1]->id()>0) {
    qorder = true;
    q1 = SpinorWaveFunction   (meMomenta()[1],mePartonData()[1],incoming);
    q2 = SpinorBarWaveFunction(meMomenta()[3],mePartonData()[3],outgoing);
  }
  else {
    qorder = false;
    q1 = SpinorWaveFunction   (meMomenta()[3],mePartonData()[3],outgoing);
    q2 = SpinorBarWaveFunction(meMomenta()[1],mePartonData()[1],incoming);
  }
  // wavefunctions for various helicities
  for(unsigned int ix=0;ix<2;++ix) {
    l1.reset(ix); f1.push_back(l1);
    l2.reset(ix); a1.push_back(l2);
    q1.reset(ix); f2.push_back(q1);
    q2.reset(ix); a2.push_back(q2);
  }
  return helicityME(f1,f2,a1,a2,lorder,qorder,false);
}

void MEChargedCurrentDIS::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // sort out the ordering
  unsigned int order[4]={0,1,2,3};
  bool lorder(true),qorder(true);
  if(abs(hard[0]->id())<6) swap(hard[0],hard[1]);
  if(abs(hard[2]->id())<6) swap(hard[2],hard[3]);
  if(hard[0]->id()<0) {
    swap(order[0],order[2]);
    lorder = false;
  }
  if(hard[1]->id()<0) {
    swap(order[1],order[3]);
    qorder = false;
  }
  vector<SpinorWaveFunction>    f1,f2;
  vector<SpinorBarWaveFunction> a1,a2;
  SpinorWaveFunction   (f1,hard[order[0]], lorder ? incoming : outgoing, !lorder,true);
  SpinorWaveFunction   (f2,hard[order[1]], qorder ? incoming : outgoing, !qorder,true);
  SpinorBarWaveFunction(a1,hard[order[2]], lorder ? outgoing : incoming,  lorder,true);
  SpinorBarWaveFunction(a2,hard[order[3]], qorder ? outgoing : incoming,  qorder,true);
  helicityME(f1,f2,a1,a2,lorder,qorder,true);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) {
    tSpinPtr spin = hard[ix]->spinInfo();
    if(ix<2) {
      tcPolarizedBeamPDPtr beam = 
	dynamic_ptr_cast<tcPolarizedBeamPDPtr>(hard[ix]->dataPtr());
      if(beam) spin->rhoMatrix() = beam->rhoMatrix();
    }
    spin->productionVertex(hardvertex);
  }
}

double MEChargedCurrentDIS::A(tcPDPtr lin, tcPDPtr,
			      tcPDPtr qin, tcPDPtr, Energy2) const {
  double output = 2.;
  if(qin->id()<0) output *= -1.;
  if(lin->id()<0) output *= -1;
  return output;
}
