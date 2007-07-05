// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2gZ2ll class.
//

#include "MEee2gZ2ll.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "HardVertex.h"

using namespace Herwig;

void MEee2gZ2ll::getDiagrams() const {
  // specific the diagrams
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  tcPDPtr em = getParticleData(ParticleID::eminus);
  // setup the processes
  for( int i =11;i<=16;++i) {
    if(_allowed==0 || (_allowed==1 && i%2==1) || (_allowed==2&&i==11)
       || (_allowed==3&&i==13) || (_allowed==4&&i==15)) {
      tcPDPtr lm = getParticleData(i);
      tcPDPtr lp = lm->CC();
      add(new_ptr((Tree2toNDiagram(2), em, ep, 1, _gamma, 3, lm, 3, lp, -1)));
      add(new_ptr((Tree2toNDiagram(2), em, ep, 1, _Z0, 3, lm, 3, lp, -2)));
    }
  }
}

Energy2 MEee2gZ2ll::scale() const {
  return sHat();
}

unsigned int MEee2gZ2ll::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2gZ2ll::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEee2gZ2ll::diagrams(const DiagramVector & diags) const {
  double lastCont(0.5),lastBW(0.5);
  if ( lastXCombPtr() ) {
    lastCont = meInfo()[0];
    lastBW = meInfo()[1];
  }
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(lastCont, i);
    else if ( diags[i]->id() == -2 ) sel.insert(lastBW, i);
  }
  return sel;
}

Selector<const ColourLines *>
MEee2gZ2ll::colourGeometries(tcDiagPtr) const {
  static ColourLines ctST(" ");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &ctST);
  return sel;
}

void MEee2gZ2ll::persistentOutput(PersistentOStream & os) const {
  os << _allowed << _theFFZVertex << _theFFPVertex << _gamma << _Z0; 
}

void MEee2gZ2ll::persistentInput(PersistentIStream & is, int) {
  is >> _allowed >> _theFFZVertex >> _theFFPVertex >> _gamma >> _Z0; 
}

ClassDescription<MEee2gZ2ll> MEee2gZ2ll::initMEee2gZ2ll;
// Definition of the static class description member.

void MEee2gZ2ll::Init() {

  static Switch<MEee2gZ2ll,int> interfaceallowed
    ("Allowed",
     "Allowed outgoing leptons",
     &MEee2gZ2ll::_allowed, 0, false, false);
  static SwitchOption interfaceallowedAll
    (interfaceallowed,
     "All",
     "Allow all leptons as outgoing particles",
     0);
  static SwitchOption interfaceallowedCharged
    (interfaceallowed,
     "Charged",
     "Only charged leptons as outgoing particles",
     1);
  static SwitchOption interfaceallowedElectron 
    (interfaceallowed, 
     "Electron",
     "Only the electron and positron as outgoing leptons",
     2);
  static SwitchOption interfaceallowedMuon 
    (interfaceallowed, 
     "Muon", 
     "Only muons as outgoing particles",
     3);
  static SwitchOption interfaceallowedTau
    (interfaceallowed,
     "Tau",
     "Only taus as outgoing particles",
     4);

  static ClassDocumentation<MEee2gZ2ll> documentation
    ("The MEee2gZ2ll class implements the matrix element for"
     "e+e- to leptons via Z and photon exchange using helicity amplitude"
     "techniques");
}

double MEee2gZ2ll::me2() const {
  int ie(0),ip(1),ilp(2),ilm(3);
  // get the order right
  if(mePartonData()[0]->id()!=11) swap(ie,ip);
  if(mePartonData()[2]->id()<0)   swap(ilm,ilp);
  vector<SpinorWaveFunction> fin,aout;
  vector<SpinorBarWaveFunction>  ain,fout;
  SpinorWaveFunction    ein  (meMomenta()[ie ],mePartonData()[ie ],incoming);
  SpinorBarWaveFunction pin  (meMomenta()[ip ],mePartonData()[ip ],incoming);
  SpinorBarWaveFunction ilmout(meMomenta()[ilm],mePartonData()[ilm],outgoing);
  SpinorWaveFunction    ilpout(meMomenta()[ilp],mePartonData()[ilp],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    ein.reset(ix)  ;fin.push_back( ein  );
    pin.reset(ix)  ;ain.push_back( pin  );
    ilmout.reset(ix);fout.push_back(ilmout);
    ilpout.reset(ix);aout.push_back(ilpout);
  }
  // compute the matrix element
  double me,lastCont,lastBW;
  HelicityME(fin,ain,fout,aout,me,lastCont,lastBW);
  // save the components
  DVector save;
  save.push_back(lastCont);
  save.push_back(lastBW);
  meInfo(save);
  // return the answer
  return me;
}

ProductionMatrixElement MEee2gZ2ll::HelicityME(vector<SpinorWaveFunction>    & fin,
					       vector<SpinorBarWaveFunction> & ain,
					       vector<SpinorBarWaveFunction> & fout,
					       vector<SpinorWaveFunction>    & aout,
					       double & me,double & cont,
					       double & BW ) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)
  // 1 incoming antifermion (vbar spinor)
  // for the outgoing       
  // 0 outgoing fermion     (ubar spinor)
  // 1 outgoing antifermion (v    spinor)
  // me to be returned
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  //   // wavefunctions for the intermediate particles
  VectorWaveFunction interZ,interG;
  // temporary storage of the different diagrams
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1,outhel2;
  double total[3]={0.,0.};
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      // intermediate Z
      interZ = _theFFZVertex->evaluate(sHat(),1,_Z0,fin[inhel1],ain[inhel2]);
      // intermediate photon
      interG = _theFFPVertex->evaluate(sHat(),1,_gamma,fin[inhel1],ain[inhel2]);
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {		
	  // first the Z exchange diagram
	  diag1 = _theFFZVertex->evaluate(sHat(),aout[outhel2],fout[outhel1],
					  interZ);
	  // then the photon exchange diagram
	  diag2 = _theFFPVertex->evaluate(sHat(),aout[outhel2],fout[outhel1],
					  interG);
	  // add up squares of individual terms
	  total[1] += real(diag1*conj(diag1));
	  total[2] += real(diag2*conj(diag2));
	  diag1+=diag2;
	  // the full thing including interference
	  diag1 +=diag2;
	  total[0] += real(diag1*conj(diag1));
	  output(inhel1,inhel2,outhel1,outhel2)=diag1;
	}
      }
    }
  }
  // results
  for(int ix=0;ix<3;++ix){total[ix]*=0.25;}
  cont = total[2];
  BW   = total[1];
  me   = total[0];
  return output;
}

void MEee2gZ2ll::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);hard.push_back(sub->outgoing()[1]);
  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  if(hard[2]->id()<hard[3]->id()) swap(hard[2],hard[3]);
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
  SpinorBarWaveFunction(fout,hard[2],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[3],outgoing,true ,true);
  // calculate the matrix element
  double me,cont,BW;
  ProductionMatrixElement prodme=HelicityME(fin,ain,fout,aout,me,cont,BW);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(prodme);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) {
    dynamic_ptr_cast<SpinfoPtr>(hard[ix]->spinInfo())->setProductionVertex(hardvertex);
  }
}
