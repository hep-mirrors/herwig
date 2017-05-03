// -*- C++ -*-
//
// MEPP2GammaGamma.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2GammaGamma class.
//

#include "MEPP2GammaGamma.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

using namespace Herwig;

IBPtr MEPP2GammaGamma::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2GammaGamma::fullclone() const {
  return new_ptr(*this);
}

unsigned int MEPP2GammaGamma::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2GammaGamma::orderInAlphaEW() const {
  return 2;
}

void MEPP2GammaGamma::doinit() {
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm)
    {_photonvertex = hwsm->vertexFFP();}
  else throw InitException() << "Wrong type of StandardModel object in "
			     << "MEPP2GammaGamma::doinit() the Herwig"
			     << " version must be used" 
			     << Exception::runerror;
  // call the base class
  HwMEBase::doinit();
}

void MEPP2GammaGamma::getDiagrams() const {
  // diagrams for q qbar to gamma gamma
  tcPDPtr p = getParticleData(ParticleID::gamma);
  if(_process==0||_process==1) {
    for ( int i = 1; i <= 5; ++i ) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
      // t channel
      add(new_ptr((Tree2toNDiagram(3), q, qb, qb, 1, p, 2, p, -1)));
      // u channel
      add(new_ptr((Tree2toNDiagram(3), q, qb, qb, 2, p, 1, p, -2)));
    }
  }
  // diagrams for g g to gamma gamma (this is garbage)
  tcPDPtr g = getParticleData(ParticleID::g);
  if(_process==0||_process==2)
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, p, 3, p, 3, p, -3)));
}

Energy2 MEPP2GammaGamma::scale() const {
  Energy2 s(sHat()),u(uHat()),t(tHat());
  return scalePreFactor_*2.*s*t*u/(s*s+t*t+u*u);
}

Selector<MEBase::DiagramIndex>
MEPP2GammaGamma::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(_diagwgt[0], i);
    else if ( diags[i]->id() == -2 )  sel.insert(_diagwgt[1], i);
    else if ( diags[i]->id() == -3 )  sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *>
MEPP2GammaGamma::colourGeometries(tcDiagPtr diag) const {
  // q qbar colour lines
  static const ColourLines cqqbar("1 -2 -3");
  // g g colour lines
  static const ColourLines cgluon("1 -2,-1 2");
  // selector
  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 || diag->id() == -2 ) sel.insert(1.0, &cqqbar);
  else                                        sel.insert(1.0, &cgluon);
  return sel;
}

void MEPP2GammaGamma::persistentOutput(PersistentOStream & os) const {
  os << _photonvertex << _maxflavour << _process << scalePreFactor_;
}

void MEPP2GammaGamma::persistentInput(PersistentIStream & is, int) {
  is >> _photonvertex >> _maxflavour >> _process >> scalePreFactor_;
}

ClassDescription<MEPP2GammaGamma> MEPP2GammaGamma::initMEPP2GammaGamma;
// Definition of the static class description member.

void MEPP2GammaGamma::Init() {

  static ClassDocumentation<MEPP2GammaGamma> documentation
    ("The MEPP2GammaGamma class implements the matrix element for photon pair"
     " production in hadron collisions.");

  static Parameter<MEPP2GammaGamma,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks in the process",
     &MEPP2GammaGamma::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Switch<MEPP2GammaGamma,unsigned int> interfaceProcess
    ("Process",
     "Subprocesses to include",
     &MEPP2GammaGamma::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all the subprocesses",
     0);
  static SwitchOption interfaceProcessqqbar
    (interfaceProcess,
     "qqbar",
     "Only include the incoming q qbar subproces",
     1);
  static SwitchOption interfaceProcessgg
    (interfaceProcess,
     "gg",
     "Only include the incoming gg subprocess",
     2);

  static Parameter<MEPP2GammaGamma,double> interfaceScalePreFactor
    ("ScalePreFactor",
     "Prefactor for the scale",
     &MEPP2GammaGamma::scalePreFactor_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

}

double MEPP2GammaGamma::me2() const {
  // total matrix element
  double me(0.);
  // g g to gamma gamma
  if(mePartonData()[0]->id()==ParticleID::g) {
    VectorWaveFunction    g1in(meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction    g2in(meMomenta()[1],mePartonData()[1],incoming);
    VectorWaveFunction   p1out(meMomenta()[2],mePartonData()[2],outgoing);
    VectorWaveFunction   p2out(meMomenta()[3],mePartonData()[3],outgoing);
    vector<VectorWaveFunction> g1,g2,p1,p2;
    for(unsigned int ix=0;ix<2;++ix) {
      g1in.reset(2*ix) ;g1.push_back( g1in);
      g2in.reset(2*ix) ;g2.push_back( g2in);
      p1out.reset(2*ix);p1.push_back(p1out);
      p2out.reset(2*ix);p2.push_back(p2out);
    }
    // calculate the matrix element
    me = ggME(g1,g2,p1,p2,false);
  }
  // q qbar to gamma gamma
  else {
    unsigned int iq(1),iqb(0);
    if(mePartonData()[0]->id()>0){iq=0;iqb=1;}
    SpinorWaveFunction    qin (meMomenta()[iq ],mePartonData()[iq ],incoming);
    SpinorBarWaveFunction qbin(meMomenta()[iqb],mePartonData()[iqb],incoming);
    VectorWaveFunction   p1out(meMomenta()[ 2 ],mePartonData()[ 2 ],outgoing);
    VectorWaveFunction   p2out(meMomenta()[ 3 ],mePartonData()[ 3 ],outgoing);
    vector<SpinorWaveFunction> fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> p1,p2;
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)    ;fin.push_back( qin );
      qbin.reset(ix)   ;ain.push_back( qbin);
      p1out.reset(2*ix); p1.push_back(p1out);
      p2out.reset(2*ix); p2.push_back(p2out);
    }
    // calculate the matrix element
    me= qqbarME(fin,ain,p1,p2,false);
  }
  return me;
}

double MEPP2GammaGamma::qqbarME(vector<SpinorWaveFunction>    & fin,
				vector<SpinorBarWaveFunction> & ain,
				vector<VectorWaveFunction>    & p1,
				vector<VectorWaveFunction>    & p2,
				bool calc) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)
  // 1 incoming antifermion (vbar spinor)
  // for the outgoing       
  // 0 first  outgoing photon
  // 1 second outgoing photon
  // me to be returned
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1,PDT::Spin1);
  // wavefunction for the intermediate particles
  SpinorWaveFunction inter;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Complex diag[3];
  double me(0.),diag1(0.),diag2(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  // first diagram
	  inter = _photonvertex->evaluate(ZERO,5,fin[inhel1].particle()->CC(),
					  fin[inhel1],p1[outhel1]);
	  diag[0] = _photonvertex->evaluate(ZERO,inter,ain[inhel2],p2[outhel2]);
	  // second diagram
	  inter = _photonvertex->evaluate(ZERO,5,fin[inhel1].particle()->CC(),
					  fin[inhel1],p2[outhel2]);
	  diag[1] = _photonvertex->evaluate(ZERO,inter,ain[inhel2],p1[outhel1]);
	  // compute the running totals
	  diag[2]=diag[0]+diag[1];
	  diag1 += norm(diag[0]);
	  diag2 += norm(diag[1]);
	  me    += norm(diag[2]);
	  // matrix element
	  if(calc) newme(inhel1,inhel2,2*outhel1,2*outhel2)=diag[2];
	}
      }		
    }
  }
  // save the info on the diagrams
  if(!calc) {
    _diagwgt[0]=diag1;
    _diagwgt[1]=diag2;
  }
  // check versus analytic result
//   Energy2 s(sHat()),u(uHat()),t(tHat());
//   double test = 2./3.*sqr(4.*Constants::pi*SM().alphaEM(ZERO))*(t/u+u/t)*
//     pow(double(mePartonData()[0]->iCharge())/3.,4);
//   cerr << "testing me " << 12./me*test << endl;
  // return the answer (including colour and spin factor)
  if(calc) _me.reset(newme);
  // this is 1/3 colour average, 1/4 spin aver, 1/2 identical particles
  return me/24.;
}

double MEPP2GammaGamma::ggME(vector<VectorWaveFunction>    &,
			     vector<VectorWaveFunction>    &,
			     vector<VectorWaveFunction>    &,
			     vector<VectorWaveFunction>    &,
			     bool calc) const {
  // we probably need some basis rotation here ?????
  // get the scales
  Energy2 s(sHat()),u(uHat()),t(tHat());
  Complex me[2][2][2][2];
  double charge(11./9.);
  // ++++
  me[1][1][1][1] = charge*ggme(s,t,u);
  // +++-
  me[1][1][1][0] =-charge;
  // ++-+
  me[1][1][0][1] =-charge;
  // ++--
  me[1][1][0][0] =-charge;
  // +-++
  me[1][0][1][1] =-charge;
  // +-+-
  me[1][0][1][0] = charge*ggme(u,t,s);
  // +--+
  me[1][0][0][1] = charge*ggme(t,s,u);
  // +---
  me[1][0][0][0] = charge;
  // -+++
  me[0][1][1][1] =-charge;
  // -++-
  me[0][1][1][0] =-me[1][0][0][1];
  // -+-+
  me[0][1][0][1] =-me[1][0][1][0];
  // -+--
  me[0][1][0][0] = charge;
  // --++
  me[0][0][1][1] = charge;
  // --+-
  me[0][0][1][0] = charge;
  // ---+
  me[0][0][0][1] = charge;
  // ----
  me[0][0][0][0] =-me[1][1][1][1];
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,
				PDT::Spin1,PDT::Spin1);
  unsigned int inhel1,inhel2,outhel1,outhel2;
  double sum(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  sum+=real(     me[inhel1][inhel2][outhel1][outhel2]*
			 conj(me[inhel1][inhel2][outhel1][outhel2]));
	  // matrix element
	  if(calc) newme(2*inhel1,2*inhel2,
			 2*outhel1,2*outhel2)=me[inhel1][inhel2][outhel1][outhel2];
	}
      }		
    }
  }
  //    double pi2(sqr(pi));
  //    Energy2 s2(sqr(s)),t2(sqr(t)),u2(sqr(u));
  //    double alntu=log(t/u);
  //    double alnst=log(-s/t);
  //    double alnsu=alnst+alntu;
  //    double test=5.*4.
  //      +sqr((2.*s2+2.*(u2-t2)*alntu+(t2+u2)*(sqr(alntu)+pi2))/s2)
  //      +sqr((2.*u2+2.*(t2-s2)*alnst+(t2+s2)* sqr(alnst)     )/u2)
  //      +sqr((2.*t2+2.*(u2-s2)*alnsu+(u2+s2)* sqr(alnsu)     )/t2)
  //      +4.*pi2*(sqr((t2-s2+(t2+s2)*alnst)/u2)+sqr((u2-s2+(u2+s2)*alnsu)/t2));
  //    cerr << "testing ratio " << sum/test/sqr(charge)*2. << endl;
  // final factors
  if(calc) _me.reset(newme);
  return 0.25*sum*sqr(SM().alphaS(scale())*SM().alphaEM(ZERO));
}

void MEPP2GammaGamma::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  // identify the process and calculate matrix element
  vector<VectorWaveFunction> g1,g2,p1,p2;
  if(hard[0]->id()==ParticleID::g) {
    VectorWaveFunction   (g1,hard[order[0]],incoming,false,true,true);
    VectorWaveFunction   (g2,hard[order[1]],incoming,false,true,true);
    VectorWaveFunction   (p1,hard[order[2]],outgoing,true ,true,true);
    VectorWaveFunction   (p2,hard[order[3]],outgoing,true ,true,true);
    g1[1]=g1[2];g2[1]=g2[2];p1[1]=p1[2];p2[1]=p2[2];
    ggME(g1,g2,p1,p2,true);
  }
  else {
    if(hard[order[0]]->id()<0) swap(order[0],order[1]);
    vector<SpinorWaveFunction> q;
    vector<SpinorBarWaveFunction>  qb;
    SpinorWaveFunction   (q ,hard[order[0]],incoming,false,true);
    SpinorBarWaveFunction(qb,hard[order[1]],incoming,false,true);
    VectorWaveFunction   (p1,hard[order[2]],outgoing,true ,true,true);
    VectorWaveFunction   (p2,hard[order[3]],outgoing,true ,true,true);
    p1[1]=p1[2];p2[1]=p2[2];
    qqbarME(q,qb,p1,p2,true);
  }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix)
    hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
}
