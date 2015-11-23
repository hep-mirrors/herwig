// -*- C++ -*-
//
// MEPP2GammaJet.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2GammaJet class.
//

#include "MEPP2GammaJet.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

using namespace Herwig;

MEPP2GammaJet::MEPP2GammaJet() : _maxflavour(5), _processopt(0), scalePreFactor_(1.) {
  massOption(vector<unsigned int>(2,0));
}

void MEPP2GammaJet::rebind(const TranslationMap & trans)
  {
  // dummy = trans.translate(dummy);
  HwMEBase::rebind(trans);
  _gluonvertex =trans.translate(_gluonvertex );
  _photonvertex=trans.translate(_photonvertex);
}

IVector MEPP2GammaJet::getReferences() {
  IVector ret = HwMEBase::getReferences();
  ret.push_back(_gluonvertex);
  ret.push_back(_photonvertex);
  return ret;
}

void MEPP2GammaJet::doinit() {
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _gluonvertex  = hwsm->vertexFFG();
    _photonvertex = hwsm->vertexFFP();
  }
  else throw InitException() << "Wrong type of StandardModel object in "
			     << "MEPP2GammaJet::doinit() the Herwig"
			     << " version must be used" 
			     << Exception::runerror;
  // call the base class
  HwMEBase::doinit();
}

void MEPP2GammaJet::getDiagrams() const {
  // need the gluon and the photon in all processes
  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr p = getParticleData(ParticleID::gamma);
  // for each quark species there are three subprocesses
  for ( int iq=1; iq<=_maxflavour; ++iq ) {
    tcPDPtr q = getParticleData(iq);
    tcPDPtr qb = q->CC();
    // q qbar to gamma gluon (two diagrams)
    if(_processopt==0||_processopt==1) {
      add(new_ptr((Tree2toNDiagram(3), q, qb, qb, 1, p, 2, g, -1)));
      add(new_ptr((Tree2toNDiagram(3), q,  q, qb, 2, p, 1, g, -2)));
    }
    // q gluon to gamma q (two diagrams)
    if(_processopt==0||_processopt==2) {
      add(new_ptr((Tree2toNDiagram(3), q, q, g, 1, p, 2, q, -3)));
      add(new_ptr((Tree2toNDiagram(2), q, g, 1, q , 3, p, 3, q, -4)));
    }
    // qbar gluon to gamma qbar (two diagrams)
    if(_processopt==0||_processopt==3) {
      add(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, p, 2, qb, -5)));
      add(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb , 3, p, 3, qb, -6)));
    }
  }
}

unsigned int MEPP2GammaJet::orderInAlphaS() const {
  return 1;
}

unsigned int MEPP2GammaJet::orderInAlphaEW() const {
  return 1;
}

Energy2 MEPP2GammaJet::scale() const {
  Energy2 s(sHat()),u(uHat()),t(tHat());
  return scalePreFactor_*2.*s*t*u/(s*s+t*t+u*u);
}

Selector<MEBase::DiagramIndex>
MEPP2GammaJet::diagrams(const DiagramVector & diags) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.
  double diag1(0.5),diag2(0.5);
  diag1 = meInfo()[0];
  diag2 = meInfo()[1];
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if      ( abs(diags[i]->id())%2 == 1 ) sel.insert(diag1, i);
    else                                   sel.insert(diag2, i);
  }
  return sel;
}

void MEPP2GammaJet::persistentOutput(PersistentOStream & os) const {
  os << _gluonvertex << _photonvertex << _maxflavour << _processopt << scalePreFactor_;
}

void MEPP2GammaJet::persistentInput(PersistentIStream & is, int) {
  is >> _gluonvertex >> _photonvertex >> _maxflavour >> _processopt >> scalePreFactor_;
}

ClassDescription<MEPP2GammaJet> MEPP2GammaJet::initMEPP2GammaJet;
// Definition of the static class description member.

void MEPP2GammaJet::Init() {

  static ClassDocumentation<MEPP2GammaJet> documentation
    ("The MEPP2GammaJet class implements the matrix element for"
     " hadron-hadron to photon+jet");

  static Parameter<MEPP2GammaJet,int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks in the process",
     &MEPP2GammaJet::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Switch<MEPP2GammaJet,unsigned int> interfaceProcesses
    ("Process",
     "Subprocesses to include",
     &MEPP2GammaJet::_processopt, 0, false, false);
  static SwitchOption interfaceProcessesAll
    (interfaceProcesses,
     "All",
     "Include all the subprocesses",
     0);
  static SwitchOption interfaceProcessesqqbar
    (interfaceProcesses,
     "qqbar",
     "Only include the incoming q qbar subprocess",
     1);
  static SwitchOption interfaceProcessesqg
    (interfaceProcesses,
     "qg",
     "Only include the incoming q g subprocess",
     2);
  static SwitchOption interfaceProcessesqbarg
    (interfaceProcesses,
     "qbarg",
     "Only include the incoming qbar g subprocess",
     3);

  static Parameter<MEPP2GammaJet,double> interfaceScalePreFactor
    ("ScalePreFactor",
     "Prefactor for the scale",
     &MEPP2GammaJet::scalePreFactor_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
}

Selector<const ColourLines *>
MEPP2GammaJet::colourGeometries(tcDiagPtr diag) const {
  // q qbar to gamma gluon colour lines
  static const ColourLines qqbar1("1 5, -5 -2 -3");
  static const ColourLines qqbar2("1 2 5, -5 -3");
  // q gluon to gamma q colour lines
  static const ColourLines qg1("1 2 -3, 3 5");
  static const ColourLines qg2("1 -2, 2 3 5");
  // qbar gluon to gamma qbar lines
  static const ColourLines qbarg1("-1 -2 3, -3 -5");
  static const ColourLines qbarg2("-1 2, -2 -3 -5");
  // only one flow per diagram so insert the right one
  Selector<const ColourLines *> sel;  
  switch (diag->id()) {
  case -1 :
    sel.insert(1.0, &qqbar1);
    break;
  case -2 :
    sel.insert(1.0, &qqbar2);
    break;
  case -3 :
    sel.insert(1.0, &qg1);
    break;
  case -4 :
    sel.insert(1.0, &qg2);
    break;
  case -5 :
    sel.insert(1.0, &qbarg1);
    break;
  case -6 :
    sel.insert(1.0, &qbarg2);
    break;
  }
  return sel;
}

double MEPP2GammaJet::me2() const {
  // total matrix element and the various components
  double me(0.);
  // first case, q qbar to gluon photon
  if(mePartonData()[0]->id()==-mePartonData()[1]->id()) {
    // order of the particles
    unsigned int iq(1),iqb(0),ip(3),ig(2);
    if(mePartonData()[0]->id()>0)              swap(iq,iqb);
    if(mePartonData()[3]->id()==ParticleID::g) swap(ig, ip);
    // calculate the spinors and polarization vectors
    vector<SpinorWaveFunction> fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> pout,gout;
    SpinorWaveFunction    qin (meMomenta()[iq ],mePartonData()[iq ],incoming);
    SpinorBarWaveFunction qbin(meMomenta()[iqb],mePartonData()[iqb],incoming);
    VectorWaveFunction   glout(meMomenta()[ig ],mePartonData()[ig ],outgoing);
    VectorWaveFunction   phout(meMomenta()[ip ],mePartonData()[ip ],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)    ; fin.push_back( qin );
      qbin.reset(ix)   ; ain.push_back( qbin);
      glout.reset(2*ix);gout.push_back(glout);
      phout.reset(2*ix);pout.push_back(phout);
    }
    // calculate the matrix element
    me = qqbarME(fin,ain,gout,pout,false)/9.;
//       Energy2 mt(scale());
//     double coupling=sqr(4.*Constants::pi)*SM().alphaEM(ZERO)*SM().alphaS(mt)*
//       sqr(mePartonData()[0]->iCharge()/3.);
//       Energy2 t(tHat()),u(uHat());
//       double me2=8./9./u/t*(t*t+u*u)*coupling;
//       cerr << "testing matrix element A" 
//  	   << me << "  " 
//  	   << me2 << " " << me/me2
//  	   << endl;
  }
  else if(mePartonData()[0]->id()>0&&mePartonData()[1]->id()) {
    // order of the particles
    unsigned int iqin(0),iqout(2),ip(3),ig(1);
    if(mePartonData()[0]->id()==ParticleID::g    ) swap(iqin,ig);
    if(mePartonData()[2]->id()==ParticleID::gamma) swap(ip,iqout);
    // calculate the spinors and polarization vectors
    vector<SpinorWaveFunction> fin;
    vector<SpinorBarWaveFunction>  fout;
    vector<VectorWaveFunction> pout,gin;
    SpinorWaveFunction     qin (meMomenta()[iqin ],mePartonData()[iqin ],incoming);
    SpinorBarWaveFunction  qout(meMomenta()[iqout],mePartonData()[iqout],outgoing);
    VectorWaveFunction     glin(meMomenta()[ig   ],mePartonData()[ig   ],incoming);  
    VectorWaveFunction    phout(meMomenta()[ip   ],mePartonData()[ip   ],outgoing); 
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)    ;fin.push_back(  qin );
      qout.reset(ix)   ;fout.push_back( qout);
      glin.reset(2*ix) ;gin.push_back(  glin);
      phout.reset(2*ix);pout.push_back(phout);
    }
    // calculate the matrix element
    me = qgME(fin,gin,pout,fout,false)/24.;
//       Energy2 mt(scale());
//     double coupling=sqr(4.*Constants::pi)*SM().alphaEM(ZERO)*SM().alphaS(mt);
//       Energy2 s(sHat()),t(tHat()),u(uHat());
//       double me2=-1./3./s/t*(s*s+t*t+2.*u*(s+t+u))*coupling*
//       sqr(mePartonData()[0]->iCharge()/3.);
//       cerr << "testing matrix element B" 
// 	   << me << "  " 
// 	   << me2 << " " << me/me2
// 	   << endl; 
  }
  else {
    // order of the particles
    unsigned int iqin(0),iqout(2),ip(3),ig(1);
    if(mePartonData()[0]->id()==ParticleID::g    ) swap(iqin,ig);
    if(mePartonData()[2]->id()==ParticleID::gamma) swap(ip,iqout);
    // calculate the spinors and polarization vectors
    vector<SpinorBarWaveFunction> ain;
    vector<SpinorWaveFunction>  aout;
    vector<VectorWaveFunction> pout,gin;
    SpinorBarWaveFunction  qin (meMomenta()[iqin ],mePartonData()[iqin ],incoming);
    SpinorWaveFunction     qout(meMomenta()[iqout],mePartonData()[iqout],outgoing);
    VectorWaveFunction     glin(meMomenta()[ig   ],mePartonData()[ig   ],incoming);  
    VectorWaveFunction    phout(meMomenta()[ip   ],mePartonData()[ip   ],outgoing); 
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)    ;ain.push_back(  qin );
      qout.reset(ix)   ;aout.push_back( qout);
      glin.reset(2*ix) ;gin.push_back(  glin);
      phout.reset(2*ix);pout.push_back(phout);
    }
    // calculate the matrix element
    me=qbargME(ain,gin,pout,aout,false)/24.;
//       Energy2 mt(scale());
//     double coupling=sqr(4.*Constants::pi)*SM().alphaEM(ZERO)*SM().alphaS(mt);
//       Energy2 s(sHat()),t(tHat()),u(uHat());
//       double me2=-1./3./s/t*(s*s+t*t+2.*u*(s+t+u))*coupling*
//       sqr(mePartonData()[0]->iCharge()/3.);
//       cerr << "testing matrix element C" 
// 	   << me << "  " 
// 	   << me2 << " " << me/me2
// 	   << endl; 
  }
  return me;
}

double MEPP2GammaJet::qqbarME(vector<SpinorWaveFunction>    & fin,
			      vector<SpinorBarWaveFunction> & ain,
			      vector<VectorWaveFunction>    & gout,
			      vector<VectorWaveFunction>    & pout,
			      bool calc) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)
  // 1 incoming antifermion (vbar spinor)
  // for the outgoing       
  // 0 outgoing gluon
  // 1 outgoing photon
  // me to be returned
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1,PDT::Spin1);
  // wavefunction for the intermediate particles
  SpinorWaveFunction inter;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Energy2 mt(scale());
  Complex diag[3];
  double me(0.),diag1(0.),diag2(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  // first diagram
	  inter = _gluonvertex->evaluate(mt,5,fin[inhel1].particle()->CC(),
					 fin[inhel1],gout[outhel1]);
	  diag[0] = _photonvertex->evaluate(ZERO,inter,ain[inhel2],pout[outhel2]);
	  // second diagram
	  inter = _photonvertex->evaluate(ZERO,5,fin[inhel1].particle()->CC(),
					  fin[inhel1],pout[outhel2]);
	  diag[1] = _gluonvertex->evaluate(mt,inter,ain[inhel2],gout[outhel1]);
	  // compute the running totals
	  diag[2]=diag[0]+diag[1];
	  diag1 +=norm(diag[0]);
	  diag2 +=norm(diag[1]);
	  me    +=norm(diag[2]);
	  // matrix element
	  if(calc) newme(inhel1,inhel2,2*outhel1,2*outhel2)=diag[2];
	}
      }		
    }
  }
  // save the info on the diagrams
  if(!calc) {
    DVector save;
    save.push_back(diag1);
    save.push_back(diag2);
    meInfo(save);
  }
  // return the answer
  if(calc) _me.reset(newme);
  return me;
}

double MEPP2GammaJet::qgME(vector<SpinorWaveFunction>    & fin,
			   vector<VectorWaveFunction>    & gin,
			   vector<VectorWaveFunction>    & pout,
			   vector<SpinorBarWaveFunction> & fout,
			   bool calc) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)       
  // 1 incoming gluon
  // for the outgoing
  // 0 outgoing photon
  // 1 outgoing fermion     (ubar spinor)
  // me to be returned
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1,
				PDT::Spin1,PDT::Spin1Half);
  // wavefunction for the intermediate particles
  SpinorWaveFunction inter;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Energy2 mt(scale());
  Complex diag[3];
  double me(0.),diag1(0.),diag2(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  // first diagram
	  inter = _photonvertex->evaluate(ZERO,5,fin[inhel1].particle()->CC(),
					  fin[inhel1],pout[outhel1]);
	  diag[0]=_gluonvertex->evaluate(mt,inter,fout[outhel2],gin[inhel2]);
	  // second diagram
	  inter = _gluonvertex->evaluate(mt,5,fin[inhel1].particle()->CC(),
					 fin[inhel1],gin[inhel2]);
	  diag[1]=_photonvertex->evaluate(ZERO,inter,fout[outhel2],pout[outhel1]);
	  // compute the running totals
	  diag[2]=diag[0]+diag[1];
	  diag1 +=norm(diag[0]);
	  diag2 +=norm(diag[1]);
	  me    +=norm(diag[2]);
	  // matrix element
	  if(calc) newme(inhel1,2*inhel2,2*outhel1,outhel2)=diag[2];
	}
      }		
    }
  }
  // save the info on the diagrams
  if(!calc) {
    DVector save;
    save.push_back(diag1);
    save.push_back(diag2);
    meInfo(save);
  }
  // return the answer
  if(calc) _me.reset(newme);
  return me;
} 

double MEPP2GammaJet::qbargME(vector<SpinorBarWaveFunction> & ain,
			      vector<VectorWaveFunction>    & gin,
			      vector<VectorWaveFunction>    & pout,
			      vector<SpinorWaveFunction>    & aout,
			      bool calc) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (vbar spinor)       
  // 1 incoming gluon
  // for the outgoing
  // 0 outgoing photon
  // 1 outgoing fermion     (v    spinor)
  //me to be returned
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1,
				PDT::Spin1,PDT::Spin1Half);
  // wavefunction for the intermediate particles
  SpinorBarWaveFunction inter;
  SpinorWaveFunction interb;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Energy2 mt(scale());
  Complex diag[3];
  double me(0.),diag1(0.),diag2(0.);
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {
	  // first diagram
	  inter = _photonvertex->evaluate(ZERO,5,ain[inhel1].particle()->CC(),
					  ain[inhel1],pout[outhel1]);
	  diag[0]=_gluonvertex->evaluate(mt,aout[outhel2],inter,gin[inhel2]);
	  // second diagram
	  inter = _gluonvertex->evaluate(mt,5,ain[inhel1].particle()->CC(),
					 ain[inhel1],gin[inhel2]);
	  diag[1]=_photonvertex->evaluate(ZERO,aout[outhel2],inter,pout[outhel1]);
	  // compute the running totals
	  diag[2]=diag[0]+diag[1];
	  diag1 +=norm(diag[0]);
	  diag2 +=norm(diag[1]);
	  me    +=norm(diag[2]);
	  // matrix element
	  if(calc) newme(inhel1,2*inhel2,2*outhel1,outhel2)=diag[2];
	}
      }		
    }
  }
  // save the info on the diagrams
  if(!calc) {
    DVector save;
    save.push_back(diag1);
    save.push_back(diag2);
    meInfo(save);
  }
  // return the answer
  if(calc) _me.reset(newme);
  return me;
}

void MEPP2GammaJet::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  // identify the process and calculate matrix element
  if(hard[0]->id()==ParticleID::g||hard[1]->id()==ParticleID::g) {
    if(hard[0]->id()==ParticleID::g    ) swap(order[0],order[1]);
    if(hard[3]->id()==ParticleID::gamma) swap(order[2],order[3]);
    if(hard[order[0]]->id()>0) {
      vector<SpinorWaveFunction> q;
      vector<SpinorBarWaveFunction>  qb;
      vector<VectorWaveFunction> p,g;
      SpinorWaveFunction   (q ,hard[order[0]],incoming,false,true);
      VectorWaveFunction   (g ,hard[order[1]],incoming,false,true,true);
      VectorWaveFunction   (p ,hard[order[2]],outgoing,true ,true,true);
      SpinorBarWaveFunction(qb,hard[order[3]],outgoing,true,true);
      qgME(q,g,p,qb,true);
    }
    else {
      vector<SpinorWaveFunction> q;
      vector<SpinorBarWaveFunction>  qb;
      vector<VectorWaveFunction> p,g;
      SpinorBarWaveFunction(qb,hard[order[0]],incoming,false,true);
      VectorWaveFunction   (g ,hard[order[1]],incoming,false,true,true);
      VectorWaveFunction   (p ,hard[order[2]],outgoing,true ,true,true);
      SpinorWaveFunction   (q ,hard[order[3]],outgoing,true,true);
      qbargME(qb,g,p,q,true);
    }
  }
  else {
    if(hard[0]->id()<0)                  swap(order[0],order[1]);
    if(hard[2]->id()==ParticleID::gamma) swap(order[2],order[3]);
    vector<SpinorWaveFunction> q;
    vector<SpinorBarWaveFunction>  qb;
    vector<VectorWaveFunction> p,g;
    SpinorWaveFunction   (q ,hard[order[0]],incoming,false,true);
    SpinorBarWaveFunction(qb,hard[order[1]],incoming,false,true);
    VectorWaveFunction   (g ,hard[order[2]],outgoing,true ,true,true);
    VectorWaveFunction   (p ,hard[order[3]],outgoing,true ,true,true);
    p[1]=p[2];g[1]=g[2];
    qqbarME(q,qb,g,p,true);
  }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) 
    hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
}
