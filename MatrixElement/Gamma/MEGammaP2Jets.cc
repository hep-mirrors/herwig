// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGammaP2Jets class.
//

#include "MEGammaP2Jets.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

using namespace Herwig;

MEGammaP2Jets::MEGammaP2Jets() : _process(0), _minflavour(1) ,_maxflavour(5) {
  massOption(vector<unsigned int>(2,0));
}

unsigned int MEGammaP2Jets::orderInAlphaS() const {
  return 1;
}

unsigned int MEGammaP2Jets::orderInAlphaEW() const {
  return 1;
}

Energy2 MEGammaP2Jets::scale() const {
  return 2.*sHat()*tHat()*uHat()/
    (sqr(sHat())+sqr(tHat())+sqr(uHat()));
}

void MEGammaP2Jets::persistentOutput(PersistentOStream & os) const {
  os << _gluonvertex << _photonvertex 
     << _process << _minflavour << _maxflavour;
}

void MEGammaP2Jets::persistentInput(PersistentIStream & is, int) {
  is >> _gluonvertex >> _photonvertex 
     >> _process >> _minflavour >> _maxflavour;
}

ClassDescription<MEGammaP2Jets> MEGammaP2Jets::initMEGammaP2Jets;
// Definition of the static class description member.

void MEGammaP2Jets::Init() {

  static ClassDocumentation<MEGammaP2Jets> documentation
    ("The MEGammaP2Jets class implements the matrix elements"
     " for pointlike photon-hadron to jets.");

  static Parameter<MEGammaP2Jets,int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The minimum flavour of the quarks",
     &MEGammaP2Jets::_minflavour, 1, 1, 5,
     false, false, Interface::limited);

  static Parameter<MEGammaP2Jets,int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks",
     &MEGammaP2Jets::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Switch<MEGammaP2Jets,unsigned int> interfaceProcess
    ("Process",
     "The allowed partonic subprocesses",
     &MEGammaP2Jets::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all the subprocesses",
     0);
  static SwitchOption interfaceProcessGluon
    (interfaceProcess,
     "Gluon",
     "Only include the gamma g -> q qbar processes",
     1);
  static SwitchOption interfaceProcessQuark
    (interfaceProcess,
     "Quark",
     "Only include the gamma q -> gluon q processes",
     2);
  static SwitchOption interfaceProcessAntiQuark
    (interfaceProcess,
     "AntiQuark",
     "Only include the gamma qbar -> gluon qbar processes",
     3);

}

void MEGammaP2Jets::getDiagrams() const {
  tcPDPtr g  = getParticleData(ParticleID::g);
  tcPDPtr gm = getParticleData(ParticleID::gamma);
  for ( int i = _minflavour; i <= _maxflavour; ++i ) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();
    // gamma g -> q qbar
    if(_process==0||_process==1) {
      add(new_ptr((Tree2toNDiagram(3), gm, q , g, 1, q, 2, qb, -1)));
      add(new_ptr((Tree2toNDiagram(3), gm, qb, g, 2, q, 1, qb, -2)));      
    }
    // gamma q -> g q
    if(_process==0||_process==2) {
      add(new_ptr((Tree2toNDiagram(3), gm, q , q, 2, g, 1, q, -3)));
      add(new_ptr((Tree2toNDiagram(2), gm, q , 1, q, 3, g, 3, q, -4)));
    }
    // gamma qbar -> g qbar
    if(_process==0||_process==3) {
      add(new_ptr((Tree2toNDiagram(3), gm, qb, qb, 2, g, 1, qb, -5)));
      add(new_ptr((Tree2toNDiagram(2), gm, qb, 1, qb, 3, g, 3, qb, -6)));
    }
  }
}

Selector<MEBase::DiagramIndex>
MEGammaP2Jets::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if      ( diags[i]->id() == -1 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
    else if ( diags[i]->id() == -3 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -4 ) sel.insert(meInfo()[1], i);
    else if ( diags[i]->id() == -5 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -6 ) sel.insert(meInfo()[1], i);
  }
  return sel;
}

IBPtr MEGammaP2Jets::clone() const {
  return new_ptr(*this);
}

IBPtr MEGammaP2Jets::fullclone() const {
  return new_ptr(*this);
}

void MEGammaP2Jets::doinit() {
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _gluonvertex  = hwsm->vertexFFG();
    _photonvertex = hwsm->vertexFFP();
  }
  else throw InitException() << "Wrong type of StandardModel object in "
			     << "MEGammaP2Jets::doinit() the Herwig"
			     << " version must be used" 
			     << Exception::runerror;
  // call the base class
  HwMEBase::doinit();
}

Selector<const ColourLines *>
MEGammaP2Jets::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c1("3 2 4,-3 -5");
  static ColourLines c2("3 4,-3 -2 -5");
  static ColourLines c3("3 4,-4 2 5");
  static ColourLines c4("2 3 4,-4 5");
  static ColourLines c5("-3 -4,4 -2 -5");
  static ColourLines c6("-2 -3 -4,4 -5");
  Selector<const ColourLines *> sel;
  if      ( diag->id() == -1 ) sel.insert(1.0, &c1);
  else if ( diag->id() == -2 ) sel.insert(1.0, &c2);
  else if ( diag->id() == -3 ) sel.insert(1.0, &c3);
  else if ( diag->id() == -4 ) sel.insert(1.0, &c4);
  else if ( diag->id() == -5 ) sel.insert(1.0, &c5);
  else if ( diag->id() == -6 ) sel.insert(1.0, &c6);
  return sel;
}

double MEGammaP2Jets::me2() const {
  // total matrix element and the various components
  double me(0.);
  // gamma g -> q qbar
  if(mePartonData()[1]->id()==ParticleID::g) {
    VectorWaveFunction    gmin (meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction    glin (meMomenta()[1],mePartonData()[1],incoming);
    SpinorBarWaveFunction qout (meMomenta()[2],mePartonData()[2],outgoing);
    SpinorWaveFunction    qbout(meMomenta()[3],mePartonData()[3],outgoing);
    vector<VectorWaveFunction> v1,v2;
    vector<SpinorBarWaveFunction> a3;
    vector<SpinorWaveFunction> f4;
    for(unsigned int ix=0;ix<2;++ix) {
      gmin.reset(2*ix);
      v1.push_back(gmin);
      glin.reset(2*ix);
      v2.push_back(glin);
      qout.reset(ix);
      a3.push_back(qout);
      qbout.reset(ix);
      f4.push_back(qbout);
    }
    // calculate the matrix element
    me = gammagluonME(v1,v2,a3,f4,false);
  }
  // gamma q -> g q
  else if(mePartonData()[1]->id()>0) {
    VectorWaveFunction    gmin (meMomenta()[0],mePartonData()[0],incoming);
    SpinorWaveFunction    qin  (meMomenta()[1],mePartonData()[1],incoming);
    VectorWaveFunction    glout(meMomenta()[2],mePartonData()[2],outgoing);
    SpinorBarWaveFunction qout (meMomenta()[3],mePartonData()[3],outgoing);
    vector<VectorWaveFunction> v1,v3;
    vector<SpinorWaveFunction> f2;
    vector<SpinorBarWaveFunction> f4;
    for(unsigned int ix=0;ix<2;++ix) {
      gmin.reset(2*ix);
      v1.push_back(gmin);
      qin.reset(ix);
      f2.push_back(qin);
      glout.reset(2*ix);
      v3.push_back(glout);
      qout.reset(ix);
      f4.push_back(qout);
    }
    me = gammaquarkME(v1,f2,v3,f4,false);
  }
  // gamma qbar -> g qbar
  else {
    VectorWaveFunction    gmin (meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction qin  (meMomenta()[1],mePartonData()[1],incoming);
    VectorWaveFunction    glout(meMomenta()[2],mePartonData()[2],outgoing);
    SpinorWaveFunction    qout (meMomenta()[3],mePartonData()[3],outgoing);
    vector<VectorWaveFunction> v1,v3;
    vector<SpinorBarWaveFunction> a2;
    vector<SpinorWaveFunction> a4;
    for(unsigned int ix=0;ix<2;++ix) {
      gmin.reset(2*ix);
      v1.push_back(gmin);
      qin.reset(ix);
      a2.push_back(qin);
      glout.reset(2*ix);
      v3.push_back(glout);
      qout.reset(ix);
      a4.push_back(qout);
    }
    me = gammaantiquarkME(v1,a2,v3,a4,false);
  }
  return me;
}

double MEGammaP2Jets::gammagluonME(vector<VectorWaveFunction> & gmin,
				   vector<VectorWaveFunction> & glin,
				   vector<SpinorBarWaveFunction> & fout, 
				   vector<SpinorWaveFunction> & aout,
				   bool calc) const {
  // the scale
  Energy2 mt(scale());
  // storage of the helicity me if needed
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,
				PDT::Spin1Half,PDT::Spin1Half);
  // loop over the helicities
  SpinorWaveFunction inter;
  vector<double> me(3,0.);
  vector<Complex> diag(2);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // first diagram
	  inter = _gluonvertex->evaluate(mt,5,aout[ohel2].particle(),
					 aout[ohel2],glin[ihel2]);
	  diag[0] = _photonvertex->evaluate(0.*GeV2,inter,fout[ohel1],
					    gmin[ihel1]);
	  // second diagram
	  inter = _photonvertex->evaluate(0.*GeV2,5,aout[ohel2].particle(),
					  aout[ohel2],gmin[ihel1]);
	  diag[1] = _gluonvertex->evaluate(mt,inter,fout[ohel1],glin[ihel2]);
	  for(unsigned int ix=0;ix<2;++ix) me[ix] += norm(diag[ix]);
	  diag[0]+=diag[1];
	  me[2] += norm(diag[0]);
	  // matrix element if needed
	  if(calc) newme(2*ihel1,2*ihel2,ohel1,ohel2)=diag[0];
	}
      }
    }
  }
  // save the info on the diagrams
  if(!calc) {
    DVector save;
    save.push_back(me[0]);
    save.push_back(me[1]);
    meInfo(save);
  }
  // return the answer with colour and spin factors
  if(calc) _me.reset(newme);
  return 0.125*me[2];
}

double MEGammaP2Jets::gammaquarkME(vector<VectorWaveFunction> & gmin,
				   vector<SpinorWaveFunction> & fin,
				   vector<VectorWaveFunction> & gout,
				   vector<SpinorBarWaveFunction> & fout,
				   bool calc) const {
  // the scale
  Energy2 mt(scale());
  // storage of the helicity me if needed
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1Half,
				PDT::Spin1,PDT::Spin1Half);
  // loop over the helicities
  SpinorWaveFunction inter;
  vector<double> me(3,0.);
  vector<Complex> diag(2);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // first diagram
	  inter = _gluonvertex->evaluate(mt,5,fin[ihel2].particle(),
					 fin[ihel2],gout[ohel1]);
	  diag[0] = _photonvertex->evaluate(0.*GeV2,inter,fout[ohel2],gmin[ihel1]);
	  // second diagram
	  inter = _photonvertex->evaluate(0.*GeV2,5,fin[ihel2].particle(),
					  fin[ihel2],gmin[ihel1]);
	  diag[1] = _gluonvertex->evaluate(mt,inter,fout[ohel2],gout[ohel1]);
	  for(unsigned int ix=0;ix<2;++ix) me[ix] += norm(diag[ix]);
	  diag[0]+=diag[1];
	  me[2] += norm(diag[0]);
	  // matrix element if needed
	  if(calc) newme(2*ihel1,ihel2,2*ohel1,ohel2)=diag[0];
	}
      }
    }
  }
  // save the info on the diagrams
  if(!calc) {
    DVector save;
    save.push_back(me[0]);
    save.push_back(me[1]);
    meInfo(save);
  }
  // return the answer with colour and spin factors
  if(calc) _me.reset(newme);
//   // test vs the analytic result
//   double test = -8./3.*norm(_photonvertex->getNorm())*norm(_gluonvertex->getNorm())
//     *sqr(int(mePartonData()[1]->iCharge()))/9.
//     *(uHat()/sHat()+sHat()/uHat());
//   cerr << "testing the matrix element " << test << " " << me[2]/3. 
//        << " " << test/me[2]*3. << "\n";
  return me[2]/3.;
}

double MEGammaP2Jets::gammaantiquarkME(vector<VectorWaveFunction> & gmin,
				       vector<SpinorBarWaveFunction> & fin,
				       vector<VectorWaveFunction> & gout,
				       vector<SpinorWaveFunction> & fout,
				       bool calc) const {
  // the scale
  Energy2 mt(scale());
  // storage of the helicity me if needed
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1Half,
				PDT::Spin1,PDT::Spin1Half);
  // loop over the helicities
  SpinorBarWaveFunction inter;
  vector<double> me(3,0.);
  vector<Complex> diag(2);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
 	  // first diagram
 	  inter = _gluonvertex->evaluate(mt,5,fin[ihel2].particle(),
 					 fin[ihel2],gout[ohel1]);
 	  diag[0] = _photonvertex->evaluate(0.*GeV2,fout[ohel2],inter,gmin[ihel1]);
 	  // second diagram
 	  inter = _photonvertex->evaluate(0.*GeV2,5,fin[ihel2].particle(),
 					  fin[ihel2],gmin[ihel1]);
 	  diag[1] = _gluonvertex->evaluate(mt,fout[ohel2],inter,gout[ohel1]);
 	  for(unsigned int ix=0;ix<2;++ix) me[ix] += norm(diag[ix]);
 	  diag[0]+=diag[1];
 	  me[2] += norm(diag[0]);
 	  // matrix element if needed
 	  if(calc) newme(2*ihel1,ihel2,2*ohel1,ohel2)=diag[0];
	}
      }
    }
  }
  // save the info on the diagrams
  if(!calc) {
    DVector save;
    save.push_back(me[0]);
    save.push_back(me[1]);
    meInfo(save);
  }
  // return the answer with colour and spin factors
  if(calc) _me.reset(newme);
  // test vs the analytic result
//   double test = -8./3.*norm(_photonvertex->getNorm())*norm(_gluonvertex->getNorm())
//     *sqr(int(mePartonData()[1]->iCharge()))/9.
//     *(uHat()/sHat()+sHat()/uHat());
//   cerr << "testing the matrix element " << test << " " << me[2]/3. 
//        << " " << test/me[2]*3. << "\n";
  return me[2]/3.;
}
