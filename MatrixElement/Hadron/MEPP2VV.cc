// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2VV class.
//

#include "MEPP2VV.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"

using namespace Herwig;

MEPP2VV::MEPP2VV() : _process(0), _maxflavour(5) 
{}

unsigned int MEPP2VV::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2VV::orderInAlphaEW() const {
  return 2;
}

ClassDescription<MEPP2VV> MEPP2VV::initMEPP2VV;
// Definition of the static class description member.

void MEPP2VV::Init() {

  static ClassDocumentation<MEPP2VV> documentation
    ("The MEPP2VV class simulates the production of W+W-, "
     "W+/-Z0 and Z0Z0 in hadron-hadron collisions using the 2->2"
     " matrix elements");

  static Switch<MEPP2VV,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &MEPP2VV::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all the processes",
     0);
  static SwitchOption interfaceProcessWW
    (interfaceProcess,
     "WW",
     "Only include W+W-",
     1);
  static SwitchOption interfaceProcessWZ
    (interfaceProcess,
     "WZ",
     "Only include W+/-Z",
     2);
  static SwitchOption interfaceProcessZZ
    (interfaceProcess,
     "ZZ",
     "Only include ZZ",
     3);

  static Parameter<MEPP2VV,int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour allowed for the incoming quarks",
     &MEPP2VV::_maxflavour, 5, 2, 5,
     false, false, Interface::limited);

}

void MEPP2VV::persistentOutput(PersistentOStream & os) const {
  os << _vertexFFP << _vertexFFW << _vertexFFZ << _vertexWWW << _process;
}

void MEPP2VV::persistentInput(PersistentIStream & is, int) {
  is >> _vertexFFP >> _vertexFFW >> _vertexFFZ >> _vertexWWW >> _process;
}

Energy2 MEPP2VV::scale() const {
  return sqr(mePartonData()[2]->mass())+sqr(mePartonData()[3]->mass());
}

IBPtr MEPP2VV::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2VV::fullclone() const {
  return new_ptr(*this);
}

void MEPP2VV::doinit() throw(InitException) {
  HwME2to2Base::doinit();
  // get the vertices we need
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " MEPP2VV::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  _vertexFFZ = hwsm->vertexFFZ();
  _vertexFFP = hwsm->vertexFFP();
  _vertexWWW = hwsm->vertexWWW();
  _vertexFFW = hwsm->vertexFFW();
}

double MEPP2VV::getCosTheta(double ctmin, double ctmax, const double * r) {
  double rand = *r;
  Energy2 m12 = sqr(meMomenta()[2].mass());
  Energy2 m22 = sqr(meMomenta()[3].mass());
  Energy2 D1 = sHat()-m12-m22;
  Energy4 lambda = sqr(D1) - 4*m12*m22;
  double D =  D1 / sqrt(lambda);
  if(mePartonData()[2]->id()==ParticleID::Z0&&
     mePartonData()[3]->id()==ParticleID::Z0) {
    double prob = 0.5;
    double costh;
    double fraction1 = (D-ctmax)/(D-ctmin);
    double fraction2 = (D+ctmin)/(D+ctmax);
    if(rand<=prob) {
      rand /=prob;
      costh = D - (D - ctmin) * pow(fraction1, rand);
    }
    else {
      rand = (rand-prob)/(1.-prob);
      costh =-D + (D + ctmax) * pow(fraction2, rand);
    }
    jacobian(1./(prob     /((costh - D) * log(fraction1))-
		 (1.-prob)/((costh + D) * log(fraction2))));
    return costh;
  }
  else {
    double fraction = (D-ctmax)/(D-ctmin);
    double costh = D - (D - ctmin) * pow(fraction, rand);
    jacobian((costh - D) * log(fraction));
    return costh;
  }
}

Selector<const ColourLines *>
MEPP2VV::colourGeometries(tcDiagPtr diag) const {
  static ColourLines cs("1 -2");
  static ColourLines ct("1 2 -3");
  Selector<const ColourLines *> sel; 
  if(diag->id()>0) sel.insert(1.0, &cs);
  else             sel.insert(1.0, &ct);
  return sel;
}

void MEPP2VV::getDiagrams() const {
  typedef std::vector<pair<tcPDPtr,tcPDPtr> > Pairvector;
  tcPDPtr wPlus  = getParticleData(ParticleID::Wplus );
  tcPDPtr wMinus = getParticleData(ParticleID::Wminus); 
  tcPDPtr z0     = getParticleData(ParticleID::Z0    );
  tcPDPtr gamma  = getParticleData(ParticleID::gamma);
  // W+ W-
  if(_process==0||_process==1) {
    for(int ix=1;ix<=_maxflavour;++ix) {
      tcPDPtr qk = getParticleData(ix);
      tcPDPtr w1 = ix%2==0 ? wPlus : wMinus;
      tcPDPtr w2 = ix%2!=0 ? wPlus : wMinus;
      for(int iy=1;iy<=_maxflavour;++iy) {
	if(abs(ix-iy)%2!=0) continue;
	tcPDPtr qb = getParticleData(-iy);
	// s channel photon
	add(new_ptr((Tree2toNDiagram(2), qk, qb, 1, gamma, 3, w1, 3, w1,  4)));
	// s-channel Z
	add(new_ptr((Tree2toNDiagram(2), qk, qb, 1,    z0, 3, w1, 3, w2,  5)));
	// t-channel
	if(ix%2==0) {
	  int idiag=0;
	  for(int iz=1;iz<=5;iz+=2) {
	    --idiag;
	    tcPDPtr tc = getParticleData(iz);
	    add(new_ptr((Tree2toNDiagram(3), qk, tc, qb, 1, w1, 2, w2, idiag)));
	  }
	}
	else {
	  int idiag=0;
	  for(int iz=2;iz<=6;iz+=2) {
	    --idiag;
	    tcPDPtr tc = getParticleData(iz);
	    add(new_ptr((Tree2toNDiagram(3), qk, tc, qb, 1, w1, 2, w2, idiag)));
	  }
	}
      }
    }
  }
  // W+/- Z
  if(_process==0||_process==2) {
    // possible parents
    Pairvector parentpair;
    parentpair.reserve(6);
    // don't even think of putting 'break' in here!
    switch(_maxflavour) {
    case 5:
//       parentpair.push_back(make_pair(getParticleData(ParticleID::b),
// 				     getParticleData(ParticleID::cbar)));
//       parentpair.push_back(make_pair(getParticleData(ParticleID::b), 
// 				     getParticleData(ParticleID::ubar)));
    case 4:
      parentpair.push_back(make_pair(getParticleData(ParticleID::s),
				     getParticleData(ParticleID::cbar)));
//       parentpair.push_back(make_pair(getParticleData(ParticleID::d),
// 				     getParticleData(ParticleID::cbar)));
    case 3:
//       parentpair.push_back(make_pair(getParticleData(ParticleID::s),
// 				     getParticleData(ParticleID::ubar)));
    case 2:
      parentpair.push_back(make_pair(getParticleData(ParticleID::d),
				     getParticleData(ParticleID::ubar)));
    default:
      ;
    }
    // W+ Z
    for(unsigned int ix=0;ix<parentpair.size();++ix) {
      add(new_ptr((Tree2toNDiagram(3), parentpair[ix].second->CC(), 
		   parentpair[ix].first, parentpair[ix].first->CC(),
		   1, wPlus, 2, z0, -1)));
      add(new_ptr((Tree2toNDiagram(3), parentpair[ix].second->CC(), 
		   parentpair[ix].second->CC() , parentpair[ix].first->CC(),
		   2, wPlus, 1, z0, -2)));
      add(new_ptr((Tree2toNDiagram(2), parentpair[ix].second->CC(),
		   parentpair[ix].first->CC(), 1, wPlus, 3, wPlus, 3, z0,  3)));
    }
    // W- Z
    for(unsigned int ix=0;ix<parentpair.size();++ix) {
      add(new_ptr((Tree2toNDiagram(3), parentpair[ix].first, 
		   parentpair[ix].second->CC(),
		   parentpair[ix].second, 1, wMinus, 2, z0, -1)));
      add(new_ptr((Tree2toNDiagram(3), parentpair[ix].first, 
		   parentpair[ix].first , parentpair[ix].second, 2, wMinus, 1, z0, -2)));
      add(new_ptr((Tree2toNDiagram(2), parentpair[ix].first,
		   parentpair[ix].second, 1, wMinus, 3, wMinus, 3, z0,  3))); 
    }
  }
  // Z Z
  if(_process==0||_process==3) {
    for(int ix=1;ix<_maxflavour;++ix) {
      tcPDPtr qk = getParticleData(ix);
      tcPDPtr qb = qk->CC();
      add(new_ptr((Tree2toNDiagram(3), qk, qk, qb, 1, z0, 2, z0, -1)));
      add(new_ptr((Tree2toNDiagram(3), qk, qk, qb, 2, z0, 1, z0, -2)));
    }
  }
}

Selector<MEBase::DiagramIndex>
MEPP2VV::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    sel.insert(meInfo()[abs(diags[i]->id())], i);
  return sel;
}

double MEPP2VV::me2() const {
  // setup momenta and particle data for the external wavefunctions 
  // incoming
  SpinorWaveFunction    em_in( meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction ep_in( meMomenta()[1],mePartonData()[1],incoming);
  // outgoing
  VectorWaveFunction v1_out(meMomenta()[2],mePartonData()[2],outgoing);
  VectorWaveFunction v2_out(meMomenta()[3],mePartonData()[3],outgoing);
  vector<SpinorWaveFunction> f1;
  vector<SpinorBarWaveFunction> a1;
  vector<VectorWaveFunction> v1,v2;
  // calculate the wavefunctions
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<2) {
      em_in.reset(ix);
      f1.push_back(em_in);
      ep_in.reset(ix);
      a1.push_back(ep_in);
    }
    v1_out.reset(ix);
    v1.push_back(v1_out);
    v2_out.reset(ix);
    v2.push_back(v2_out);
  }
  return helicityME(f1,a1,v1,v2,false);
}

double MEPP2VV::helicityME(vector<SpinorWaveFunction>    & f1,
			   vector<SpinorBarWaveFunction> & a1,
			   vector<VectorWaveFunction>    & v1,
			   vector<VectorWaveFunction>    & v2,
			   bool calc) const {
//   cerr << "testing process" 
//        << mePartonData()[0]->PDGName() << " "
//        << mePartonData()[1]->PDGName() << " -> "
//        << mePartonData()[2]->PDGName() << " "
//        << mePartonData()[3]->PDGName() << "\n";
  double output(0.);
  vector<double> me(5,0.0);
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1,PDT::Spin1);
  // q qbar -> Z Z
  if(mePartonData()[2]->id()==ParticleID::Z0&&
     mePartonData()[3]->id()==ParticleID::Z0) {
    vector<Complex> diag(2,0.0);
    SpinorWaveFunction inter;
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	for(unsigned int ohel1=0;ohel1<3;++ohel1) {
	  for(unsigned int ohel2=0;ohel2<3;++ohel2) {
	    inter   = _vertexFFZ->evaluate(scale(),1,f1[ihel1].getParticle(),
					   f1[ihel1],v1[ohel1]);
 	    diag[0] = _vertexFFZ->evaluate(scale(),inter,a1[ihel2],v2[ohel2]);
 	    inter   = _vertexFFZ->evaluate(scale(),1,f1[ihel1].getParticle(),
					   f1[ihel1] ,v2[ohel2]);
	    diag[1] = _vertexFFZ->evaluate(scale(),inter,a1[ihel2],v1[ohel1]);
	    // individual diagrams
	    for (size_t ii=0; ii<2; ++ii) me[ii] += std::norm(diag[ii]);
	    // full matrix element
	    diag[0] += diag[1];
	    output += std::norm(diag[0]);
	    // storage of the matrix element for spin correlations
	    if(calc) newme(ihel1,ihel2,ohel1,ohel2) = diag[0];
	  }
	}
      }
    }
    // identical particle factor
    output /= 2.;
  }
  // q qbar -> W W
  else if(abs(mePartonData()[2]->id())==ParticleID::Wplus&&
	  abs(mePartonData()[3]->id())==ParticleID::Wplus) {
    // particle data for the t-channel intermediate
    tcPDPtr tc[3];
    if(f1[0].getParticle()->id()%2==0) {
      for(unsigned int ix=0;ix<3;++ix) tc[ix] = getParticleData(1+2*ix);
    }
    else {
      for(unsigned int ix=0;ix<3;++ix) tc[ix] = getParticleData(2+2*ix);
    }
    tcPDPtr gamma = getParticleData(ParticleID::gamma);
    tcPDPtr z0    = getParticleData(ParticleID::Z0);
    vector<Complex> diag(5,0.0);
    VectorWaveFunction interP,interZ;
    bool sChannel = 
      f1[0].getParticle()->id() == -a1[0].getParticle()->id();
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	if(sChannel) {
	  interP = _vertexFFP->evaluate(scale(),1,gamma,f1[ihel1],a1[ihel2]);
	  interZ = _vertexFFZ->evaluate(scale(),1,z0   ,f1[ihel1],a1[ihel2]);
	}
	for(unsigned int ohel1=0;ohel1<3;++ohel1) {
	  for(unsigned int ohel2=0;ohel2<3;++ohel2) {
	    // s-channel photon
	    diag[3] = sChannel ? 
	      _vertexWWW->evaluate(scale(),interP,v2[ohel2],v1[ohel1]) : 0.;
	    // s-channel Z0
	    diag[4] = sChannel ? 
	      _vertexWWW->evaluate(scale(),interZ,v2[ohel2],v1[ohel1]) : 0.;
	    // t-channel
	    for(unsigned int ix=0;ix<3;++ix) {
// 	      cerr << "testing in loop " << f1[ihel1].getParticle()->PDGName()
// 		   << " " << tc[ix]->PDGName() << " " 
// 		   << a1[ihel2].getParticle()->PDGName()
// 		   << "\n";
	      SpinorWaveFunction inter = 
		_vertexFFW->evaluate(scale(),1,tc[ix],f1[ihel1],v1[ohel1]);
	      diag[ix] = 
		_vertexFFW->evaluate(scale(),inter,a1[ihel2],v2[ohel2]);
	    }
	    // individual diagrams
	    for (size_t ii=0; ii<5; ++ii) me[ii] += std::norm(diag[ii]);
	    // full matrix element
// 	    cerr << "testing diagrams " 
// 		 << diag[0] << " " << diag[1] << " "
// 		 << diag[2] << " " << diag[3] << " "
// 		 << diag[4] << "\n"; 
	    diag[0] += diag[1]+diag[2]+diag[3]+diag[4];
	    output += std::norm(diag[0]);
	    // storage of the matrix element for spin correlations
	    if(calc) newme(ihel1,ihel2,ohel1,ohel2) = diag[0];
	  }
	}
      }
    }
  }
  // q qbar -> W Z
  else {
    vector<Complex> diag(3,0.0);
    SpinorWaveFunction inter;
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	VectorWaveFunction interW =
	  _vertexFFW->evaluate(scale(),1,v1[0].getParticle()->CC(),
			       f1[ihel1],a1[ihel2]);
	for(unsigned int ohel1=0;ohel1<3;++ohel1) {
	  for(unsigned int ohel2=0;ohel2<3;++ohel2) {
	    // t-channel diagrams
	    inter   = _vertexFFW->evaluate(scale(),1,a1[ihel1].getParticle(),
					   f1[ihel1],v1[ohel1]);
 	    diag[0] = _vertexFFZ->evaluate(scale(),inter,a1[ihel2],v2[ohel2]);
 	    inter   = _vertexFFZ->evaluate(scale(),1,f1[ihel1].getParticle(),
					   f1[ihel1] ,v2[ohel2]);
	    diag[1] = _vertexFFW->evaluate(scale(),inter,a1[ihel2],v1[ohel1]);
	    // s-channel diagram
	    diag[2] = _vertexWWW->evaluate(scale(),interW,v1[ohel1],v2[ohel2]);
	    // individual diagrams
	    for (size_t ii=0; ii<3; ++ii) me[ii] += std::norm(diag[ii]);
	    // full matrix element
	    diag[0] += diag[1]+diag[2];
	    output += std::norm(diag[0]);
	    // storage of the matrix element for spin correlations
	    if(calc) newme(ihel1,ihel2,ohel1,ohel2) = diag[0];
	  }
	}
      }
    }
  }
  DVector save(5);
  for (size_t i = 0; i < 5; ++i) {
    save[i] = 0.25 * me[i];
  }
  meInfo(save);
  return 0.25*output/3.;
}
