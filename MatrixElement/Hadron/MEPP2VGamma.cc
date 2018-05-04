// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2VGamma class.
//

#include "MEPP2VGamma.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

MEPP2VGamma::MEPP2VGamma()  : process_(0), maxflavour_(5), massOption_(2) 
{}

unsigned int MEPP2VGamma::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2VGamma::orderInAlphaEW() const {
  return 2;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2VGamma,HwMEBase>
describeHerwigMEPP2VGamma("Herwig::MEPP2VGamma", "HwMEHadron.so");

void MEPP2VGamma::Init() {

  static ClassDocumentation<MEPP2VGamma> documentation
    ("The MEPP2VGamma class simulates the production of"
     " W+/-gamma and Z0gamma in hadron-hadron collisions"
     " using the 2->2 matrix elements");

  static Switch<MEPP2VGamma,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &MEPP2VGamma::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all the processes",
     0);
  static SwitchOption interfaceProcessWGamma
    (interfaceProcess,
     "WGamma",
     "Only include W+/-gamma",
     1);
  static SwitchOption interfaceProcessZGamma
    (interfaceProcess,
     "ZGamma",
     "Only include ZGamma",
     2);

  static Parameter<MEPP2VGamma,int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour allowed for the incoming quarks",
     &MEPP2VGamma::maxflavour_, 5, 2, 5,
     false, false, Interface::limited);

  static Switch<MEPP2VGamma,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the boson masses",
     &MEPP2VGamma::massOption_, 1, false, false);
  static SwitchOption interfaceMassOptionOnMassShell
    (interfaceMassOption,
     "OnMassShell",
     "The boson is produced on its mass shell",
     1);
  static SwitchOption interfaceMassOption2
    (interfaceMassOption,
     "OffShell",
     "The bosons are generated off-shell using the mass and width generator.",
     2);

}

void MEPP2VGamma::persistentOutput(PersistentOStream & os) const {
  os << FFPvertex_ << FFWvertex_ << FFZvertex_ << WWWvertex_ 
     << process_ << massOption_;
}

void MEPP2VGamma::persistentInput(PersistentIStream & is, int) {
  is >> FFPvertex_ >> FFWvertex_ >> FFZvertex_ >> WWWvertex_ 
     >> process_ >> massOption_;
}

Energy2 MEPP2VGamma::scale() const {
  return sHat();
}

IBPtr MEPP2VGamma::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2VGamma::fullclone() const {
  return new_ptr(*this);
}

void MEPP2VGamma::doinit() {
  HwMEBase::doinit();
  // mass option
  vector<unsigned int> mopt(2,1);
  mopt[0]=massOption_;
  massOption(mopt);
  rescalingOption(2);
  // get the vertices we need
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " MEPP2VGamma::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  FFZvertex_ = hwsm->vertexFFZ();
  FFPvertex_ = hwsm->vertexFFP();
  WWWvertex_ = hwsm->vertexWWW();
  FFWvertex_ = hwsm->vertexFFW();
}

Selector<const ColourLines *>
MEPP2VGamma::colourGeometries(tcDiagPtr diag) const {
  static ColourLines cs("1 -2");
  static ColourLines ct("1 2 -3");
  Selector<const ColourLines *> sel; 
  if(diag->id()<-2) sel.insert(1.0, &cs);
  else              sel.insert(1.0, &ct);
  return sel;
}

void MEPP2VGamma::getDiagrams() const {
  typedef std::vector<pair<tcPDPtr,tcPDPtr> > Pairvector;
  tcPDPtr wPlus  = getParticleData(ParticleID::Wplus );
  tcPDPtr wMinus = getParticleData(ParticleID::Wminus); 
  tcPDPtr z0     = getParticleData(ParticleID::Z0    );
  tcPDPtr gamma  = getParticleData(ParticleID::gamma);
  // W+/- gamma
  if(process_==0||process_==1) {
    // possible parents
    Pairvector parentpair;
    parentpair.reserve(6);
    // don't even think of putting 'break' in here!
    switch(maxflavour_) {
    case 5:
      parentpair.push_back(make_pair(getParticleData(ParticleID::b),
				     getParticleData(ParticleID::cbar)));
      parentpair.push_back(make_pair(getParticleData(ParticleID::b), 
				     getParticleData(ParticleID::ubar)));
      [[fallthrough]];
    case 4:
      parentpair.push_back(make_pair(getParticleData(ParticleID::s),
				     getParticleData(ParticleID::cbar)));
      parentpair.push_back(make_pair(getParticleData(ParticleID::d),
				     getParticleData(ParticleID::cbar)));
      [[fallthrough]];
    case 3:
      parentpair.push_back(make_pair(getParticleData(ParticleID::s),
				     getParticleData(ParticleID::ubar)));
      [[fallthrough]];
    case 2:
      parentpair.push_back(make_pair(getParticleData(ParticleID::d),
				     getParticleData(ParticleID::ubar)));
      [[fallthrough]];
    default:
      ;
    }
    // W+ gamma
    for(unsigned int ix=0;ix<parentpair.size();++ix) {
      add(new_ptr((Tree2toNDiagram(3), parentpair[ix].second->CC(), 
		   parentpair[ix].first, parentpair[ix].first->CC(),
		   1, wPlus, 2, gamma, -1)));
      add(new_ptr((Tree2toNDiagram(3), parentpair[ix].second->CC(), 
		   parentpair[ix].second->CC() , parentpair[ix].first->CC(),
		   2, wPlus, 1, gamma, -2)));
      add(new_ptr((Tree2toNDiagram(2), parentpair[ix].second->CC(),
		   parentpair[ix].first->CC(), 1, wPlus, 3, wPlus, 3, gamma,  -3)));
    }
    // W- gamma
    for(unsigned int ix=0;ix<parentpair.size();++ix) {
      add(new_ptr((Tree2toNDiagram(3), parentpair[ix].first, 
		   parentpair[ix].second->CC(),
		   parentpair[ix].second, 1, wMinus, 2, gamma, -1)));
      add(new_ptr((Tree2toNDiagram(3), parentpair[ix].first, parentpair[ix].first ,
		   parentpair[ix].second, 2, wMinus, 1, gamma, -2)));
      add(new_ptr((Tree2toNDiagram(2), parentpair[ix].first,
		   parentpair[ix].second, 1, wMinus, 3, wMinus, 3, gamma,  -3))); 
    }
  }
  if(process_==0||process_==2) {
    for(int ix=1;ix<=maxflavour_;++ix) {
      tcPDPtr qk = getParticleData(ix);
      tcPDPtr qb = qk->CC();
      add(new_ptr((Tree2toNDiagram(3), qk, qk, qb, 1, z0, 2, gamma, -1)));
      add(new_ptr((Tree2toNDiagram(3), qk, qk, qb, 2, z0, 1, gamma, -2)));
    }
  }
}

Selector<MEBase::DiagramIndex>
MEPP2VGamma::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    sel.insert(meInfo()[abs(diags[i]->id()) - 1], i);
  return sel;
}

double MEPP2VGamma::me2() const {
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
    if(ix!=1) {
      v2_out.reset(ix);
      v2.push_back(v2_out);
    }
  }
  if(mePartonData()[2]->id()==ParticleID::Z0) {
    return ZGammaME(f1,a1,v1,v2,false);
  }
  else {
    return WGammaME(f1,a1,v1,v2,false);
  }
}

double MEPP2VGamma::ZGammaME(vector<SpinorWaveFunction>    & f1,
			     vector<SpinorBarWaveFunction> & a1,
			     vector<VectorWaveFunction>    & v1,
			     vector<VectorWaveFunction>    & v2,
			     bool calc) const {
  double output(0.);
  vector<double> me(3,0.0);
  if(calc) me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1,PDT::Spin1));
  vector<Complex> diag(2,0.0);
  SpinorWaveFunction inter;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<3;++ohel1) {
 	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
 	  inter   = FFZvertex_->evaluate(scale(),5,f1[ihel1].particle()->CC(),
 					 f1[ihel1],v1[ohel1]);
 	  diag[0] = FFPvertex_->evaluate(scale(),inter,a1[ihel2],v2[ohel2]);
 	  inter   = FFPvertex_->evaluate(scale(),5,f1[ihel1].particle()->CC(),
 					 f1[ihel1] ,v2[ohel2]);
 	  diag[1] = FFZvertex_->evaluate(scale(),inter,a1[ihel2],v1[ohel1]);
 	  // individual diagrams
	  for (size_t ii=0; ii<2; ++ii) me[ii] += std::norm(diag[ii]);
	  // full matrix element
	  diag[0] += diag[1];
	  output += std::norm(diag[0]);
	  // storage of the matrix element for spin correlations
	  if(calc) me_(ihel1,ihel2,ohel1,ohel2) = diag[0];
	}
      }
    }
  }
  DVector save(3);
  for (size_t i = 0; i < 3; ++i) {
    save[i] = 0.25 * me[i];
  }
  meInfo(save);
  return 0.25*output/3.;
}

double MEPP2VGamma::WGammaME(vector<SpinorWaveFunction>    & f1,
			     vector<SpinorBarWaveFunction> & a1,
			     vector<VectorWaveFunction>    & v1,
			     vector<VectorWaveFunction>    & v2,
			     bool calc) const {
  double output(0.);
  vector<double> me(3,0.0);
  if(calc) me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1,PDT::Spin1));
  vector<Complex> diag(3,0.0);
  SpinorWaveFunction inter;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      VectorWaveFunction interW =
 	FFWvertex_->evaluate(scale(),3,v1[0].particle(),
 			     f1[ihel1],a1[ihel2]);
      for(unsigned int ohel1=0;ohel1<3;++ohel1) {
 	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // t-channel diagrams
	  inter   = FFWvertex_->evaluate(scale(),5,a1[ihel1].particle(),
					 f1[ihel1],v1[ohel1]);
	  diag[0] = FFPvertex_->evaluate(scale(),inter,a1[ihel2],v2[ohel2]);
	  inter   = FFPvertex_->evaluate(scale(),5,f1[ihel1].particle()->CC(),
					 f1[ihel1] ,v2[ohel2]);
	  diag[1] = FFWvertex_->evaluate(scale(),inter,a1[ihel2],v1[ohel1]);
	  // s-channel diagram
	  diag[2] = WWWvertex_->evaluate(scale(),interW,v1[ohel1],v2[ohel2]);
	  // individual diagrams
	  for (size_t ii=0; ii<3; ++ii) me[ii] += std::norm(diag[ii]);
 	  // full matrix element
 	  diag[0] += diag[1]+diag[2];
 	  output += std::norm(diag[0]);
 	  // storage of the matrix element for spin correlations
	  if(calc) me_(ihel1,ihel2,ohel1,ohel2) = diag[0];
	}
      }
    }
  }
  DVector save(3);
  for (size_t i = 0; i < 3; ++i) save[i] = 0.25 * me[i];
  meInfo(save);
  // spin and colour factors
  output *= 0.25/3.;
  // testing code
//   int iu = abs(mePartonData()[0]->id());
//   int id = abs(mePartonData()[1]->id());
//   if(iu%2!=0) swap(iu,id);
//   iu = (iu-2)/2;
//   id = (id-1)/2;
//   double ckm = SM().CKM(iu,id);
//   InvEnergy4 dsigdt =  Constants::pi*sqr(SM().alphaEM(scale()))
//     /6./sqr(sHat())/SM().sin2ThetaW()*sqr(1./(1.+tHat()/uHat())-1./3.)*
//     (sqr(tHat())+sqr(uHat())+2.*sqr(getParticleData(ParticleID::Wplus)->mass())*sHat())/
//     tHat()/uHat();
//   double test = 16.*ckm*Constants::pi*sqr(sHat())*dsigdt;
//   cerr << "testing W gamma " << test << " " << output << " " 
//        << (test-output)/(test+output) << "\n";
  return output;
}

void MEPP2VGamma::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  if(hard[order[0]]->id()<0) swap(order[0],order[1]);
  vector<SpinorWaveFunction> q;
  vector<SpinorBarWaveFunction>  qb;
  SpinorWaveFunction   (q ,hard[order[0]],incoming,false);
  SpinorBarWaveFunction(qb,hard[order[1]],incoming,false);
  vector<VectorWaveFunction> w1,w2;
  if(hard[order[2]]->id()==ParticleID::gamma)
    swap(order[2],order[3]);
  VectorWaveFunction   (w1,hard[order[2]],outgoing,true ,false);
  VectorWaveFunction   (w2,hard[order[3]],outgoing,true ,true );
  w2[1]=w2[2];
  // q qbar -> Z gamma
  if(hard[order[2]]->id()==ParticleID::Z0) {
    ZGammaME(q,qb,w1,w2,true);
  }
  // q qbar -> W gamma
  else {
    WGammaME(q,qb,w1,w2,true);
  }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix)
    hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
}
