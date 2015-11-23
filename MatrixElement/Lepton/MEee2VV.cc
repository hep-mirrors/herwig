// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2VV class.
//

#include "MEee2VV.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"

using namespace Herwig;

MEee2VV::MEee2VV() : process_(0), massOption_(2)  {}

void MEee2VV::doinit() {
  HwMEBase::doinit();
  massOption(vector<unsigned int>(2,massOption_));
  rescalingOption(2);
  // get the vertices we need
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " MEee2VV::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  FFZvertex_ = hwsm->vertexFFZ();
  FFPvertex_ = hwsm->vertexFFP();
  WWWvertex_ = hwsm->vertexWWW();
  FFWvertex_ = hwsm->vertexFFW();
}

void MEee2VV::getDiagrams() const {
  // get the particle data objects we need
  tcPDPtr wPlus  = getParticleData(ParticleID::Wplus );
  tcPDPtr wMinus = getParticleData(ParticleID::Wminus); 
  tcPDPtr z0     = getParticleData(ParticleID::Z0    );
  tcPDPtr gamma  = getParticleData(ParticleID::gamma);
  tcPDPtr em = getParticleData(ParticleID::eminus);
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  tcPDPtr nu_e = getParticleData(ParticleID::nu_e);
  if(process_==0||process_==1) {
    // s-channel Z0 for W+W- production
    add(new_ptr((Tree2toNDiagram(2), em, ep, 1,    z0, 3, wMinus, 3, wPlus, -2)));
    // s-channel photon for W+W- production
    add(new_ptr((Tree2toNDiagram(2), em, ep, 1, gamma, 3, wMinus, 3, wPlus, -1)));
    // t channel for W+W- production
    add(new_ptr((Tree2toNDiagram(3), em, nu_e, ep, 1, wMinus, 2, wPlus, -3)));
  }
  if(process_==0||process_==2) {
    add(new_ptr((Tree2toNDiagram(3), em, em, ep, 1, z0, 2, z0, -1)));
    add(new_ptr((Tree2toNDiagram(3), em, em, ep, 2, z0, 1, z0, -2)));
  }
}

Energy2 MEee2VV::scale() const {
  return sHat();
}

unsigned int MEee2VV::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2VV::orderInAlphaEW() const {
  return 2;
}

Selector<const ColourLines *>
MEee2VV::colourGeometries(tcDiagPtr ) const {
  static ColourLines cl("");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &cl);
  return sel;
}

IBPtr MEee2VV::clone() const {
  return new_ptr(*this);
}

IBPtr MEee2VV::fullclone() const {
  return new_ptr(*this);
}

ClassDescription<MEee2VV> MEee2VV::initMEee2VV;
// Definition of the static class description member.

void MEee2VV::Init() {

  static ClassDocumentation<MEee2VV> documentation
    ("The MEee2VV class simulates the processes e+e->W+W-"
     " and e+e-->Z0Z0 using a 2->2 matrix element");

  static Switch<MEee2VV,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &MEee2VV::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include WW and ZZ",
     0);
  static SwitchOption interfaceProcessWW
    (interfaceProcess,
     "WW",
     "Only include WW",
     1);
  static SwitchOption interfaceProcessZZ
    (interfaceProcess,
     "ZZ",
     "Only include ZZ",
     2);

  static Switch<MEee2VV,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the W/Z mass",
     &MEee2VV::massOption_, 1, false, false);
  static SwitchOption interfaceMassOptionOnMassShell
    (interfaceMassOption,
     "OnMassShell",
     "The W/Z is produced on its mass shell",
     1);
  static SwitchOption interfaceMassOption2
    (interfaceMassOption,
     "OffShell",
     "The W/Z is generated off-shell using the mass and width generator.",
     2);

}

void MEee2VV::persistentOutput(PersistentOStream & os) const {
  os << process_ << massOption_
     << FFPvertex_ << FFWvertex_ << FFZvertex_ << WWWvertex_; 
}

void MEee2VV::persistentInput(PersistentIStream & is, int) {
  is >> process_ >> massOption_
     >> FFPvertex_ >> FFWvertex_ >> FFZvertex_ >> WWWvertex_; 
}

double MEee2VV::me2() const {
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

  // e+e- > Z Z 
  if(v1[0].particle()->id()==ParticleID::Z0) {
    return ZZME(f1,a1,v1,v2);
  }
  // e+e- > W+W-
  else {
    return WWME(f1,a1,v1,v2);
  }
  
}

double MEee2VV::WWME(vector<SpinorWaveFunction>    & f1,
		     vector<SpinorBarWaveFunction> & a1,
		     vector<VectorWaveFunction>    & v1,
		     vector<VectorWaveFunction>    & v2) const {
  double output(0.);
  vector<double> me(3,0.0);
  me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
				    PDT::Spin1,PDT::Spin1));
  ProductionMatrixElement hme[3]={ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
							  PDT::Spin1,PDT::Spin1),
				  ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
							  PDT::Spin1,PDT::Spin1),
				  ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
							  PDT::Spin1,PDT::Spin1)};
  // particle data for the t-channel intermediate
  tcPDPtr nu_e  = getParticleData(ParticleID::nu_e);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr z0    = getParticleData(ParticleID::Z0);
  vector<Complex> diag(3,0.0);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      VectorWaveFunction interP =
	FFPvertex_->evaluate(scale(),3,gamma,f1[ihel1],a1[ihel2]);
      VectorWaveFunction interZ =
	FFZvertex_->evaluate(scale(),3,z0   ,f1[ihel1],a1[ihel2]);
      for(unsigned int ohel1=0;ohel1<3;++ohel1) {
	for(unsigned int ohel2=0;ohel2<3;++ohel2) {
	  diag[0] = WWWvertex_->evaluate(scale(),interP,v2[ohel2],v1[ohel1]);
	  // s-channel Z0
	  diag[1] = WWWvertex_->evaluate(scale(),interZ,v2[ohel2],v1[ohel1]);
	  // t-channel neutrino
	  SpinorWaveFunction inter_nu_e = 
	    FFWvertex_->evaluate(scale(),1,nu_e,f1[ihel1],v1[ohel1]);
	  diag[2] = 
	    FFWvertex_->evaluate(scale(),inter_nu_e,a1[ihel2],v2[ohel2]);
	  // individual diagrams
	  for (size_t ii=0; ii<3; ++ii) {
	    me[ii] += std::norm(diag[ii]);
	    hme[ii](ihel1,ihel2,ohel1,ohel2) = diag[ii];
	  }
	  // full matrix element
	  diag[0] += diag[1]+diag[2];
	  output += std::norm(diag[0]);
	  // storage of the matrix element for spin correlations
	  me_(ihel1,ihel2,ohel1,ohel2) = diag[0];
	}
      }
    }
  }
  DVector save(3);
  for (size_t i = 0; i < 3; ++i) save[i] = 0.25 * me[i];
  output *= 0.25;
  // polarization stuff
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
			 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    for(unsigned int i=0;i<3;++i) me[i] = hme[i].average(rho[0],rho[1]);
    output = me_.average(rho[0],rho[1]);
  }
  meInfo(save);
  // testing code
//   double xW = SM().sin2ThetaW();
//   double Q=-1.;
//   double l = 2.*(-0.5-Q*xW);
//   double r =-2.*Q*xW;
//   Energy2 mW2 = sqr(getParticleData(ParticleID::Wplus)->mass());
//   Energy2 mZ2 = sqr(getParticleData(ParticleID::Z0)->mass());
//   Energy2 sh=sHat(),th=tHat(),uh=uHat();
//   double A = (th*uh/sqr(mW2)-1.)*(0.25-mW2/sh+3.*sqr(mW2/sh))
//     +sh/mW2-4;
//   double bracket = 
//     A*(sqr(Q+0.25*(l+r)/xW*sh/(sh-mZ2))+
//        sqr(  0.25*(l-r)/xW*sh/(sh-mZ2)));
//   if(Q>0) swap(uh,th);
//   double I = (th*uh/sqr(mW2)-1.)*(0.25-0.5*mW2/sh-sqr(mW2)/sh/th)
//     +sh/mW2-2.+2.*mW2/th;

//   double E = (th*uh/sqr(mW2)-1.)*(0.25+sqr(mW2/th))+sh/mW2;
//   if(Q<0.)
//     bracket += 0.5/xW*(Q+0.5*l/xW*sh/(sh-mZ2))*I+0.125/sqr(xW)*E;
//   else
//     bracket +=-0.5/xW*(Q+0.5*l/xW*sh/(sh-mZ2))*I+0.125/sqr(xW)*E;
//   InvEnergy4 dsigdt = 2.*Constants::pi*sqr(SM().alphaEM(scale()))/sqr(sh)*bracket;
//   double test = 16.*Constants::pi*sqr(sHat())*dsigdt;
//   cerr << "testing " << test << " " << output << " " << test/output << "\n";
  return output;
}

double MEee2VV::ZZME(vector<SpinorWaveFunction>    & f1,
		     vector<SpinorBarWaveFunction> & a1,
		     vector<VectorWaveFunction>    & v1,
		     vector<VectorWaveFunction>    & v2) const {
  double output(0.);
  vector<double> me(3,0.0);
  me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
				    PDT::Spin1,PDT::Spin1));
  ProductionMatrixElement hme[2]={ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
							  PDT::Spin1,PDT::Spin1),
				  ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
							  PDT::Spin1,PDT::Spin1)};
  tcPDPtr em  = getParticleData(ParticleID::eminus);
  vector<Complex> diag(2,0.0);
  SpinorWaveFunction inter;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<3;++ohel1) {
	for(unsigned int ohel2=0;ohel2<3;++ohel2) {
	  inter   = FFZvertex_->evaluate(scale(),1,em,f1[ihel1] ,v1[ohel1]);
	  diag[0] = FFZvertex_->evaluate(scale(),inter,a1[ihel2],v2[ohel2]);
	  inter   = FFZvertex_->evaluate(scale(),1,em,f1[ihel1] ,v2[ohel2]);
	  diag[1] = FFZvertex_->evaluate(scale(),inter,a1[ihel2],v1[ohel1]);
	  // individual diagrams
	  for (size_t ii=0; ii<2; ++ii) {
	    me[ii] += std::norm(diag[ii]);
	    hme[ii](ihel1,ihel2,ohel1,ohel2) = diag[ii];
	  }
	  // full matrix element
	  diag[0] += diag[1];
	  output += std::norm(diag[0]);
	  // storage of the matrix element for spin correlations
	  me_(ihel1,ihel2,ohel1,ohel2) = diag[0];
	}
      }
    }
  }
  DVector save(3);
  for (size_t i = 0; i < 3; ++i) save[i] = 0.25 * me[i];
  meInfo(save);
  // spin average
  output *= 0.25;
  // polarization stuff
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
			 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    for(unsigned int i=0;i<2;++i) me[i] = hme[i].average(rho[0],rho[1]);
    output = me_.average(rho[0],rho[1]);
  }
  // identical particle factor
  output /= 2.;
  // testing code
//   double xW = SM().sin2ThetaW();
//   double Q=-1.;
//   double l = 2.*(-0.5-Q*xW);
//   double r =-2.*Q*xW;

//   Energy2 mZ2 = sqr(getParticleData(ParticleID::Z0)->mass());
//   Energy2 sh=sHat(),th=tHat(),uh=uHat();
//   InvEnergy4 dsigdt = Constants::pi*sqr(SM().alphaEM(scale()))/32.*
//     (pow(l,4)+pow(r,4))/sqr(xW)/sqr(1.-xW)/sqr(sh)
//     *(th/uh+uh/th+4.*mZ2*sh/th/uh-sqr(mZ2)*(1./sqr(th)+1./sqr(uh)));
//   double test = 16.*Constants::pi*sqr(sHat())*dsigdt;
//   cerr << "testing " << (test-output)/(test+output) << "\n";
  return output;
}

Selector<MEBase::DiagramIndex>
MEee2VV::diagrams(const DiagramVector & diags) const {
  vector<double> last(3);
  if ( lastXCombPtr() ) {
    for(unsigned int ix=0;ix<3;++ix) last[ix] = meInfo()[ix];
  }
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if(diags[i]->id() >= -3 ) sel.insert(last[-diags[i]->id() - 1],i);
  }
  return sel;
}

double MEee2VV::getCosTheta(double ctmin, double ctmax, const double r) {
  double rand = r;
  Energy2 m12 = sqr(meMomenta()[2].mass());
  Energy2 m22 = sqr(meMomenta()[3].mass());
  Energy2 D1 = sHat()-m12-m22;
  Energy4 lambda = sqr(D1) - 4*m12*m22;
  double D =  D1 / sqrt(lambda);
  if(abs(mePartonData()[2]->id())==ParticleID::Wplus) {
    double fraction = (D-ctmax)/(D-ctmin);
    double costh = D - (D - ctmin) * pow(fraction, rand);
    jacobian((costh - D) * log(fraction));
    return costh;
  }
  else {
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
}


void MEee2VV::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  if(hard[order[0]]->id()<0) swap(order[0],order[1]);
  if(hard[order[3]]->id()<0) swap(order[2],order[3]);
  vector<SpinorWaveFunction> f;
  vector<SpinorBarWaveFunction>  fbar;
  SpinorWaveFunction   (f   ,hard[order[0]],incoming,false);
  SpinorBarWaveFunction(fbar,hard[order[1]],incoming,false);
  vector<VectorWaveFunction> w1,w2;
  VectorWaveFunction   (w1,hard[order[2]],outgoing,true ,false);
  VectorWaveFunction   (w2,hard[order[3]],outgoing,true ,false);
  if(hard[order[2]]->id()==ParticleID::Z0) {
    ZZME(f,fbar,w1,w2);
  }
  else {
    WWME(f,fbar,w1,w2);
  }  
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) {
    tcSpinPtr spin = hard[order[ix]]->spinInfo();
    if(ix<2) {
      tcPolarizedBeamPDPtr beam = 
	dynamic_ptr_cast<tcPolarizedBeamPDPtr>(hard[ix]->dataPtr());
      if(beam) spin->rhoMatrix() = beam->rhoMatrix();
    }
    spin->productionVertex(hardvertex);
  }
}
