// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PP2VGammaJet class.
//

#include "MEPP2VGammaJet.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Cuts/Cuts.h"
#include <numeric>
#include "Herwig++/Models/StandardModel/StandardModel.h"

using namespace Herwig;

Energy2 MEPP2VGammaJet::scale() const {
  return sqr(mePartonData()[2]->mass());
}

int MEPP2VGammaJet::nDim() const {
  return 5;
}

unsigned int MEPP2VGammaJet::orderInAlphaS() const {
  return 1;
}

unsigned int MEPP2VGammaJet::orderInAlphaEW() const {
  return 2;
}

void MEPP2VGammaJet::doinit() {
  HwMEBase::doinit();
  vector<unsigned int> mopt(3,1);
  massOption(mopt);
  // get the vertices we need
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " MEPP2VGamma::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  FFZvertex_ = dynamic_ptr_cast<FFVVertexPtr>(hwsm->vertexFFZ());
  FFPvertex_ = hwsm->vertexFFP();
  WWWvertex_ = hwsm->vertexWWW();
  FFWvertex_ = dynamic_ptr_cast<FFVVertexPtr>(hwsm->vertexFFW());
  FFGvertex_ = hwsm->vertexFFG();
//   // test of the phase-space volume
//   Energy pTgammaMin = ZERO;
//   Energy pTjetMin   = ZERO;
//   Energy rootS = 1000.*GeV;
//   Energy2 sHat = sqr(rootS);
//   Energy mV = getParticleData(ParticleID::Z0)->mass();
//   unsigned int npoint = 100000000;
//   double wgtsum(0.),wgtsq(0.);
//   for(unsigned int ix=0;ix<npoint;++ix) {
//     double wgt = 1.;
//     Energy eMax = 0.5*(sHat-sqr(mV))/rootS;
//     if(eMax<pTgammaMin) continue;
//     Energy eGamma = pTgammaMin+UseRandom::rnd()*(eMax-pTgammaMin);
//     wgt *= (eMax-pTgammaMin)/rootS;
//     // generate the jet energy
//     Energy2 m122 = sHat-2.*rootS*eGamma;
//     Energy m12 = sqrt(m122);
//     Energy E2star = 0.5*(m122+sqr(mV))/m12;
//     Energy E3star = 0.5*(sHat-m122)/m12;
//     eMax = 
//       0.5/rootS*(sHat-sqr(mV)-2.*E3star*(E2star-sqrt(sqr(E2star)-sqr(mV))));
//     Energy eMin = 
//       0.5/rootS*(sHat-sqr(mV)-2.*E3star*(E2star+sqrt(sqr(E2star)-sqr(mV))));
//     eMin = max(eMin,pTjetMin);
//     if(eMax<eMin) continue;
//     Energy eJet = eMin+UseRandom::rnd()*(eMax-eMin);
//     wgt *= (eMax-eMin)/rootS;
//     // construct the momenta
//     Energy eV = rootS-eGamma-eJet;
//     if(eV<mV) continue;
//     Energy pV = sqrt(sqr(eV)-sqr(mV));
//     double cos3 = 0.5*(sqr(eJet)+sqr(eGamma)-sqr(pV))/eJet/eGamma;
//     double sin3 = sqrt(1.-sqr(cos3));
//     // finally rotate
//     wgt *= 8.*sqr(Constants::pi);
//     wgt *= 1./pow(Constants::twopi,5)/16.;
//     wgtsum += wgt;
//     wgtsq += sqr(wgt);
//   }
//   Energy2 mV2 = sqr(mV);
//   double volume;
//   if(mV>ZERO) 
//     volume = 1./64./pow(Constants::twopi,3)/sqr(sHat)*(sqr(sHat)-sqr(mV2)
// 						       -2.*sHat*mV2*log(sHat/mV2));
//   else
//     volume = 1./64./pow(Constants::twopi,3)/sqr(sHat)*sqr(sHat);
//   wgtsum /= double(npoint);
//   wgtsq  /= double(npoint);
//   wgtsq = max(wgtsq-sqr(wgtsum),0.)/double(npoint);
//   wgtsq = sqrt(wgtsq);
//   cerr << "testing volume A " << wgtsum << " " << volume << "\n";
//   cerr << "testing volume B " << wgtsum/volume << " +/- " << wgtsq/volume << "\n"; 
//   wgtsum = wgtsq = 0.;
//   for(unsigned int ix=0;ix<npoint;++ix) {
//     Energy m12 = mV+UseRandom::rnd()*(rootS-mV);
//     Energy2 m122 = sqr(m12);
//     Energy p1 = 0.5*(m122-mV2)/m12;
//     Energy p3 = 0.5*(sHat-m122)/rootS;
//     double wgt = 1./pow(Constants::twopi,5)/16./sHat/rootS*p1*p3*sqr(4.*Constants::pi)*(rootS-mV);
//     wgtsum += wgt;
//     wgtsq += sqr(wgt);
//   }
//   wgtsum /= double(npoint);
//   wgtsq  /= double(npoint);
//   wgtsq = max(wgtsq-sqr(wgtsum),0.)/double(npoint);
//   wgtsq = sqrt(wgtsq);
//   cerr << "testing volume C " << wgtsum << " " << volume << "\n";
//   cerr << "testing volume D " << wgtsum/volume << " +/- " << wgtsq/volume << "\n";
  
}

MEPP2VGammaJet::MEPP2VGammaJet()
  : process_(0),maxFlavour_(5),boson_(0),
    phase_(0),alphaS_(0.118),fixedAlphaS_(false)
{}


void MEPP2VGammaJet::getDiagrams() const {
  typedef std::vector<pair<tcPDPtr,tcPDPtr> > Pairvector;
  tcPDPtr wPlus  = getParticleData(ParticleID::Wplus );
  tcPDPtr wMinus = getParticleData(ParticleID::Wminus); 
  tcPDPtr z0     = getParticleData(ParticleID::Z0    );
  tcPDPtr gamma  = getParticleData(ParticleID::gamma);
  tcPDPtr g      = getParticleData(ParticleID::g);
  // W+/- gamma
  if(boson_==0||boson_==1) {
    // possible parents
    Pairvector parentpair;
    parentpair.reserve(6);
    // don't even think of putting 'break' in here!
    switch(maxFlavour_) {
    case 5:
      parentpair.push_back(make_pair(getParticleData(ParticleID::b),
				     getParticleData(ParticleID::cbar)));
      parentpair.push_back(make_pair(getParticleData(ParticleID::b), 
				     getParticleData(ParticleID::ubar)));
    case 4:
      parentpair.push_back(make_pair(getParticleData(ParticleID::s),
				     getParticleData(ParticleID::cbar)));
      parentpair.push_back(make_pair(getParticleData(ParticleID::d),
				     getParticleData(ParticleID::cbar)));
    case 3:
      parentpair.push_back(make_pair(getParticleData(ParticleID::s),
				     getParticleData(ParticleID::ubar)));
    case 2:
      parentpair.push_back(make_pair(getParticleData(ParticleID::d),
				     getParticleData(ParticleID::ubar)));
    default:
      ;
    }
    Pairvector::const_iterator parent = parentpair.begin();
    for (; parent != parentpair.end(); ++parent) {
      tcPDPtr qNeg1 = parent->first ;
      tcPDPtr qNeg2 = parent->second;
      tcPDPtr qPos1 = qNeg2->CC();
      tcPDPtr qPos2 = qNeg1->CC();
      if(process_==0 || process_ == 1) {
	assert(false);
      }
      if(process_==0 || process_ == 2) {
	assert(false);
      }
      if(process_==0 || process_ == 3) {
	assert(false);
      }
    }
  }
  // Z processes
  if(boson_==0||boson_==2) {
    for(int ix=1;ix<=int(maxFlavour_);++ix) {
      tcPDPtr qk = getParticleData(ix);
      tcPDPtr qb = qk->CC();
      if(process_== 0 || process_ == 1) {
	// Z + g + gamma
	add(new_ptr((Tree2toNDiagram(4), qk, qk, qk, qb, 1,    z0,
		     2, gamma, 3, g, -10)));
	add(new_ptr((Tree2toNDiagram(4), qk, qk, qk, qb, 1,    z0,
		     3, gamma, 2, g, -11)));
	add(new_ptr((Tree2toNDiagram(4), qk, qk, qk, qb, 3,    z0,
		     2, gamma, 1, g, -12)));
	add(new_ptr((Tree2toNDiagram(4), qk, qk, qk, qb, 3,    z0,
		     1, gamma, 2, g, -13)));
      }
      if(process_== 0 || process_ == 2) {
	// Z + q + gamma
	add(new_ptr((Tree2toNDiagram(4),qk,qk,qk,g,1,z0,2,gamma,3,qk, -20)));
	add(new_ptr((Tree2toNDiagram(4),qk,qk,qk,g,2,z0,1,gamma,3,qk, -21)));
	add(new_ptr((Tree2toNDiagram(3),qk,qk,g,1,z0,2,qk,5,gamma,5,qk,-22)));
      }
      // Z + qbar + gamma
      if(process_== 0 || process_ == 3) {
	add(new_ptr((Tree2toNDiagram(4),g,qb,qb,qb,3,z0,2,gamma,1,qb     ,-30)));
	add(new_ptr((Tree2toNDiagram(4),g,qb,qb,qb,2,z0,3,gamma,1,qb     ,-31)));
	add(new_ptr((Tree2toNDiagram(3),g,qb,qb   ,2,z0,1,qb,5,gamma,5,qb,-32)));
      }
    }
  }
}

bool MEPP2VGammaJet::generateKinematics(const double * r) {
  Energy pTgammaMin = lastCuts().minKT(mePartonData()[3]);
  Energy pTjetMin   = lastCuts().minKT(mePartonData()[4]);
  Energy rootS = sqrt(sHat());
  Energy mV = mePartonData()[2]->mass();
  jacobian(1.0);
  // standard 3 body phase-space
  if(phase_==0) {
    double r0 = r[0];
    if(r0<0.5) {
      r0 *=2.;
      // generate the photon energy
      Energy eMax = 0.5*(sHat()-sqr(mV))/rootS;
      if(eMax<pTgammaMin) return false;
      Energy eGamma = pTgammaMin+r0*(eMax-pTgammaMin);
      jacobian(jacobian()*(eMax-pTgammaMin)/rootS);
      // generate the jet energy
      Energy2 m122 = sHat()-2.*rootS*eGamma;
      Energy m12 = sqrt(m122);
      Energy E2star = 0.5*(m122+sqr(mV))/m12;
      Energy E3star = 0.5*(sHat()-m122)/m12;
      eMax = 
	0.5/rootS*(sHat()-sqr(mV)-2.*E3star*(E2star-sqrt(sqr(E2star)-sqr(mV))));
      Energy eMin = 
	0.5/rootS*(sHat()-sqr(mV)-2.*E3star*(E2star+sqrt(sqr(E2star)-sqr(mV))));
      eMin = max(eMin,pTjetMin);
      if(eMax<eMin) return false;
      Energy eJet = eMin+r[1]*(eMax-eMin);
      jacobian((eMax-eMin)/rootS*jacobian());
      // construct the momenta
      Energy eV = rootS-eGamma-eJet;
      if(eV<mV) return false;
      Energy pV = sqrt(sqr(eV)-sqr(mV));
      double cos3 = 0.5*(sqr(eJet)+sqr(eGamma)-sqr(pV))/eJet/eGamma;
      double sin3 = sqrt(1.-sqr(cos3));
      meMomenta()[2] = Lorentz5Momentum( eGamma*sin3,ZERO,-eJet+eGamma*cos3,eV    ,mV  );
      meMomenta()[3] = Lorentz5Momentum(-eGamma*sin3,ZERO,-eGamma*cos3     ,eGamma,ZERO);
      meMomenta()[4] = Lorentz5Momentum(ZERO,ZERO,eJet,eJet,ZERO);
      // finally rotate
      LorentzRotation rot;
      rot.rotateZ(r[2]*Constants::twopi);
      double ctmin = sqrt(1.-sqr(pTjetMin/eJet));
      rot.rotateX(acos(-ctmin+2.*ctmin*r[3]));
      rot.rotateZ(r[4]*Constants::twopi);
      for(unsigned int ix=2;ix<meMomenta().size();++ix)
	meMomenta()[ix] *= rot;
      jacobian(jacobian()*8.*sqr(Constants::pi)*ctmin);
    }
    else {
      r0 = 2.*(r0-0.5);
      // generate the jet energy
      Energy eMax = 0.5*(sHat()-sqr(mV))/rootS;
      if(eMax<pTjetMin) return false;
      Energy eJet = pTjetMin+r0*(eMax-pTjetMin);
      jacobian(jacobian()*(eMax-pTjetMin)/rootS);
      // generate the gamma energy
      Energy2 m122 = sHat()-2.*rootS*eJet;
      Energy m12 = sqrt(m122);
      Energy E2star = 0.5*(m122+sqr(mV))/m12;
      Energy E3star = 0.5*(sHat()-m122)/m12;
      eMax = 
	0.5/rootS*(sHat()-sqr(mV)-2.*E3star*(E2star-sqrt(sqr(E2star)-sqr(mV))));
      Energy eMin = 
	0.5/rootS*(sHat()-sqr(mV)-2.*E3star*(E2star+sqrt(sqr(E2star)-sqr(mV))));
      eMin = max(eMin,pTgammaMin);
      if(eMax<eMin) return false;
      Energy eGamma = eMin+r[1]*(eMax-eMin);
      jacobian((eMax-eMin)/rootS*jacobian());
      // construct the momenta
      Energy eV = rootS-eGamma-eJet;
      if(eV<mV) return false;
      Energy pV = sqrt(sqr(eV)-sqr(mV));
      double cos3 = 0.5*(sqr(eJet)+sqr(eGamma)-sqr(pV))/eJet/eGamma;
      double sin3 = sqrt(1.-sqr(cos3));
      meMomenta()[2] = Lorentz5Momentum( eJet*sin3,ZERO,-eGamma+eJet*cos3,eV    ,mV  );
      meMomenta()[4] = Lorentz5Momentum(-eJet*sin3,ZERO,-eJet*cos3     ,eJet,ZERO);
      meMomenta()[3] = Lorentz5Momentum(ZERO,ZERO,eGamma,eGamma,ZERO);
      // finally rotate
      LorentzRotation rot;
      rot.rotateZ(r[2]*Constants::twopi);
      double ctmin = sqrt(1.-sqr(pTgammaMin/eGamma));
      rot.rotateX(acos(-ctmin+2.*ctmin*r[3]));
      rot.rotateZ(r[4]*Constants::twopi);
      for(unsigned int ix=2;ix<meMomenta().size();++ix)
	meMomenta()[ix] *= rot;
      jacobian(jacobian()*8.*sqr(Constants::pi)*ctmin);
    }
  }
  else if(phase_==1) {
  }
  else
    assert(false);
  // check the cuts
  vector<LorentzMomentum> out(3);
  tcPDVector tout(3);
  for(unsigned int ix=0;ix<3;++ix) {
    out[ix]  = meMomenta()[2+ix];
    tout[ix] = mePartonData()[2+ix];
  }
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

double MEPP2VGammaJet::me2() const {
  alphaEM_ = SM().alphaEM();
  if(!fixedAlphaS_) alphaS_ = SM().alphaEM(scale());
  return alphaS_*alphaEM_*realME(mePartonData(),meMomenta())*16.*sqr(Constants::pi);
}

CrossSection MEPP2VGammaJet::dSigHatDR() const {
  return 1./pow(Constants::twopi,5)/16.*sqr(hbarc)*me2()*jacobian()/sHat();
}

Selector<MEBase::DiagramIndex>
MEPP2VGammaJet::diagrams(const DiagramVector & diags) const {
    Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) { 
    sel.insert(1., i);
  }
  return sel;
}

Selector<const ColourLines *>
MEPP2VGammaJet::colourGeometries(tcDiagPtr diag) const {
  // colour lines for q qbar -> V gamma g
  static const ColourLines cqqbarg[4]={ColourLines("1 2 3 7,-4 -7"),
				       ColourLines("1 2 7,-4 3 -7"),
				       ColourLines("1 7,-4 3 2 -7"),
				       ColourLines("1 2 7,-4 3 -7")};
  // colour lines for q g -> V gamma q
  static const ColourLines cqgq  [3]={ColourLines("1 2 3 -4,4 7"),
				      ColourLines("1 2 3 -4,4 7"),
				      ColourLines("1 2 -3,3 5 7")};
  // colour lines for gbar -> V gamma qbar
  static const ColourLines cqbargqbar[3]={ColourLines("1 -2 -3 -4,-1 -7"),
					  ColourLines("1 -2 -3 -4,-1 -7"),
					  ColourLines("1 -2 -3,-1 -5 -7")};
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
  case 10: case 11: case 12: case 13:
    sel.insert(1.0, &cqqbarg[abs(diag->id())-10]);
    break;
  case 20: case 21: case 22:
    sel.insert(1.0, &cqgq[abs(diag->id())-20]);
    break;
  case 30: case 31: case 32:
    sel.insert(1.0, &cqbargqbar[abs(diag->id())-30]);
    break;
  default:
    assert(false);
  }
  return sel;
}


IBPtr MEPP2VGammaJet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2VGammaJet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2VGammaJet::persistentOutput(PersistentOStream & os) const {
  os << process_ << maxFlavour_ << phase_ << alphaS_ << fixedAlphaS_ << alphaEM_
     << FFPvertex_ << FFWvertex_ << FFGvertex_ << FFZvertex_ << WWWvertex_;
}

void MEPP2VGammaJet::persistentInput(PersistentIStream & is, int) {
  is >> process_ >> maxFlavour_ >> phase_ >> alphaS_ >> fixedAlphaS_ >> alphaEM_
     >> FFPvertex_ >> FFWvertex_ >> FFGvertex_ >> FFZvertex_ >> WWWvertex_;
}

ClassDescription<MEPP2VGammaJet> MEPP2VGammaJet::initMEPP2VGammaJet;
// Definition of the static class description member.

void MEPP2VGammaJet::Init() {

  static ClassDocumentation<MEPP2VGammaJet> documentation
    ("There is no documentation for the MEPP2VGammaJet class");

  static Switch<MEPP2VGammaJet,unsigned int> interfaceBoson
    ("Boson",
     "Which electroweak vector bosons to include",
     &MEPP2VGammaJet::boson_, 0, false, false);
  static SwitchOption interfaceBosonAll
    (interfaceBoson,
     "All",
     "Include both W and Z",
     0);
  static SwitchOption interfaceBosonW
    (interfaceBoson,
     "W",
     "Only include W bosons",
     1);
  static SwitchOption interfaceBosonZ
    (interfaceBoson,
     "Z",
     "Only include Z",
     2);

  static Switch<MEPP2VGammaJet,unsigned int> interfaceProcess
    ("Process",
     "Processes to include",
     &MEPP2VGammaJet::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all processes",
     0);
  static SwitchOption interfaceProcessqqbar
    (interfaceProcess,
     "qqbar",
     "Only include q qbar -> V gamma g",
     1);
  static SwitchOption interfaceProcessqg
    (interfaceProcess,
     "qg",
     "Only include q g -> V gamma q",
     2);
  static SwitchOption interfaceProcessgqbar
    (interfaceProcess,
     "gqbar",
     "Only include g qbar -> V gamma qbar",
     3);

  static Parameter<MEPP2VGammaJet,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour allowed for the incoming quarks",
     &MEPP2VGammaJet::maxFlavour_, 5, 1, 5,
     false, false, Interface::limited);

  static Parameter<MEPP2VGammaJet,double> interfacealphaS
    ("AlphaS",
     "The value of alphaS to use if using a fixed alphaS",
     &MEPP2VGammaJet::alphaS_, 0.118, 0.0, 0.2,
     false, false, Interface::limited);

  static Switch<MEPP2VGammaJet,bool> interfaceFixedAlphaS
    ("FixedAlphaS",
     "Use a fixed value of alphaS",
     &MEPP2VGammaJet::fixedAlphaS_, false, false, false);
  static SwitchOption interfaceFixedAlphaSYes
    (interfaceFixedAlphaS,
     "Yes",
     "Use a fixed alphaS",
     true);
  static SwitchOption interfaceFixedAlphaSNo
    (interfaceFixedAlphaS,
     "No",
     "Use a running alphaS",
     false);

}

double MEPP2VGammaJet::realME(const cPDVector & particles,
			      const vector<Lorentz5Momentum> & momenta) const {
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<VectorWaveFunction> wout,pout,gout;
  SpinorWaveFunction    q_in;  
  SpinorBarWaveFunction qbar_in;
  VectorWaveFunction    g_out;  
  VectorWaveFunction    v_out  (momenta[2],particles[2],outgoing);
  VectorWaveFunction    p_out  (momenta[3],particles[3],outgoing);
  // q qbar -> V gamma g
  if(particles[4]->id()==ParticleID::g) {
    q_in    = SpinorWaveFunction    (momenta[0],particles[0],incoming);
    qbar_in = SpinorBarWaveFunction (momenta[1],particles[1],incoming);
    g_out   = VectorWaveFunction    (momenta[4],particles[4],outgoing);
  }
  // q g -> V gamma q
  else if(particles[4]->id()>0) {
    q_in    = SpinorWaveFunction    (momenta[0],particles[0],incoming);
    qbar_in = SpinorBarWaveFunction (momenta[4],particles[4],outgoing);
    g_out   = VectorWaveFunction    (momenta[1],particles[1],incoming);
  }
  else if(particles[4]->id()<0) {
    q_in    = SpinorWaveFunction    (momenta[4],particles[4],outgoing);
    qbar_in = SpinorBarWaveFunction (momenta[1],particles[1],incoming);
    g_out   = VectorWaveFunction    (momenta[0],particles[0],incoming);
  }
  else assert(false);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<2) {
      q_in.reset(ix);
      qin.push_back(q_in);
      qbar_in.reset(ix);
      qbarin.push_back(qbar_in);
      g_out.reset(2*ix);
      gout.push_back(g_out);
      p_out.reset(2*ix);
      pout.push_back(p_out);
    }
    v_out.reset(ix);
    wout.push_back(v_out);
  }
  vector<Complex> diag(particles[2]->id()==ParticleID::Z0 ? 6 : 8, 0.);
  AbstractFFVVertexPtr vertex = 
    particles[2]->id()==ParticleID::Z0 ? FFZvertex_ : FFWvertex_;
  Energy2 mu2 = scale();
  double sum(0.);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int whel=0;whel<3;++whel) {
	for(unsigned int phel=0;phel<2;++phel) {
	  for(unsigned int ghel=0;ghel<2;++ghel) {
	    //  Diagrams which are the same for W/Z gamma
	    // first diagram
	    SpinorWaveFunction inters1 = 
	      FFPvertex_->evaluate(ZERO,5,qin[ihel1].particle(),qin[ihel1],pout[phel]);
	    SpinorBarWaveFunction inters2 = 
	      vertex->evaluate(mu2,5,qin[ihel1].particle()->CC(),
			       qbarin[ihel2],wout[whel]);
	    diag[0] = FFGvertex_->evaluate(mu2,inters1,inters2,gout[ghel]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      FFGvertex_->evaluate(mu2,5,qin[ihel1].particle(),qin[ihel1],gout[ghel]);
	    SpinorBarWaveFunction inters4 = 
	      FFPvertex_->evaluate(ZERO,5,qbarin[ihel2].particle(),
				   qbarin[ihel2],pout[phel]);
	    diag[1] = vertex->evaluate(mu2,inters3,inters4,wout[whel]);
	    // fourth diagram
	    diag[2] = FFPvertex_->evaluate(ZERO,inters3,inters2,pout[phel]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      FFGvertex_->evaluate(mu2,5,qbarin[ihel2].particle(),
				   qbarin[ihel2],gout[ghel]);
	    diag[3] = 
	      vertex->evaluate(mu2,inters1,inters5,wout[whel]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      vertex->evaluate(mu2,5,qbarin[ihel2].particle()->CC(),
			       qin[ihel1],wout[whel]);
	    diag[4] = FFGvertex_->evaluate(mu2,inters6,inters4,gout[ghel]);
	    // eighth diagram
	    diag[5] = FFPvertex_->evaluate(ZERO,inters6,inters5,pout[phel]);
 	    //  Diagrams only for W gamma
	    if(particles[2]->id()!=ParticleID::Z0) {
	      // third diagram
	      VectorWaveFunction interv = 
		WWWvertex_->evaluate(mu2,3,particles[2]->CC(),pout[phel],wout[whel]);
	      diag[6] = vertex->evaluate(mu2,inters3,qbarin[ihel2],interv);
	      // seventh diagram
	      diag[7] = vertex->evaluate(mu2,qin[ihel1],inters5,interv);
	    }
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // divide out the em and strong couplings
  sum /= norm(FFGvertex_->norm()*FFPvertex_->norm());
  // final spin and colour factors spin = 1/4 colour = 4/9
  if(particles[4]->id()==ParticleID::g) sum /= 9.;
  // final spin and colour factors spin = 1/4 colour = 4/(3*8)
  else                                  sum /= 24.;
  return sum*UnitRemoval::InvE2*sHat();
}
