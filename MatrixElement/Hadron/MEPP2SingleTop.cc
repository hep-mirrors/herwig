// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2SingleTop class.
//

#include "MEPP2SingleTop.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

MEPP2SingleTop::MEPP2SingleTop() : process_(1), maxflavour_(5), topOption_(1),
				   wOption_(1)
{}

IBPtr MEPP2SingleTop::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2SingleTop::fullclone() const {
  return new_ptr(*this);
}

void MEPP2SingleTop::doinit() {
  HwMEBase::doinit();
  vector<unsigned int> massopt(2,0);
  massopt[0] = topOption_;
  if(process_==1) {
    massopt[1]= 0;
  }
  else if(process_==2) {
    massopt[1]= 1;
  }
  else if(process_==3) {
    massopt[1]= wOption_;
  }
  massOption(massopt);
  // mass option
  rescalingOption(2);
  // get the vertices we need
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " MEPP2SingleTop::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  FFWvertex_ = hwsm->vertexFFW();
  FFGvertex_ = hwsm->vertexFFG();
}

void MEPP2SingleTop::persistentOutput(PersistentOStream & os) const {
  os << process_ << FFWvertex_ << FFGvertex_ << maxflavour_
     << topOption_ << wOption_;
}

void MEPP2SingleTop::persistentInput(PersistentIStream & is, int) {
  is >> process_ >> FFWvertex_ >> FFGvertex_ >> maxflavour_
     >> topOption_ >> wOption_;
}

ClassDescription<MEPP2SingleTop> MEPP2SingleTop::initMEPP2SingleTop;
// Definition of the static class description member.

void MEPP2SingleTop::Init() {

  static ClassDocumentation<MEPP2SingleTop> documentation
    ("The MEPP2SingleTop class implements the production of "
     "a single top quark.");

  static Switch<MEPP2SingleTop,unsigned int> interfaceProcess
    ("Process",
     "The process to generate",
     &MEPP2SingleTop::process_, 1, false, false);
  static SwitchOption interfaceProcesstChannel
    (interfaceProcess,
     "tChannel",
     "Generate t-channel W exchange processes",
     1);
  static SwitchOption interfaceProcesssChannel
    (interfaceProcess,
     "sChannel",
     "Generate s-channel W exchange processes",
     2);
  static SwitchOption interfaceProcesstW
    (interfaceProcess,
     "tW",
     "Generate top production in association with a W boson",
     3);

  static Parameter<MEPP2SingleTop,int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the non-top quarks",
     &MEPP2SingleTop::maxflavour_, 5, 1, 5,
     false, false, Interface::limited);

  static Switch<MEPP2SingleTop,int> interfaceTopMassOption
    ("TopMassOption",
     "Option for the treatment of the top quark masses",
     &MEPP2SingleTop::topOption_, 1, false, false);
  static SwitchOption interfaceTopMassOptionOnMassShell
    (interfaceTopMassOption,
     "OnMassShell",
     "The top is produced on its mass shell",
     1);
  static SwitchOption interfaceTopMassOptionOffShell
    (interfaceTopMassOption,
     "OffShell",
     "The top is generated off-shell using the mass and width generator.",
     2);

  static Switch<MEPP2SingleTop,int> interfaceWMassOption
    ("WMassOption",
     "Option for the treatment of the top quark masses",
     &MEPP2SingleTop::wOption_, 1, false, false);
  static SwitchOption interfaceWMassOptionOnMassShell
    (interfaceWMassOption,
     "OnMassShell",
     "The W is produced on its mass shell",
     1);
  static SwitchOption interfaceWMassOptionOffShell
    (interfaceWMassOption,
     "OffShell",
     "The W is generated off-shell using the mass and width generator.",
     2);

}

Energy2 MEPP2SingleTop::scale() const{
  if(process_==2)
    return sHat();
  Energy2 s(0.5*sHat()),u(0.5*(uHat()-meMomenta()[3].mass2())),
    t(0.5*(tHat()-meMomenta()[2].mass2()));
  return 4.*s*t*u/(s*s+t*t+u*u);
}

unsigned int MEPP2SingleTop::orderInAlphaS() const {
  return process_!=3 ? 0 : 1;
}

unsigned int MEPP2SingleTop::orderInAlphaEW() const {
  return process_!=3 ? 2 : 1;
}

void MEPP2SingleTop::getDiagrams() const {
  tcPDPtr g  = getParticleData(ParticleID::g);
  tcPDPtr Wp = getParticleData(ParticleID::Wplus );
  tcPDPtr Wm = getParticleData(ParticleID::Wminus);
  tcPDPtr tp = getParticleData(ParticleID::t);
  tcPDPtr tb = tp->CC();
  // light particles
  typedef std::vector<pair<tcPDPtr,tcPDPtr> > Pairvector;
  Pairvector lightPair;
  switch(maxflavour_) {
  case 5:
    lightPair.push_back(make_pair(getParticleData(ParticleID::b),
				  getParticleData(ParticleID::cbar)));
    lightPair.push_back(make_pair(getParticleData(ParticleID::b),
				  getParticleData(ParticleID::ubar)));
    [[fallthrough]];
  case 4:
    lightPair.push_back(make_pair(getParticleData(ParticleID::s), 
				  getParticleData(ParticleID::cbar)));
    lightPair.push_back(make_pair(getParticleData(ParticleID::d), 
				  getParticleData(ParticleID::cbar)));
    [[fallthrough]];
  case 3:
    lightPair.push_back(make_pair(getParticleData(ParticleID::s),
				  getParticleData(ParticleID::ubar)));
    [[fallthrough]];
  case 2:
    lightPair.push_back(make_pair(getParticleData(ParticleID::d),
				  getParticleData(ParticleID::ubar)));
    [[fallthrough]];
  default:
    ;
  }
  // top partner
  vector<tcPDPtr> topPartner;
  for(long ix=1;ix<=maxflavour_;ix+=2)
    topPartner.push_back(getParticleData(ix));
  // t-channel
  if(process_==1) {
    for (Pairvector::const_iterator light = lightPair.begin();
	 light != lightPair.end(); ++light) {
      // lights
      tcPDPtr qNeg1 = light->first;
      tcPDPtr qNeg2 = light->second->CC();
      tcPDPtr qPos1 = light->first ->CC();
      tcPDPtr qPos2 = light->second;
      if(SM().CKM(*qNeg2,*qNeg1)<1e-30) continue;
      for(unsigned int ix=0;ix<topPartner.size();++ix) {
	if(SM().CKM(*tp,*topPartner[ix])<1e-30) continue;
	// diagrams
	add(new_ptr((Tree2toNDiagram(3), topPartner[ix]      ,
		     Wp, qNeg2, 1, tp, 2, qNeg1, -1)));
	add(new_ptr((Tree2toNDiagram(3), topPartner[ix]      ,
		     Wp, qPos1, 1, tp, 2, qPos2, -2)));
	add(new_ptr((Tree2toNDiagram(3), topPartner[ix]->CC(),
		     Wm, qNeg1, 1, tb, 2, qNeg2, -3)));
	add(new_ptr((Tree2toNDiagram(3), topPartner[ix]->CC(),
		     Wm, qPos2, 1, tb, 2, qPos1, -4)));
      }
    }
  }
  else if(process_==2) {
    Pairvector::const_iterator light = lightPair.begin();
    for (; light != lightPair.end(); ++light) {
      // lights
      tcPDPtr qNeg1 = light->first ;
      tcPDPtr qNeg2 = light->second;
      tcPDPtr qPos1 = qNeg2->CC();
      tcPDPtr qPos2 = qNeg1->CC();
      if(SM().CKM(*qPos1,*qNeg1)<1e-30) continue;
      for(unsigned int ix=0;ix<topPartner.size();++ix) {
	if(SM().CKM(*tp,*topPartner[ix])<1e-30) continue;
	// diagrams
	add(new_ptr((Tree2toNDiagram(2), qNeg1, qNeg2, 1, Wm, 
		     3, tb, 3, topPartner[ix]      , -11)));
	add(new_ptr((Tree2toNDiagram(2), qPos1, qPos2, 1, Wp, 
		     3, tp, 3, topPartner[ix]->CC(), -12)));
      }
    }
  }
  else if(process_==3) {
    for(unsigned int ix=0;ix<topPartner.size();++ix) {
      if(SM().CKM(*tp,*topPartner[ix])<1e-30) continue;
      add(new_ptr((Tree2toNDiagram(2), topPartner[ix], g,
		   1, topPartner[ix], 3, tp, 3, Wm, -21)));
      add(new_ptr((Tree2toNDiagram(3), topPartner[ix], tp, g,
		   3, tp, 1, Wm, -22)));
      add(new_ptr((Tree2toNDiagram(2), topPartner[ix]->CC(), g,
		   1, topPartner[ix]->CC(), 3, tb, 3, Wp, -23)));
      add(new_ptr((Tree2toNDiagram(3), topPartner[ix]->CC(), tb, g,
		   3, tb, 1, Wp, -24)));
    }
  }
  else
    assert(false);
}

Selector<MEBase::DiagramIndex>
MEPP2SingleTop::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    // only one diagram for t- and s-channel processes
    if      ( abs(diags[i]->id()) <20 )
      sel.insert(1.0, i);
    else if ( abs(diags[i]->id()) <23 )
      sel.insert(meInfo()[abs(diags[i]->id())-21],i);
    else 
      sel.insert(meInfo()[abs(diags[i]->id())-23],i);
  }
  return sel;
}

Selector<const ColourLines *>
MEPP2SingleTop::colourGeometries(tcDiagPtr diag) const {
  // t-channel
  static ColourLines ct[4]={ColourLines(" 1  4,  3  5"),
			    ColourLines(" 1  4, -3 -5"),
			    ColourLines("-1 -4,  3  5"),
			    ColourLines("-1 -4, -3 -5")};
  // s-channel
  static ColourLines cs[2]={ColourLines("1 -2, -4  5"),
			    ColourLines("1 -2,  4 -5")};
  // tW
  static ColourLines ctw[4]={ColourLines(" 1  -2,  2   3   4"),
			     ColourLines(" 1   2  -3,  3   4"),
			     ColourLines("-1   2, -2  -3  -4"),
			     ColourLines("-1  -2   3, -3  -4")};
  // select the right one
  Selector<const ColourLines *> sel;
  int id = abs(diag->id()); 
  switch (id) {
  case 1: case 2: case 3: case 4:
    sel.insert(1.0, &ct[id-1]);
    break;
  case 11: case 12:
    sel.insert(1.0, &cs[id-11]);
    break;
  case 21: case 22: case 23: case 24:
    sel.insert(1.0, &ctw[id-21]);
    break;
  default:
    assert(false);
  };
  return sel;
}

double MEPP2SingleTop::me2() const {
  if(process_==1) {
    vector<SpinorWaveFunction>    f1,f2;
    vector<SpinorBarWaveFunction> a1,a2;
    SpinorWaveFunction    sf1,sf2;
    SpinorBarWaveFunction sa1,sa2;
    if(mePartonData()[0]->id()>0) {
      sf1 = SpinorWaveFunction   (meMomenta()[0],mePartonData()[0],incoming);
      sa1 = SpinorBarWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
    }
    else {
      sa1 = SpinorBarWaveFunction(meMomenta()[0],mePartonData()[0],incoming);
      sf1 = SpinorWaveFunction   (meMomenta()[2],mePartonData()[2],outgoing);
    }
    if(mePartonData()[1]->id()>0) {
      sf2 = SpinorWaveFunction   (meMomenta()[1],mePartonData()[1],incoming);
      sa2 = SpinorBarWaveFunction(meMomenta()[3],mePartonData()[3],outgoing);
    }
    else {
      sa2 = SpinorBarWaveFunction(meMomenta()[1],mePartonData()[1],incoming);
      sf2 = SpinorWaveFunction   (meMomenta()[3],mePartonData()[3],outgoing);
    }
    for(unsigned int ix=0;ix<2;++ix) {
      sf1.reset(ix); f1.push_back(sf1);
      sf2.reset(ix); f2.push_back(sf2);
      sa1.reset(ix); a1.push_back(sa1);
      sa2.reset(ix); a2.push_back(sa2);
    }
    return tChannelME(f1,a1,f2,a2,false);
  }
  else if(process_==2) {
    vector<SpinorWaveFunction>    fin,aout;
    vector<SpinorBarWaveFunction> ain,fout;
    SpinorWaveFunction    q   (meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction qbar(meMomenta()[1],mePartonData()[1],incoming);
    SpinorBarWaveFunction f;
    SpinorWaveFunction    fbar;
    if(mePartonData()[2]->id()==ParticleID::t) {
      f    = SpinorBarWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
      fbar = SpinorWaveFunction   (meMomenta()[3],mePartonData()[3],outgoing);
    }
    else {
      f    = SpinorBarWaveFunction(meMomenta()[3],mePartonData()[3],outgoing);
      fbar = SpinorWaveFunction   (meMomenta()[2],mePartonData()[2],outgoing);
    }
    for(unsigned int ix=0;ix<2;++ix) {
      q.reset(ix)   ; fin.push_back(q);
      qbar.reset(ix); ain.push_back(qbar);
      f.reset(ix)   ;fout.push_back(f);
      fbar.reset(ix);aout.push_back(fbar);
    }
    return sChannelME(fin,ain,fout,aout,false);
  }
  else {
    SpinorWaveFunction    sf;
    SpinorBarWaveFunction sa;
    vector<SpinorWaveFunction>    f1;
    vector<SpinorBarWaveFunction> a1;
    vector<VectorWaveFunction> gl,vec;
    if(mePartonData()[0]->id()>0) {
      sf = SpinorWaveFunction   (meMomenta()[0],mePartonData()[0],incoming);
      sa = SpinorBarWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
    }
    else {
      sa = SpinorBarWaveFunction(meMomenta()[0],mePartonData()[0],incoming);
      sf = SpinorWaveFunction   (meMomenta()[2],mePartonData()[2],outgoing);
    }
    VectorWaveFunction gluon(meMomenta()[1],mePartonData()[1],incoming);
    VectorWaveFunction wbos (meMomenta()[3],mePartonData()[3],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      sf.reset(ix); f1.push_back(sf);
      sa.reset(ix); a1.push_back(sa);
      gluon.reset(2*ix); gl.push_back(gluon);
      wbos.reset(ix); vec.push_back(wbos);
    }
    wbos.reset(2);
    vec.push_back(wbos);
    return tWME(f1,gl,a1,vec,false);
  }
}

double MEPP2SingleTop::sChannelME(vector<SpinorWaveFunction>    & fin ,
				  vector<SpinorBarWaveFunction> & ain ,
				  vector<SpinorBarWaveFunction> & fout,
				  vector<SpinorWaveFunction>    & aout,
				  bool calc) const {
  // matrix element to be stored
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // positive or negative W boson
  bool positive = mePartonData()[0]->iCharge() + mePartonData()[1]->iCharge() > 0;
  tcPDPtr wBoson = positive ? 
    getParticleData(ParticleID::Wplus) : getParticleData(ParticleID::Wminus);
  // sum over helicities to get the matrix element
  double me = 0.;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      VectorWaveFunction inter = FFWvertex_->evaluate(scale(),1,wBoson,fin[ihel1],ain[ihel2]);
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  Complex  diag = FFWvertex_->evaluate(scale(),aout[ohel2],fout[ohel1],inter);
	  // sum over helicities
	  me += norm(diag);
	  if(calc) {
	    if(fout[ohel1].id()==ParticleID::t)
	      newme(ihel1,ihel2,ohel1,ohel2) = diag;
	    else
	      newme(ihel1,ihel2,ohel2,ohel1) = diag;
	  }
	}
      }
    }
  }
  // results
  // spin and colour factor
  double colspin=1./12.;
  if(abs(fout[0].id())<=6) colspin*=3.;
  if(calc) me_.reset(newme);
  return me*colspin;
}

double MEPP2SingleTop::tChannelME(vector<SpinorWaveFunction>    & f1,
				  vector<SpinorBarWaveFunction> & a1,
				  vector<SpinorWaveFunction>    & f2,
				  vector<SpinorBarWaveFunction> & a2,
				  bool calc) const {
  // matrix element to be stored
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  tcPDPtr wBoson = getParticleData(ParticleID::Wplus);
  // sum over helicities to get the matrix element
  double me = 0.;
  unsigned int ihel[4];
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ohel1=0;ohel1<2;++ohel1) {
      VectorWaveFunction inter = FFWvertex_->evaluate(scale(),3,wBoson,
						      f1[ihel1],a1[ohel1]);
      ihel[0] = ihel1;
      ihel[2] = ohel1;
      if(f1[ihel1].direction()==outgoing) swap(ihel[0],ihel[2]);
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  Complex  diag = FFWvertex_->evaluate(scale(),f2[ihel2],a2[ohel2],inter);
	  // sum over helicities
	  me += norm(diag);
	  ihel[1] = ihel2;
	  ihel[3] = ohel2;
	  if(f2[ihel2].direction()==outgoing) swap(ihel[1],ihel[3]);
	  if(calc) 
	    newme(ihel[0],ihel[1],ihel[2],ihel[3]) = diag;
	}
      }
    }
  }
  // results
  // spin and colour factor
  double colspin = 1./4.;
  if(calc) me_.reset(newme);
  // Energy mt = getParticleData(ParticleID::t)->mass();
  // Energy mw = getParticleData(ParticleID::Wplus)->mass();
  // double test1 
  //   = 0.25*sqr(4.*Constants::pi*SM().alphaEM(scale())/SM().sin2ThetaW())*
  //   sHat()*(sHat()-sqr(mt))/sqr(tHat()-sqr(mw));
  // double test2 
  //   = 0.25*sqr(4.*Constants::pi*SM().alphaEM(scale())/SM().sin2ThetaW())*
  //   uHat()*(uHat()-sqr(mt))/sqr(tHat()-sqr(mw));
  return me*colspin;
}

double MEPP2SingleTop::tWME(vector<SpinorWaveFunction> & fin,
			   vector<VectorWaveFunction> & gin,
			   vector<SpinorBarWaveFunction> & fout,
			   vector<VectorWaveFunction> & Wout,
			   bool calc) const {
  if(calc) me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
					     PDT::Spin1Half,PDT::Spin1));
  double me[3]={0.,0.,0.};
  Complex diag[2];
  Energy2 mt=scale();
  for(unsigned int ihel=0;ihel<2;++ihel) {
    for(unsigned int ghel=0;ghel<2;++ghel) {
      for(unsigned int ohel=0;ohel<2;++ohel) {
	for(unsigned int whel=0;whel<3;++whel) {
	  if(fin[ihel].direction()==incoming) {
	    // 1st diagram
	    SpinorWaveFunction inters = 
	      FFGvertex_->evaluate(mt,5,mePartonData()[0],
				   fin[ihel],gin[ghel]);
	    diag[0]= FFWvertex_->evaluate(mt,inters,fout[ohel],Wout[whel]);
	    // 2nd diagram
	    SpinorBarWaveFunction interb = 
	      FFGvertex_->evaluate(mt,3,mePartonData()[2],
				   fout[ohel],gin[ghel]);
	    diag[1]= FFWvertex_->evaluate(mt,fin[ihel],interb,Wout[whel]);
	  }
	  else {
	    // 1st diagram
	    SpinorBarWaveFunction interb =
	      FFGvertex_->evaluate(mt,5,mePartonData()[0],
				   fout[ihel],gin[ghel]);
 	    diag[0]= FFWvertex_->evaluate(mt,fin[ohel],interb,Wout[whel]);
	    // 2nd diagram
	    SpinorWaveFunction inters =
	      FFGvertex_->evaluate(mt,3,mePartonData()[2],
				   fin[ohel],gin[ghel]);
	    diag[1]= FFWvertex_->evaluate(mt,inters,fout[ihel],Wout[whel]);
	  }
	  // diagram contributions
	  me[1] += norm(diag[0]);
	  me[2] += norm(diag[1]);
	  // total
	  diag[0] += diag[1];
	  me[0]   += norm(diag[0]);
	  if(calc) 
	    me_(ihel,ghel,ohel,whel) = diag[0];
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin=1./24./4.;
  // and C_F N_c from matrix element
  colspin *=4.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix] *= colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0];
}

void MEPP2SingleTop::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  if(mePartonData()[0]->id()!=hard[0]->id()) {
    swap(hard[0],hard[1]);
  }
  if(mePartonData()[2]->id()!=hard[2]->id()) {
    swap(hard[2],hard[3]);
  }
  // on-shell for matrix element
  vector<Lorentz5Momentum> momenta;
  cPDVector data;
  for(unsigned int ix=0;ix<4;++ix) {
    momenta.push_back(hard[ix]->momentum());
    data   .push_back(hard[ix]->dataPtr());
  }
  rescaleMomenta(momenta,data);
  if(process_==1) {
    vector<SpinorWaveFunction>    f1,f2;
    vector<SpinorBarWaveFunction> a1,a2;
    // off-shell wavefunctions for the spin correlations (and on shell for matrix elements)
    SpinorWaveFunction    sf1,sf2;
    SpinorBarWaveFunction sa1,sa2;
    if(hard[0]->id()>0) {
      SpinorWaveFunction   (f1,hard[0],incoming,false,true);
      SpinorBarWaveFunction(a1,hard[2],outgoing,true,true);
      sf1 = SpinorWaveFunction   (rescaledMomenta()[0],data[0],incoming);
      sa1 = SpinorBarWaveFunction(rescaledMomenta()[2],data[2],outgoing);
    }
    else {
      SpinorBarWaveFunction(a1,hard[0],incoming,false,true);
      SpinorWaveFunction   (f1,hard[2],outgoing,true,true);
      sa1 = SpinorBarWaveFunction(rescaledMomenta()[0],data[0],incoming);
      sf1 = SpinorWaveFunction   (rescaledMomenta()[2],data[2],outgoing);
    }
    if(hard[1]->id()>0) {
      SpinorWaveFunction   (f2,hard[1],incoming,false,true);
      SpinorBarWaveFunction(a2,hard[3],outgoing,true,true);
      sf2 = SpinorWaveFunction   (rescaledMomenta()[1],data[1],incoming);
      sa2 = SpinorBarWaveFunction(rescaledMomenta()[3],data[3],outgoing);
    }
    else {
      SpinorBarWaveFunction(a2,hard[1],incoming,false,true);
      SpinorWaveFunction   (f2,hard[3],outgoing,true,true);
      sa2 = SpinorBarWaveFunction(rescaledMomenta()[1],data[1],incoming);
      sf2 = SpinorWaveFunction   (rescaledMomenta()[3],data[3],outgoing);
    }
    for(unsigned int ix=0;ix<2;++ix) {
      sf1.reset(ix); f1[ix] = sf1;
      sf2.reset(ix); f2[ix] = sf2;
      sa1.reset(ix); a1[ix] = sa1;
      sa2.reset(ix); a2[ix] = sa2;
    }
    tChannelME(f1,a1,f2,a2,true);
  }
  else if(process_==2) {
    vector<SpinorWaveFunction>    fin,aout;
    vector<SpinorBarWaveFunction> ain,fout;
    SpinorWaveFunction   (fin,hard[0],incoming,false,true);
    SpinorBarWaveFunction(ain,hard[1],incoming,false,true);
    SpinorWaveFunction    q   (rescaledMomenta()[0],data[0],incoming);
    SpinorBarWaveFunction qbar(rescaledMomenta()[1],data[1],incoming);
    SpinorBarWaveFunction f;
    SpinorWaveFunction    fbar;
    if(hard[2]->id()==ParticleID::t) {
      SpinorBarWaveFunction(fout,hard[2],outgoing,true,true);
      SpinorWaveFunction   (aout,hard[3],outgoing,true,true);
      f    = SpinorBarWaveFunction(rescaledMomenta()[2],data[2],outgoing);
      fbar = SpinorWaveFunction   (rescaledMomenta()[3],data[3],outgoing);
    }
    else {
      SpinorBarWaveFunction(fout,hard[3],outgoing,true,true);
      SpinorWaveFunction   (aout,hard[2],outgoing,true,true);
      f    = SpinorBarWaveFunction(rescaledMomenta()[3],data[3],outgoing);
      fbar = SpinorWaveFunction   (rescaledMomenta()[2],data[2],outgoing);
    }
    for(unsigned int ix=0;ix<2;++ix) {
      q.reset(ix)   ; fin[ix] = q;
      qbar.reset(ix); ain[ix] = qbar;
      f.reset(ix)   ;fout[ix] = f;
      fbar.reset(ix);aout[ix] = fbar;
    }
    sChannelME(fin,ain,fout,aout,true);
  }
  else {
    SpinorWaveFunction    sf;
    SpinorBarWaveFunction sa;
    vector<SpinorWaveFunction>    f1;
    vector<SpinorBarWaveFunction> a1;
    vector<VectorWaveFunction> gl,vec;
    if(hard[0]->id()>0) {
      SpinorWaveFunction   (f1,hard[0],incoming,false,true);
      SpinorBarWaveFunction(a1,hard[2],outgoing,true ,true);
      sf = SpinorWaveFunction   (rescaledMomenta()[0],data[0],incoming);
      sa = SpinorBarWaveFunction(rescaledMomenta()[2],data[2],outgoing);
    }
    else {
      SpinorBarWaveFunction(a1,hard[0],incoming,false,true);
      SpinorWaveFunction   (f1,hard[2],outgoing,true ,true);
      sa = SpinorBarWaveFunction(rescaledMomenta()[0],data[0],incoming);
      sf = SpinorWaveFunction   (rescaledMomenta()[2],data[2],outgoing);
    }
    VectorWaveFunction(gl,hard[1],incoming,false,true ,true);
    VectorWaveFunction(vec,hard[3],outgoing,true,false,true);
    VectorWaveFunction gluon(rescaledMomenta()[1],data[1],incoming);
    VectorWaveFunction wbos (rescaledMomenta()[3],data[3],outgoing);
    gl[1]=gl[2];
    for(unsigned int ix=0;ix<2;++ix) {
      sf.reset(ix); f1[ix] = sf;
      sa.reset(ix); a1[ix] = sa;
      gluon.reset(2*ix); gl[ix] = gluon;
      wbos.reset(ix); vec[ix] = wbos;
    }
    wbos.reset(2);
    vec[2] = wbos;
    tWME(f1,gl,a1,vec,true);
  }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) 
    hard[ix]->spinInfo()->productionVertex(hardvertex);
}
