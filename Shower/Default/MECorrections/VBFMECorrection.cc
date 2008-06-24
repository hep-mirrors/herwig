// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VBFMECorrection class.
//

#include "VBFMECorrection.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <numeric>
#include "Herwig++/Utilities/Maths.h"

using namespace Herwig;

void VBFMECorrection::persistentOutput(PersistentOStream & os) const {
  os << _procprob << _comptonint;
}

void VBFMECorrection::persistentInput(PersistentIStream & is, int) {
  is >> _procprob >> _comptonint;
}

ClassDescription<VBFMECorrection> VBFMECorrection::initVBFMECorrection;
// Definition of the static class description member.

void VBFMECorrection::Init() {

  static ClassDocumentation<VBFMECorrection> documentation
    ("The VBFMECorrection class implements the matrix element"
     " correction for VBF Higgs production");

  static Parameter<VBFMECorrection,double> interfaceProcessProbability
    ("ProcessProbability",
     "The probabilty of the QCD compton process for the process selection",
     &VBFMECorrection::_procprob, 0.3, 0.0, 1.,
     false, false, Interface::limited);

}

bool VBFMECorrection::canHandle(ShowerTreePtr tree,double & initial,
				double & final,EvolverPtr evolver) {
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // three outgoing particles
  if(tree->outgoingLines().size()!=3) return false;
  // extract the incoming quarks
  vector<PPtr> incoming;
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) 
      incoming.push_back(cit->first->progenitor());
  }
  if(incoming.size()!=2) return false;
  // extract the outgoing quarks and the higgs
  bool higgs = false;
  vector<PPtr> outgoing;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()==ParticleID::h0) higgs = true;
    else if(QuarkMatcher::Check(cjt->first->progenitor()->data()))
      outgoing.push_back(cjt->first->progenitor());
  }
  if(!higgs||outgoing.size()!=2) return false;
  // check things match up
  unsigned int nmatch(0);
  for(unsigned int ix=0;ix<incoming.size();++ix) {
    if(incoming[ix]->colourLine()) {
      for(unsigned int iy=0;iy<outgoing.size();++iy) {
	if(outgoing[iy]->colourLine()==incoming[ix]->colourLine()) {
	  ++nmatch;
	  break;
	}
      }
    }
    else {
      for(unsigned int iy=0;iy<outgoing.size();++iy) {
	if(outgoing[iy]->antiColourLine()==incoming[ix]->antiColourLine()) {
	  ++nmatch;
	  break;
	}
      }
    }
  }
  if(nmatch!=2) return false;
  // can handle it
  initial = 1.;
  final   = 1.;
  return true;
}

void VBFMECorrection::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  ++_ntry;
  //cerr << "testing in apply hard\n";
  //cerr << *generator()->currentEvent() << "\n";
  vector<tChannelPair> systems;
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) {
      systems.push_back(tChannelPair());
      systems.back().hadron   = cit->first->original()->parents()[0];
      systems.back().beam     = cit->first->beam();
      systems.back().incoming = cit->first->progenitor();
      systems.back().pdf      = systems.back().beam->pdf();
    }
  }
  assert(systems.size()==2);
  // extract the outgoing quarks and the higgs
  PPtr higgs;
  vector<ShowerParticlePtr> outgoing;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()==ParticleID::h0) 
      higgs = cjt->first->progenitor();
    else if(QuarkMatcher::Check(cjt->first->progenitor()->data()))
      outgoing.push_back(cjt->first->progenitor());
  }
  assert(outgoing.size()==2&&higgs);
  // match up the quarks
  for(unsigned int ix=0;ix<systems.size();++ix) {
    if(systems[ix].incoming->colourLine()) {
      for(unsigned int iy=0;iy<outgoing.size();++iy) {
	if(outgoing[iy]->colourLine()==systems[ix].incoming->colourLine()) {
	  systems[ix].outgoing=outgoing[iy];
	  break;
	}
      }
    }
    else {
      for(unsigned int iy=0;iy<outgoing.size();++iy) {
	if(outgoing[iy]->antiColourLine()==systems[ix].incoming->antiColourLine()) {
	  systems[ix].outgoing=outgoing[iy];
	  break;
	}
      }
    }
  }
  assert(systems[0].outgoing&&systems[1].outgoing);
//   for(unsigned int ix=0;ix<systems.size();++ix) {
//     cerr << *systems[ix].hadron;
//     cerr << *systems[ix].incoming << "\n" << *systems[ix].outgoing << "\n";
//     cerr << "testing scale " << sqrt(-(systems[ix].incoming->momentum()-systems[ix].outgoing->momentum()).m2()/GeV2) << "\n";
//   }
  // select emitting line
  if(UseRandom::rndbool()) swap(systems[0],systems[1]);
  // extract the born variables
  Lorentz5Momentum q = systems[0].outgoing->momentum()-
                       systems[0].incoming->momentum();
  // extract born variables
  Energy2 q2 = -q.m2();
  Energy Q = sqrt(q2);
  double xB = systems[0].incoming->x();
//   cerr << "testing q2 " << q2/GeV2 << " " << Q/GeV << " " << xB << "\n";
  // construct lorentz transform from lab to breit frame
  Lorentz5Momentum phadron =  systems[0].hadron->momentum();
  phadron.setMass(0.*GeV);
  phadron.rescaleEnergy();
  Lorentz5Momentum pcmf = phadron+0.5/xB*q;
  pcmf.rescaleMass();
  LorentzRotation rot(-pcmf.boostVector());
  Lorentz5Momentum pbeam = rot*phadron;
  Axis axis(pbeam.vect().unit());
  double sinth(sqrt(1.-sqr(axis.z())));
  rot.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  Lorentz5Momentum pout = rot*(systems[1].outgoing->momentum()+higgs->momentum());
  rot.rotateZ(-atan2(pout.y(),pout.x()));
//   cerr << "testing momenta in the Breit frame\n";
//   cerr << rot*systems[0].incoming->momentum()/GeV << "\n";
//   cerr << rot*systems[0].outgoing->momentum()/GeV << "\n";
//   cerr << rot*systems[1].incoming->momentum()/GeV << "\n";
//   cerr << rot*systems[1].outgoing->momentum()/GeV << "\n";
//   cerr << rot*higgs->momentum()/GeV << "\n";
//   cerr << rot*(systems[1].outgoing->momentum()+higgs->momentum())/GeV << "\n";
  // calculate the A coefficient for the correlations
  double acoeff = A(systems[0].incoming->dataPtr(),systems[0].outgoing->dataPtr(),
		    systems[1].incoming->dataPtr(),systems[1].outgoing->dataPtr());
  vector<double> azicoeff;
  // select the type of process
  bool BGF = UseRandom::rnd()>_procprob;
  tcPDPtr gluon = getParticleData(ParticleID::g);
  double wgt,xp,zp,x1,x2,x3,xperp;
  LorentzVector<double> l(2.*(rot*systems[1].incoming->momentum())/Q);
  LorentzVector<double> m(2.*(rot*systems[1].outgoing->momentum())/Q);
  //cerr << "testing scale " << Q/2/GeV << "\n";
  //cerr << "testing other momentum " << rot*systems[1].incoming->momentum()/GeV << "\n";
  //cerr << "testing other momentum " << rot*systems[1].outgoing->momentum()/GeV << "\n";
  //cerr << "testing l " << l << "\n";
  //cerr << "testing m " << m << "\n";
  if(!BGF) {
    wgt = generateComptonPoint(xp,zp);
    if(xp<1e-6) return;
    // common pieces
    Energy2 scale = q2*((1.-xp)*(1-zp)*zp/xp+1);
    //cerr << "testing weight A1 " << wgt << "\n";
    wgt *= 2./3./Constants::pi*coupling()->value(scale)/_procprob;
    //cerr << "testing weight B1 " << wgt << "\n";
    // PDF piece
    wgt *= systems[0].pdf->xfx(systems[0].beam,systems[0].incoming->dataPtr(),scale,xB/xp)/
           systems[0].pdf->xfx(systems[0].beam,systems[0].incoming->dataPtr(),q2   ,xB);
    //cerr << "testing weight C1 " << wgt << "\n";
    // numerator factors
    wgt /= (1.-xp)*(1.-zp);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = ComptonME(xp,x2,xperp,acoeff,l,m);



//     Energy2 dp1p2=systems[0].outgoing->momentum()*systems[1].outgoing->momentum();
//     Energy2 dk1k2=systems[0].incoming->momentum()*systems[1].incoming->momentum();
//     Energy2 dk1p2=systems[0].incoming->momentum()*systems[1].outgoing->momentum();
//     Energy2 dk2p1=systems[0].outgoing->momentum()*systems[1].incoming->momentum();
//     Energy4 meb = 8*(dp1p2*dk1k2*(1.+0.5*acoeff)+dk1p2*dk2p1*(1.-0.5*acoeff));
//     cerr << "testing acoeff " << acoeff << "\n";
//     double denom = -l.z()*m.z()+l.t()*m.t()+0.5*acoeff*(l.t()*m.z()-l.z()*m.t());
//     cerr << "testing ME " << meb/pow<4,1>(Q) << " " << denom << " " 
// 	 << meb/pow<4,1>(Q)/denom << "\n";
  }
  else {
    wgt = generateBGFPoint(xp,zp);
    if(xp<1e-6) return;
    // common pieces 
    Energy2 scale = q2*((1.-xp)*(1-zp)*zp/xp+1);
    wgt *= 0.25/Constants::pi*coupling()->value(scale)/(1.-_procprob);
    // PDF piece
    wgt *= systems[0].pdf->xfx(systems[0].beam,gluon                         ,scale,xB/xp)/
           systems[0].pdf->xfx(systems[0].beam,systems[0].incoming->dataPtr(),q2   ,xB);
    // numerator factors
    wgt /= (1.-zp);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = BGFME(xp,x2,x3,xperp,acoeff,l,m);
  }
  // compute the azimuthal average of the weight
  wgt *= azicoeff[0]+0.5*azicoeff[2]+0.5*azicoeff[4];
  // finally factor as picked one line
  wgt *= 2.;
  // decide whether or not to accept the weight
  if(UseRandom::rnd()>wgt) return;
  // if accepted generate generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi)),sphi(sin(phi));
    phiwgt =  azicoeff[0]+azicoeff[5]*sphi*cphi
      +azicoeff[1]*cphi+azicoeff[2]*sqr(cphi)
      +azicoeff[3]*sphi+azicoeff[4]*sqr(sphi);
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in VBFMECorrection"
				  << "::applyHardMatrixElementCorrection() to"
				  << " generate phi" << Exception::eventerror;
  // compute the new incoming and outgoing momenta
  Lorentz5Momentum p1 = Lorentz5Momentum( 0.5*Q*xperp*cos(phi), 0.5*Q*xperp*sin(phi),
					  -0.5*Q*x2,0.*GeV,0.*GeV);
  p1.rescaleEnergy();
  Lorentz5Momentum p2 = Lorentz5Momentum(-0.5*Q*xperp*cos(phi),-0.5*Q*xperp*sin(phi),
					 -0.5*Q*x3,0.*GeV,0.*GeV);
  p2.rescaleEnergy();
  Lorentz5Momentum pin(0.*GeV,0.*GeV,-0.5*x1*Q,-0.5*x1*Q,0.*GeV);
  // we need inverse of the rotation, i.e back to lab from breit
  rot.invert();
  // transform the momenta to lab frame
  pin *= rot;
  p1  *= rot;
  p2  *= rot;
  // test to ensure outgoing particles can be put on-shell
  if(!BGF) {
    if(p1.e()<systems[0].outgoing->dataPtr()->constituentMass()) return;
    if(p2.e()<gluon              ->constituentMass()) return;
  }
  else {
    if(p1.e()<systems[0].outgoing->dataPtr()      ->constituentMass()) return;
    if(p2.e()<systems[0].incoming->dataPtr()->CC()->constituentMass()) return;
  }
  // stats for weights > 1
  if(wgt>1.) {
    ++_nover;
    if(!BGF) {
      _maxwgt.first  = max(_maxwgt.first ,wgt);
      _comptonover.push_back(make_pair(xp,zp));
    }
    else {
      _maxwgt.second = max(_maxwgt.second,wgt);
      _bgfover.push_back(make_pair(xp,zp));
    }
  }
  // points into hist
  ++_ngen;
  if(!BGF) _compton.push_back(make_pair(xp,zp));
  else     _bgf    .push_back(make_pair(xp,zp));
  // create the new particles and add to ShowerTree
  bool isquark = systems[0].incoming->colourLine();
  if(!BGF) {
    PPtr newin  = new_ptr(Particle(*systems[0].incoming));
    newin->set5Momentum(pin);
    PPtr newg   = gluon                         ->produceParticle(p2 );
    PPtr newout = systems[0].outgoing->dataPtr()->produceParticle(p1 ); 
    ColinePtr col=isquark ? 
      systems[0].incoming->colourLine() : systems[0].incoming->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine());
    // final-state emission
    if(xp>zp) {
      col->removeColoured(newout,!isquark);
      col->addColoured(newin,!isquark);
      col->addColoured(newg,!isquark);
      newline->addColoured(newg,isquark);
      newline->addColoured(newout,!isquark);
    }
    // initial-state emission
    else {
      col->removeColoured(newin ,!isquark);
      col->addColoured(newout,!isquark);
      col->addColoured(newg,isquark);
      newline->addColoured(newg,!isquark);
      newline->addColoured(newin,!isquark);
    }
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=systems[0].incoming) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB/xp);
      cit->first->perturbative(xp>zp);
      if(xp<=zp) orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=systems[0].outgoing) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(xp<=zp);
      if(xp>zp) orig=cit->first->original();
    }
    assert(orig);
    // add the gluon
    ShowerParticlePtr sg=new_ptr(ShowerParticle(*newg,1,true));
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(orig,newg,sg));
    gluon->perturbative(false);
    tree->outgoingLines().insert(make_pair(gluon,sg));
    tree->hardMatrixElementCorrection(true);
    //cerr << *newin << "\n";
    //cerr << *newout << "\n";
    //cerr << *sg << "\n";
  }
  else {
    PPtr newin   = gluon                    ->produceParticle(pin);
    PPtr newqbar = systems[0].incoming->dataPtr()->CC()->produceParticle(p2 );
    PPtr newout  = systems[0].outgoing->dataPtr()      ->produceParticle(p1 );
    ColinePtr col=isquark ? systems[0].incoming->colourLine() : systems[0].incoming->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine()); 
    col    ->addColoured(newin  ,!isquark);
    newline->addColoured(newin  , isquark);
    col    ->addColoured(newout ,!isquark);
    newline->addColoured(newqbar, isquark);
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=systems[0].incoming) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB/xp);
      cit->first->perturbative(false);
      orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=systems[0].outgoing) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(true);
    }
    assert(orig);
    // add the (anti)quark
    ShowerParticlePtr sqbar=new_ptr(ShowerParticle(*newqbar,1,true));
    ShowerProgenitorPtr qbar=new_ptr(ShowerProgenitor(orig,newqbar,sqbar));
    qbar->perturbative(false);
    tree->outgoingLines().insert(make_pair(qbar,sqbar));
    tree->hardMatrixElementCorrection(true);
    //cerr << *newin  << "\n";
    //cerr << *newout << "\n";
    //cerr << *sqbar  << "\n";
  }
}

bool VBFMECorrection::softMatrixElementVeto(ShowerProgenitorPtr initial,
					    ShowerParticlePtr parent,Branching br) {
  return false;
}

void VBFMECorrection::dofinish() {
  MECorrectionBase::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  outfile << "SET FONT DUPLEX\n";
  outfile << "SET ORDER X Y\n";
  outfile << "SET WINDOW X 2 9 Y 2 9\n";
  outfile << "TITLE BOTTOM \"x0p1\"\n";
  outfile << "CASE         \" X X\"\n";
  outfile << "TITLE LEFT \"z0p1\"\n";
  outfile << "CASE       \" X X\"\n";
  outfile << "SET LIMITS X 0 1 Y 0 1\n";
  for(unsigned int ix=0;ix<_compton.size();++ix) {
    outfile << _compton[ix].first << "\t" << _compton[ix].second << "\n";
    if(ix%50000==0&&ix>0) outfile << "PLOT\n";
  }
  if(!_compton.empty()) outfile << "PLOT\n";
  for(double xp=0;xp<=1.001;xp+=0.01) {
    outfile << xp << "\t" << 1./(1.+xp-sqr(xp)) << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double z=0;z<=1.001;z+=0.01) {
    outfile << 1./(1+z*(1.-z)) << "\t" << z << "\n";
  }
  outfile << "JOIN BLUE\n";;
  outfile << "NEW FRAME\n";
  outfile << "SET FONT DUPLEX\n";
  outfile << "SET ORDER X Y\n";
  outfile << "SET WINDOW X 2 9 Y 2 9\n";
  outfile << "TITLE BOTTOM \"x0p1\"\n";
  outfile << "CASE         \" X X\"\n";
  outfile << "TITLE LEFT \"z0p1\"\n";
  outfile << "CASE       \" X X\"\n";
  outfile << "SET LIMITS X 0 1 Y 0 1\n";
  for(unsigned int ix=0;ix<_comptonover.size();++ix) {
    outfile << _comptonover[ix].first << "\t" << _comptonover[ix].second << "\n";
    if(ix%50000==0&&ix>0) outfile << "PLOT RED\n";
  }
  if(!_comptonover.empty()) outfile << "PLOT RED\n";
  for(double xp=0;xp<=1.001;xp+=0.01) {
    outfile << xp << "\t" << 1./(1.+xp-sqr(xp)) << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double z=0;z<=1.001;z+=0.01) {
    outfile << 1./(1+z*(1.-z)) << "\t" << z << "\n";
  }
  outfile << "JOIN BLUE\n";
  outfile << "NEW FRAME\n";
  outfile << "SET FONT DUPLEX\n";
  outfile << "SET ORDER X Y\n";
  outfile << "SET WINDOW X 2 9 Y 2 9\n";
  outfile << "TITLE BOTTOM \"x0p1\"\n";
  outfile << "CASE         \" X X\"\n";
  outfile << "TITLE LEFT \"z0p1\"\n";
  outfile << "CASE       \" X X\"\n";
  outfile << "SET LIMITS X 0 1 Y 0 1\n";
  for(unsigned int ix=0;ix<_bgf.size();++ix) {
    outfile << _bgf[ix].first << "\t" << _bgf[ix].second << "\n";
    if(ix%50000==0&&ix>0) outfile << "PLOT\n";
  }
  if(!_bgf.empty()) outfile << "PLOT\n";
  for(double z=0;z<=1.001;z+=0.01) {
    double xp = 2.*z/(1.+(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    double zp = 0.5* (1.-(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    outfile << xp << "\t" << zp << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double z=0;z<=1.001;z+=0.01) {
    double xp = 2.*z/(1.+(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    double zp = 0.5* (1.-(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    outfile << xp << "\t" << 1.-zp << "\n";
  }
  outfile << "JOIN BLUE\n";;
  outfile << "NEW FRAME\n";
  outfile << "SET FONT DUPLEX\n";
  outfile << "SET ORDER X Y\n";
  outfile << "SET WINDOW X 2 9 Y 2 9\n";
  outfile << "TITLE BOTTOM \"x0p1\"\n";
  outfile << "CASE         \" X X\"\n";
  outfile << "TITLE LEFT \"z0p1\"\n";
  outfile << "CASE       \" X X\"\n";
  outfile << "SET LIMITS X 0 1 Y 0 1\n";
  for(unsigned int ix=0;ix<_bgfover.size();++ix) {
    outfile << _bgfover[ix].first << "\t" << _bgfover[ix].second << "\n";
    if(ix%50000==0&&ix>0) outfile << "PLOT RED\n";
  }
  if(!_bgfover.empty()) outfile << "PLOT RED\n";
  for(double z=0;z<=1.001;z+=0.01) {
    double xp = 2.*z/(1.+(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    double zp = 0.5* (1.-(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    outfile << xp << "\t" << zp << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double z=0;z<=1.001;z+=0.01) {
    double xp = 2.*z/(1.+(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    double zp = 0.5* (1.-(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    outfile << xp << "\t" << 1.-zp << "\n";
  }
  outfile << "JOIN BLUE\n";;
  outfile.close();
  if(_ntry==0) return;
  generator()->log() << "VBFMECorrection when applying the hard correction "
		     << "generated " << _ntry << " trial emissions of which "
		     << _ngen << " were accepted\n";
  if(_nover==0) return;
  generator()->log() << "VBFMECorrection when applying the hard correction " 
		     << _nover << " weights larger than one were generated of which"
		     << " the largest was " << _maxwgt.first << " for the QCD compton"
		     << " processes and " << _maxwgt.second << " for the BGF process\n";
}

void VBFMECorrection::doinit() throw(InitException) {
  MECorrectionBase::doinit();
  // integrals of me over phase space
  double r5=sqrt(5.),darg((r5-1.)/(r5+1.)),ath(0.5*log((1.+1./r5)/(1.-1./r5)));
  _comptonint = 2.*(-21./20.-6.*r5/25.*ath+sqr(Constants::pi)/3.
		    -2.*Math::ReLi2(1.-darg)-2.*Math::ReLi2(1.-1./darg));
}
