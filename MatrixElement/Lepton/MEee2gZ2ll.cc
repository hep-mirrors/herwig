// -*- C++ -*-
//
// MEee2gZ2ll.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2gZ2ll class.
//

#include "MEee2gZ2ll.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <numeric>
#include "Herwig/Shower/RealEmissionProcess.h"

using namespace Herwig;

void MEee2gZ2ll::getDiagrams() const {
  // specific the diagrams
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  tcPDPtr em = getParticleData(ParticleID::eminus);
  // setup the processes
  for( int i =11;i<=16;++i) {
    if(allowed_==0 || (allowed_==1 && i%2==1) || (allowed_==2&&i==11)
       || (allowed_==3&&i==13) || (allowed_==4&&i==15)) {
      tcPDPtr lm = getParticleData(i);
      tcPDPtr lp = lm->CC();
      add(new_ptr((Tree2toNDiagram(2), em, ep, 1, gamma_, 3, lm, 3, lp, -1)));
      add(new_ptr((Tree2toNDiagram(2), em, ep, 1, Z0_, 3, lm, 3, lp, -2)));
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
  os << allowed_ << FFZVertex_ << FFPVertex_ << gamma_ << Z0_ 
     << alphaQED_ << ounit(pTmin_,GeV) << preFactor_; 
}

void MEee2gZ2ll::persistentInput(PersistentIStream & is, int) {
  is >> allowed_ >> FFZVertex_ >> FFPVertex_ >> gamma_ >> Z0_
     >> alphaQED_ >> iunit(pTmin_,GeV) >> preFactor_; 
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MEee2gZ2ll,HwMEBase>
describeMEee2gZ2ll("Herwig::MEee2gZ2ll", "HwMELepton.so");

void MEee2gZ2ll::Init() {

  static ClassDocumentation<MEee2gZ2ll> documentation
    ("The MEee2gZ2ll class implements the matrix element for"
     "e+e- to leptons via Z and photon exchange using helicity amplitude"
     "techniques");

  static Switch<MEee2gZ2ll,int> interfaceallowed
    ("Allowed",
     "Allowed outgoing leptons",
     &MEee2gZ2ll::allowed_, 0, false, false);
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

  static Parameter<MEee2gZ2ll,Energy> interfacepTMin
    ("pTMin",
     "Minimum pT for hard radiation",
     &MEee2gZ2ll::pTmin_, GeV, 1.0*GeV, 0.001*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEee2gZ2ll,double> interfacePrefactor
    ("Prefactor",
     "Prefactor for the overestimate of the emission probability",
     &MEee2gZ2ll::preFactor_, 6.0, 1.0, 100.0,
     false, false, Interface::limited);

  static Reference<MEee2gZ2ll,ShowerAlpha> interfaceEMCoupling
    ("AlphaQED",
     "Pointer to the object to calculate the EM coupling for the correction",
     &MEee2gZ2ll::alphaQED_, false, false, true, false, false);
}

double MEee2gZ2ll::me2() const {
  return loME(mePartonData(),rescaledMomenta(),true);
}

double MEee2gZ2ll::loME(const vector<cPDPtr> & partons, 
			const vector<Lorentz5Momentum> & momenta,
			bool first) const {
  vector<SpinorWaveFunction> fin,aout;
  vector<SpinorBarWaveFunction>  ain,fout;
  SpinorWaveFunction    ein   (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction pin   (momenta[1],partons[1],incoming);
  SpinorBarWaveFunction ilmout(momenta[2],partons[2],outgoing);
  SpinorWaveFunction    ilpout(momenta[3],partons[3],outgoing);
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
  if(first) {
    DVector save;
    save.push_back(lastCont);
    save.push_back(lastBW);
    meInfo(save);
  }
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
  ProductionMatrixElement gamma (PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement Zboson(PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  //   // wavefunctions for the intermediate particles
  VectorWaveFunction interZ,interG;
  // temporary storage of the different diagrams
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1,outhel2;
  double total[3]={0.,0.,0.};
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      // intermediate Z
      interZ = FFZVertex_->evaluate(sHat(),1,Z0_,fin[inhel1],ain[inhel2]);
      // intermediate photon
      interG = FFPVertex_->evaluate(sHat(),1,gamma_,fin[inhel1],ain[inhel2]);
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {		
	  // first the Z exchange diagram
	  diag1 = FFZVertex_->evaluate(sHat(),aout[outhel2],fout[outhel1],
					  interZ);
	  // then the photon exchange diagram
	  diag2 = FFPVertex_->evaluate(sHat(),aout[outhel2],fout[outhel1],
					  interG);
	  // add up squares of individual terms
	  total[1] += norm(diag1);
	  Zboson(inhel1,inhel2,outhel1,outhel2) = diag1;
	  total[2] += norm(diag2);
	  gamma (inhel1,inhel2,outhel1,outhel2) = diag2;
	  // the full thing including interference
	  diag1 += diag2;
	  total[0] += norm(diag1);
	  output(inhel1,inhel2,outhel1,outhel2) = diag1;
	}
      }
    }
  }
  // results
  for(int ix=0;ix<3;++ix) total[ix] *= 0.25;
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
			 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    total[0] = output.average(rho[0],rho[1]);
    total[1] = Zboson.average(rho[0],rho[1]);
    total[2] = gamma .average(rho[0],rho[1]);
  }
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
    tSpinPtr spin = hard[ix]->spinInfo();
    if(ix<2) {
      tcPolarizedBeamPDPtr beam = 
	dynamic_ptr_cast<tcPolarizedBeamPDPtr>(hard[ix]->dataPtr());
      if(beam) spin->rhoMatrix() = beam->rhoMatrix();
    }
    spin->productionVertex(hardvertex);
  }
}

void MEee2gZ2ll::doinit() {
  HwMEBase::doinit();
  // set the particle data objects
  Z0_=getParticleData(ThePEG::ParticleID::Z0);
  gamma_=getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm) throw InitException() << "Wrong type of StandardModel object in "
				  << "MEee2gZ2ll::doinit() the Herwig"
				  << " version must be used"
				  << Exception::runerror;
  FFZVertex_ = hwsm->vertexFFZ();
  FFPVertex_ = hwsm->vertexFFP();
}

void MEee2gZ2ll::rebind(const TranslationMap & trans) {
  FFZVertex_ = trans.translate(FFZVertex_);
  FFPVertex_ = trans.translate(FFPVertex_);
  Z0_        = trans.translate(Z0_);
  gamma_     = trans.translate(gamma_);
  HwMEBase::rebind(trans);
}

IVector MEee2gZ2ll::getReferences() {
  IVector ret = HwMEBase::getReferences();
  ret.push_back(FFZVertex_);
  ret.push_back(FFPVertex_);
  ret.push_back(Z0_       );
  ret.push_back(gamma_    );
  return ret;
}

RealEmissionProcessPtr MEee2gZ2ll::generateHardest(RealEmissionProcessPtr born,
						   ShowerInteraction inter) {
  // check if QED switched on
  if(inter==ShowerInteraction::QCD) return RealEmissionProcessPtr();
  // generate the momenta for the hard emission
  vector<Lorentz5Momentum> emission;
  unsigned int iemit,ispect;
  Energy pTveto = generateHard(born,emission,iemit,ispect,false);
  // check if charged
  if(!partons_[2]->charged()) return RealEmissionProcessPtr();
  // maximum pT of emission
  if(pTveto<=ZERO) {
    born->pT()[ShowerInteraction::QED] = pTmin_;
    return born;
  }
  else {
    born->pT()[ShowerInteraction::QED] = pTveto;
  }
  born->interaction(ShowerInteraction::QED);
  // get the quark and antiquark
  ParticleVector qq;
  for(unsigned int ix=0;ix<2;++ix) qq.push_back(born->bornOutgoing()[ix]);
  bool order = qq[0]->id()>0;
  if(!order) swap(qq[0],qq[1]);
  // create the new quark, antiquark and gluon
  PPtr newq = qq[0]->dataPtr()->produceParticle(emission[2]);
  PPtr newa = qq[1]->dataPtr()->produceParticle(emission[3]);
  PPtr newg = gamma_->produceParticle(emission[4]);
  // create the output real emission process
  for(unsigned int ix=0;ix<born->bornIncoming().size();++ix) {
    born->incoming().push_back(born->bornIncoming()[ix]->dataPtr()->
			       produceParticle(born->bornIncoming()[ix]->momentum()));
  }
  if(order) {
    born->outgoing().push_back(newq);
    born->outgoing().push_back(newa);
  }
  else {
    born->outgoing().push_back(newa);
    born->outgoing().push_back(newq);
    swap(iemit,ispect);
  }
  born->outgoing().push_back(newg);
  // set emitter and spectator
  born->emitter   (iemit);
  born->spectator(ispect);
  born->emitted(4);
  return born;
}

double MEee2gZ2ll::meRatio(vector<cPDPtr> partons, 
			   vector<Lorentz5Momentum> momenta,
			   unsigned int iemitter, bool subtract) const {
  Lorentz5Momentum q = momenta[2]+momenta[3]+momenta[4];
  Energy2 Q2=q.m2();
  Energy2 lambda = sqrt((Q2-sqr(momenta[2].mass()+momenta[3].mass()))*
			(Q2-sqr(momenta[2].mass()-momenta[3].mass())));
  InvEnergy2 D[2];
  double lome[2];
  for(unsigned int iemit=0;iemit<2;++iemit) {
    unsigned int ispect = iemit==0 ? 1 : 0;
    Energy2 pipj = momenta[4      ] * momenta[2+iemit ];
    Energy2 pipk = momenta[4      ] * momenta[2+ispect];
    Energy2 pjpk = momenta[2+iemit] * momenta[2+ispect];
    double y = pipj/(pipj+pipk+pjpk);
    double z = pipk/(     pipk+pjpk);
    Energy mij = sqrt(2.*pipj+sqr(momenta[2+iemit].mass()));
    Energy2 lamB = sqrt((Q2-sqr(mij+momenta[2+ispect].mass()))*
 			(Q2-sqr(mij-momenta[2+ispect].mass())));
    Energy2 Qpk = q*momenta[2+ispect];
    Lorentz5Momentum pkt = 
      lambda/lamB*(momenta[2+ispect]-Qpk/Q2*q)
      +0.5/Q2*(Q2+sqr(momenta[2+ispect].mass())-sqr(momenta[2+ispect].mass()))*q;
    Lorentz5Momentum pijt = 
      q-pkt;
    double muj = momenta[2+iemit ].mass()/sqrt(Q2);
    double muk = momenta[2+ispect].mass()/sqrt(Q2);
    double vt = sqrt((1.-sqr(muj+muk))*(1.-sqr(muj-muk)))/(1.-sqr(muj)-sqr(muk));
    double v  = sqrt(sqr(2.*sqr(muk)+(1.-sqr(muj)-sqr(muk))*(1.-y))-4.*sqr(muk))
      /(1.-y)/(1.-sqr(muj)-sqr(muk));
    // dipole term
    D[iemit] = 0.5/pipj*(2./(1.-(1.-z)*(1.-y))
 			 -vt/v*(2.-z+sqr(momenta[2+iemit].mass())/pipj));
    // matrix element
    vector<Lorentz5Momentum> lomom(4);
    lomom[0] = momenta[0];
    lomom[1] = momenta[1];
    if(iemit==0) {
      lomom[2] = pijt;
      lomom[3] = pkt ;
    }
    else {
      lomom[3] = pijt;
      lomom[2] = pkt ;
    }
    lome[iemit]  = loME(partons,lomom,false);
  }
  InvEnergy2 ratio = realME(partons,momenta)
    *abs(D[iemitter])/(abs(D[0]*lome[0])+abs(D[1]*lome[1]));
  double output = Q2*ratio;
  if(subtract) output -= 2.*Q2*D[iemitter];
  return output;
}

InvEnergy2 MEee2gZ2ll::realME(const vector<cPDPtr> & partons, 
			      const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<SpinorWaveFunction> fin,aout;
  vector<SpinorBarWaveFunction>  ain,fout;
  vector<VectorWaveFunction> gout;
  SpinorWaveFunction    ein  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction pin  (momenta[1],partons[1],incoming);
  SpinorBarWaveFunction qkout(momenta[2],partons[2],outgoing);
  SpinorWaveFunction    qbout(momenta[3],partons[3],outgoing);
  VectorWaveFunction    photon(momenta[4],partons[4],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    ein.reset(ix)  ;
    fin.push_back( ein  );
    pin.reset(ix)  ;
    ain.push_back( pin  );
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
    photon.reset(2*ix);
    gout.push_back(photon);
  }
  vector<Complex> diag(4,0.);
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1);
  double total(0.);
  for(unsigned int inhel1=0;inhel1<2;++inhel1) {
    for(unsigned int inhel2=0;inhel2<2;++inhel2) {
      // intermediate Z
      VectorWaveFunction interZ = 
	FFZVertex_->evaluate(scale(),1,Z0_,fin[inhel1],ain[inhel2]);
      // intermediate photon
      VectorWaveFunction interG = 
	FFPVertex_->evaluate(scale(),1,gamma_,fin[inhel1],ain[inhel2]);
      for(unsigned int outhel1=0;outhel1<2;++outhel1) {
	for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	  for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	    SpinorBarWaveFunction off1 =
	      FFPVertex_->evaluate(scale(),3,partons[2],fout[outhel1],gout[outhel3]);
	    diag[0] = FFZVertex_->evaluate(scale(),aout[outhel2],off1,interZ);
	    diag[1] = FFPVertex_->evaluate(scale(),aout[outhel2],off1,interG);
	    SpinorWaveFunction    off2 = 
	      FFPVertex_->evaluate(scale(),3,partons[3],aout[outhel2],gout[outhel3]);
	    diag[2] = FFZVertex_->evaluate(scale(),off2,fout[outhel1],interZ);
	    diag[3] = FFPVertex_->evaluate(scale(),off2,fout[outhel1],interG);
	    // sum of diagrams
	    Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    // matrix element
	    output(inhel1,inhel2,outhel1,outhel2,outhel3)=sum;
	    // me2
	    total += norm(sum);
	  }
	}
      }		
    }
  }
  // spin average
  total *= 0.25;
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(partons[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(partons[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = 
      {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
       beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    total = output.average(rho[0],rho[1]);
  }
  // divide out the coupling and charge
  total /= norm(FFPVertex_->norm())*
    sqr(double(mePartonData()[2]->iCharge())/3.);
  // return the total
  return total*UnitRemoval::InvE2;
}

Energy MEee2gZ2ll::generateHard(RealEmissionProcessPtr born, 
				vector<Lorentz5Momentum> & emmision,
				unsigned int & iemit, unsigned int & ispect,
				bool applyVeto) {
  // get the momenta of the incoming and outgoing particles 
  // incoming
  tPPtr em = born->bornIncoming()[0];
  tPPtr ep = born->bornIncoming()[1];
  if(em->id()<0) swap(em,ep);
  // outgoing
  tPPtr qk = born->bornOutgoing()[0];
  tPPtr qb = born->bornOutgoing()[1];
  if(qk->id()<0) swap(qk,qb);
  // extract the momenta 
  loMomenta_.clear();
  loMomenta_.push_back(em->momentum());
  loMomenta_.push_back(ep->momentum());
  loMomenta_.push_back(qk->momentum());
  loMomenta_.push_back(qb->momentum());
  // and ParticleData objects
  partons_.resize(5);
  partons_[0]=em->dataPtr();
  partons_[1]=ep->dataPtr();
  partons_[2]=qk->dataPtr();
  partons_[3]=qb->dataPtr();
  partons_[4]=gamma_;
  // boost from lab to CMS frame with outgoing particles
  // along the z axis
  LorentzRotation eventFrame( ( loMomenta_[2] + loMomenta_[3] ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*loMomenta_[2];
  eventFrame.rotateZ( -spectator.phi() );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  // mass of the final-state system
  Energy2 M2 = (loMomenta_[2]+loMomenta_[3]).m2();
  Energy  M  = sqrt(M2);
  double mu1 = loMomenta_[2].mass()/M;
  double mu2 = loMomenta_[3].mass()/M;
  double mu12 = sqr(mu1), mu22 = sqr(mu2);
  double lambda = sqrt(1.+sqr(mu12)+sqr(mu22)-2.*mu12-2.*mu22-2.*mu12*mu22);
  // max pT
  Energy pTmax = 0.5*sqrt(M2)*
    (1.-sqr(loMomenta_[2].mass()+loMomenta_[3].mass())/M2);
  // max y
  double ymax = acosh(pTmax/pTmin_);
  double a = alphaQED_->overestimateValue()/Constants::twopi*
    2.*ymax*preFactor_*sqr(double(mePartonData()[2]->iCharge())/3.);
  // variables for the emission
  Energy pT[2];
  double  y[2],phi[2],x3[2],x1[2][2],x2[2][2];
  double contrib[2][2];
  // storage of the real emission momenta
  vector<Lorentz5Momentum> realMomenta[2][2]=
    {{vector<Lorentz5Momentum>(5),vector<Lorentz5Momentum>(5)},
     {vector<Lorentz5Momentum>(5),vector<Lorentz5Momentum>(5)}};
  for(unsigned int ix=0;ix<2;++ix)
    for(unsigned int iy=0;iy<2;++iy)
      for(unsigned int iz=0;iz<2;++iz)
	realMomenta[ix][iy][iz] = loMomenta_[iz];
  // generate the emission
  for(unsigned int ix=0;ix<2;++ix) {
    if(ix==1) {
      swap(mu1 ,mu2 );
      swap(mu12,mu22);
    }
    pT[ix] = pTmax;
    y [ix] = 0.;
    bool reject = true;
    do {
      // generate pT
      pT[ix] *= pow(UseRandom::rnd(),1./a);
      if(pT[ix]<pTmin_) {
	pT[ix] = -GeV;
	break;
      }
      // generate y
      y[ix] = -ymax+2.*UseRandom::rnd()*ymax;
      // generate phi
      phi[ix] = UseRandom::rnd()*Constants::twopi;
      // calculate x3 and check in allowed region
      x3[ix] = 2.*pT[ix]*cosh(y[ix])/M;
      if(x3[ix] < 0. || x3[ix] > 1. -sqr( mu1 + mu2 ) ) continue;
      // find the possible solutions for x1
      double xT2 = sqr(2./M*pT[ix]);
      double root = (-sqr(x3[ix])+xT2)*
	(xT2*mu22+2.*x3[ix]-sqr(mu12)+2.*mu22+2.*mu12-sqr(x3[ix])-1.
	 +2.*mu12*mu22-sqr(mu22)-2.*mu22*x3[ix]-2.*mu12*x3[ix]);
      double c1=2.*sqr(x3[ix])-4.*mu22-6.*x3[ix]+4.*mu12-xT2*x3[ix]
	+2.*xT2-2.*mu12*x3[ix]+2.*mu22*x3[ix]+4.;
      if(root<0.) continue;
      x1[ix][0] = 1./(4.-4.*x3[ix]+xT2)*(c1-2.*sqrt(root));
      x1[ix][1] = 1./(4.-4.*x3[ix]+xT2)*(c1+2.*sqrt(root));
      // change sign of y if 2nd particle emits
      if(ix==1) y[ix] *=-1.;
      // loop over the solutions
      for(unsigned int iy=0;iy<2;++iy) {
	contrib[ix][iy]=0.;
	// check x1 value allowed
	if(x1[ix][iy]<2.*mu1||x1[ix][iy]>1.+mu12-mu22) continue;
	// calculate x2 value and check allowed
	x2[ix][iy] = 2.-x3[ix]-x1[ix][iy];
	double root = max(0.,sqr(x1[ix][iy])-4.*mu12);
	root = sqrt(root);
	double x2min = 1.+mu22-mu12
	  -0.5*(1.-x1[ix][iy]+mu12-mu22)/(1.-x1[ix][iy]+mu12)*(x1[ix][iy]-2.*mu12+root);
	double x2max = 1.+mu22-mu12
	  -0.5*(1.-x1[ix][iy]+mu12-mu22)/(1.-x1[ix][iy]+mu12)*(x1[ix][iy]-2.*mu12-root);
	if(x2[ix][iy]<x2min||x2[ix][iy]>x2max) continue;
	// check the z components
	double z1 =  sqrt(sqr(x1[ix][iy])-4.*mu12-xT2);
	double z2 = -sqrt(sqr(x2[ix][iy])-4.*mu22);
	double z3 =  pT[ix]*sinh(y[ix])*2./M;
	if(ix==1) z3 *=-1.;
	if(abs(-z1+z2+z3)<1e-9) z1 *= -1.;
	if(abs(z1+z2+z3)>1e-5) continue;
	// if using as an ME correction the veto
	if(applyVeto) {
// 	  double xb = x1[ix][iy], xc = x2[ix][iy];
// 	  double b  = mu12, c = mu22;
// 	  double r = 0.5*(1.+b/(1.+c-xc));
// 	  double z1  = r + (xb-(2.-xc)*r)/sqrt(sqr(xc)-4.*c);
// 	  double kt1 = (1.-b+c-xc)/z1/(1.-z1);
// 	  r = 0.5*(1.+c/(1.+b-xb));
// 	  double z2  = r + (xc-(2.-xb)*r)/sqrt(sqr(xb)-4.*b);
// 	  double kt2 = (1.-c+b-xb)/z2/(1.-z2);
// 	  if(ix==1) {
// 	    swap(z1 ,z2);
// 	    swap(kt1,kt2);
// 	  }
// 	  // veto the shower region
// 	  if( kt1 < d_kt1_ || kt2 < d_kt2_ ) continue;
	}
	// construct the momenta
	realMomenta[ix][iy][4] =
	  Lorentz5Momentum(pT[ix]*cos(phi[ix]),pT[ix]*sin(phi[ix]),
			   pT[ix]*sinh(y[ix]) ,pT[ix]*cosh(y[ix]),ZERO);
	if(ix==0) {
	  realMomenta[ix][iy][2] =
	    Lorentz5Momentum(-pT[ix]*cos(phi[ix]),-pT[ix]*sin(phi[ix]),
			     z1*0.5*M,x1[ix][iy]*0.5*M,M*mu1);
	  realMomenta[ix][iy][3] =
	    Lorentz5Momentum(ZERO,ZERO, z2*0.5*M,x2[ix][iy]*0.5*M,M*mu2);
	}
	else {
	  realMomenta[ix][iy][2] =
	    Lorentz5Momentum(ZERO,ZERO,-z2*0.5*M,x2[ix][iy]*0.5*M,M*mu2);
	  realMomenta[ix][iy][3] =
	    Lorentz5Momentum(-pT[ix]*cos(phi[ix]),-pT[ix]*sin(phi[ix]),
			     -z1*0.5*M,x1[ix][iy]*0.5*M,M*mu1);
	}
	// boost the momenta back to the lab
	for(unsigned int iz=2;iz<5;++iz)
	  realMomenta[ix][iy][iz] *= eventFrame;
	// jacobian and prefactors for the weight
	Energy J = M/sqrt(xT2)*abs(-x1[ix][iy]*x2[ix][iy]+2.*mu22*x1[ix][iy]
				   +x2[ix][iy]+x2[ix][iy]*mu12+mu22*x2[ix][iy]
				   -sqr(x2[ix][iy]))
	  /pow(sqr(x2[ix][iy])-4.*mu22,1.5);
	// prefactors etc
	contrib[ix][iy] = 0.5*pT[ix]/J/preFactor_/lambda;
	// matrix element piece
	contrib[ix][iy] *= meRatio(partons_,realMomenta[ix][iy],ix,false);
	// coupling piece
	contrib[ix][iy] *= alphaQED_->ratio(sqr(pT[ix]));
      }
      if(contrib[ix][0]+contrib[ix][1]>1.) {
	ostringstream s;
	s << "MEee2gZ2qq::generateHardest weight for channel " << ix
	  << "is " << contrib[ix][0]+contrib[ix][1] 
	  << " which is greater than 1";
	generator()->logWarning( Exception(s.str(), Exception::warning) );
      }
      reject =  UseRandom::rnd() > contrib[ix][0] + contrib[ix][1];
    }
    while (reject);
    if(pT[ix]<pTmin_)
      pT[ix] = -GeV;
  }
  // now pick the emmision with highest pT
  Energy pTemit(ZERO);
  // no emission
  if(pT[0]<ZERO&&pT[1]<ZERO) return -GeV;
  // which one emitted
  if(pT[0]>pT[1]) {
    iemit  = 2;
    ispect = 3;
    pTemit = pT[0];
    if(UseRandom::rnd()<contrib[0][0]/(contrib[0][0]+contrib[0][1]))
      emmision = realMomenta[0][0];
    else
      emmision = realMomenta[0][1];
  }
  else {
    iemit  = 3;
    ispect = 2;
    pTemit = pT[1];
    if(UseRandom::rnd()<contrib[1][0]/(contrib[1][0]+contrib[1][1]))
      emmision = realMomenta[1][0];
    else
      emmision = realMomenta[1][1];
  }
  // return pT of emmision
  return pTemit;
}






