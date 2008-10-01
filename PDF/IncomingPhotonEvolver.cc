// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IncomingPhotonEvolver class.
//

#include "IncomingPhotonEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/PDF/PartonExtractor.h"

using namespace Herwig;

IncomingPhotonEvolver::IncomingPhotonEvolver() 
  : PDFMax_(50.), PDFPower_(1.0)
{}

void IncomingPhotonEvolver::
handle(EventHandler & eh, const tPVector & ,
       const Hint & ) throw(Veto, Stop, Exception) {
  // extract the incoming partons from the hard process
  PPair incomingPartons = 
    eh.currentEvent()->primarySubProcess()->incoming();
  // check one and only one is a photon
  if(!( (incomingPartons.first->id()  == ParticleID::gamma && 
	 incomingPartons.second->id() != ParticleID::gamma) ||
	(incomingPartons.first->id()  != ParticleID::gamma && 
	 incomingPartons.second->id() == ParticleID::gamma))) return;
  // extract the incoming hadrons from the event
  PPair incomingHadrons = 
    eh.currentEvent()->incoming();
  // calculate the value of x
  double x[2] = 
    {incomingPartons.first ->momentum().t()/incomingHadrons.first ->momentum().t(),
     incomingPartons.second->momentum().t()/incomingHadrons.second->momentum().t()};
  // swap order so photon first
  if(incomingPartons.first->id()!=ParticleID::gamma) {
    swap(incomingPartons.first,incomingPartons.second);
    swap(incomingHadrons.first,incomingHadrons.second);
    swap(x[0],x[1]);
  }
  // calculate CMF momentum get the pt scale for the process
  Energy pt(generator()->maximumCMEnergy());
  Energy ptmin=2.*GeV;
  Lorentz5Momentum pcmf;
  for(ParticleVector::const_iterator 
	cit=eh.currentEvent()->primarySubProcess()->outgoing().begin(),
	end=eh.currentEvent()->primarySubProcess()->outgoing().end();
      cit!=end;++cit) {
    pcmf += (**cit).momentum();
    Energy pttest = (**cit).momentum().perp();
    if(pttest<pt&&pt>2.*ptmin) pt=pttest;
  }
  // limits for the z integrand
  double lower = 1./PDFPower_;
  double upper = lower/pow(x[0],PDFPower_);
  // extract the PDF for the beam particle
  Ptr<BeamParticleData>::transient_const_pointer beam = 
    dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (incomingHadrons.first->dataPtr());
  assert(beam);
  // get the PDF 
  tcPDFPtr pdf = beam->pdf();
  assert(pdf);
  // generate the values of pt and z for the branching
  double wgt = PDFMax_*SM().alphaEM()/Constants::pi/PDFPower_*
    (1./pow(x[0],PDFPower_)-1.);
  double z,pdftotal,rwgt;
  unsigned int ncount(0);
  Energy scale;
  do {
    scale = max(pt,ptmin);
    // new value of pT
    pt *= pow(UseRandom::rnd(),0.5/wgt);
    // new value of z
    z = lower+UseRandom::rnd()*(upper-lower);
    z = pow(1./z/PDFPower_,1./PDFPower_);
    // the weight
    rwgt = 0.5*(1.+sqr(1.-z));
    // denominator of the pdf bit
    rwgt /= pdf->xfx(beam,photon_,sqr(scale),x[0]);
    // numerator of the pdf bit
    pdftotal = 0.;
    for(unsigned int ix=0;ix<partons_.size();++ix) { 
      double pdfval =  pdf->xfx(beam,partons_[ix],sqr(scale),x[0]/z);
      if(pdfval>0.) pdftotal += pdfval*sqr(double(partons_[ix]->iCharge())/3.);
    }
    rwgt *= pdftotal;
    // finally divide by the overestimate of the PDF bit
    rwgt /= PDFMax_/pow(z,PDFPower_);
    if(rwgt>1.) cerr << "testing weight problem " << rwgt << " " << PDFMax_ << "\n";
    if(UseRandom::rnd()>rwgt) {
      ++ncount;
      continue;
    }
    break;
  }
  while(true&&ncount<5000);
  if(ncount==5000) {
    throw Exception() << "Too many attempts to generate backward "
		      << "photon evolution in IncomingPhotonEvolver::handle()"
		      << Exception::eventerror; 
  }
  // now select the flavour of the emitted parton
  pdftotal *= UseRandom::rnd();
  tcPDPtr quark;
  for(unsigned int ix=0;ix<partons_.size();++ix) {
    double pdfval =  pdf->xfx(beam,partons_[ix],sqr(scale),x[0]/z);
    if(pdfval<=0.) continue;
    pdfval *= sqr(double(partons_[ix]->iCharge())/3.);
    if(pdftotal<pdfval) {
      quark = partons_[ix];
      break;
    }
    pdftotal -= pdfval;
  }
  assert(quark);
  // construct the kinematics
  // p and n vectors
  Energy mag = incomingHadrons.first->momentum().t();
  Lorentz5Momentum p(0.*GeV,0.*GeV, mag,mag);
  Lorentz5Momentum n(0.*GeV,0.*GeV,-mag,mag);
  if(incomingHadrons.first->momentum().z()<0.*GeV) swap(p,n);
  Energy2 pdotn = p*n;
  // calculate the momenta of the partons involved in the branching
  double betaq = 0.5*z*sqr(pt)/(x[0]*(1.-z)*pdotn);
  double phi = Constants::twopi*UseRandom::rnd();
  Lorentz5Momentum pin    = x[0]/z*p;
  Lorentz5Momentum pout   = (1.-z)*x[0]/z*p+betaq*n+
    Lorentz5Momentum(pt*cos(phi),pt*sin(phi),0.*GeV,0.*GeV);
  Lorentz5Momentum pgamma = pin-pout;
  pout  .rescaleMass();
  pgamma.rescaleMass();
  // compute the boosts for momentum conservation
  Energy2 shat = (incomingPartons.first ->momentum()+
		  incomingPartons.second->momentum()).m2();
  Energy2 S = (p+n).m2();
  Lorentz5Momentum pother = incomingPartons.second->momentum();
  // find alphas and betas in terms of desired basis
  double a[2] = {pgamma*n/pdotn,pother*n/pdotn};
  double b[2] = {pgamma*p/pdotn,pother*p/pdotn};
  Lorentz5Momentum p1p = pgamma - a[0]*p - b[0]*n;
  Lorentz5Momentum p2p = pother - a[1]*p - b[1]*n;
  // compute kappa
  Energy2 A = a[0]*b[1]*S;
  Energy2 B = shat - (a[0]*b[0]+a[1]*b[1])*S - (p1p+p2p).mag2();
  Energy2 C = a[1]*b[0]*S; 
  double rad = 1.-4.*A*C/sqr(B);
  if(rad < 0.) throw Exception() << "Can't generate backward evolution of the photon"
				 << " in IncomingPhotonEvolver::handle()"
				 << Exception::eventerror;
  double kp = B/(2.*A)*(1.+sqrt(rad));
  // now compute k1, k2
  rad = kp*(b[0]+kp*b[1])/(kp*a[0]+a[1])*(x[0]/x[1]);  
  if(rad <= 0.) throw Exception() << "Can't generate backward evolution of the photon"
				  << " in IncomingPhotonEvolver::handle()"
				  << Exception::eventerror;
  double k1 = sqrt(rad);
  double k2 = kp/k1;
  double beta[2] = 
    {getBeta((a[0]+b[0]), (a[0]-b[0]), (k1*a[0]+b[0]/k1), (k1*a[0]-b[0]/k1)),
     getBeta((a[1]+b[1]), (a[1]-b[1]), (a[1]/k2+k2*b[1]), (a[1]/k2-k2*b[1]))};
  if (p.z() > 0*MeV) {
    beta[0] = -beta[0]; 
    beta[1] = -beta[1];
  }
  // apply the boosts
  Boost betaboost(0, 0, beta[0]);
  pin   .boost(betaboost);
  pgamma.boost(betaboost);
  pout  .boost(betaboost);
  betaboost = Boost(0, 0, beta[1]);
  pother.boost(betaboost);
  Lorentz5Momentum newcmf = pother+pgamma;
  pcmf.rescaleMass();
  if(pin.e()/p.e()>1.||pin.z()/p.z()>1.)  
    throw Exception() << "Can't generate backward evolution of the photon"
		      << " in IncomingPhotonEvolver::handle()"
		      << Exception::eventerror;
  if(pother.e()/n.e()>1.||pother.z()/n.z()>1.)
    throw Exception() << "Can't generate backward evolution of the photon"
		      << " in IncomingPhotonEvolver::handle()"
		      << Exception::eventerror;
  if(newcmf.m()<0.*GeV||newcmf.e()<0.*GeV)  
    throw Exception() << "Can't generate backward evolution of the photon"
		      << " in IncomingPhotonEvolver::handle()"
		      << Exception::eventerror;
  Boost toRest   = pcmf.findBoostToCM();
  Boost fromRest = newcmf.boostVector();
  // apply the boosts to the outgoing particles
  for(ParticleVector::const_iterator 
	cit=eh.currentEvent()->primarySubProcess()->outgoing().begin(),
	end=eh.currentEvent()->primarySubProcess()->outgoing().end();
      cit!=end;++cit) {
    (**cit).deepBoost(  toRest);
    (**cit).deepBoost(fromRest);
    newcmf -= (**cit).momentum();
  }
  // now sort out the event record
  // make the new outgoing parton
  PPtr newOutgoing = quark->produceParticle(pout);
  // make the new incoming parton
  PPtr newIncoming = quark->produceParticle(pin);
  // and new other
  PPtr newOther    = incomingPartons.second->dataPtr()->
    produceParticle(pother);
  // colour connections
  if(newIncoming->id()>0) {
    ColinePtr newline = ColourLine::create(newIncoming);
    newline->addColoured(newOutgoing);
  }
  else {
    ColinePtr newline = ColourLine::create(newIncoming,true);
    newline->addAntiColoured(newOutgoing);
  }
  if(incomingPartons.second->colourLine())
    incomingPartons.second->colourLine()->addColoured(newOther);
  if(incomingPartons.second->antiColourLine())
    incomingPartons.second->antiColourLine()->addAntiColoured(newOther);
  // sort out the remnants
  // get the parton extractor
  PartonExtractor & pex = *eh.lastExtractor();
  // get the new partons
  tPPair newp = make_pair(newIncoming,newOther);
  // Creates the new remnants and returns the new PartonBinInstances
  PBIPair newbins = pex.newRemnants(incomingPartons, newp, eh.currentStep());
  // reset the momenta of the old incoming partons
  incomingPartons.first ->set5Momentum(pgamma);
  incomingPartons.second->set5Momentum(pother);  
  // sort out the mother/child stuff
  incomingHadrons.first ->abandonChild(incomingPartons.first);
  incomingHadrons.second->abandonChild(incomingPartons.second);
  newIncoming->addChild(incomingPartons.first);
  newOther->addChild(incomingPartons.second);
  newIncoming->addChild(newOutgoing);
  // add the new particles
  eh.currentEvent()->primarySubProcess()->addOutgoing(newOutgoing,false);
  eh.currentStep()->addParticle(newOutgoing);
  eh.currentEvent()->primarySubProcess()->changeIncoming(newIncoming,
  							 incomingPartons.first);
  eh.currentEvent()->primarySubProcess()->changeIncoming(newOther,
  							 incomingPartons.second);
  eh.currentStep()->addIntermediate(newIncoming);
  eh.currentStep()->addIntermediate(newOther   );
  // clean up the remnants
  ParticleVector children = incomingHadrons.first->children();
  for(unsigned int ix=0;ix<children.size();++ix) {
    if(children[ix]==newIncoming) continue;
    PPtr temp = children[ix];
    incomingHadrons.first->abandonChild(temp);
    ParticleVector childrenB = temp->children();
    for(unsigned int iy=0;iy<childrenB.size();++iy) {
      temp->abandonChild(childrenB[ix]);
      incomingHadrons.first->addChild(childrenB[ix]);
    }
    eh.currentStep()->removeParticle(temp);
  }
  children = incomingHadrons.second->children();
  for(unsigned int ix=0;ix<children.size();++ix) {
    if(children[ix]==newOther) continue;
    PPtr temp = children[ix];
    incomingHadrons.second->abandonChild(temp);
    ParticleVector childrenB = temp->children();
    for(unsigned int iy=0;iy<childrenB.size();++iy) {
      temp->abandonChild(childrenB[ix]);
      incomingHadrons.second->addChild(childrenB[ix]);
    }
    eh.currentStep()->removeParticle(temp);
  }
}

IBPtr IncomingPhotonEvolver::clone() const {
  return new_ptr(*this);
}

IBPtr IncomingPhotonEvolver::fullclone() const {
  return new_ptr(*this);
}

void IncomingPhotonEvolver::persistentOutput(PersistentOStream & os) const {
  os << PDFMax_ << PDFPower_ << photon_ << partons_;
}

void IncomingPhotonEvolver::persistentInput(PersistentIStream & is, int) {
  is >> PDFMax_ >> PDFPower_ >> photon_ >> partons_;
}

ClassDescription<IncomingPhotonEvolver> IncomingPhotonEvolver::initIncomingPhotonEvolver;
// Definition of the static class description member.

void IncomingPhotonEvolver::Init() {

  static ClassDocumentation<IncomingPhotonEvolver> documentation
    ("There is no documentation for the IncomingPhotonEvolver class");

}

void IncomingPhotonEvolver::doinit() throw(InitException) {
  StepHandler::doinit();
  photon_ = getParticleData(ParticleID::gamma);
  for(int ix=1;ix<=5;++ix) {
    partons_.push_back(getParticleData(ix));
    partons_.push_back(partons_.back()->CC());
  }
  // test of the pdf
//   ofstream output("pdf.top");
//   tcPDPtr proton = getParticleData(ParticleID::pplus);
//   Ptr<BeamParticleData>::transient_const_pointer beam = 
//     dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
//     (proton);
//   assert(beam);
//   // get the PDF 
//   tcPDFPtr pdf = beam->pdf();
//   for(int ix=-5;ix<=22;++ix) {
//     if(ix==0) ++ix;
//     if(ix==6) ix=21;
//     tcPDPtr quark = getParticleData(ix);
//     Energy2 scale = 20*GeV2;
//     output << "NEW FRAME\n";
//     output << "SET FONT DUPLEX\n";
//     output << "SET SCALE Y LOG\n";
//     string title = "distribution for" + quark->PDGName();
//     output << "TITLE TOP \""<< title << "\"\n";
//     for(double x=1e-3;x<1.;x*=1.1) {
//       double pdfval =  pdf->xfx(proton,quark,scale,x);
//       output << x << "\t" << pdfval << "\n";
//     }
//     output << "JOIN\n";
//   }
//   output.close();
}
