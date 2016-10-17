
//-*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWFermionsPOWHEGDecayer class.
//

#include "SMWFermionsPOWHEGDecayer.h"
#include <numeric>
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

using namespace Herwig;

SMWFermionsPOWHEGDecayer::SMWFermionsPOWHEGDecayer() 
  : CF_(4./3.), pTmin_(1.*GeV)
{  }

IBPtr SMWFermionsPOWHEGDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SMWFermionsPOWHEGDecayer::fullclone() const {
  return new_ptr(*this);
}

void SMWFermionsPOWHEGDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFGVertex_ << FFWVertex_ << gluon_ << ounit( pTmin_, GeV );
}

void SMWFermionsPOWHEGDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFGVertex_ >> FFWVertex_ >> gluon_ >> iunit( pTmin_, GeV );
}

ClassDescription<SMWFermionsPOWHEGDecayer> 
SMWFermionsPOWHEGDecayer::initSMWFermionsPOWHEGDecayer;
// Definition of the static class description member.

void SMWFermionsPOWHEGDecayer::Init() {

  static ClassDocumentation<SMWFermionsPOWHEGDecayer> documentation
    ("There is no documentation for the SMWFermionsPOWHEGDecayer class");

  static Parameter<SMWFermionsPOWHEGDecayer, Energy> interfacePtMin
    ("minpT",
     "The pt cut on hardest emision generation",
     &SMWFermionsPOWHEGDecayer::pTmin_, GeV, 1.*GeV, 0*GeV, 100000.0*GeV,
     false, false, Interface::limited);
}

RealEmissionProcessPtr SMWFermionsPOWHEGDecayer::
generateHardest(RealEmissionProcessPtr born) {
  assert(born->bornOutgoing().size()==2);
  // check coloured
  if(!born->bornOutgoing()[0]->dataPtr()->coloured()) return RealEmissionProcessPtr();
  // extract required info
  partons_.resize(2);
  quark_.resize(2);
  vector<PPtr> hardProcess;
  wboson_ = born->bornIncoming()[0];
  hardProcess.push_back(wboson_);
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    partons_[ix] = born->bornOutgoing()[ix]->dataPtr();
    quark_[ix]   = born->bornOutgoing()[ix]->momentum();
    quark_[ix].setMass(partons_[ix]->mass());
    hardProcess.push_back(born->bornOutgoing()[ix]);
  }
  bool order = partons_[0]->id()<0;
  if(order) {
    swap(partons_[0]   ,partons_[1]   );
    swap(quark_[0]     ,quark_[1]     );
    swap(hardProcess[1],hardProcess[2]);
  }
  gauge_.setMass(0.*MeV);
  // Get the W boson mass.
  mw2_ = (quark_[0] + quark_[1]).m2();
  // Generate emission and set _quark[0,1] and _gauge to be the 
  // momenta of q, qbar and g after the hardest emission:
  if(!getEvent(hardProcess)) {
    born->pT()[ShowerInteraction::QCD] = pTmin_;
    return born;
  }
  // Ensure the energies are greater than the constituent masses:
  for (int i=0; i<2; i++) {
    if (quark_[i].e() < partons_[i]->constituentMass()) return RealEmissionProcessPtr();
  }
  if (gauge_.e()    < gluon_     ->constituentMass()) return RealEmissionProcessPtr();
  // set masses
  quark_[0].setMass( partons_[0]->mass() );
  quark_[1].setMass( partons_[1]->mass() );
  gauge_   .setMass( ZERO );
  // // assign the emitter based on evolution scales
  unsigned int iemitter   = quark_[0]*gauge_ > quark_[1]*gauge_ ? 2 : 1;
  unsigned int ispectator = iemitter==1                         ? 1 : 2;
  // create new partices and insert
  PPtr wboson = wboson_->dataPtr()->produceParticle(wboson_->momentum());
  born->incoming().push_back(wboson);
  PPtr newq = partons_[0]->produceParticle(quark_[0]);
  PPtr newa = partons_[1]->produceParticle(quark_[1]);
  PPtr newg = gluon_->produceParticle(gauge_);
  // make colour connections
  newg->colourNeighbour(newq);
  newa->colourNeighbour(newg);
  // insert in output structure
  if(!order) {
    born->outgoing().push_back(newq);
    born->outgoing().push_back(newa);
  }
  else {
    born->outgoing().push_back(newa);
    born->outgoing().push_back(newq);
    swap(iemitter,ispectator);
  }
  born->outgoing().push_back(newg);
  born->emitter  (iemitter  );
  born->spectator(ispectator);
  born->emitted  (3);
  born->pT()[ShowerInteraction::QCD] = pT_;
  // return process
  born->interaction(ShowerInteraction::QCD);
  return born;
}

double SMWFermionsPOWHEGDecayer::
me2(const int ichan, const Particle & part,
    const ParticleVector & decay, MEOption meopt) const {
  // leading-order result
  double output = SMWDecayer::me2(ichan,part,decay,meopt);
  // check decay products coloured, otherwise return
  if(!decay[0]->dataPtr()->coloured()) return output;
  // inital masses, couplings  etc
  // W mass
  mW_ = part.mass();
  // strong coupling
  aS_ = SM().alphaS(sqr(mW_));
  // reduced mass
  double mu1_  = (decay[0]->dataPtr()->mass())/mW_;
  double mu2_  = (decay[1]->dataPtr()->mass())/mW_;
  // scale
  scale_ = sqr(mW_);
  // now for the nlo loop correction
  double virt = CF_*aS_/Constants::pi;
  // now for the real correction
  double realFact=0.;
  for(int iemit=0;iemit<2;++iemit) {
    double phi  = UseRandom::rnd()*Constants::twopi;
    // set the emitter and the spectator
    double muj  = iemit==0 ? mu1_ : mu2_;
    double muk  = iemit==0 ? mu2_ : mu1_;
    double muj2 = sqr(muj);
    double muk2 = sqr(muk);
    // calculate y
    double yminus = 0.; 
    double yplus  = 1.-2.*muk*(1.-muk)/(1.-muj2-muk2);
    double y = yminus + UseRandom::rnd()*(yplus-yminus);
    double v = sqrt(sqr(2.*muk2 + (1.-muj2-muk2)*(1.-y))-4.*muk2)
      /(1.-muj2-muk2)/(1.-y);
    double zplus  = (1.+v)*(1.-muj2-muk2)*y/2./(muj2+(1.-muj2-muk2)*y);
    double zminus = (1.-v)*(1.-muj2-muk2)*y/2./(muj2+(1.-muj2-muk2)*y);
    double z = zminus + UseRandom::rnd()*(zplus-zminus);
    double jac = (1.-y)*(yplus-yminus)*(zplus-zminus);
    // calculate x1,x2,x3,xT
    double x2 = 1.-y*(1.-muj2-muk2)-muj2+muk2;
    double x1 = 1.+muj2-muk2-z*(x2-2.*muk2);
    // copy the particle objects over for calculateRealEmission
    vector<PPtr> hardProcess(3);
    hardProcess[0] = const_ptr_cast<PPtr>(&part);
    hardProcess[1] = decay[0];
    hardProcess[2] = decay[1];
    realFact = 0.25*jac*sqr(1.-muj2-muk2)/
      sqrt((1.-sqr(muj-muk))*(1.-sqr(muj+muk)))/Constants::twopi
      *2.*CF_*aS_*calculateRealEmission(x1, x2, hardProcess, phi, 
					muj, muk, iemit, true);
  }
  // the born + virtual + real
  output = output*(1. + virt + realFact);
  return output;
}

double SMWFermionsPOWHEGDecayer::meRatio(vector<cPDPtr> partons, 
					 vector<Lorentz5Momentum> momenta,
				 	 unsigned int iemitter, bool subtract) const {
  Lorentz5Momentum q = momenta[1]+momenta[2]+momenta[3];
  Energy2 Q2=q.m2();
  Energy2 lambda = sqrt((Q2-sqr(momenta[1].mass()+momenta[2].mass()))*
			(Q2-sqr(momenta[1].mass()-momenta[2].mass())));
  InvEnergy2 D[2];
  double lome[2];
  for(unsigned int iemit=0;iemit<2;++iemit) {
    unsigned int ispect = iemit==0 ? 1 : 0;    
    Energy2 pipj = momenta[3      ] * momenta[1+iemit ];
    Energy2 pipk = momenta[3      ] * momenta[1+ispect];
    Energy2 pjpk = momenta[1+iemit] * momenta[1+ispect];
    double y = pipj/(pipj+pipk+pjpk); 
    double z = pipk/(     pipk+pjpk);
    Energy mij = sqrt(2.*pipj+sqr(momenta[1+iemit].mass()));
    Energy2 lamB = sqrt((Q2-sqr(mij+momenta[1+ispect].mass()))*
			(Q2-sqr(mij-momenta[1+ispect].mass())));
    Energy2 Qpk = q*momenta[1+ispect];
    Lorentz5Momentum pkt = 
      lambda/lamB*(momenta[1+ispect]-Qpk/Q2*q)
      +0.5/Q2*(Q2+sqr(momenta[1+ispect].mass())-sqr(momenta[1+ispect].mass()))*q;
    Lorentz5Momentum pijt = 
      q-pkt;
    double muj = momenta[1+iemit ].mass()/sqrt(Q2);
    double muk = momenta[1+ispect].mass()/sqrt(Q2);
    double vt = sqrt((1.-sqr(muj+muk))*(1.-sqr(muj-muk)))/(1.-sqr(muj)-sqr(muk));
    double v  = sqrt(sqr(2.*sqr(muk)+(1.-sqr(muj)-sqr(muk))*(1.-y))-4.*sqr(muk))
      /(1.-y)/(1.-sqr(muj)-sqr(muk));
    // dipole term
    D[iemit] = 0.5/pipj*(2./(1.-(1.-z)*(1.-y))
			 -vt/v*(2.-z+sqr(momenta[1+iemit].mass())/pipj));
    // matrix element
    vector<Lorentz5Momentum> lomom(3);
    lomom[0] = momenta[0];
    if(iemit==0) {
      lomom[1] = pijt;
      lomom[2] = pkt ;
    }
    else {
      lomom[2] = pijt;
      lomom[1] = pkt ;
    }
    lome[iemit]  = loME(partons,lomom);
  }
  InvEnergy2 ratio = realME(partons,momenta)*abs(D[iemitter])
    /(abs(D[0]*lome[0])+abs(D[1]*lome[1]));
  if(subtract)
    return Q2*(ratio-2.*D[iemitter]);
  else
    return Q2*ratio;
}

double SMWFermionsPOWHEGDecayer::loME(const vector<cPDPtr> & partons, 
			const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<VectorWaveFunction>    vin;
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  VectorWaveFunction    win  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  for(unsigned int ix=0;ix<3;++ix){
    win.reset(ix);
    vin.push_back(win);
  }
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  double total(0.);
  for(unsigned int inhel=0;inhel<3;++inhel) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	Complex diag1 = FFWVertex()->evaluate(scale_,aout[outhel2],fout[outhel1],vin[inhel]);
	total += norm(diag1);
      }
    }
  }
  // return the answer
  return total;
}
 
InvEnergy2 SMWFermionsPOWHEGDecayer::realME(const vector<cPDPtr> & partons, 
			      const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<VectorWaveFunction>     vin;
  vector<SpinorWaveFunction>     aout;
  vector<SpinorBarWaveFunction>  fout;
  vector<VectorWaveFunction>     gout;
  VectorWaveFunction    win  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  VectorWaveFunction    gluon(momenta[3],partons[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
    gluon.reset(2*ix);
    gout.push_back(gluon);
  }
  for(unsigned int ix=0;ix<3;++ix){
    win.reset(ix);
    vin.push_back(win);
  }
  vector<Complex> diag(2,0.);

  double total(0.);
  for(unsigned int inhel1=0;inhel1<3;++inhel1) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	  SpinorBarWaveFunction off1 =
	    FFGVertex()->evaluate(scale_,3,partons[1],fout[outhel1],gout[outhel3]);
	  diag[0] = FFWVertex()->evaluate(scale_,aout[outhel2],off1,vin[inhel1]);

	  SpinorWaveFunction off2 = 
	    FFGVertex()->evaluate(scale_,3,partons[2],aout[outhel2],gout[outhel3]);
	  diag[1] = FFWVertex()->evaluate(scale_,off2,fout[outhel1],vin[inhel1]);

	  // sum of diagrams
	  Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	  // me2
	  total += norm(sum);
	}
      }
    }
  }

  // divide out the coupling
  total /= norm(FFGVertex()->norm());
  // return the total
  return total*UnitRemoval::InvE2;
}

void SMWFermionsPOWHEGDecayer::doinit() {
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm) throw InitException() 
	      << "Wrong type of StandardModel object in "
	      << "SMWFermionsPOWHEGDecayer::doinit() "
	      << "the Herwig version must be used." 
	      << Exception::runerror;
  // cast the vertices
  FFWVertex_ = hwsm->vertexFFW();
  FFWVertex_->init();
  FFGVertex_ = hwsm->vertexFFG();
  FFGVertex_->init();
  SMWDecayer::doinit();
  gluon_ = getParticleData(ParticleID::g);
}

bool SMWFermionsPOWHEGDecayer::getEvent(vector<PPtr> hardProcess) {
  vector<Energy> particleMass;
  for(unsigned int ix=0;ix<hardProcess.size();++ix) {
    if(abs(hardProcess[ix]->id())==ParticleID::Wplus) {
      mW_ = hardProcess[ix]->mass();
    }
    else {
      particleMass.push_back(hardProcess[ix]->mass());
    }
  }
  if (particleMass.size()!=2)  {
    throw Exception()
      << "Number of outgoing particles is not equal to 2 in "
      << "SMWFermionPOWHEGDecayer::getEvent()" 
      << Exception::runerror;
  }
  // reduced mass
  mu1_ = particleMass[0]/mW_;
  mu2_ = particleMass[1]/mW_;
  // scale
  scale_ = sqr(mW_);
  // max pT
  Energy pTmax = 0.5*sqrt(mw2_);
  if(pTmax<pTmin_) return false;
  // Define over valued y_max & y_min according to the associated pt_min cut.
  double ymax  =  acosh(pTmax/pTmin_);
  double ymin  = -ymax;
  // pt of the emmission
  pT_ = pTmax;
  // prefactor
  double overEst = 4.;
  double prefactor = overEst*alphaS()->overestimateValue()*CF_*
    (ymax-ymin)/Constants::twopi;
  // loop to generate the pt and rapidity
  bool reject;  
  //arrays to hold the temporary  probabilities whilst the for loop progresses
  double probTemp[2][2]={{0.,0.},{0.,0.}};
  probTemp[0][0]=probTemp[0][1]=probTemp[1][0]=probTemp[1][1]=0.;
  double x1Solution[2][2] = {{0.,0.},{0.,0.}};
  double x2Solution[2][2] = {{0.,0.},{0.,0.}};
  double x3Solution[2]    = {0.,0.};
  Energy pT[2]            = {pTmax,pTmax};
  double yTemp[2]         = {0.,0.};
  double phi              = 0.;
  // do the competition
  for(int i=0; i<2; i++) {
    // set the emitter and the spectator
    double muj  = i==0 ? mu1_ : mu2_;
    double muk  = i==0 ? mu2_ : mu1_;
    double muj2 = sqr(muj);
    double muk2 = sqr(muk);
    do {
      // generation of phi
      phi = UseRandom::rnd() * Constants::twopi;
      // reject the emission
      reject = true; 
      // generate pt
      pT[i] *= pow(UseRandom::rnd(),1./prefactor);
      if(pT[i]<pTmin_) {
        pT[i] = -GeV;
        break;
      }
      // generate xT
      double xT2 = sqr(2./mW_*pT[i]);
      // generate y
      yTemp[i] = ymin + UseRandom::rnd()*(ymax-ymin);
      // generate x3 & x1 from pT & y
      double x1Plus  = 1-muk2+muj2;
      double x1Minus = 2.*muj;
      x3Solution[i]  = 2.*pT[i]*cosh(yTemp[i])/mW_;
      // prefactor
      double weightPrefactor = 0.5/sqrt((1.-sqr(muj-muk))*(1.-sqr(muj+muk)))/overEst;
      // calculate x1 & x2 solutions
      double discrim2 = (-sqr(x3Solution[i])+xT2)*
	(xT2*muk2+2.*x3Solution[i]-sqr(muj2)+2.*muk2+2.*muj2-sqr(x3Solution[i])-1.
	 +2.*muj2*muk2-sqr(muk2)-2.*muk2*x3Solution[i]-2.*muj2*x3Solution[i]);
      // check discrim2 is > 0
      if( discrim2 < ZERO) continue;
      double fact1 =2.*sqr(x3Solution[i])-4.*muk2-6.*x3Solution[i]+4.*muj2-xT2*x3Solution[i]
	+2.*xT2-2.*muj2*x3Solution[i]+2.*muk2*x3Solution[i]+4.;
      double fact2 = (4.-4.*x3Solution[i]+xT2);
      double discriminant = sqrt(discrim2);
      // two solns for x1
      x1Solution[i][0] = (fact1 + 2.*discriminant)/fact2;
      x1Solution[i][1] = (fact1 - 2.*discriminant)/fact2;
      bool found = false;
      for(unsigned int j=0;j<2;++j) {
	// calculate x2
	x2Solution[i][j] = 2.-x3Solution[i]-x1Solution[i][j];
	// set limits on x2
	double root = max(0.,sqr(x1Solution[i][j])-4.*muj2);
	root = sqrt(root);
	double x2Plus  = 1.+muk2-muj2
	  -0.5*(1.-x1Solution[i][j]+muj2-muk2)/(1.-x1Solution[i][j]+muj2)
	  *(x1Solution[i][j]-2.*muj2-root);
	double x2Minus = 1.+muk2-muj2
	  -0.5*(1.-x1Solution[i][j]+muj2-muk2)/(1.-x1Solution[i][j]+muj2)
	  *(x1Solution[i][j]-2.*muj2+root);

        if(x1Solution[i][j]>=x1Minus && x1Solution[i][j]<=x1Plus &&
	   x2Solution[i][j]>=x2Minus && x2Solution[i][j]<=x2Plus &&
           checkZMomenta(x1Solution[i][j], x2Solution[i][j], x3Solution[i], yTemp[i], pT[i], 
			 muj, muk)) {
          probTemp[i][j] = weightPrefactor*pT[i]*
            calculateJacobian(x1Solution[i][j], x2Solution[i][j], pT[i], muj, muk)*
	    calculateRealEmission(x1Solution[i][j], x2Solution[i][j], 
				  hardProcess, phi, muj, muk, i, false);
          found = true;
        }
        else {
          probTemp[i][j] = 0.;
        }
      }
      if(!found) continue;
      // alpha S piece
      double wgt = (probTemp[i][0]+probTemp[i][1])*alphaS()->ratio(sqr(pT[i]));
      // matrix element weight
      reject = UseRandom::rnd()>wgt;
    }
    while(reject);
  } // end of emitter for loop
  // no emission
  if(pT[0]<ZERO&&pT[1]<ZERO) return false;
  //pick the spectator and x1 x2 values
  double x1,x2,y;
  // particle 1 emits, particle 2 spectates
  unsigned int iemit=0;
  if(pT[0]>pT[1]){ 
    pT_ = pT[0];
    y=yTemp[0];
    if(probTemp[0][0]>UseRandom::rnd()*(probTemp[0][0]+probTemp[0][1])) {
      x1 = x1Solution[0][0];
      x2 = x2Solution[0][0];
    }
    else {
      x1 = x1Solution[0][1];
      x2 = x2Solution[0][1];
    }
  }
  // particle 2 emits, particle 1 spectates
  else {
    iemit=1;
    pT_ = pT[1];
    y=yTemp[1];
    if(probTemp[1][0]>UseRandom::rnd()*(probTemp[1][0]+probTemp[1][1])) {
      x1 = x1Solution[1][0];
      x2 = x2Solution[1][0];
    }
    else {
      x1 = x1Solution[1][1];
      x2 = x2Solution[1][1];
    }
  }
  // find spectator
  unsigned int ispect = iemit == 0 ? 1 : 0;
  double muk = iemit == 0 ? mu2_ : mu1_;
  double muk2 = sqr(muk);
  double muj = iemit == 0 ? mu1_ : mu2_;
  double muj2 = sqr(muj);
  double xT2 = sqr(2./mW_*pT_);
  // Find the boost from the lab to the c.o.m with the spectator 
  // along the -z axis, and then invert it.
  LorentzRotation eventFrame( ( quark_[0] + quark_[1] ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*quark_[ispect];
  eventFrame.rotateZ( -spectator.phi() );
  eventFrame.rotateY( -spectator.theta() - Constants::pi );
  eventFrame.invert();
  // spectator
  quark_[ispect].setT( 0.5*x2*mW_ );
  quark_[ispect].setX( ZERO );
  quark_[ispect].setY( ZERO );
  quark_[ispect].setZ( -sqrt(0.25*mw2_*x2*x2-mw2_*muk2) );
  // gluon
  gauge_.setT( pT_*cosh(y)  );
  gauge_.setX( pT_*cos(phi) );
  gauge_.setY( pT_*sin(phi) );
  gauge_.setZ( pT_*sinh(y)  );
  gauge_.setMass(ZERO);
  // emitter
  quark_[iemit].setX( -pT_*cos(phi) );
  quark_[iemit].setY( -pT_*sin(phi) );
  quark_[iemit].setZ(  0.5*mW_*sqrt(sqr(x1)-xT2-4.*muj2) );
  if(sqrt(0.25*mw2_*x2*x2-mw2_*muk2)-pT_*sinh(y)<ZERO)
    quark_[iemit].setZ(-quark_[iemit].z());
  quark_[iemit].setT( 0.5*mW_*x1 );
  // boost constructed vectors into the event frame
  quark_[0] = eventFrame * quark_[0];
  quark_[1] = eventFrame * quark_[1];
  gauge_    = eventFrame * gauge_;
  // need to reset masses because for whatever reason the boost  
  // touches the mass component of the five-vector and can make  
  // zero mass objects acquire a floating point negative mass(!).
  gauge_.setMass( ZERO );
  quark_[iemit] .setMass(partons_[iemit ]->mass());
  quark_[ispect].setMass(partons_[ispect]->mass());

  return true;
}

InvEnergy SMWFermionsPOWHEGDecayer::calculateJacobian(double x1, double x2, Energy pT, 
						      double muj, double muk) const{
  double xPerp = abs(2.*pT/mW_);
  Energy jac = mW_/xPerp*fabs((x2*sqr(muj)+2.*sqr(muk)*x1
			       +sqr(muk)*x2-x1*x2-sqr(x2)+x2)/pow((sqr(x2)-4.*sqr(muk)),1.5));
  
  return 1./jac; //jacobian as defined is dptdy=jac*dx1dx2, therefore we have to divide by it
}

bool SMWFermionsPOWHEGDecayer::checkZMomenta(double x1, double x2, double x3, 
					     double y, Energy pT, double muj,
					     double muk) const {
  double xPerp2 = 4.*pT*pT/mW_/mW_;
  double root1 = sqrt(max(0.,sqr(x2)-4.*sqr(muk)));
  double root2 = sqrt(max(0.,sqr(x1)-xPerp2 - 4.*sqr(muj)));
  static double tolerance = 1e-6; 
  bool isMomentaReconstructed = false;  

  if(pT*sinh(y) > ZERO) {
    if(abs(-root1 + sqrt(sqr(x3)-xPerp2)  + root2) <= tolerance ||
       abs(-root1 + sqrt(sqr(x3)-xPerp2)  - root2)  <= tolerance)
      isMomentaReconstructed=true;
  }
  else if(pT*sinh(y) < ZERO){
    if(abs(-root1 - sqrt(sqr(x3)-xPerp2)  + root2) <= tolerance ||
       abs(-root1 - sqrt(sqr(x3)-xPerp2)  - root2)  <= tolerance)
	isMomentaReconstructed=true;
  }
  else 
    if(abs(-root1+ sqrt(sqr(x1)-xPerp2 - 4.*(muj))) <= tolerance)
      isMomentaReconstructed=true;
      
  return isMomentaReconstructed;
}

double SMWFermionsPOWHEGDecayer::calculateRealEmission(double x1, double x2, 
						       vector<PPtr> hardProcess,
						       double phi, double muj,
						       double muk, int iemit, 
						       bool subtract) const {
  // make partons data object for meRatio
  vector<cPDPtr> partons (3);
  for(int ix=0; ix<3; ++ix)
    partons[ix] = hardProcess[ix]->dataPtr();
  partons.push_back(gluon_);
  // calculate x3
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3)-0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1)-4.*sqr(muk)+4.*sqr(muj))
		       /(sqr(x2)-4.*sqr(muk))));
  // calculate the momenta
  Energy M = mW_;
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*sqr(muk),0.)),
			  0.5*M*x2,M*muk); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*sqr(muj),0.)),
			  0.5*M*x1,M*muj);
  Lorentz5Momentum pgluon(0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // loop over the possible emitting partons
  double realwgt(0.);

  // boost and rotate momenta
  LorentzRotation eventFrame( ( hardProcess[1]->momentum() +
				hardProcess[2]->momentum() ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*hardProcess[iemit+1]->momentum();
  eventFrame.rotateZ( -spectator.phi()    );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  vector<Lorentz5Momentum> momenta(3);
  momenta[0]   = hardProcess[0]->momentum();
  if(iemit==0) {
    momenta[2] = eventFrame*pspect;
    momenta[1] = eventFrame*pemit ;
  }
  else {
    momenta[1] = eventFrame*pspect;
    momenta[2] = eventFrame*pemit ;
  }
  momenta.push_back(eventFrame*pgluon);
  // calculate the weight
  if(1.-x1>1e-5 && 1.-x2>1e-5) 
    realwgt += meRatio(partons,momenta,iemit,subtract);
  
  return realwgt;
}
