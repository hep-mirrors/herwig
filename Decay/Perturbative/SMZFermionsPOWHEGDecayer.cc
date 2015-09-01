//-*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMZFermionsPOWHEGDecayer class.
//

#include "SMZFermionsPOWHEGDecayer.h"
#include <numeric>
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Shower/Base/HardTree.h"
#include "Herwig/Shower/Base/ShowerTree.h"
#include "Herwig/Shower/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "Herwig/Shower/Base/Branching.h"

using namespace Herwig;

SMZFermionsPOWHEGDecayer::SMZFermionsPOWHEGDecayer() 
  : CF_(4./3.), pTmin_(1.*GeV)
{  }

IBPtr SMZFermionsPOWHEGDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SMZFermionsPOWHEGDecayer::fullclone() const {
  return new_ptr(*this);
}

void SMZFermionsPOWHEGDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFGVertex_ << FFZVertex_ << gluon_ << ounit( pTmin_, GeV );
}

void SMZFermionsPOWHEGDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFGVertex_ >> FFZVertex_ >> gluon_ >> iunit( pTmin_, GeV );
}

ClassDescription<SMZFermionsPOWHEGDecayer> 
SMZFermionsPOWHEGDecayer::initSMZFermionsPOWHEGDecayer;
// Definition of the static class description member.

void SMZFermionsPOWHEGDecayer::Init() {

  static ClassDocumentation<SMZFermionsPOWHEGDecayer> documentation
    ("There is no documentation for the SMZFermionsPOWHEGDecayer class");

  static Parameter<SMZFermionsPOWHEGDecayer, Energy> interfacePtMin
    ("minpT",
     "The pt cut on hardest emision generation",
     &SMZFermionsPOWHEGDecayer::pTmin_, GeV, 1.*GeV, 0*GeV, 100000.0*GeV,
     false, false, Interface::limited);
}

HardTreePtr SMZFermionsPOWHEGDecayer::
generateHardest(ShowerTreePtr tree) {
  // Get the progenitors: Q and Qbar.
  ShowerProgenitorPtr 
    QProgenitor    = tree->outgoingLines().begin()->first,
    QbarProgenitor = tree->outgoingLines().rbegin()->first;
  if(QProgenitor->id()<0) swap( QProgenitor, QbarProgenitor );
  partons_.resize(2);
  partons_[0] = QProgenitor->progenitor()   ->dataPtr();
  partons_[1] = QbarProgenitor->progenitor()->dataPtr();
  if(!partons_[0]->coloured()) return HardTreePtr();
  // momentum of the partons
  quark_.resize(2);
  quark_[0] = QProgenitor   ->copy()->momentum();
  quark_[1] = QbarProgenitor->copy()->momentum();
  // Set the existing mass entries in partons 5 vectors with the
  // once and for all.
  quark_[0].setMass(partons_[0]->mass());
  quark_[1].setMass(partons_[1]->mass());
  gauge_.setMass(0.*MeV);
  // Get the Z boson.
  zboson_ = tree->incomingLines().begin()->first->copy();
  // copy the particle objects
  vector<PPtr> hardProcess (3);
  hardProcess[0] = zboson_;
  hardProcess[1] = QbarProgenitor->copy();
  hardProcess[2] = QProgenitor   ->copy();
  // Get the Z boson mass.
  mz2_ = (quark_[0] + quark_[1]).m2();
  // Generate emission and set _quark[0,1] and _gauge to be the 
  // momenta of q, qbar and g after the hardest emission:
  if(!getEvent(hardProcess)) {
    QProgenitor   ->maximumpT(pTmin_,ShowerInteraction::QCD);
    QbarProgenitor->maximumpT(pTmin_,ShowerInteraction::QCD);
    return HardTreePtr();
  }
  // Ensure the energies are greater than the constituent masses:
  for (int i=0; i<2; i++) {
    if (quark_[i].e() < partons_[i]->constituentMass()) return HardTreePtr();
    if (gauge_.e()    < gluon_     ->constituentMass()) return HardTreePtr();
  }
  // set masses
  quark_[0].setMass( partons_[0]->mass() );
  quark_[1].setMass( partons_[1]->mass() );
  gauge_   .setMass( ZERO );
  // assign the emitter based on evolution scales
  unsigned int iemitter   = quark_[0]*gauge_ > quark_[1]*gauge_ ? 1 : 0;
  unsigned int ispectator = iemitter==1                         ? 0 : 1;
  // Make the particles for the HardTree:
  ShowerParticlePtr emitter  (new_ptr(ShowerParticle(partons_[iemitter  ],true)));
  ShowerParticlePtr spectator(new_ptr(ShowerParticle(partons_[ispectator],true)));
  ShowerParticlePtr gauge    (new_ptr(ShowerParticle(gluon_              ,true)));
  ShowerParticlePtr zboson   (new_ptr(ShowerParticle(zboson_->dataPtr() ,false)));
  ShowerParticlePtr parent   (new_ptr(ShowerParticle(partons_[iemitter  ],true)));
  emitter  ->set5Momentum(quark_[iemitter  ] ); 
  spectator->set5Momentum(quark_[ispectator] );  
  gauge    ->set5Momentum(gauge_             ); 
  zboson   ->set5Momentum(zboson_->momentum());  
  Lorentz5Momentum parentMomentum(quark_[iemitter]+gauge_);
  parentMomentum.rescaleMass();
  parent->set5Momentum(parentMomentum);
  // Create the vectors of HardBranchings to create the HardTree:
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  // Incoming boson:
  spaceBranchings.push_back(new_ptr(HardBranching(zboson,SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  // Outgoing particles from hard emission:
  HardBranchingPtr spectatorBranch(new_ptr(HardBranching(spectator,SudakovPtr(),
							 HardBranchingPtr(),
							 HardBranching::Outgoing)));
  HardBranchingPtr emitterBranch(new_ptr(HardBranching(parent,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(emitter,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(gauge,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  emitterBranch->type(emitterBranch->branchingParticle()->id()>0 ? 
		      ShowerPartnerType::QCDColourLine : 
		      ShowerPartnerType::QCDAntiColourLine);
  allBranchings.push_back(emitterBranch);
  allBranchings.push_back(spectatorBranch);
  if(iemitter==1) swap(allBranchings[0],allBranchings[1]);
  // Add incoming boson to allBranchings
  allBranchings.push_back( spaceBranchings.back() );
  // Make the HardTree from the HardBranching vectors.
  HardTreePtr hardtree = new_ptr(HardTree(allBranchings,spaceBranchings,
					   ShowerInteraction::QCD));
  // Set the maximum pt for all other emissions
  Energy ptveto(pT_);
  QProgenitor   ->maximumpT(ptveto,ShowerInteraction::QCD);
  QbarProgenitor->maximumpT(ptveto,ShowerInteraction::QCD);
  // Connect the particles with the branchings in the HardTree
  hardtree->connect( QProgenitor->progenitor()   , allBranchings[0] );
  hardtree->connect( QbarProgenitor->progenitor(), allBranchings[1] );
  // colour flow
  ColinePtr newline=new_ptr(ColourLine());
  for(set<HardBranchingPtr>::const_iterator cit=hardtree->branchings().begin();
      cit!=hardtree->branchings().end();++cit) {
    if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3)
      newline->addColoured((**cit).branchingParticle());
    else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar)
      newline->addAntiColoured((**cit).branchingParticle());
  }
  ColinePtr newLine2=new_ptr(ColourLine());
  if(emitterBranch->branchingParticle()->dataPtr()->iColour()==PDT::Colour3) {
    emitterBranch->branchingParticle()->colourLine()->addColoured(gauge);
    newLine2->addColoured(emitter);
    newLine2->addAntiColoured(gauge);
  }
  else {
    emitterBranch->branchingParticle()->antiColourLine()->addAntiColoured(gauge);
    newLine2->addAntiColoured(emitter);
    newLine2->addColoured(gauge);
  }
  // return the tree
  return hardtree;
}

double SMZFermionsPOWHEGDecayer::
me2(const int ichan, const Particle & part,
    const ParticleVector & decay, MEOption meopt) const {
  // leading-order result
  double output = SMZDecayer::me2(ichan,part,decay,meopt);
  // check decay products coloured, otherwise return
  if(!decay[0]->dataPtr()->coloured()) return output;
  // inital masses, couplings  etc
  // fermion mass
  Energy particleMass = decay[0]->dataPtr()->mass();
  // Z mass
  mZ_ = part.mass();
  // strong coupling
  aS_ = SM().alphaS(sqr(mZ_));
  // reduced mass
  mu_  = particleMass/mZ_;
  mu2_ = sqr(mu_);
  // scale
  scale_ = sqr(mZ_);
  // cast the vertices
  tcFFVVertexPtr Zvertex = dynamic_ptr_cast<tcFFVVertexPtr>(FFZVertex());
  // compute the spinors
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  vector<VectorWaveFunction>    vin;
  SpinorBarWaveFunction qkout(decay[0]->momentum(),decay[0]->dataPtr(),outgoing);
  SpinorWaveFunction    qbout(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  VectorWaveFunction    zin  (part.momentum()     ,part.dataPtr()     ,incoming);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  for(unsigned int ix=0;ix<3;++ix){
    zin.reset(ix);
    vin.push_back(zin);
  }
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  double total=0.;
  if(mu_!=0.) {
    LorentzPolarizationVector momDiff = 
      (decay[0]->momentum()-decay[1]->momentum())/2./
      (decay[0]->momentum().mass()+decay[1]->momentum().mass());
    // scalars
    Complex scalar1 = zin.wave().dot(momDiff);
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {		
	for(unsigned int inhel=0;inhel<3;++inhel) {
	  // first the LO bit
	  Complex diag1 = FFZVertex()->evaluate(scale_,aout[outhel2],fout[outhel1],vin[inhel]);
	  // extra stuff for NLO 
	  LorentzPolarizationVector left  = 
	    aout[outhel2].wave().leftCurrent(fout[outhel1].wave());
	  LorentzPolarizationVector right = 
	    aout[outhel2].wave().rightCurrent(fout[outhel1].wave());
	  Complex scalar = 
	    aout[outhel2].wave().scalar(fout[outhel1].wave());
	  // nlo specific pieces
	  Complex diag3 =
	    Complex(0.,1.)*Zvertex->norm()*
	    (Zvertex->right()*( left.dot(zin.wave())) +
	     Zvertex-> left()*(right.dot(zin.wave())) -
	     ( Zvertex-> left()+Zvertex->right())*scalar1*scalar);
	  // nlo piece
	  total += real(diag1*conj(diag3) + diag3*conj(diag1));
	}
      }
    }
    // rescale
    total *= UnitRemoval::E2/scale_;
  }
  else {
    total = ZERO;
  }
  // now for the NLO bit
  double mu4 = sqr(mu2_);
  double lmu = mu_!=0. ? log(mu_) : 0.;
  double v = sqrt(1.-4.*mu2_),v2(sqr(v));
  double omv = 4.*mu2_/(1.+v);
  double f1,f2,fNS,VNS;
  double r = omv/(1.+v);
  double lr = mu_!=0. ? log(r) : 0.;
  // normal form
  if(mu_>1e-4) {
    f1 = CF_*aS_/Constants::pi*
      ( +1. + 3.*log(0.5*(1.+v)) - 1.5*log(0.5*(1.+v2)) + sqr(Constants::pi)/6.
	- 0.5*sqr(lr) - (1.+v2)/v*(lr*log(1.+v2) + sqr(Constants::pi)/12. 
				       -0.5*log(4.*mu2_)*lr + 0.25*sqr(lr)));
    fNS = -0.5*(1.+2.*v2)*lr/v + 1.5*lr - 2./3.*sqr(Constants::pi) + 0.5*sqr(lr)
      + (1.+v2)/v*(Herwig::Math::ReLi2(r) + sqr(Constants::pi)/3. - 0.25*sqr(lr) 
		   + lr*log((2.*v/ (1.+v))));
    VNS = 1.5*log(0.5*(1.+v2)) 
      + 0.5*(1.+v2)/v*( 2.*lr*log(2.*(1.+v2)/sqr(1.+v))  
			+ 2.*Herwig::Math::ReLi2(sqr(r)) 
			- 2.*Herwig::Math::ReLi2(2.*v/(1.+v)) - sqr(Constants::pi)/6.)
      + log(1.-mu_) - 2.*log(1.-2.*mu_) - 4.*mu2_/(1.+v2)*log(mu_/(1.-mu_)) 
      - mu_/(1.-mu_)
      + 4.*(2.*mu2_-mu_)/(1.+v2) + 0.5*sqr(Constants::pi); 
    f2 = CF_*aS_/Constants::pi*mu2_*lr/v;
  }
  // small mass limit
  else {
    f1 = -CF_*aS_/Constants::pi/6.*
      ( - 6. - 24.*lmu*mu2_ - 15.*mu4 - 12.*mu4*lmu - 24.*mu4*sqr(lmu) 
	+ 2.*mu4*sqr(Constants::pi) - 12.*mu2_*mu4 - 96.*mu2_*mu4*sqr(lmu) 
	+ 8.*mu2_*mu4*sqr(Constants::pi) - 80.*mu2_*mu4*lmu);
    fNS = - mu2_/18.*( + 36.*lmu - 36. - 45.*mu2_ + 216.*lmu*mu2_ - 24.*mu2_*sqr(Constants::pi) 
		      + 72.*mu2_*sqr(lmu) - 22.*mu4 + 1032.*mu4 * lmu
		      - 96.*mu4*sqr(Constants::pi) + 288.*mu4*sqr(lmu));
    VNS = - mu2_/1260.*(-6930. + 7560.*lmu + 2520.*mu_ - 16695.*mu2_ 
		       + 1260.*mu2_*sqr(Constants::pi) 
		       + 12600.*lmu*mu2_ + 1344.*mu_*mu2_ - 52780.*mu4 + 36960.*mu4*lmu 
		       + 5040.*mu4*sqr(Constants::pi) - 12216.*mu_*mu4);
    f2 = CF_*aS_*mu2_/Constants::pi*( 2.*lmu + 4.*mu2_*lmu + 2.*mu2_ + 12.*mu4*lmu + 7.*mu4);
  }
  // add up bits for f1
  f1 += CF_*aS_/Constants::pi*(fNS+VNS);
  // now for the real correction
  double phi  = UseRandom::rnd()*Constants::twopi;
  // calculate y
  double yminus = 0.; 
  double yplus  = 1.-2.*mu_*(1.-mu_)/(1.-2*mu2_);
  double y = yminus + UseRandom::rnd()*(yplus-yminus);
  // calculate z
  double v1  = sqrt(sqr(2.*mu2_+(1.-2.*mu2_)*(1.-y))-4.*mu2_)/(1.-2.*mu2_)/(1.-y);
  double zplus  = (1.+v1)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
  double zminus = (1.-v1)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
  double z = zminus + UseRandom::rnd()*(zplus-zminus);
  double jac = (1.-y)*(yplus-yminus)*(zplus-zminus);
  // calculate x1,x2
  double x2 = 1. - y*(1.-2.*mu2_);
  double x1 = 1. - z*(x2-2.*mu2_);
  // copy the particle objects over for calculateRealEmission
  vector<PPtr> hardProcess(3);
  hardProcess[0] = const_ptr_cast<PPtr>(&part);
  hardProcess[1] = decay[0];
  hardProcess[2] = decay[1];
  // total real emission contribution
  double realFact = 0.25*jac*sqr(1.-2.*mu2_)/
    sqrt(1.-4.*mu2_)/Constants::twopi
    *2.*CF_*aS_*calculateRealEmission(x1, x2, hardProcess, phi, true);
  // the born + virtual + real
  output = output*(1. + f1 + realFact) + f2*total;
  return output;
}

double SMZFermionsPOWHEGDecayer::meRatio(vector<cPDPtr> partons, 
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

double SMZFermionsPOWHEGDecayer::loME(const vector<cPDPtr> & partons, 
			const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<VectorWaveFunction>    vin;
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  VectorWaveFunction    zin  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  for(unsigned int ix=0;ix<3;++ix){
    zin.reset(ix);
    vin.push_back(zin);
  }
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  double total(0.);
  for(unsigned int inhel=0;inhel<3;++inhel) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	Complex diag1 = FFZVertex()->evaluate(scale_,aout[outhel2],fout[outhel1],vin[inhel]);
	total += norm(diag1);
      }
    }
  }
  // return the answer
  return total;
}
 
InvEnergy2 SMZFermionsPOWHEGDecayer::realME(const vector<cPDPtr> & partons, 
			      const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<VectorWaveFunction>     vin;
  vector<SpinorWaveFunction>     aout;
  vector<SpinorBarWaveFunction>  fout;
  vector<VectorWaveFunction>     gout;
  VectorWaveFunction    zin  (momenta[0],partons[0],incoming);
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
    zin.reset(ix);
    vin.push_back(zin);
  }
  vector<Complex> diag(2,0.);

  double total(0.);
  for(unsigned int inhel1=0;inhel1<3;++inhel1) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	  SpinorBarWaveFunction off1 =
	    FFGVertex()->evaluate(scale_,3,partons[1],fout[outhel1],gout[outhel3]);
	  diag[0] = FFZVertex()->evaluate(scale_,aout[outhel2],off1,vin[inhel1]);

	  SpinorWaveFunction off2 = 
	    FFGVertex()->evaluate(scale_,3,partons[2],aout[outhel2],gout[outhel3]);
	  diag[1] = FFZVertex()->evaluate(scale_,off2,fout[outhel1],vin[inhel1]);

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

void SMZFermionsPOWHEGDecayer::doinit() {
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm) throw InitException() 
	      << "Wrong type of StandardModel object in "
	      << "SMZFermionsPOWHEGDecayer::doinit() "
	      << "the Herwig version must be used." 
	      << Exception::runerror;
  // cast the vertices
  FFZVertex_ = hwsm->vertexFFZ();
  FFZVertex_->init();
  FFGVertex_ = hwsm->vertexFFG();
  FFGVertex_->init();
  SMZDecayer::doinit();
  gluon_ = getParticleData(ParticleID::g);
}

bool SMZFermionsPOWHEGDecayer::getEvent(vector<PPtr> hardProcess) {
  Energy particleMass = ZERO;
  for(unsigned int ix=0;ix<hardProcess.size();++ix) {
    if(hardProcess[ix]->id()==ParticleID::Z0) {
      mZ_ = hardProcess[ix]->mass();
    }
    else {
      particleMass =  hardProcess[ix]->mass();
    }
  }
  // reduced mass
  mu_  = particleMass/mZ_;
  mu2_ = sqr(mu_);
  // scale
  scale_ = sqr(mZ_);
  // max pT
  Energy pTmax = 0.5*sqrt(mz2_);
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
    do {
      //generation of phi
      phi = UseRandom::rnd() * Constants::twopi;
      // reject the emission
      reject = true; 
      // generate pt
      pT[i] *= pow(UseRandom::rnd(),1./prefactor);
      Energy2 pT2 = sqr(pT[i]);
      if(pT[i]<pTmin_) {
        pT[i] = -GeV;
        break;
      }
      // generate y
      yTemp[i] = ymin + UseRandom::rnd()*(ymax-ymin);
      //generate x3 & x1 from pT & y
      double x1Plus  = 1.;
      double x1Minus = 2.*mu_;
      x3Solution[i] = 2.*pT[i]*cosh(yTemp[i])/mZ_;
      // prefactor
      double weightPrefactor = 0.5/sqrt(1.-4.*mu2_)/overEst;
      // calculate x1 & x2 solutions
      Energy4 discrim2 = (sqr(x3Solution[i]*mZ_) - 4.*pT2)*
        (mz2_*(x3Solution[i]-1.)*(4.*mu2_+x3Solution[i]-1.)-4.*mu2_*pT2);
      //check discriminant2 is > 0
      if( discrim2 < ZERO) continue;
      Energy2 discriminant = sqrt(discrim2);
      Energy2 fact1 = 3.*mz2_*x3Solution[i]-2.*mz2_+2.*pT2*x3Solution[i]
	-4.*pT2-mz2_*sqr(x3Solution[i]);
      Energy2 fact2 = 2.*mz2_*(x3Solution[i]-1.)-2.*pT2;
      // two solns for x1
      x1Solution[i][0] = (fact1 + discriminant)/fact2;
      x1Solution[i][1] = (fact1  - discriminant)/fact2;

      bool found = false;
      for(unsigned int j=0;j<2;++j) {
	x2Solution[i][0] = 2.-x3Solution[i]-x1Solution[i][0];
	x2Solution[i][1] = 2.-x3Solution[i]-x1Solution[i][1];
	// set limits on x2
	double root = max(0.,sqr(x1Solution[i][j])-4.*mu2_);
	root = sqrt(root);
	double x2Plus  = 1.-0.5*(1.-x1Solution[i][j])/(1.-x1Solution[i][j]+mu2_)
	  *(x1Solution[i][j]-2.*mu2_-root);
	double x2Minus = 1.-0.5*(1.-x1Solution[i][j])/(1.-x1Solution[i][j]+mu2_)
	  *(x1Solution[i][j]-2.*mu2_+root);
	if(x1Solution[i][j]>=x1Minus && x1Solution[i][j]<=x1Plus &&
	   x2Solution[i][j]>=x2Minus && x2Solution[i][j]<=x2Plus &&
           checkZMomenta(x1Solution[i][j], x2Solution[i][j], x3Solution[i], yTemp[i], pT[i])) {
          probTemp[i][j] = weightPrefactor*pT[i]*
            calculateJacobian(x1Solution[i][j], x2Solution[i][j], pT[i])*
	    calculateRealEmission(x1Solution[i][j], x2Solution[i][j], 
				  hardProcess, phi, false, i);
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
  } //end of emitter for loop
  // no emission
  if(pT[0]<ZERO&&pT[1]<ZERO) return false;
  //pick the spectator and x1 x2 values
  double x1,x2,y;
  //particle 1 emits, particle 2 spectates
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
  // Find the boost from the lab to the c.o.m with the spectator 
  // along the -z axis, and then invert it.
  LorentzRotation eventFrame( ( quark_[0] + quark_[1] ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*quark_[ispect];
  eventFrame.rotateZ( -spectator.phi() );
  eventFrame.rotateY( -spectator.theta() - Constants::pi );
  eventFrame.invert();
  // spectator
  quark_[ispect].setT( 0.5*x2*mZ_ );
  quark_[ispect].setX( ZERO );
  quark_[ispect].setY( ZERO );
  quark_[ispect].setZ( -sqrt(0.25*mz2_*x2*x2-mz2_*mu2_) );
  // gluon
  gauge_.setT( pT_*cosh(y)  );
  gauge_.setX( pT_*cos(phi) );
  gauge_.setY( pT_*sin(phi) );
  gauge_.setZ( pT_*sinh(y)  );
  gauge_.setMass(ZERO);
  // emitter reconstructed from gluon & spectator
  quark_[iemit] = - gauge_ - quark_[ispect];
  quark_[iemit].setT( 0.5*mZ_*x1 );
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

InvEnergy SMZFermionsPOWHEGDecayer::calculateJacobian(double x1, double x2, Energy pT) const{
  double xPerp = abs(2.*pT/mZ_);
  Energy jac = mZ_*fabs((x1*x2-2.*mu2_*(x1+x2)+sqr(x2)-x2)/xPerp/pow(sqr(x2)-4.*mu2_,1.5));   
  return 1./jac; //jacobian as defined is dptdy=jac*dx1dx2, therefore we have to divide by it
}

bool SMZFermionsPOWHEGDecayer::checkZMomenta(double x1, double x2, double x3, double y, Energy pT) const {
  double xPerp2 = 4.*pT*pT/mZ_/mZ_;
  static double tolerance = 1e-6; 
  bool isMomentaReconstructed = false;  

  if(pT*sinh(y)>ZERO) {
    if(abs(-sqrt(sqr(x2)-4.*mu2_)+sqrt(sqr(x3)-xPerp2) + sqrt(sqr(x1)-xPerp2 - 4.*mu2_)) <= tolerance ||
       abs(-sqrt(sqr(x2)-4.*mu2_)+sqrt(sqr(x3)-xPerp2)  - sqrt(sqr(x1)-xPerp2 - 4.*mu2_))  <= tolerance) isMomentaReconstructed=true;
  }
  else if(pT*sinh(y) < ZERO){
      if(abs(-sqrt(sqr(x2)-4.*mu2_)-sqrt(sqr(x3)-xPerp2) + sqrt(sqr(x1)-xPerp2 - 4.*mu2_)) <= tolerance ||
         abs(-sqrt(sqr(x2)-4.*mu2_)-sqrt(sqr(x3)-xPerp2)  - sqrt(sqr(x1)-xPerp2 - 4.*mu2_))  <= tolerance) isMomentaReconstructed=true;
  }
  else 
    if(abs(-sqrt(sqr(x2)-4.*mu2_)+ sqrt(sqr(x1)-xPerp2 - 4.*mu2_)) <= tolerance) isMomentaReconstructed=true;
      
  return isMomentaReconstructed;
}

double SMZFermionsPOWHEGDecayer::calculateRealEmission(double x1, double x2, 
						       vector<PPtr> hardProcess,
						       double phi,
						       bool subtract) const {
  // make partons data object for meRatio
  vector<cPDPtr> partons (3);
  for(int ix=0; ix<3; ++ix)
    partons[ix] = hardProcess[ix]->dataPtr();
  partons.push_back(gluon_);
  // calculate x3
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3) -0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1))/(sqr(x2)-4.*mu2_)));
  // calculate the momenta
  Energy M = mZ_;
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*mu2_,0.)),0.5*M*x2,M*mu_); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
  			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*mu2_,0.)),0.5*M*x1,M*mu_);
  Lorentz5Momentum pgluon(0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // loop over the possible emitting partons
  double realwgt(0.);
  for(unsigned int iemit=0;iemit<2;++iemit) {
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
  }
  
  // total real emission contribution
  return realwgt;
}

double SMZFermionsPOWHEGDecayer::calculateRealEmission(double x1, double x2, 
						       vector<PPtr> hardProcess,
						       double phi,
						       bool subtract,
						       int emitter) const {
  // make partons data object for meRatio
  vector<cPDPtr> partons (3);
  for(int ix=0; ix<3; ++ix)
    partons[ix] = hardProcess[ix]->dataPtr();
  partons.push_back(gluon_);
  // calculate x3
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3) -0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1))/(sqr(x2)-4.*mu2_)));
  // calculate the momenta
  Energy M = mZ_;
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*mu2_,0.)),0.5*M*x2,M*mu_); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
  			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*mu2_,0.)),0.5*M*x1,M*mu_);
  Lorentz5Momentum pgluon( 0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
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
  Lorentz5Momentum spectator = eventFrame*hardProcess[emitter+1]->momentum();
  eventFrame.rotateZ( -spectator.phi()    );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  vector<Lorentz5Momentum> momenta(3);
  momenta[0]   = hardProcess[0]->momentum();
  if(emitter==0) {
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
    realwgt += meRatio(partons,momenta,emitter,subtract);  
  // total real emission contribution
  return realwgt;
}
