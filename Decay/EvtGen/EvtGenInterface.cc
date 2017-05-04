// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EvtGenInterface class.
//

#include "EvtGenInterface.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenExternal/EvtExternalGenList.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtRaritaSchwingerParticle.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtStringParticle.hh"
#include "EvtGenBase/EvtHighSpinParticle.hh"
#include "EvtGenBase/EvtDecayTable.hh"

#ifndef EVTGEN_PREFIX
#error Makefile.am needs to define EVTGEN_PREFIX
#endif

using namespace Herwig;

namespace {

const string prefix=EVTGEN_PREFIX "";
const string p8data=PYTHIA8DATA "";

}

EvtGenInterface::EvtGenInterface() : decayName_(prefix+"/share/DECAY_2010.DEC"),
				     pdtName_(prefix+"/share/evt.pdl"),
				     reDirect_(true), checkConv_(false),
				     p8Data_(p8data)
{}

EvtGenInterface::EvtGenInterface(const EvtGenInterface & x)
  : Interfaced(x), decayName_(x.decayName_), pdtName_(x.pdtName_),
    userDecays_(x.userDecays_), reDirect_(x.reDirect_),
    checkConv_(x.checkConv_), convID_(x.convID_),
    p8Data_(x.p8Data_), evtrnd_(x.evtrnd_),evtgen_(x.evtgen_)
{}

EvtGenInterface::~EvtGenInterface() {}

IBPtr EvtGenInterface::clone() const {
  return new_ptr(*this);
}

IBPtr EvtGenInterface::fullclone() const {
  return new_ptr(*this);
}

void EvtGenInterface::persistentOutput(PersistentOStream & os) const {
  os << decayName_ << pdtName_ << reDirect_
     << userDecays_ << checkConv_ << convID_ << p8Data_;
}

void EvtGenInterface::persistentInput(PersistentIStream & is, int) {
  is >> decayName_ >> pdtName_ >> reDirect_
     >> userDecays_ >> checkConv_ >> convID_ >> p8Data_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EvtGenInterface,Interfaced>
describeHerwigEvtGenInterface("Herwig::EvtGenInterface",
			      "HwEvtGenInterface.so");

void EvtGenInterface::Init() {

  static ClassDocumentation<EvtGenInterface> documentation
    ("The EvtGenInterface class is the main class for the use of the EvtGen "
     "decay package with Herwig");

  static Parameter<EvtGenInterface,string> interfaceDecay_File
    ("Decay_File",
     "The name of the file for the EvtGen decays.",
     &EvtGenInterface::decayName_, prefix+"/share/DECAY_2010.DEC",
     false, false);

  static Parameter<EvtGenInterface,string> interfacePDT_File
    ("PDTFile",
     "The name of the file for the EvtGen particle data.",
     &EvtGenInterface::pdtName_, prefix+"/share/evt.pdl",
     false, false);

  static Switch<EvtGenInterface,bool> interfaceRedirect
    ("Redirect",
     "By default cerr and cout are redirected when EvtGen is running to"
     " allow us to catch errors, due to EvtGen's poor internal error "
     "handling. This can be a problem for debugging and so can be switched off.",
     &EvtGenInterface::reDirect_, true, false, false);
  static SwitchOption interfaceRedirectYes
    (interfaceRedirect,
     "Yes",
     "Redirect the output",
     true);
  static SwitchOption interfaceRedirectNo
    (interfaceRedirect,
     "No",
     "Don't redirect the output",
     false);

  static Switch<EvtGenInterface,bool> interfaceCheckConversion
    ("CheckConversion",
     "Check the conversion of particles to and from EvtGen",
     &EvtGenInterface::checkConv_, false, false, false);
  static SwitchOption interfaceCheckConversionfalse
    (interfaceCheckConversion,
     "No",
     "Don't check the conversion",
     false);
  static SwitchOption interfaceCheckConversionCheck
    (interfaceCheckConversion,
     "Yes",
     "Check the conversion",
     true);

  static ParVector<EvtGenInterface,long> interfaceOutputModes
    ("OutputModes",
     "Particles for which to output the EvtGen decay modes"
     " so they can be read into Herwig",
     &EvtGenInterface::convID_, -1, long(0), 0, 0,
     false, false, Interface::nolimits);

  static ParVector<EvtGenInterface,string> interfaceUserDecays
    ("UserDecays",
     "List of user decay files to be loaded",
     &EvtGenInterface::userDecays_, -1, "", "", "",
     false, false, Interface::nolimits);

  static Parameter<EvtGenInterface,string> interfacePythia8Data
    ("Pythia8Data",
     "Location of the Pythia8 data directory",
     &EvtGenInterface::p8Data_, p8data,
     false, false);

}

void EvtGenInterface::doinitrun() {
  Interfaced::doinitrun();
  // redirect cerr and cout if needed
  std::streambuf *temp[2]={cout.rdbuf(),cerr.rdbuf()};
  if(reDirect_ && ! generator()->useStdOut() ) {
    const string f = generator()->filename() + "-EvtGen.log";
    logFile_.open(f);
    std::streambuf *psbuf = logFile_.rdbuf();
    cout.rdbuf(psbuf);
    cerr.rdbuf(psbuf);
  }
  // output the EvtGen initialization info to the log file
  cout << "Initializing EvtGen \n";
  // set up the random number generator for EvtGen
  cout << "Setting EvtGen random number generator"
  		     << " to the Herwig one.\n";
  evtrnd_ = new EvtGenRandom(const_ptr_cast<Ptr<RandomGenerator>::pointer>
  			     (&(UseRandom::current())));
  EvtRandom::setRandomEngine(evtrnd_);
  // PHOTOS STUFF
  EvtAbsRadCorr* radCorrEngine = 0;
  std::list<EvtDecayBase*> extraModels;
  EvtExternalGenList genList(true,p8Data_,"gamma",true);
  radCorrEngine = genList.getPhotosModel();
  extraModels = genList.getListOfModels();
  // Initialize EvtGen
  evtgen_ = new EvtGen(decayName_.c_str(),pdtName_.c_str(),
  		       evtrnd_, radCorrEngine, &extraModels);
  // additional user decays
  for(unsigned int ix=0;ix<userDecays_.size();++ix)
    evtgen_->readUDecay(userDecays_[ix].c_str());
  // check the conversion of particles if needed
  if(checkConv_) checkConversion();
  // convert any particle decay modes which are needed
  for(unsigned int ix=0;ix<convID_.size();++ix) {
    outputEvtGenDecays(convID_[ix]);
  }
  // that's the lot
  cout << "Finished initialisation of EvtGen \n";
  if(reDirect_  && ! generator()->useStdOut() ) {
    cout.rdbuf(temp[0]);
    cerr.rdbuf(temp[1]);
  }
}

void EvtGenInterface::dofinish() {
  Interfaced::dofinish();
  if(logFile_.is_open())
    logFile_.close();
}

ParticleVector EvtGenInterface::decay(const Particle &parent,
				      bool recursive, const DecayMode & dm) const {
  // redirect cout to the log file
  ostringstream stemp;
  std::streambuf *temp[2]={cout.rdbuf(),cerr.rdbuf()};
  if(reDirect_ && ! generator()->useStdOut() ) {
    std::streambuf *psbuf[2] ={logFile_.rdbuf(),stemp.rdbuf()};
    cout.rdbuf(psbuf[0]);
    cerr.rdbuf(psbuf[1]);
  }
  // generate the decay
  ParticleVector output;
  EvtId parID(EvtGenID(parent.id()));
  try {
    EvtDecayBase *decayer = NULL;
    // create evtgen particle for the parent
    EvtParticle* particle = EvtGenParticle(parent);
    EvtParticle* evtParent = NULL;
    // special if parent is the same to avoid doing mixing twice
    if(parent.parents().size()==1 && abs(parent.parents()[0]->id())==abs(parent.id())) {
      evtParent = EvtGenParticle(*parent.parents()[0]);
      particle->addDaug(evtParent);
    }
    // simplest case, evtgen selects mode and recursively decays, use their standard mechanism
    if(dm.wildProductMatcher() && recursive) {
      evtgen_->generateDecay(particle);
    }
    // otherwise we're in control
    else {
      // select the EvtGen decayer to use
      if(dm.wildProductMatcher()) {
  	decayer = EvtDecayTable::getInstance()->getDecayFunc(particle);
      }
      // otherwise we should pick one
      else {
  	int imode(EvtGenChannel(dm));
  	decayer = EvtDecayTable::getInstance()->getDecay(parID.getAlias(),imode);
  	particle->setChannel(imode);
      }
      // must be a decayer
      if(!decayer) throw Exception() << "Could find EvtGen decayer in EvtGen::decay()" 
  				     << Exception::runerror;
      // masses of the children
      bool massTreeOK(true);
      if ( particle->getNDaug() == 0 ) 
  	massTreeOK = particle->generateMassTree();
      if ( ! massTreeOK ) {
  	// delete the EvtGen particle
  	particle->deleteDaughters();
  	delete particle;
  	particle = 0;
	if(evtParent) {
	  delete evtParent;
	}
  	throw Exception() << "EvtGen could not decay " << EvtPDL::name(particle->getId())
  			  <<" with mass "<< particle->mass()
  			  <<" to decay channel number "<< particle->getChannel() << "\n"
  			  << Exception::eventerror;
      }
      assert(decayer !=0 );
      // perform decay
      decayer->makeDecay(particle,recursive);
    }
    // cast the decayer
    EvtDecayAmp *damp  = dynamic_cast<EvtDecayAmp*>(decayer);
    // translate the decay products
    output=decayProducts(particle);
    // set spin information if needed
    tSpinPtr pspin(getSpinInfo(parent));
    // if has spin information translate it
    if(damp) {
      ParticleVector products;
      unsigned int nbeforerad=decayer->getNDaug();
      for(unsigned int ix=0;ix<nbeforerad;++ix) products.push_back(output[ix]);
      constructVertex(parent,products,damp);
    }
    // otherwise
    else {
      pspin->develop();
      RhoDMatrix rhotemp(pspin->iSpin());
      pspin->DMatrix()=rhotemp;
    }
    pspin->decay();
    if(recursive) {
      pspin->developed();
      pspin->DMatrix()=ThePEGSpinDensity(particle->getSpinDensityBackward(),
    					 ThePEGID(particle->getId()));
    }
    // delete the EvtGen particle
    particle->deleteDaughters();
    delete particle; 
    if(evtParent) {
      delete evtParent;
    }
    particle = 0;
  }
  catch ( ... ) {
    // change stream back
    if(reDirect_ && ! generator()->useStdOut() ) {
      cout.rdbuf(temp[0]);
      cerr.rdbuf(temp[1]);
    }
    throw;
  }
  // change stream back
  if(reDirect_ && ! generator()->useStdOut() ) {
    cout.rdbuf(temp[0]);
    cerr.rdbuf(temp[1]);
    string stemp2=stemp.str();
    if(stemp2.length()>0)
      throw Exception() << "EvtGen report error in EvtGen::decay "
    			<< "killing event\n"
    			<< "Error was " << stemp2
    			<< Exception::eventerror;
  }
  // return the decay products
  return output;
}

EvtParticle * EvtGenInterface::EvtGenParticle(const Particle & part) const {
  // convert the momentum
  Lorentz5Momentum inmom(part.momentum());
  EvtVector4R p4(EvtGenMomentum(inmom));
  EvtParticle *evtpart;
  // EvtGen ID and spin type
  EvtId id = EvtGenID(part.id());
  EvtSpinType::spintype thisSpin=EvtPDL::getSpinType(id);
  // boost to the rest frame to get the basis states right
  PPtr decay(const_ptr_cast<PPtr>(new_ptr(part)));
  decay->transform(LorentzRotation(-inmom.boostVector()));
  // get spin information so can transfer to EvtGen
  tcSpinPtr spin(dynamic_ptr_cast<tcSpinPtr>(part.spinInfo()));
  if(spin) spin->decay();
  // don't convert neutrinos, aren't decays so waste of time
  if         ( thisSpin == EvtSpinType::NEUTRINO) {
    throw Exception() << "Tried to convert a neutrino to EvtGen in "
  		      << "EvtGen::EvtParticle. This should not be needed as"
  		      << "EvtGen does not decay neutrinos" 
  		      << Exception::eventerror;
  }
  // don't convert photons, aren't decays so waste of time
  else if ( thisSpin == EvtSpinType::PHOTON ) {
    throw Exception() << "Tried to convert a photon to EvtGen in "
  		      << "EvtGen::EvtParticle. This should not be needed as"
  		      << "EvtGen does not decay photons" 
  		      << Exception::eventerror;
  }
  // scalar particles
  else if ( thisSpin == EvtSpinType::SCALAR ) {
    // create particle
    EvtScalarParticle * myPart = new EvtScalarParticle;
    myPart->init(id, p4);
    // set rho matrix if needed
    tcScalarSpinPtr sp(dynamic_ptr_cast<tcScalarSpinPtr>(spin));
    if(sp) {
      myPart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
    }
    else {
      EvtSpinDensity rho;
      rho.setDiag(EvtSpinType::getSpinStates(thisSpin));
      myPart->setSpinDensityForward(rho);
    }
    evtpart=myPart;
  }
  else if ( thisSpin == EvtSpinType::DIRAC ) {
    EvtDiracParticle *myPart =new EvtDiracParticle;
    // get the spin info
    tcFermionSpinPtr sp(dynamic_ptr_cast<tcFermionSpinPtr>(spin));
    // has spin info transfer
    if(sp) {
      vector<EvtDiracSpinor> prod,decay;
      // transfer spinors
      for(unsigned int ix=0;ix<2;++ix) {
  	prod .push_back(EvtGenSpinor(sp->getProductionBasisState(ix)));
  	decay.push_back(EvtGenSpinor(sp->getDecayBasisState     (ix)));
      }
      myPart->init(id, p4,prod[0],prod[1],decay[0],decay[1]);
      // set the density matrix
      myPart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
    }
    // no spin info
    else {
      myPart->init(id, p4);
      EvtSpinDensity rho;
      rho.setDiag(EvtSpinType::getSpinStates(thisSpin));
      myPart->setSpinDensityForward(rho);
    }
    evtpart=myPart;
  }
  // vector particles
  else if ( thisSpin == EvtSpinType::VECTOR ) {
    EvtVectorParticle *myPart=new EvtVectorParticle;
    // get the spin info
    tcVectorSpinPtr sp(dynamic_ptr_cast<tcVectorSpinPtr>(spin));
    // has spin info transfer
    if(sp) {
      // transfer polarization vectors
      vector<EvtVector4C> eps;
      for(unsigned int ix=0;ix<3;++ix) {
      	eps.push_back(EvtGenPolarization(sp->getDecayBasisState(ix)));
      }
      myPart->init(id, p4,eps[0],eps[1],eps[2]);
      // set spin density matrix
      myPart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
    }
    // no spin info
    else {
      myPart->init(id, p4);
      EvtSpinDensity rho;
      rho.setDiag(EvtSpinType::getSpinStates(thisSpin));
      myPart->setSpinDensityForward(rho);
    }
    evtpart=myPart;
  }    
  // spin 3/2 particles
  else if ( thisSpin == EvtSpinType::RARITASCHWINGER ) {
    EvtRaritaSchwingerParticle *myPart=new EvtRaritaSchwingerParticle;
    // get the spin info
    tcRSFermionSpinPtr sp(dynamic_ptr_cast<tcRSFermionSpinPtr>(spin));
    // has spin info transfer
    if(sp) {
      // transfer spinors
      vector<EvtRaritaSchwinger> prod,decay;
      for(unsigned int ix=0;ix<4;++ix) {
   	prod .push_back(EvtGenRSSpinor(sp->getProductionBasisState(ix)));
   	decay.push_back(EvtGenRSSpinor(sp->getDecayBasisState(ix)     ));
      }
      myPart->init(id, p4,
   		   prod[0] ,prod[1] ,prod[2] , prod[3],
   		   decay[0],decay[1],decay[2],decay[3]);
      // set spin density matrix
      myPart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
    }
    // no spin info
    else {
      myPart->init(id, p4);
      EvtSpinDensity rho;
      rho.setDiag(EvtSpinType::getSpinStates(thisSpin));
      myPart->setSpinDensityForward(rho);
    }
    evtpart=myPart;
  }
  // spin-2 particles
  else if ( thisSpin == EvtSpinType::TENSOR ) {
    EvtTensorParticle *myPart =new EvtTensorParticle;
    tcTensorSpinPtr sp(dynamic_ptr_cast<tcTensorSpinPtr>(spin));
    // has spin info transfer
    if(sp) {
      // transfer the polarization tensors
      vector<EvtTensor4C> eps;
      for(unsigned int ix=0;ix<5;++ix) {
  	eps.push_back(EvtGenTensor(sp->getDecayBasisState(ix)));
      }
      myPart->init(id, p4,eps[0],eps[1],eps[2],eps[3],eps[4]);
      // set spin density matrix
      myPart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
    }
    // no spin info
    else {
      myPart->init(id, p4);
      EvtSpinDensity rho;
      rho.setDiag(EvtSpinType::getSpinStates(thisSpin));
      myPart->setSpinDensityForward(rho);
    }
    evtpart=myPart;
  }
  // shouldn't be doing this but here for safety
  else if ( thisSpin == EvtSpinType::STRING ) {
    EvtStringParticle *myPart;
    myPart=new EvtStringParticle;
    myPart->init(id, p4);
    EvtSpinDensity rho;
    rho.setDiag(EvtSpinType::getSpinStates(thisSpin));
    myPart->setSpinDensityForward(rho);
    evtpart=myPart;
  }
  // high spin particles
  else if ( thisSpin == EvtSpinType::SPIN3 || thisSpin == EvtSpinType::SPIN5HALF ||
   	    thisSpin == EvtSpinType::SPIN4 || thisSpin == EvtSpinType::SPIN7HALF) {
    EvtHighSpinParticle *myPart;
    myPart=new EvtHighSpinParticle;
    myPart->init(id, p4);
    EvtSpinDensity rho;
    rho.setDiag(EvtSpinType::getSpinStates(thisSpin));
    myPart->setSpinDensityForward(rho);
    evtpart=myPart;
  }
  else {
    throw Exception() << "Can't convert particle to EvtGen " 
  		      << part << Exception::eventerror;
  }
  // boost particle back and return the particle
  decay->boost(inmom.boostVector());
  return evtpart;
}

// convert from ThePEG to EvtGen
EvtId EvtGenInterface::EvtGenID(int id, bool exception) const {
  EvtId output;
  int absid=abs(id),ispin(absid%10),isgn(id/absid);
  // handle the easy cases
  if(
     //quarks(+excited)
     absid<=8|| 
     // leptons(+exicted)
     (absid>=11&&absid<=18)|| 
     // SM(+extensions) gauge bosons and higgs
     (absid>=21&&absid<=25)||(absid>=32&&absid<=37)||
     // 1 1S0, 1 3S1 and 1 3P2 1 3D3 mesons are the same 
     (absid>100&&absid<600&&ispin%2==1&&ispin<=9)|| 
     // 2 1S0, 2 3S1  are the same 
     (absid>100100&&absid<100600&&(ispin==1||ispin==3))||
     // 1 3d1 goes to 3s in evtgen
     (absid>30100&&absid<30600&&ispin==3) ||
     // mixed kaons and diffractive states
     (absid>=100&&absid<=3000&&ispin==0)) {
    output = EvtPDL::evtIdFromStdHep(id);
  }
  // lowest baryon multiplets and diquarks are almost the same
  else if(absid>1000&&absid<6000&&ispin>=1&&ispin<=4) {
    output = EvtPDL::evtIdFromStdHep(id);
  }
  // 1 1P1 mesons are the same apart from D_s1
  else if(absid>10100&&absid<10600&&(ispin==3)) {
    if(absid==10433) output = EvtPDL::evtIdFromStdHep(isgn*20433);
    else             output = EvtPDL::evtIdFromStdHep(id);
  }
  // 1 3p0 same apart from f'0
  else if (absid>10100&&absid<10600&&ispin==1) {
    if(absid==10221) output = EvtPDL::evtIdFromStdHep(isgn*30221);
    else             output = EvtPDL::evtIdFromStdHep(id);
  }
  // 1 3P1 mesons are the same apart from D_s1
  else if(absid>20100&&absid<20600&&(ispin==3)) {
    if(absid==20433) output = EvtPDL::evtIdFromStdHep(isgn*10433);
    else             output = EvtPDL::evtIdFromStdHep(id);
  }
  // excited baryons
  else if(((absid>10000&&absid<16000)||(absid>20000&&absid<26000)||
	   (absid>30000&&absid<36000))&&(ispin==2||ispin==4)) {
    output=EvtPDL::evtIdFromStdHep(id);
  }
  // special for any bottomonium states
  else if((absid%1000)/10==55) {
    output=EvtPDL::evtIdFromStdHep(id);
  }
  // special meson codes
  else if (absid>9000000) {
    // f_0/a_0(980) goes into f_0/a_0(980)
    if(absid==9000111||absid==9000211||absid==9010221) {
      output=EvtPDL::evtIdFromStdHep(id);
    }
    // sigma
    else if(absid==9000221) {
      output=EvtPDL::evtIdFromStdHep(id);
    }
    // psi(4040)
    else if(absid==9000443) {
      output=EvtPDL::evtIdFromStdHep(50443);
    }
    // f_0(1500)
    else if(absid==9030221) {
      output=EvtPDL::evtIdFromStdHep(9020221);
    }
  }
  // check its O.K.
  if(output.getAlias()==-1&&exception) {
    throw Exception() << "Can't find the EvtGen Id for particle " 
		      << getParticleData(id)->PDGName()
		      << " with PDG code = " 
		      << id << " in EvtGen::EvtGenID" 
		      << Exception::eventerror;
  }
  return output;
}

ParticleVector EvtGenInterface::decayProducts(EvtParticle *part, bool boost) const {
  ParticleVector output,temp;
  for(unsigned int ix=0,N=part->getNDaug();ix<N;++ix) {
    // get the daughter particle
    EvtParticle * daug = part->getDaug(ix);
    // may just have been used for mass generation, so check has valid momentum
    if(!daug->hasValidP4()) continue;
    // get the Herwig ParticleData object
    int id=ThePEGID(daug->getId());
    tcPDPtr pd=getParticleData(id);
    // if easy to convert do it
    if(pd) output.push_back(ThePEGParticle(daug,pd));
    // special for EvtGen particles like string we don't need to include
    else if(id==90) {
      // check if needs to be decayed by EvtGen
      if(EvtPDL::getStdHep(daug->getId())==92||daug->isDecayed()) {
   	// add the particles
	ParticleVector temp(decayProducts(daug,EvtPDL::getStdHep(daug->getId())!=92));
   	if(EvtPDL::getStdHep(daug->getId())==92) {
   	  Lorentz5Momentum pstring(ThePEGMomentum(daug->getP4(),daug->mass()));
   	  Boost bv(-pstring.x()/pstring.e(),-pstring.y()/pstring.e(),-pstring.z()/pstring.e());
   	  for(unsigned int ix=0;ix<temp.size();++ix) {
   	    temp[ix]->deepBoost(bv);
   	  }
   	}
   	for(unsigned int iy=0;iy<temp.size();++iy) output.push_back(temp[iy]);
      }
      else if(!daug->isDecayed()) {
   	EvtDecayBase *decayer = EvtDecayTable::getInstance()->getDecayFunc(daug);
   	// must be a decayer
   	if(!decayer) throw Exception() << "Could find EvtGen decayer in "
  				       << "EvtGenInterface::decayProducts()" 
  				       << Exception::runerror;
  	// If there are already daughters, then this step is already done!
  	// figure out the masses
  	if ( daug->getNDaug() == 0 ) daug->generateMassTree();
  	// perform decay
  	decayer->makeDecay(daug,false);
  	// add the particles
  	ParticleVector temp(decayProducts(daug));
  	for(unsigned int iy=0;iy<temp.size();++iy) output.push_back(temp[iy]);
      }
    }
    else {
      throw Exception() << "Tried to convert an unknown particle, PDG code = " << id 
  			<< " from EvtGen in EvtGen::decayproducts which is "
  			<< " unknown to Herwig. Killing event." 
  			<< Exception::eventerror;
    }
  }
  // boost to lab
  if(output.size()>0) {
    if(boost) {
      Boost bv=ThePEGMomentum(part->getP4(),part->mass()).boostVector();
      for(unsigned int ix=0; ix<output.size(); ++ix) output[ix]->deepBoost(bv); 
    }
  }
  return output;
}

void EvtGenInterface::ThePEGSpin(PPtr peg,EvtParticle *evt) const {
  // if no corresponding ThePEG spin type return
  if(evt->getSpinType()==EvtSpinType::STRING||
     evt->getSpinType()==EvtSpinType::SPIN3||evt->getSpinType()==EvtSpinType::SPIN4||
     evt->getSpinType()==EvtSpinType::SPIN5HALF||
     evt->getSpinType()==EvtSpinType::SPIN7HALF) return;
  // now deal with the other cases
  // scalars
  unsigned int ix;
  SpinPtr spin;
  if(evt->getSpinType()==EvtSpinType::SCALAR) {
    ScalarSpinPtr sca(new_ptr(ScalarSpinInfo(peg->momentum(),true)));
    spin=sca;
  }
  // massless spin 1/2
  else if(evt->getSpinType()==EvtSpinType::NEUTRINO) {
    FermionSpinPtr fer(new_ptr(FermionSpinInfo(peg->momentum(),true)));
    spin=fer;
    ix=peg->id()<0;
    fer->setBasisState(ix,ThePEGSpinor(evt->spParentNeutrino()));
    ix=peg->id()>0;
    fer->setBasisState(ix,LorentzSpinor<SqrtEnergy>());
  }
  // massive spin-1/2
  else if(evt->getSpinType()==EvtSpinType::DIRAC) {
    FermionSpinPtr fer(new_ptr(FermionSpinInfo(peg->momentum(),true)));
    spin=fer;
    for(ix=0;ix<2;++ix) {
      fer->setBasisState(ix,ThePEGSpinor(evt->spParent(ix)));
    }
  }
  // massless vectors
  else if(evt->getSpinType()==EvtSpinType::PHOTON) {
    VectorSpinPtr vec(new_ptr(VectorSpinInfo(peg->momentum(),true)));
    spin=vec;
    for(ix=0;ix<1;++ix) {
      vec->setBasisState(2*ix,ThePEGPolarization(evt->epsParentPhoton(ix)));}
    vec->setBasisState(1,LorentzVector<Complex>());
  }
  // massive vectors
  else if(evt->getSpinType()==EvtSpinType::VECTOR) {
    VectorSpinPtr vec(new_ptr(VectorSpinInfo(peg->momentum(),true)));
    spin=vec;
    for(ix=0;ix<3;++ix) {
      vec->setBasisState(ix,ThePEGPolarization(evt->epsParent(ix)));
    }
  }
  // massive spin 3/2
  else if(evt->getSpinType()==EvtSpinType::RARITASCHWINGER) {
    RSFermionSpinPtr rs(new_ptr(RSFermionSpinInfo(peg->momentum(),true)));
    spin=rs;
    for(ix=0;ix<4;++ix) {
      rs->setBasisState(ix,ThePEGRSSpinor(evt->spRSParent(ix)));}
  }
  // tensors
  else if(evt->getSpinType()==EvtSpinType::TENSOR) {
    TensorSpinPtr ten(new_ptr(TensorSpinInfo(peg->momentum(),true)));
    spin=ten;
    for(ix=0;ix<5;++ix) {
      ten->setBasisState(ix,ThePEGTensor(evt->epsTensorParent(ix)));
    }
  }
  else {
    throw Exception() << "Unknown EvtGen spin type in EvtGen::ThePEGSpin() " 
		      << Exception::runerror;
  }
  // set the spininfo
  peg->spinInfo(spin);
  // final set up if the particle has decayed
  if(evt->getSpinDensityForward().getDim()!=0) {
    spin->rhoMatrix()=ThePEGSpinDensity(evt->getSpinDensityForward(),peg->id());
  }
  if(evt->isDecayed()) {
    spin->decayed(true);
    spin->develop();
    if(evt->getSpinDensityBackward().getDim()!=0) {
      spin->DMatrix()=ThePEGSpinDensity(evt->getSpinDensityBackward(),peg->id());
    }
  }
}

// convert from EvtGen to ThePEG
int EvtGenInterface::ThePEGID(EvtId eid,bool exception) const {
  int output(0);
  int id(EvtPDL::getStdHep(eid)),absid(abs(id)),ispin(absid%10),isgn(id/absid);
  // handle the easy cases
  // special particles like string which we delete and include decay products
  if(absid==92||absid==41||absid==42||absid==43||absid==44||
     absid==30343||absid==30353||
     absid==30363||absid==30373||absid==30383) {
    output=90;
  }
  //quarks(+excited)
  else if(absid<=8|| 
     // leptons(+excited)
     (absid>=11&&absid<=18)|| 
     // SM gauge bosons and higgs
     (absid>=21&&absid<=25)||(absid>=32&&absid<=37)|| 
     // 1 1S0, 1 3S1 and 1 3P2 mesons are the same
     (absid>100&&absid<600&&(ispin%2==1&&ispin<9))|| 
     // 2 1S0 and 2 3S1 mesons are the same
     (absid>100100&&absid<100600&&(ispin==1||ispin==3))||
     // 1 3p0 are the same 
     (absid>10100&&absid<10600&&ispin==1) ||
     // 3d1 are the same
     (absid>30100&&absid<30600&&ispin==3) ||
     // lowest baryon multiplets and diquarks
     (absid>1000&&absid<6000&&ispin>=1&&ispin<=4)|| 
     // mixed kaons and diffractive states
     (absid>=100&&absid<=3000&&ispin==0)) {
    output=id;
  }
  // 1 1P1 mesons are the same apart from D_s1
  else if (absid>10100&&absid<10600&&(ispin==3)) {
    if(absid==10433) output=isgn*20433;
    else             output=id;
  }
  // 1 3P1 mesons are the same apart from D_s1
  else if(absid>20100&&absid<20600&&(ispin==3)) {
    if(absid==20433) output=isgn*10433;
    else             output=id;
  }
  // 1 3P0 special for f'0
  else if(absid==30221)
    output = 10221;
  // special for the virtual W (we don't distinguish so return W)
  else if(absid==89) return id/absid*24;
  // excited baryons
  else if(((absid>10000&&absid<16000)||(absid>20000&&absid<26000)||
	   (absid>30000&&absid<36000))&&(ispin==2||ispin==4)) {
    output=id;
  }
  // bottomium
  else if((absid%1000)/10==55) {
    output = id;
  }
  // special mesons
  else if(absid>9000000) {
    // f_0/a_0(980) goes into f_0/a_0(980)
    if(absid==9000111||absid==9000211||absid==9010221) {
      output=id;
    }
    // sigma
    else if(absid==9000221) {
      output=id;
    }
    // psi(4040)
    else if(absid==9000443) {
      output=id;
    }
    // f_0(1500)
    else if(absid==9020221) {
      output=9030221;
    }
  }
  // check its O.K.
  if(output==0&&exception) {
    throw Exception() << "Can't find the ThePEG id for particle " 
		      << EvtPDL::name(eid)
		      << " with PDG code = " 
		      << id << " in EvtGen::ThePEGID" 
		      << Exception::eventerror;
  }
  return output;
}

RhoDMatrix EvtGenInterface::ThePEGSpinDensity(const EvtSpinDensity & rho, int id) const {
  RhoDMatrix output;
  unsigned int ix,iy;
  // special for neutrinos
  if(abs(id)%2==0&&abs(id)>=12&&abs(id)<=16) {
    output = RhoDMatrix(PDT::Spin1Half, false);
    if(id>0) output(0,0)=ThePEGComplex(rho.get(0,0));
    else     output(1,1)=ThePEGComplex(rho.get(0,0));
  }
  // special for photons
  else if(id==ParticleID::gamma) {
    output = RhoDMatrix(PDT::Spin1, false);
    for(ix=0;ix<2;++ix) {
      for(iy=0;iy<2;++iy) output(2*ix,2*iy)=ThePEGComplex(rho.get(ix,iy));
    }
  }
  // normal behaviour
  else {
    unsigned int ndim(abs(rho.getDim()));
    if(ndim!=0) {
      output = RhoDMatrix(PDT::Spin(ndim));
      for(ix=0;ix<ndim;++ix) {
	for(iy=0;iy<ndim;++iy) output(ix,iy)=ThePEGComplex(rho.get(ix,iy));
      }
    }
    else if(ndim==0) {
      output = RhoDMatrix(PDT::Spin0);
    }
  }
  return output;
}

void EvtGenInterface::checkConversion() const {
  // check the translation of particles from ThePEG to EvtGen.
  ParticleMap::const_iterator pit  = generator()->particles().begin();
  ParticleMap::const_iterator pend = generator()->particles().end();
  cout << "Testing conversion of particles from ThePEG to EvtGen\n";
  EvtId etemp;
  for(;pit!=pend;++pit) {
    cout << pit->first << "     ";
    etemp=EvtGenID(pit->first,false);
    Energy mass = pit->second->mass();
    Energy width = pit->second->width();
    if(etemp.getAlias()>=0) {
      cout << pit->second->PDGName() << "\t becomes " 
			 << EvtPDL::name(etemp) << "\t " 
			 << EvtPDL::getStdHep(etemp);
      if(ThePEGID(etemp,false)-pit->first!=0) {
	cout << " and converting back to ThePEG fails";
      }
      double mass2  = EvtPDL::getMeanMass(etemp);
      double width2 = EvtPDL::getWidth(etemp);
      if(mass!=ZERO) {
	double delta = (mass-mass2*GeV)/mass;
	if(abs(delta)>1e-6)
	  cout << " Mass Difference " << mass/GeV-mass2;
      }
      if(width>ZERO) {
	double delta = (width-width2*GeV)/width;
	if(abs(delta)>1e-6)
	  cout << " Width Difference " << width/GeV-width2;
      }
      cout << "\n";
    }
    else {
      cout << pit->second->PDGName()  
			 << " has no match in EvtGen \n";
    }
  }
  // test the conversion the other way
  cout << "Testing conversion of particles from EvtGen to ThePEG\n";
  int idtemp;
  tcPDPtr ptemp;
  for(unsigned int ix=0;ix<EvtPDL::entries();++ix) {
    cout << EvtPDL::getStdHep(EvtPDL::getEntry(ix)) << "     "
		       << EvtPDL::name(EvtPDL::getEntry(ix));
    idtemp=ThePEGID(EvtPDL::getEntry(ix),false);
    ptemp=getParticleData(idtemp);
    if(ptemp) {
      cout << " becomes " << ptemp->PDGName()  << "   "
			 << ptemp->id();
      if(EvtGenID(ptemp->id(),false)!=EvtPDL::getEntry(ix)) {
	cout << " and converting back to EvtGEN fails ";
      }
      cout << "\n";
    }
    else {
      cout << " has no match in ThePEG \n";
    }
  }
}

void EvtGenInterface::outputEvtGenDecays(long parentid) const {
  // output the decay modes from EvtGen so we can read them in
  EvtId parent(EvtGenID(parentid));
  // get evtids for children
  cout << "Outputting decays for " 
		     << getParticleData(parentid)->PDGName() << "\n";
  for(int ix=0;ix<EvtDecayTable::getInstance()->getNMode(parent.getAlias());++ix) {
    EvtDecayBase* decayer=EvtDecayTable::getInstance()->getDecay(parent.getAlias(),ix);
    vector<long> pegid;
    for(int iy=0;iy<decayer->getNDaug();++iy) {
      pegid.push_back(ThePEGID(decayer->getDaug(iy),false));
    }
    for(unsigned int iy=pegid.size();iy<7;++iy) pegid.push_back(0);
    double br=decayer->getBranchingFraction();
    cout 
      << "insert into decay_modes (incomingID,BR,decayon,star,outgoingID1,"
      << "outgoingID2,outgoingID3,outgoingID4,outgoingID5,outgoingID6,outgoingID7,"
      << "description,decayer) values (" << parentid << "," << br << ",'on','*',";
    for(int iy=0;iy<7;++iy) {
      cout << pegid[iy] << ",";
    }
    cout << "'Decay of %name% with branching "
	     << "ratio taken from EvtGen.',3);\n";
  }
}

int EvtGenInterface::EvtGenChannel(const DecayMode & dm) const {
  // get evtid for parent
  EvtId parent(EvtGenID(dm.parent()->id()));
  // get evtids for children
  EvtId daugs[20];
  unsigned int ndaug(0);
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  EvtId etemp;
  for( ;pit!=pend&&ndaug<20;++pit) {
    etemp=EvtGenID((**pit).id());
    daugs[ndaug]=etemp;
    ndaug++;
  }
  if(ndaug>=20) throw Exception() << "Too many decay products "
				  << "in EvtGen::EvtGenChannel"
				  << Exception::eventerror;
  // find the channel
  int output=EvtDecayTable::getInstance()->inChannelList(parent,ndaug,daugs);
  if(output<0) {
    string dmode;
    dmode = "decay of " + dm.parent()->PDGName() + " -> ";
    pit  = dm.products().begin();
    pend = dm.products().end();
    for( ;pit!=pend&&ndaug<20;++pit) dmode += (**pit).PDGName()+" ";
    throw Exception() << "Can't find EVtGen decay mode in EvtGenChannel for " 
		      << dmode << " in EvtGen::EvtGenChannel" 
		      << Exception::runerror;
  }
  return output;
}

// translate the matrix element to our form
void EvtGenInterface::constructVertex(const Particle & parent,
				      ParticleVector products,
				      EvtDecayAmp* damp) const {
  unsigned int ix,iy;
  vector <PDT::Spin> outspin;
  for(ix=0;ix<products.size();++ix) {
    outspin.push_back(products[ix]->dataPtr()->iSpin());
  }
  vector<int> constants(products.size()+2,0);
  unsigned int nstate(1),idtemp;
  for(ix=outspin.size();ix>0;--ix) {
    idtemp=abs(products[ix-1]->id());
    if(idtemp==ParticleID::gamma) {
      nstate*=2;
      constants[ix]=nstate;
    }
    else if(idtemp%2==0&&idtemp>=12&&idtemp<=16) {
      constants[ix]=nstate;
    }
    else {
      nstate*=outspin[ix-1];
      constants[ix]=nstate;
    }
  }
  PDT::Spin inspin(parent.dataPtr()->iSpin());
  nstate*=inspin;
  constants[0]=nstate;
  constants[outspin.size()+1]=1;
  int eind[10]={0,0,0,0,0,0,0,0,0,0};
  unsigned int iloc;
  vector<unsigned int> hind(outspin.size()+1,0);
  DecayMEPtr newME = new_ptr(GeneralDecayMatrixElement(inspin,outspin));
  for(ix=0;ix<nstate;++ix) {
    iloc=0;
    hind[0]=ix/constants[1];
    if(inspin!=PDT::Spin0) {
      eind[iloc]=hind[0];
      ++iloc;
    }
    for(iy=0;iy<outspin.size();++iy) {
      if(outspin[iy]==PDT::Spin0) {
	hind[iy+1]=0;
      }
      else if(products[iy]->id()==ParticleID::gamma) {
	hind[iy+1]=2*(ix%constants[iy+1])/constants[iy+2];
	eind[iloc]=hind[iy+1]/2;
	++iloc;
      }
      else if(constants[iy+1]==1) {
	hind[iy+1]=products[iy]->id()<0;
      }
      else {
	hind[iy+1]=(ix%constants[iy+1])/constants[iy+2];
	eind[iloc]=hind[iy+1];++iloc;
      }
    }
    if(iloc==0) eind[0]=0;
    (*newME)(hind)=ThePEGComplex(damp->amplitude().getAmp(eind));
  }
  // create the decay vertex
  VertexPtr vertex(new_ptr(DecayVertex()));
  DVertexPtr Dvertex(dynamic_ptr_cast<DVertexPtr>(vertex));
  // set the incoming particle for the decay vertex
  dynamic_ptr_cast<tcSpinPtr>(parent.spinInfo())->decayVertex(vertex);
  for(ix=0;ix<products.size();++ix) {
    dynamic_ptr_cast<tcSpinPtr>(products[ix]->spinInfo())
      ->productionVertex(vertex);
  }
  // set the matrix element
  Dvertex->ME(newME);
}
