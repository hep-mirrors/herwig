// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayIntegrator2 class.
//

#include "DecayIntegrator2.h"
#include "PhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void DecayIntegrator2::persistentOutput(PersistentOStream & os) const {
  os << modes_ << nIter_ << nPoint_ << nTry_
     << photonGen_ << generateInter_ << ounit(eps_,GeV);
}

void DecayIntegrator2::persistentInput(PersistentIStream & is, int) {
  is >> modes_ >> nIter_ >> nPoint_ >> nTry_
     >> photonGen_ >> generateInter_ >> iunit(eps_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<DecayIntegrator2,HwDecayerBase>
describeHerwigDecayIntegrator2("Herwig::DecayIntegrator2", "Herwig.so");

void DecayIntegrator2::Init() {

  static ClassDocumentation<DecayIntegrator2> documentation
    ("The DecayIntegrator class is a base decayer class "
     "including a multi-channel integrator.");

  static Parameter<DecayIntegrator2,unsigned int> interfaceIteration
    ("Iteration",
     "Number of iterations for the initialization of the phase space",
     &DecayIntegrator2::nIter_, 10, 0, 100,
     false, false, true);  
  
  static Parameter<DecayIntegrator2,unsigned int> interfacePoints
    ("Points",
     "number of phase space points to generate in the initialisation.",
     &DecayIntegrator2::nPoint_, 10000, 1, 1000000000,
     false, false, true);
  
  static Parameter<DecayIntegrator2,unsigned int> interfaceNtry
    ("Ntry",
     "Number of attempts to generate the decay",
     &DecayIntegrator2::nTry_, 500, 0, 100000,
     false, false, true);

  static Reference<DecayIntegrator2,DecayRadiationGenerator> interfacePhotonGenerator
    ("PhotonGenerator",
     "Object responsible for generating photons in the decay.",
     &DecayIntegrator2::photonGen_, false, false, true, true, false);
 
  static Switch<DecayIntegrator2,bool> interfaceGenerateIntermediates
    ("GenerateIntermediates",
     "Whether or not to include intermediate particles in the output",
     &DecayIntegrator2::generateInter_, false, false, false);
  static SwitchOption interfaceGenerateIntermediatesNoIntermediates
    (interfaceGenerateIntermediates,
     "No",
     "Don't include the intermediates",
     false);
  static SwitchOption interfaceGenerateIntermediatesIncludeIntermediates
    (interfaceGenerateIntermediates,
     "Yes",
     "include the intermediates",
     true);
}

double DecayIntegrator2::oneLoopVirtualME(unsigned int ,
					  const Particle &, 
					  const ParticleVector &) {
  throw DecayIntegrator2Error()
    << "DecayIntegrator2::oneLoopVirtualME() called. This should"
    << " have been overidden in an inheriting class if it is used"
    << Exception::runerror;
}

InvEnergy2 DecayIntegrator2::realEmissionME(unsigned int,
					    const Particle &, 
					    ParticleVector &,
					    unsigned int,
					    double, double, 
					    const LorentzRotation &,
					    const LorentzRotation &) {
  throw DecayIntegrator2Error()
    << "DecayIntegrator2::realEmmisionME() called. This should"
    << " have been overidden in an inheriting class if it is used"
    << Exception::runerror;
}

ParticleVector DecayIntegrator2::decay(const Particle & parent,
				       const tPDVector & children) const {
  // return empty vector if products heavier than parent
  Energy mout(ZERO);
  for(tPDPtr pd : children)  mout += pd->massMin();
  if(mout>parent.mass()) return ParticleVector();
  // generate the decay
  bool cc;
  iMode_ = modeNumber(cc,parent.dataPtr(),children);
  return modes_[iMode_]->generateDecay(parent,this,generateInter_,cc);
}

void DecayIntegrator2::doinitrun() {
  HwDecayerBase::doinitrun();
  if ( initialize() && Debug::level > 1 ) 
    CurrentGenerator::current().log() << "Start of the initialisation for " 
				      << name() << "\n";
  for(unsigned int ix=0;ix<modes_.size();++ix) {
    if(!modes_[ix]) continue;
    modes_[ix]->initrun();
    iMode_=ix;
    modes_[ix]->initializePhaseSpace(initialize(),this);
  }
}

// output the information for the database
void DecayIntegrator2::dataBaseOutput(ofstream & output,bool header) const {
  // header for MySQL
  if(header) output << "update decayers set parameters=\"";
  output << "newdef " << name() << ":Iteration " << nIter_ << "\n";
  output << "newdef " << name() << ":Ntry " << nTry_ << "\n";
  output << "newdef " << name() << ":Points " << nPoint_ << "\n";
  //if(_photongen){;}
  output << "newdef " << name() << ":GenerateIntermediates " << generateInter_ << " \n";
  // footer for MySQL
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
  }
}

// set the code for the partial width
void DecayIntegrator2::setPartialWidth(const DecayMode & dm, int imode) {
  int ifound = findMode(dm);
  if(ifound>=0) modes_[ifound]->setPartialWidth(imode);
}

WidthCalculatorBasePtr 
DecayIntegrator2::threeBodyMEIntegrator(const DecayMode &) const {
  return WidthCalculatorBasePtr();
}

int DecayIntegrator2::findMode(const DecayMode & dm) {
  int imode(-1);
  vector<int> extid;
  bool found(false);
  int id;
  unsigned int ix(0),iy,N,iz,tmax,nmatched;
  if(modes_.size()==0) return -1;
  do {
    if(!modes_[ix]) {
      ++ix;
      continue;
    }
    tcPDPtr in = modes_[ix]->incoming().first;
    tcPDPtr cc = modes_[ix]->incoming().first->CC();
    tmax=1;if(!cc){++tmax;}
    for(iz=0;iz<tmax;++iz) {
      extid.clear();
      // check the parent
      if(dm.parent()!=in && dm.parent()!=cc) continue;
      if(dm.parent()->id()==in->id()&&iz==0) {
	for(iy=0,N=modes_[ix]->numberOfParticles();iy<N;++iy) {
	  extid.push_back(modes_[ix]->outgoing()[iy]->id());
	}
      }
      else if(dm.parent()->id()==in->id()&&iz==1) {
	for(iy=0,N=modes_[ix]->numberOfParticles();iy<N;++iy) {
	  tcPDPtr cc2=modes_[ix]->outgoing()[iy]->CC();
	  extid.push_back( cc2 ? cc2->id() : modes_[ix]->outgoing()[iy]->id());
	}
      }
      else if(cc&&dm.parent()->id()==cc->id()) {
	for(iy=0,N=modes_[ix]->numberOfParticles();iy<N;++iy) {
	  tcPDPtr cc2 = modes_[ix]->outgoing()[iy]->CC();
	  extid.push_back( cc2 ? cc2->id() : modes_[ix]->outgoing()[iy]->id());
	}
      }
      // if the parents match
      if(!extid.empty()) {
	vector<bool> matched(extid.size(),false);
	bool done;
	nmatched=0;
	ParticleMSet::const_iterator pit = dm.products().begin();
	do {
	  id=(**pit).id();
	  done=false;
	  iy=0;
	  do {
	    if(id==extid[iy]&&!matched[iy]) {
	      matched[iy]=true;
	      ++nmatched;
	      done=true;
	    }
	    ++iy;
	  }
	  while(iy<extid.size()&&!done);
	  ++pit;
	}
	while(pit!=dm.products().end());
	if(nmatched==extid.size()) {
	  imode=ix;
	  found=true;
	}
      }
    }
    ++ix;
  }
  while(!found&&ix<modes_.size());
  return imode;
}

// the matrix element to be integrated for the me
double DecayIntegrator2::threeBodyMatrixElement(const int,const Energy2,
					       const Energy2,
					       const Energy2,const Energy2,
					       const Energy, const Energy, 
					       const Energy) const {
  throw DecayIntegratorError() 
    << "Calling the virtual DecayIntegrator2::threeBodyMatrixElement"
    << "method. This must be overwritten in the classes "
    << "inheriting from DecayIntegrator where it is needed"
    << Exception::runerror;
}

// the differential three body decay rate with one integral performed
InvEnergy DecayIntegrator2::threeBodydGammads(const int, const Energy2,
					     const Energy2,
					     const Energy, const Energy, 
					     const Energy) const {
  throw DecayIntegratorError() 
    << "Calling the virtual DecayIntegrator2::threeBodydGammads()" 
    <<"method. This must be overwritten in the classes "
    << "inheriting from DecayIntegrator where it is needed"
    << Exception::runerror;
}

// generate the momenta for the decay
ParticleVector DecayIntegrator2::generate(bool inter,bool cc,
					 const unsigned int & imode,
					 const Particle & inpart) const {
  iMode_=imode;
  return modes_[imode]->generateDecay(inpart,this,inter,cc);
}

void  DecayIntegrator2::addMode(PhaseSpaceModePtr mode) const {
  modes_.push_back(mode);
  if(mode) mode->init();
}

ostream & Herwig::operator<<(ostream & os, const DecayIntegrator2 & decay) { 
  os << "The integrator has " << decay.modes_.size() << " modes"  << endl;
  for(unsigned int ix=0;ix<decay.modes_.size();++ix) {
    os << "Information on mode " << ix << endl;
    os << *(decay.modes_[ix]);
  }
  return os;
}
  
// reset the properities of all intermediates
void DecayIntegrator2::resetIntermediate(tcPDPtr part, Energy mass, Energy width) {
  if(!part) return;
  for(unsigned int ix=0,N=modes_.size();ix<N;++ix) {
    modes_[ix]->resetIntermediate(part,mass,width);
  }
}

Energy DecayIntegrator2::initializePhaseSpaceMode(unsigned int imode,bool init, bool onShell) const{
  tcPhaseSpaceModePtr cmodeptr=mode(imode);
  tPhaseSpaceModePtr modeptr = const_ptr_cast<tPhaseSpaceModePtr>(cmodeptr);
  modeptr->init();
  return modeptr->initializePhaseSpace(init,this,onShell);
}
