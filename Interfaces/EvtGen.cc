// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EvtGen class.
//

#include "EvtGen.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Helicity/Correlations/DecayVertex.h"
#include "Herwig++/Helicity/Correlations/DecayMatrixElement.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenModels/EvtModelReg.hh"
#include "EvtGenBase/EvtStatus.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenModels/EvtPHOTOS.hh"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void EvtGen::doinitrun() {
  Interfaced::doinitrun();
  // output the EvtGen initialization info to the log file
  generator()->log() << "Initializing EvtGen \n";
  // set up the random number generator for EvtGen
  generator()->log() << "Setting EvtGen random number generator to the Herwig++ one.\n";
  _evtrnd = new EvtGenRandom(const_ptr_cast<Ptr<RandomGenerator>::pointer>
 			     (&(UseRandom::current())));
  EvtRandom::setRandomEngine(_evtrnd);
  generator()->log() << "Storing known decay models \n";
  // tell EvtGen about all the different models
  EvtModelReg dummy;
  // setup for photos here (change this to Keith's code if possible)
  report(INFO,"EvtGen") << "Initializing RadCorr=PHOTOS \n";
  generator()->log() << "Setting PHOTOS as the Radiateive Correction"
 		     << " Engine for EvtGEN\n";
  static EvtPHOTOS defaultRadCorrEngine;
  EvtRadCorr::setRadCorrEngine(&defaultRadCorrEngine);
  // read the particle data tables
  generator()->log() << "PDT table file name  : " << _pdtname   << "\n";
  _pdl.readPDT(_pdtname);
  // read the decay file
  generator()->log() << "Main decay file name : " << _decayname << "\n";
  EvtDecayTable::readDecayFile(_decayname);
  // check the conversion of particles if needed
  if(_checkconv) checkConversion();
  // convert any particle decay modes which are needed
  for(unsigned int ix=0;ix<_convid.size();++ix) {
    outputEvtGenDecays(_convid[ix]);
  }
  // that's the lot
  generator()->log() << "Finished initialisation of EvtGen \n";
}

// persistent output
void EvtGen::persistentOutput(PersistentOStream & os) const {
  os << _decayname << _pdtname << _maxtry << _maxunwgt << _checkconv << _convid;
}

// persistent input
void EvtGen::persistentInput(PersistentIStream & is, int) { 
  is >> _decayname >> _pdtname >> _maxtry >> _maxunwgt >> _checkconv >> _convid;
}

ClassDescription<EvtGen> EvtGen::initEvtGen;
// Definition of the static class description member.

void EvtGen::Init() {

  static ClassDocumentation<EvtGen> documentation
    ("The EvtGen class is the main class for the use of the EvtGen "
     "decay package with Herwig++");

  static Parameter<EvtGen,string> interfaceDecay_File
    ("Decay_File",
     "The name of the file for the EvtGen decays.",
     &EvtGen::_decayname, "./DECAY.DEC",
     false, false);

  static Parameter<EvtGen,string> interfacePDT_File
    ("PDT_File",
     "The name of the file for the EvtGen particle data.",
     &EvtGen::_pdtname, "./evt.pdl",
     false, false);

  static Parameter<EvtGen,unsigned int> interfaceMaximumTries
    ("MaximumTries",
     "Maximum number of attempts to generate a decay using EvtGen.",
     &EvtGen::_maxtry, 1000, 1, 10000,
     false, false, Interface::limited);

  static Parameter<EvtGen,unsigned int> interfaceMaximumUnWeight
    ("MaximumUnWeight",
     "The maximum number of tries for EvtGen to unweight a decay.",
     &EvtGen::_maxunwgt, 10000, 100, 1000000,
     false, false, Interface::limited);

  static Switch<EvtGen,bool> interfaceCheckConversion
    ("CheckConversion",
     "Check the conversion of particles to and from EvtGen",
     &EvtGen::_checkconv, false, false, false);
  static SwitchOption interfaceCheckConversionfalse
    (interfaceCheckConversion,
     "NoCheck",
     "Don't check the conversion",
     false);
  static SwitchOption interfaceCheckConversionCheck
    (interfaceCheckConversion,
     "Check",
     "Check the conversion",
     true);

  static ParVector<EvtGen,long> interfaceOutputModes
    ("OutputModes",
     "Particles for which to output the EvtGen decay modes"
     " so they can be read into Herwig++",
     &EvtGen::_convid, -1, long(0), 0, 0,
     false, false, Interface::nolimits);

}

ParticleVector EvtGen::randomDecayAll(const Particle & parent) {
  // create the particle from EvtGen using the PDG code, set's id and momentum
  EvtParticle *part=EvtGenParticle(parent);
  unsigned int ntry(1);
  int ii;
  bool done(false);
  do 
    {
      EvtStatus::initRejectFlag();
      part->decay();
      done=(EvtStatus::getRejectFlag()==0);
      if(!done)
	{
	  for (ii=0;ii<part->getNDaug();ii++)
	    {
	      EvtParticle *temp=part->getDaug(ii);
	      temp->deleteTree();
	    }
	  part->resetFirstOrNot();
	  part->resetNDaug();
	}
    }
  while(ntry<_maxtry&&!done);
  if(!done)
    {throw Exception() << "Failed to generate a decay in EvtGen::randomDecayAll"
			<< Exception::eventerror;}
  // translate the decay products
  ParticleVector output(decayProducts(part,false));
  // sort out the spin info
  tSpinfoPtr pspin(getSpinInfo(parent));
  pspin->decayed(true);
  pspin->setDeveloped(true);
  pspin->DMatrix()=ThePEGSpinDensity(part->getSpinDensityBackward(),
				     ThePEGID(part->getId()));
  // delete the EvtGen particle
  part->deleteDaughters();
  delete part;
  // return the decay products
  return output;
}

// obtain and convert a particle's decay products from EvtGen
ParticleVector EvtGen::decayProducts(EvtParticle *part,bool spin)
{
  ParticleVector output,temp;
  tcPDPtr pd;
  EvtParticle * daug;
  Hep3Vector bv;
  unsigned int ix,iy;
  int id;
  for(ix=0;ix<abs(part->getNDaug());++ix)
    {
      daug=part->getDaug(ix);
      id=ThePEGID(daug->getId());
      pd=getParticleData(id);      
      if(pd)
	{output.push_back(ThePEGParticle(daug,pd,spin));}
      // special for EvtGen particles like string we don't need to include
      else if(id==90)
	{
	  // check if needs to be decayed by EvtGen
	  bool evtdec(daug->getNDaug()==0);
	  if(daug->getNDaug()!=0)
	    {if(abs(EvtPDL::getStdHep(daug->getDaug(0)->getId()))<=6){evtdec=true;}}
	  if(!evtdec)
	    {
	      temp=decayProducts(daug,spin);
	      if(temp.size()!=0)
		{for(iy=0;iy<temp.size();++iy){output.push_back(temp[iy]);}}
	      else
		{throw Exception() << "Found EvtGen special particle with no decay"
				   << " products in EvtGen::decayProducts()" 
				   << Exception::eventerror;}
	    }
	  // get EvtGen to decay the particle
	  else
	    {
	      EvtDecayAmp* damp=0;
	      EvtDecayIncoherent* dinc=0;
	      EvtDecayProb* dprob=0;
	      unsigned int nbeforerad(0);
	      evtDecay(daug,0,damp,dinc,dprob,nbeforerad);
	      // add the particles
	      ParticleVector temp(decayProducts(daug,spin));
	      for(iy=0;iy<temp.size();++iy){output.push_back(temp[iy]);}
	    }
	}
      else
	{throw Exception() << "Tried to convert an unknown particle, PDG code = " << id 
			   << " from EvtGen in EvtGen::decayproducts which is "
			   << " unknown to Herwig++. Killing event." 
			   << Exception::eventerror;}
    }
  if(output.size()>0)
    {
      bv=ThePEGMomentum(part->getP4(),part->mass()).boostVector();
      for(ix=0;ix<output.size();++ix){output[ix]->deepBoost(bv);}
    }
  return output;
}

// convert EvtGen spin density to ThePEG 
RhoDMatrix EvtGen::ThePEGSpinDensity(const EvtSpinDensity & rho, int id)
{
  RhoDMatrix output;
  unsigned int ix,iy;
  // special for neutrinos
  if(abs(id)%2==0&&abs(id)>=12&&abs(id)<=16)
    {
      output.iSpin(PDT::Spin1Half);output.zero();
      if(id>0){output(0,0)=ThePEGComplex(rho.Get(0,0));}
      else{output(1,1)=ThePEGComplex(rho.Get(0,0));}
    }
  // special for photons
  else if(id==ParticleID::gamma)
    {
      output.iSpin(PDT::Spin1);output.zero();
      for(ix=0;ix<2;++ix)
	{for(iy=0;iy<2;++iy)
	    {output(2*ix,2*iy)=ThePEGComplex(rho.Get(ix,iy));}}
    }
  // normal behaviour
  else
    {
      unsigned int ndim(abs(rho.GetDim()));
      if(ndim!=0)
	{
	  output.iSpin(PDT::Spin(ndim));
	  for(ix=0;ix<ndim;++ix)
	    {for(iy=0;iy<ndim;++iy)
		{output(ix,iy)=ThePEGComplex(rho.Get(ix,iy));}}
	}
      else if(ndim==0)
	{output.iSpin(PDT::Spin0);output.average();}
    }
  return output;
}

// location of mode in EvtGen list
int EvtGen::EvtGenChannel(const DecayMode & dm)
{
  // get evtid for parent
  EvtId parent(EvtGenID(dm.parent()->id()));
  // get evtids for children
  EvtId daugs[20];
  unsigned int ndaug(0);
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  EvtId etemp;
  for( ;pit!=pend&&ndaug<20;++pit)
    {
      etemp=EvtGenID((**pit).id());
      daugs[ndaug]=etemp;
      ndaug++;
    }
  if(ndaug>=20){throw Exception() << "Too many decay products in EvtGen::EvtGenChannel"
				  << Exception::eventerror;}
  // find the channel
  int output=EvtDecayTable::inChannelList(parent,ndaug,daugs);
  if(output<0)
    {
      string dmode;
      dmode = "decay of " + dm.parent()->PDGName() + " -> ";
      pit  = dm.products().begin();
      pend = dm.products().end();
      for( ;pit!=pend&&ndaug<20;++pit){dmode += (**pit).PDGName()+" ";}
      throw Exception() << "Can't find EVtGen decay mode in EvtGenChannel for " 
		       << dmode << " in EvtGen::EvtGenChannel" << Exception::runerror;
    }
  return output;
}

// evtgen performing all decays of particle but herwig picks inital mode
ParticleVector EvtGen::decayAll(const DecayMode & dm,const Particle & parent)
{
  // find the location of the decay mode in the EvtGen list
  int imode(EvtGenChannel(dm));
  EvtId parID(EvtGenID(parent.id()));
  EvtDecayBase *decayer(EvtDecayTable::getDecay(parID.getAlias(),imode));
  if(!decayer){throw Exception() << "Could find EvtGen decayer in EvtGen::decayAll()" 
				 << Exception::runerror;}
  // create the particle from EvtGen using the PDG code, set's id and momentum
  EvtParticle *part=EvtGenParticle(parent);
  part->setChannel(imode);
  unsigned int ntry(1);
  int ii;
  bool done(false);
  do 
    {
      EvtStatus::initRejectFlag();
      // generate masses
      if (part->getNDaug()==0){part->generateMassTree();}
      // perform the decay
      decayer->makeDecay(part);
      done=(EvtStatus::getRejectFlag()==0);
      if(!done)
	{
	  for (ii=0;ii<part->getNDaug();ii++)
	    {
	      EvtParticle *temp=part->getDaug(ii);
	      temp->deleteTree();
	    }
	  part->resetFirstOrNot();
	  part->resetNDaug();
	}
    }
  while(ntry<_maxtry&&!done);
  if(!done)
    {throw Exception() << "Failed to generate a decay in EvtGen::decayAll"
			<< Exception::eventerror;}
  // translate the decay products
  ParticleVector output(decayProducts(part,false));
  // sort out the spin info
  tSpinfoPtr pspin(getSpinInfo(parent));
  pspin->decayed(true);
  pspin->setDeveloped(true);
  pspin->DMatrix()=ThePEGSpinDensity(part->getSpinDensityBackward(),
				     ThePEGID(part->getId()));
  // delete the EvtGen particle
  part->deleteDaughters();
  delete part;
  // return the decay products
  return output;
}

ParticleVector EvtGen::randomDecay(const Particle & parent) {
  // create the particle from EvtGen using the PDG code, set it's id and momentum
  EvtParticle *part=EvtGenParticle(parent);
  EvtDecayAmp *damp=0;
  EvtDecayIncoherent* dinc=0;
  EvtDecayProb* dprob=0;
  unsigned int nbeforerad(0);
  evtDecay(part,0,damp,dinc,dprob,nbeforerad);
  // translate the decay products
  ParticleVector output(decayProducts(part,damp));
  // sort out the spin info
  tSpinfoPtr pspin(getSpinInfo(parent));
  // if using amplitudes
  if(damp)
    {
      ParticleVector products;
      for(unsigned int ix=0;ix<nbeforerad;++ix)
	{products.push_back(output[ix]);}
      constructVertex(parent,products,damp);
    }
  // otherwise
  else
    {
      pspin->setDeveloped(true);
      RhoDMatrix rhotemp(pspin->iSpin());rhotemp.average();
      pspin->DMatrix()=rhotemp;
    }
  pspin->decayed(true);
  // delete the EvtGen particle
  part->deleteDaughters();
  delete part;
  // return the decay products
  return output;
}


void EvtGen::ThePEGSpin(PPtr peg,EvtParticle *evt) {
  // if no corresponding ThePEG spin type return
  if(evt->getSpinType()==EvtSpinType::STRING||
     evt->getSpinType()==EvtSpinType::SPIN3||evt->getSpinType()==EvtSpinType::SPIN4||
     evt->getSpinType()==EvtSpinType::SPIN5HALF||
     evt->getSpinType()==EvtSpinType::SPIN7HALF) return;
  // now deal with the other cases
  // scalars
  unsigned int ix;
  SpinfoPtr spin;
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
    fer->setBasisState(ix,LorentzSpinor());
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
    vec->setBasisState(1,LorentzVector());
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
    cerr << "testing rs access function not implemented\n";
    exit(0);
//     RSFermionSpinPtr rs(new_ptr(RSFermionSpinInfo(peg->momentum(),true)));
//     spin=rs;
//     for(ix=0;ix<4;++ix) {
//       rs->setBasisState(ix,ThePEGRSSpinor(evt->spRSParent(ix)));}
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
  if(evt->getSpinDensityForward().GetDim()!=0) {
    spin->rhoMatrix()=ThePEGSpinDensity(evt->getSpinDensityForward(),peg->id());
  }
  if(evt->getNDaug()>0) {
    spin->decayed(true);
    spin->setDeveloped(true);
    if(evt->getSpinDensityBackward().GetDim()!=0) {
      spin->DMatrix()=ThePEGSpinDensity(evt->getSpinDensityBackward(),peg->id());
    }
  }
}

// translate the matrix element to our form
void EvtGen::constructVertex(const Particle & parent,ParticleVector products,
			     EvtDecayAmp* damp) {
  unsigned int ix,iy;
  vector <PDT::Spin> outspin;
  for(ix=0;ix<products.size();++ix){outspin.push_back(products[ix]->dataPtr()->iSpin());}
  vector<int> constants(products.size()+2,0);
  unsigned int nstate(1),idtemp;
  for(ix=outspin.size();ix>0;--ix) {
    idtemp=abs(products[ix-1]->id());
    if(idtemp==ParticleID::gamma){nstate*=2;constants[ix]=nstate;}
    else if(idtemp%2==0&&idtemp>=12&&idtemp<=16){constants[ix]=nstate;}
    else{nstate*=outspin[ix-1];constants[ix]=nstate;}
  }
  PDT::Spin inspin(parent.dataPtr()->iSpin());
  nstate*=inspin;constants[0]=nstate;
  constants[outspin.size()+1]=1;
  int eind[10];
  unsigned int iloc;
  vector<unsigned int> hind(outspin.size()+1,0);
  Helicity::DecayMatrixElement newME(inspin,outspin);
  for(ix=0;ix<nstate;++ix) {
    iloc=0;
    hind[0]=ix/constants[1];
    if(inspin!=PDT::Spin0) {
      eind[iloc]=hind[0];
      ++iloc;
    }
    for(iy=0;iy<outspin.size();++iy) {
      if(outspin[iy]==PDT::Spin0)
	{hind[iy+1]=0;}
      else if(products[iy]->id()==ParticleID::gamma)
	{
	  hind[iy+1]=2*(ix%constants[iy+1])/constants[iy+2];
	  eind[iloc]=hind[iy+1]/2;++iloc;
	}
      else if(constants[iy+1]==1)
	{hind[iy+1]=products[iy]->id()<0;}
      else
	{
	  hind[iy+1]=(ix%constants[iy+1])/constants[iy+2];
	  eind[iloc]=hind[iy+1];++iloc;
	}
    }
    cerr << "testing can't get amplitude\n";
    //newME(hind)=ThePEGComplex(damp->amplitude().getAmp(eind));
  }
  // create the decay vertex
  VertexPtr vertex(new_ptr(Helicity::DecayVertex()));
  Helicity::DVertexPtr Dvertex(dynamic_ptr_cast<Helicity::DVertexPtr>(vertex));
  // set the incoming particle for the decay vertex
  dynamic_ptr_cast<tcSpinfoPtr>(parent.spinInfo())->setDecayVertex(vertex);
  for(ix=0;ix<products.size();++ix)
    {dynamic_ptr_cast<tcSpinfoPtr>
	(products[ix]->spinInfo())->setProductionVertex(vertex);}
  // set the matrix element
  Dvertex->ME().reset(newME);
}

// convert from ThePEG to EvtGen
EvtId EvtGen::EvtGenID(int id, bool exception) {
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
     // 1 1P1 mesons are the same
     (absid>10100&&absid<10600&&(ispin==3))|| 
     // 1 3P1 mesons are the same
     (absid>20100&&absid<20600&&(ispin==3))|| 
     // mixed kaons and diffractive states
     (absid>=100&&absid<=3000&&ispin==0)) {
    output = EvtPDL::evtIdFromStdHep(id);
  }
  // lowest baryon multiplets and diquarks are almost the same
  else if(absid>1000&&absid<6000&&ispin>=1&&ispin<=4) {
    // special for lambda_c(2625)
    if(absid==4124) output = EvtPDL::evtIdFromStdHep(isgn*(absid+10000));
    else            output = EvtPDL::evtIdFromStdHep(id);
  } 
  // 2 1S0 multiplet is all different
  else if((absid>100100&&absid<100600)&&ispin==1) {
    output=EvtPDL::evtIdFromStdHep(isgn*(absid%1000+20000));
  }
  // 2 3S1 multiplet is all different
  else if((absid>100100&&absid<100600)&&ispin==3) {
    // special for kaons change 2s and 1d
    if((absid%1000)/100==3) output=EvtPDL::evtIdFromStdHep(id);
    else                    output=EvtPDL::evtIdFromStdHep(isgn*(absid%1000+30000));
  }
  // 1 3p0 most same but some changes
  else if (absid>10100&&absid<10600&&ispin==1) {
    // heavy and strange
    if(absid%1000>300&&absid!=10331) output=EvtPDL::evtIdFromStdHep(id);
    // light
    else if(absid==10221)            output=EvtPDL::evtIdFromStdHep(10331);
  }
  // 1 3d1 goes to 3s in evtgen
  else if(absid>30100&&absid<30600&&ispin==3) {
    unsigned iq((absid%1000)/100);
    if(iq==3)       output=EvtPDL::evtIdFromStdHep(id);
    else if (iq==5) output=EvtPDL::evtIdFromStdHep(isgn*(absid%1000+120000));
    else            output=EvtPDL::evtIdFromStdHep(isgn*(absid%1000+40000));
  }
  // excited baryons
  else if(((absid>10000&&absid<16000)||(absid>20000&&absid<26000)||
	   (absid>30000&&absid<36000))&&(ispin==2||ispin==4)) {
    output=EvtPDL::evtIdFromStdHep(id);
  }
  // special for any bottomonium states not handled yet
  else if((absid%1000)/10==55) {
    if(absid==200553)      output=EvtPDL::evtIdFromStdHep( 60553);
    else if(absid==300553) output=EvtPDL::evtIdFromStdHep( 70553);
    else if(absid==110551) output=EvtPDL::evtIdFromStdHep( 30551);
    else if(absid==120553) output=EvtPDL::evtIdFromStdHep( 50553);
    else if(absid==100555) output=EvtPDL::evtIdFromStdHep( 10555);
    else if(absid==20555)  output=EvtPDL::evtIdFromStdHep( 30555);
    else if(absid==100557) output=EvtPDL::evtIdFromStdHep( 10557);
    else if(absid==120555) output=EvtPDL::evtIdFromStdHep( 50555);
    else if(absid==130553) output=EvtPDL::evtIdFromStdHep(130553);
    else if(absid==200555) output=EvtPDL::evtIdFromStdHep( 20555);
    else if(absid==210551) output=EvtPDL::evtIdFromStdHep( 50551);
    else if(absid==220553) output=EvtPDL::evtIdFromStdHep(110553);
  }
  // special meson codes
  else if (absid>9000000) {
    // a_0(980) goes into 3p0
    if(absid==9000111||absid==9000211) {
      output=EvtPDL::evtIdFromStdHep(isgn*(absid-9000000+10000));
    }
    // f_0(980) goes into 3p0 as well
    else if(absid==9010221) {
      output=EvtPDL::evtIdFromStdHep(isgn*(absid-9000000));
    }
    // psi(4040)
    else if(absid==9000443) {
      output=EvtPDL::evtIdFromStdHep(50443);
    }
    // f_0(1500)
    else if(absid==9030221) {
      output=EvtPDL::evtIdFromStdHep(50221);
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

// convert from EvtGen to ThePEG
int EvtGen::ThePEGID(EvtId eid,bool exception) {
  int output(0);
  int id(EvtPDL::getStdHep(eid)),absid(abs(id)),ispin(absid%10),isgn(id/absid);
  // handle the easy cases
  //quarks(+excited)
  if(absid<=8|| 
     // leptons(+excited)
     (absid>=11&&absid<=18)|| 
     // SM gauge bosons and higgs
     (absid>=21&&absid<=25)||(absid>=32&&absid<=37)|| 
     // 1 1S0, 1 3S1 and 1 3P2 mesons are the same
     (absid>100&&absid<600&&(ispin%2==1&&ispin<9))|| 
     // 1 1P1 mesons are the same
     (absid>10100&&absid<10600&&(ispin==3))|| 
     // 1 3P1 mesons are the same
     (absid>20100&&absid<20600&&(ispin==3))|| 
     // lowest baryon multiplets and diquarks
     (absid>1000&&absid<6000&&ispin>=1&&ispin<=4)|| 
     // mixed kaons and diffractive states
     (absid>=100&&absid<=3000&&ispin==0)) {
    output=id;
  }
  // special particles like string which we delete and include decay products
  else if(absid==92||absid==41||absid==42||absid==30343||absid==30353||
	  absid==30363||absid==30373||absid==30383) {
    output=90;
  }
  // 2 1S0 multiplet is all different
  else if(absid>20100&&absid<20600&&ispin==1) {
    output=isgn*(absid%1000+100000);
  }
  // 2s kaons are the same
  else if(absid>100300&&absid<100400&&ispin==3) {
    output=id;
  }
  // 2 3S1 multiplet is all different
  else if(absid>30100&&absid<30600&&ispin==1) {
    // special for kaons change 2s and 1d
    unsigned iq((absid%1000)/100);
    if(iq==3)      output=id;
    else if(iq==5) output=isgn*(absid%1000+110000);
    else           output=isgn*(absid%1000+100000);
  }
  // 1 3p0 most same some different
  else if (absid>10100&&absid<10600&&ispin==1) {
    // heavy and strange
    if(absid%1000>300&&absid!=10331) {
      output=id;
    }
    // light
    else {
      if(absid==10111||absid==10211) output=isgn*(absid-10000+9000000);
      else if(absid==10221)          output=id+9000000;
      else if(absid==10331)          output=10221;
    }
  }
  // special for the virtual W (we don't distinguish so return W)
  else if(absid==89) return id/absid*24;
  // excited baryons
  else if(((absid>10000&&absid<16000)||(absid>20000&&absid<26000)||
	   (absid>30000&&absid<36000))&&(ispin==2||ispin==4)) {
    // most cases
    if(absid!=14124) output=id;
    // lambda_c(2625)
    else             output=isgn*(absid-10000);
  }
  // 2s is a problem
  else if(absid>30100&&absid<30600&&ispin==3) {
    unsigned iq((absid%1000)/100);
    if(iq!=3) output=isgn*(absid%1000+100000);
    else      output=id;
  }
  // 3s is a problem
  else if(absid>40100&&absid<40600&&ispin==3) {
    unsigned iq((absid%1000)/100);
    if(iq<3||iq==4)              output=isgn*(absid%1000+30000);
    else if(iq==5&&absid!=40553) output=isgn*(absid%1000+120000);
  }
  else if((absid%1000)/10==55) {
    if(absid==60553)       output=200553;
    else if(absid== 70553) output=300553;
    else if(absid== 50553) output=120553;
    else if(absid== 10555) output=100555;
    else if(absid==120553) output= 30553;
    else if(absid== 30555) output= 20555;
    else if(absid== 10557) output=100557;
    else if(absid== 50555) output=120555;
    else if(absid==130553) output=130553;
    else if(absid== 20555) output=200555;
    else if(absid== 50551) output=210551;
    else if(absid==110553) output=220553;
  }
  // things that came from special codes
  else {
    if(absid==50443)      output=9000443;
    else if(absid==50221) output=9030221;
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

// translate a particle to EvtGen
EvtParticle * EvtGen::EvtGenParticle(const Particle & part) {
  // create the new particle with the correct momentum
  Lorentz5Momentum inmom(part.momentum());
  EvtParticle *evtpart;
  evtpart=EvtParticleFactory::particleFactory(EvtGenID(part.id()),EvtGenMomentum(inmom));
  // boost to the rest frame to get the basis states right
  PPtr decay(const_ptr_cast<PPtr>(new_ptr(part)));
  decay->transform(LorentzRotation(-inmom.boostVector()));
  // transfer the spin information to EvtGen
  // N.B. don't handle special massless (i.e. photon and neutrino cases as shouldn't
  // be needed)
  tcSpinfoPtr spin(dynamic_ptr_cast<tcSpinfoPtr>(part.spinInfo()));
  if(spin) {
    spin->decay();
    RhoDMatrix rho;
    // scalar particles
    if(part.dataPtr()->iSpin()==PDT::Spin0) {
      tcScalarSpinPtr sp(dynamic_ptr_cast<tcScalarSpinPtr>(spin));
      if(sp) evtpart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
    }
    // spin 1/2 particles
    else if(part.dataPtr()->iSpin()==PDT::Spin1Half) {
      tcFermionSpinPtr sp(dynamic_ptr_cast<tcFermionSpinPtr>(spin));
      if(sp) {
	unsigned int ix,id(abs(part.id()));
	// special for neutrinos as only one state in EvtGen
	if(id%2==0&&id>=12&&id<=16) {
	  throw Exception() << "Tried to convert a neutrino to EvtGen in "
			    << "EvtGen::EvtParticle. This should not be needed as"
			    << "EvtGen does not decay neutrinos" 
			    << Exception::eventerror;
	}
	else {
	  evtpart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
	  cerr << "can't set fermion spin\n";
	  for(ix=0;ix<2;++ix) {
	    //evtpart->
	    //  setspParent(ix,EvtGenSpinor(sp->getProductionBasisState(ix)));
	    //evtpart->setsp(ix,EvtGenSpinor(sp->getDecayBasisState(ix)));
	  }
	}
      }
    }
    // vector particles
    else if(part.dataPtr()->iSpin()==PDT::Spin1) {
      tcVectorSpinPtr sp(dynamic_ptr_cast<tcVectorSpinPtr>(spin));
      if(sp) {
	if(part.id()==ParticleID::gamma) {
	  throw Exception() << "Tried to convert a photon to EvtGen in "
			    << "EvtGen::EvtParticle. This should not be needed as"
			    << "EvtGen does not decay photons" 
			    << Exception::eventerror;
	}
	else {
	  evtpart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
	  cerr << "can't set spin vector\n";
	  //for(unsigned int ix=0;ix<3;++ix) {
	  //  evtpart->seteps(ix,EvtGenPolarization(sp->getDecayBasisState(ix)));
	  //}
	}
      }
    }
    // spin 3/2 particles
    else if(part.dataPtr()->iSpin()==PDT::Spin3Half) {
      tcRSFermionSpinPtr sp(dynamic_ptr_cast<tcRSFermionSpinPtr>(spin));
      if(sp) {
	evtpart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
	cerr << "can't set spin 3/2\n";
// 	for(unsigned ix=0;ix<4;++ix) {
// 	  evtpart->
// 	    setspRSParent(ix,EvtGenRSSpinor(sp->getProductionBasisState(ix)));
// 	  evtpart->setspRS(ix,EvtGenRSSpinor(sp->getDecayBasisState(ix)));
// 	}
      }
    }
    // spin-2 particles
    else if(part.dataPtr()->iSpin()==PDT::Spin2) {
      tcTensorSpinPtr sp(dynamic_ptr_cast<tcTensorSpinPtr>(spin));
      if(sp) {
	evtpart->setSpinDensityForward(EvtGenSpinDensity(sp->rhoMatrix()));
	cerr << "can't set spin 2\n";
// 	for(unsigned int ix=0;ix<5;++ix) {
// 	  evtpart->setepsTensor(ix,EvtGenTensor(sp->getDecayBasisState(ix)));
// 	}
      }
    }
  }
  // boost particle back and return the particle
  decay->boost(inmom.boostVector());
  return evtpart;
}

void EvtGen::evtDecay(EvtParticle * part,EvtDecayBase* decayerin,
		      EvtDecayAmp* damp,
		      EvtDecayIncoherent* dinc, EvtDecayProb* dprob,
		      unsigned int & nbeforerad) {
  unsigned int ntry(1);
  int ii;
  bool done(false);
  nbeforerad=0;  
  EvtDecayBase *decayer;
  EvtParticle * unmix;
  // evtgen particle ids for mixing
  static EvtId BS0=EvtPDL::getId("B_s0");
  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BD0=EvtPDL::getId("B0");
  static EvtId BDB=EvtPDL::getId("anti-B0"); 
  EvtId pid(part->getId());
  do {
    bool selected(true);
    cout << "testing the particleA " << EvtPDL::name(part->getId()) << endl;
    do {
      EvtStatus::initRejectFlag();
      // get the decayer
      if(decayerin) decayer=decayerin;
      else          decayer= EvtDecayTable::GetDecayFunc(part);
      if (part->getNDaug()==0) part->generateMassTree();
      // special if mixing
      if (part->getNDaug()==1&&(pid==BS0||pid==BSB||pid==BD0||pid==BDB)) {
	// if we have a decayer don't allow mixing
	if(decayerin) {
	  for (ii=0;ii<part->getNDaug();ii++) {
	    EvtParticle *temp=part->getDaug(ii);
	    temp->deleteTree();
	  }
	  part->resetFirstOrNot();
	  part->resetNDaug();
	  selected=false;
	}
	else {
	  unmix=part;
	  part=unmix->getDaug(0);
	  selected=false;
	}
      }
      else{selected=true;}
    }
    while(!selected);
    cout << "testing the particleB " << EvtPDL::name(part->getId()) << endl;
    // work out which type it is and perform decay
    damp  = dynamic_cast<EvtDecayAmp*>(decayer);
    dinc  = dynamic_cast<EvtDecayIncoherent*>(decayer);
    dprob = dynamic_cast<EvtDecayProb*>(decayer);
    // EvtDecayAmp is base class can do the correlations
    if(damp) {
      cerr << "testing can't get amplitude\n";
      //damp->amplitude().init(part->getId(),damp->getNDaug(),damp->getDaugs());
      double prob,prob_max;
      EvtSpinDensity rho;
      bool fail;
      int ntimes(_maxunwgt);
      do {
	cerr << "testing get set daugsDecayedByParentModel\n";
	//damp->setdaugsDecayedByParentModel(false);
	damp->setWeight(1.);
	nbeforerad=damp->getNDaug();
	damp->decay(part);
	cerr << "testing can't get amplitude\n";
	//rho=damp->amplitude().getSpinDensity();
	prob=part->getSpinDensityForward().NormalizedProb(rho);
	if(prob<0.)
	  {throw Exception() << "Negative probablity in EvtGen::randomDecay() " 
			     << Exception::runerror;}
	cerr << "testing can't get weight\n";
	//prob/=damp->getWeight();
	prob_max = damp->getProbMax(prob);
	cerr << "testing can't set decay prob\n";
	//part->setDecayProb(prob/prob_max);
	--ntimes;
	fail=prob<UseRandom::rnd()*prob_max;
      }
      while(ntimes>0&&fail);
      if(ntimes==0)
	{throw Exception() << "Tried accept/reject: " << _maxunwgt 
			   << " times and rejected all the times! in"
			   << " EvtGen::randomDecay()" << Exception::eventerror;}
    } 
    // other classes have no correlations
    else if(dinc) {
      cerr << "testing can't use setdaugsDecayedByParentModel\n";
      //dinc->setdaugsDecayedByParentModel(false);
      dinc->decay(part);
      cerr << "can't set decayprob\n";
      //part->setDecayProb(1.);
    }
    else if(dprob) {
      int ntimes(_maxunwgt);
      double prob,dummy;
      do {
	dprob->setWeight(1.);
	cerr << "testing can't use setdaugsDecayedByParentModel\n";
	//dprob->setdaugsDecayedByParentModel(false);
	dprob->decay(part);
	ntimes--;
	cerr << "can't getprob or weight\n";
	//prob=dprob->getProb()/dprob->getWeight();
	dummy=dprob->getProbMax(prob)*UseRandom::rnd();
	cerr << "can't set decayprob\n";
	//part->setDecayProb(prob/dprob->getProbMax(prob));
      }
      while(ntimes&&prob<dummy);
      if(ntimes==0) {
	throw Exception() << "Tried accept/reject: " << _maxunwgt 
			  << " times and rejected all the times! in"
			  << " EvtGen::randomDecay()" << Exception::eventerror;
      }
    }
    else {
      throw Exception() << "Unknown type of EvtGen decayer in EvtGen::randomDecay()"
			<< Exception::abortnow;
    }
    // call photos if needed
    if (decayer->getPHOTOS()||EvtRadCorr::alwaysRadCorr()) {
      EvtRadCorr::doRadCorr(part);
    }
    // check everything o.k.
    done=(EvtStatus::getRejectFlag()==0);
    if(!done) {
      for (ii=0;ii<part->getNDaug();ii++) {
	EvtParticle *temp=part->getDaug(ii);
	temp->deleteTree();
      }
      part->resetFirstOrNot();
      part->resetNDaug();
    }
  }
  while(ntry<_maxtry&&!done);
  if(!done) {
    throw Exception() << "Failed to generate a decay in EvtGen::randomDecayAll"
		      << Exception::eventerror;
  }
  // evtgen will have created decay products of the unstable particles, get rid of them
  int idtemp;
  for(ii=0;ii<part->getNDaug();ii++)
    {
      EvtParticle *temp=part->getDaug(ii);
      idtemp=ThePEGID(temp->getId());
      if(idtemp!=0&&idtemp!=90){temp->deleteDaughters(false);}
    }
  part=unmix;
}

void EvtGen::checkConversion() {
  // check the translation of particles from ThePEG to EvtGen.
  ParticleMap::const_iterator pit  = generator()->particles().begin();
  ParticleMap::const_iterator pend = generator()->particles().end();
  generator()->log() << "Testing conversion of particles from ThePEG to EvtGen\n";
  EvtId etemp;
  for(;pit!=pend;++pit) {
    generator()->log() << pit->first << "     ";
    etemp=EvtGenID(pit->first,false);
    if(etemp.getAlias()>=0) {
      generator()->log() << pit->second->PDGName() << "\t becomes " 
			 << EvtPDL::name(etemp) << "\t " 
			 << EvtPDL::getStdHep(etemp);
      if(ThePEGID(etemp)-pit->first!=0) {
	generator()->log() << " and converting back to ThePEG fails";
      }
      generator()->log() << "\n";
    }
    else {
      generator()->log() << pit->second->PDGName()  
			 << " has no match in EvtGen \n";
    }
  }
  // test the conversion the other way
  generator()->log() << "Testing conversion of particles from EvtGen to ThePEG\n";
  int idtemp;
  tcPDPtr ptemp;
  for(int ix=0;ix<EvtPDL::entries();++ix) {
    generator()->log() << EvtPDL::partlist()[ix].getStdHep() << "     "
		       << EvtPDL::partlist()[ix].getName();
    idtemp=ThePEGID(EvtPDL::partlist()[ix].getId(),false);
    ptemp=getParticleData(idtemp);
    if(ptemp) {
      generator()->log() << " becomes " << ptemp->PDGName()  << "   "
			 << ptemp->id();
      if(EvtGenID(ptemp->id())!=EvtPDL::partlist()[ix].getId()) {
	generator()->log() << " and converting back to EvtGEN fails ";
      }
      generator()->log() << "\n";
    }
    else {
      generator()->log() << " has no match in ThePEG \n";
    }
  }
}

void EvtGen::outputEvtGenDecays(long parentid) {
  // output the decay modes from EvtGen so we can read them in
  EvtId parent(EvtGenID(parentid));
  // get evtids for children
  generator()->log() << "Outputting decays for " 
		     << getParticleData(parentid)->PDGName() << "\n";
  for(int ix=0;ix<EvtDecayTable::getNMode(parent.getAlias());++ix) {
    EvtDecayBase* decayer=EvtDecayTable::getDecay(parent.getAlias(),ix);
    vector<long> pegid;
    for(int iy=0;iy<decayer->getNDaug();++iy) {
      pegid.push_back(ThePEGID(decayer->getDaug(iy),false));
    }
    for(unsigned int iy=pegid.size();iy<7;++iy) pegid.push_back(0);
    double br=decayer->getBranchingFraction();
    generator()->log() 
      << "insert into decay_modes (incomingID,BR,decayon,star,outgoingID1,"
      << "outgoingID2,outgoingID3,outgoingID4,outgoingID5,outgoingID6,outgoingID7,"
      << "description,decayer) values (" << parentid << "," << br << ",'on','*',";
    for(int iy=0;iy<7;++iy) {
      generator()->log() << pegid[iy] << ",";
    }
    generator()->log() << "'Decay of %name% with branching "
		       << "ratio taken from EvtGen.',3);\n";
  }
}


















/*
class EvtGen{


  void readUDecay(const char* const udecay_name);
  
void EvtGen::readUDecay(const char* const uDecayName){

  ifstream indec;

  if ( uDecayName[0] == 0) {
    report(INFO,"EvtGen") << "Is not reading a user decay file!"<<endl;
  }
  else{  
    indec.open(uDecayName);
    if (indec) {
      EvtDecayTable::readDecayFile(uDecayName);
      
      report(INFO,"EvtGen") << "Reading "<<uDecayName
			    <<" to override decay table."<<endl;
    }    
    else{
      
      report(INFO,"EvtGen") << "Can not find UDECAY file '"
			    <<uDecayName<<"'.  Stopping"<<endl;
      ::abort();
    }
  }
  
}
 */
