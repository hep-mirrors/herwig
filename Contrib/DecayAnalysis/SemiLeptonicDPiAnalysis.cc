// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SemiLeptonicDPiAnalysis class.
//

#include "SemiLeptonicDPiAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;
using namespace ThePEG;

void SemiLeptonicDPiAnalysis::findChildren(tPPtr part,ParticleVector & prod) {
  if(part->children().empty()) {
    prod.push_back(part);
  }
  else {
    for(unsigned ix=0;ix<part->children().size();++ix) {
      findChildren(part->children()[ix],prod);
    }
  }
}


void SemiLeptonicDPiAnalysis::analyze(tEventPtr event, long, int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event); 
  tPVector particles;
  // find B mesons
  for(unsigned int ix=0; ix<event->primaryCollision()->steps().size();++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      if((**iter).children().empty()) continue;
      int id=abs((**iter).id());
      if(id==ParticleID::B0||id==ParticleID::Bplus) particles.push_back(*iter);
    }
  }
  // analyse them
  analyze(particles);
}

LorentzRotation SemiLeptonicDPiAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void SemiLeptonicDPiAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void SemiLeptonicDPiAnalysis::analyze(tPPtr part) {
  // find the stable decay products
  ParticleVector products;
  findChildren(part,products);
  // check them
  if(products.size()!=4) return;
  // find the particles
  int ipi=-1,idm=-1,il=-1,inu=-1,id;
  for(unsigned int ix=0;ix<products.size();++ix) {
    id=abs(products[ix]->id());
    if(id==ParticleID::piplus||id==ParticleID::pi0) ipi=ix;
    else if(id==ParticleID::eminus||id==ParticleID::muminus||
	    id==ParticleID::tauminus) il=ix;
    else if(id==ParticleID::nu_e||id==ParticleID::nu_mu||
	    id==ParticleID::nu_tau) inu=ix;
    else if(id==ParticleID::Dplus||id==ParticleID::D0||
	    id==ParticleID::Dstarplus||id==ParticleID::Dstar0) idm=ix;
  }
  // check we have everything
  if(ipi<0||idm<0||il<0||inu<0) return;
  tPPtr lep[2]={products[il],products[inu]};
  unsigned int lid=(abs(lep[0]->id())-9)/2;
  unsigned int ix=0; 
  bool found(false);
  while(!found&&ix<_incoming.size()) {
    if(_incoming[ix]==part->id()&&
       _outgoingD[ix]==products[idm]->id()&&
       _outgoingP[ix]==products[ipi]->id()&&
       lid==_outgoingL[ix]) {
      found=true;
      break;
    }
    else {
      ++ix;
    }
  }
  if(!found) {
    ix=_incoming.size();
    _incoming.push_back(part->id());
    _outgoingD.push_back(products[idm]->id());
    _outgoingP.push_back(products[ipi]->id());
    _outgoingL.push_back(lid);
    _energy.push_back(new_ptr(Histogram(0.0,
					(part->nominalMass()+part->dataPtr()->widthUpCut()
					-products[idm]->nominalMass()
					 +products[idm]->dataPtr()->widthLoCut())/GeV,200)));
    _scale.push_back(new_ptr(Histogram(0.0,
				       (part->nominalMass()+part->dataPtr()->widthUpCut()
				       -products[ipi]->nominalMass()
				       +products[ipi]->dataPtr()->widthLoCut())/GeV,200)));
    _mDpi.push_back(new_ptr(Histogram(4.0,25.,200)));
  }
  // add the results to the histogram
  Lorentz5Momentum ptemp;
  ptemp = lep[0]->momentum()+lep[1]->momentum();
  ptemp.rescaleMass();
  *_scale[ix]+=ptemp.mass()/GeV;
  ptemp = products[idm]->momentum()+lep[1]->momentum();
  ptemp.rescaleMass();
  Energy ee = 1./2./part->mass()*
    (part->mass()*part->mass()-ptemp.mass()*ptemp.mass()+lep[1]->mass()*lep[1]->mass());
  *_energy[ix]+=ee/GeV;
  double mdpi=(products[idm]->momentum()+
	       products[ipi]->momentum()).m2()/GeV2;
  *_mDpi[ix]+=mdpi;
}

void SemiLeptonicDPiAnalysis::persistentOutput(PersistentOStream &) const {}

void SemiLeptonicDPiAnalysis::persistentInput(PersistentIStream &, int) {}

ClassDescription<SemiLeptonicDPiAnalysis> SemiLeptonicDPiAnalysis::initSemiLeptonicDPiAnalysis;
// Definition of the static class description member.

void SemiLeptonicDPiAnalysis::Init() {

  static ClassDocumentation<SemiLeptonicDPiAnalysis> documentation
    ("There is no documentation for the SemiLeptonicDPiAnalysis class");

}

void SemiLeptonicDPiAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  string title,temp;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    title= getParticleData(_incoming[ix])->PDGName() +
      " -> " 
      + getParticleData(_outgoingD[ix])->PDGName()+  " "
      + getParticleData(_outgoingP[ix])->PDGName()+  " " +
      getParticleData(9+2*_outgoingL[ix])->PDGName() + " " +
      getParticleData(10+2*_outgoingL[ix])->PDGName();
    temp = "Mass for l nu in " +title;
    using namespace HistogramOptions;
    _scale[ix]->topdrawOutput(outfile,Frame,"BLACK",temp);
    temp = "Lepton energy for in " +title;
    _energy[ix]->topdrawOutput(outfile,Frame,"BLACK",temp);
    temp = "D(*) pi mass in " +title;
    _mDpi[ix]->topdrawOutput(outfile,Frame,"BLACK",temp);
  }
}

