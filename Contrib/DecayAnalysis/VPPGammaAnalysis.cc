// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VPPGammaAnalysis class.
//

#include "VPPGammaAnalysis.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

VPPGammaAnalysis::VPPGammaAnalysis() {
  _id.push_back(ParticleID::rho0);
  _outgoing1.push_back(ParticleID::piplus);
  _outgoing2.push_back(ParticleID::piminus);
  _id.push_back(ParticleID::phi);
  _outgoing1.push_back(ParticleID::piplus);
  _outgoing2.push_back(ParticleID::piminus);
  _id.push_back(ParticleID::Jpsi);
  _outgoing1.push_back(ParticleID::piplus);
  _outgoing2.push_back(ParticleID::piminus);
  _id.push_back(ParticleID::phi);
  _outgoing1.push_back(ParticleID::Kplus);
  _outgoing2.push_back(ParticleID::Kminus);
  _id.push_back(ParticleID::Jpsi);
  _outgoing1.push_back(ParticleID::Kplus);
  _outgoing2.push_back(ParticleID::Kminus);
  _id.push_back(ParticleID::Kstar0);
  _outgoing1.push_back(ParticleID::Kplus);
  _outgoing2.push_back(ParticleID::piminus);
}

void VPPGammaAnalysis::analyze(tEventPtr event, long, int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find all decaying particles
  tPVector particles;
  int id;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      id=abs((**iter).id());
      if((**iter).children().size()>=2&&
	 find(_id.begin(),_id.end(),id)!=_id.end())
	particles.push_back(*iter);
    }
  }
}

void VPPGammaAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void VPPGammaAnalysis::analyze(tPPtr part) {
  int imode(-1);
  long id [3]={part->id(),part->children()[0]->id(),part->children()[1]->id()};
  long idb[3];
  if(part->dataPtr()->CC()) idb[0]=part->dataPtr()->CC()->id();
  else                      idb[0]=id[0];
  for(unsigned int ix=0;ix<2;++ix) {
    if(part->children()[ix]->dataPtr()->CC()) 
      idb[ix+1]=part->children()[ix]->dataPtr()->CC()->id();
    else
      idb[ix+1]=id[ix+1];
  }
  for(unsigned int ix=0;ix<_id.size();++ix) {
    if(id[0]==_id[ix]) {
      if((_outgoing1[ix]==id[1]&&_outgoing2[ix]==id[2])||
	 (_outgoing1[ix]==id[2]&&_outgoing2[ix]==id[1])) {
	imode=ix;
      }
    }
    else if(idb[0]==_id[ix]) {
      if((_outgoing1[ix]==idb[1]&&_outgoing2[ix]==idb[2])||
	 (_outgoing1[ix]==idb[2]&&_outgoing2[ix]==idb[1])) {
	imode=ix;
      }
    }
    if(imode>0) break;
  }
  Lorentz5Momentum pphoton;
  unsigned int mult=0;
  Energy emax=ZERO;
  int imax=-1;
  Lorentz5Momentum pferm(part->children()[0]->momentum()+
			 part->children()[1]->momentum());
  pferm.boost(-part->momentum().boostVector());
  pferm.rescaleMass();
  *_masstotal[imode] += pferm.mass()/MeV;
  if(part->children().size()==2) {
    *_etotal[ imode]+=0.;
    *_esingle[imode]+=0.;
    *_eall[   imode]+=0.;
    *_nphoton[imode]+=0.;
    return;
  }
  Lorentz5Momentum ptemp;
  for(unsigned int ix=2;ix<part->children().size();++ix) {
    if(part->children()[ix]->id()!=ParticleID::gamma){return;}
    pphoton+=part->children()[ix]->momentum();
    ptemp=part->children()[ix]->momentum();
    ptemp.boost(-part->momentum().boostVector());
    *_eall[imode]+=ptemp.e()/MeV;
    if(part->children()[ix]->momentum().e()>emax) {
      emax=part->children()[ix]->momentum().e();
      imax=ix;
    }
    ++mult;
  }
  pphoton.boost(-part->momentum().boostVector());
  pphoton.rescaleMass();
  *_etotal[imode]+=pphoton.e()/MeV;
  ptemp=part->children()[imax]->momentum();
  ptemp.boost(-part->momentum().boostVector());
  *_esingle[imode]+=ptemp.e()/MeV;
  *_nphoton[imode]+=mult;
}

NoPIOClassDescription<VPPGammaAnalysis> VPPGammaAnalysis::initVPPGammaAnalysis;
// Definition of the static class description member.

void VPPGammaAnalysis::Init() {

  static ClassDocumentation<VPPGammaAnalysis> documentation
    ("There is no documentation for the VPPGammaAnalysis class");

  static ParVector<VPPGammaAnalysis,long> interfaceIncoming
    ("Incoming",
     "The PDG code of the incoming particle",
     &VPPGammaAnalysis::_id, -1, long(0), 0, 0,
     false, false, Interface::nolimits);

  static ParVector<VPPGammaAnalysis,long> interfaceOutgoing1
    ("Outgoing1",
     "The PDG code of the first outgoing particle",
     &VPPGammaAnalysis::_outgoing1, -1, long(0), 0, 0,
     false, false, Interface::nolimits);

  static ParVector<VPPGammaAnalysis,long> interfaceOutgoing2
    ("Outgoing2",
     "The PDG code of the second outgoing particle",
     &VPPGammaAnalysis::_outgoing2, -1, long(0), 0, 0,
     false, false, Interface::nolimits);

}

void VPPGammaAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  string titlea,titleb;
  for(unsigned int ix=0;ix<_id.size();++ix) {
    titlea = " for decay " + getParticleData(_id[ix])->PDGName() + " -> "
      + getParticleData(_outgoing1[ix])->PDGName() + " "
      + getParticleData(_outgoing2[ix])->PDGName();
    // number of photons
    titleb = "Photon multiplicity " +titlea;
    using namespace HistogramOptions;
    _nphoton[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",titleb,""
				"1/SdS/dN0G1",
				"  G G   XGX",
				"N0G1",
				" XGX");
    // fermion masses
    titleb = "Mass of the charged decay products "+titlea;
    _masstotal[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
				  "RED",titleb,""
				  "1/SdS/dm0l2+3l2-31/GeV2-13",
				  "  G G   X X X X XX    X  X",
				  "m0l2+3l2-31",
				  " X X X X XX");
    // total photon energy
    _masstotal[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
				  "RED",titleb,""
				  "1/SdS/dE0G1/GeV2-13",
				  "  G G   XGX    X  X",
				  "E0G1",
				  " XGX");
    titleb="Single photon energy for all events"+titlea;
    _esingle[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",titleb,""
				"1/SdS/dE0G1/GeV2-13",
				"  G G   XGX    X  X",
				"E0G1",
				" XGX");
    // all photon
    titleb="All photon energy for all events"+titlea;
    _eall[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",titleb,""
			     "1/SdS/dE0G1/GeV2-13",
			     "  G G   XGX    X  X",
			     "E0G1",
			     " XGX");
  }
}

void VPPGammaAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  for(unsigned int ix=0;ix<_id.size();++ix) {
    double mass = getParticleData(_id[ix])->mass()/MeV;
    for(unsigned iy=0;iy<3;++iy)  {
      _masstotal.push_back(new_ptr(Histogram(0.,mass,1000)));
      _etotal.push_back(new_ptr(Histogram(0.,mass,1000)));
      _eall.push_back(new_ptr(Histogram(0.,mass,1000)));
      _esingle.push_back(new_ptr(Histogram(0.,mass,1000)));
      _nphoton.push_back(new_ptr(Histogram(-0.5,20.5,21)));
    }
  }
}
