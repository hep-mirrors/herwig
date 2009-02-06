// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VffGammaAnalysis class.
//

#include "VffGammaAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

VffGammaAnalysis::VffGammaAnalysis() {
  // c cbar mesons
  _id.push_back(443);_id.push_back(100443);_id.push_back(30443);
  // b bbar mesons
  _id.push_back(553);_id.push_back(100553);
  _id.push_back(200553);_id.push_back(300553);
}

void VffGammaAnalysis::analyze(tEventPtr event, long ,
			       int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find all particles
  tPVector particles;
  unsigned int idtemp[2];
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      if((**iter).children().size()>=2) {
	idtemp[0]=abs((**iter).children()[0]->id());
	idtemp[1]=abs((**iter).children()[1]->id());
	if(idtemp[0]>10&&idtemp[0]<17&&idtemp[0]%2==1&&
	   idtemp[1]>10&&idtemp[1]<17&&idtemp[1]%2==1) {
	  for(unsigned int iy=0;iy<_id.size();++iy) {
	    if((**iter).id()==_id[iy]){particles.push_back(*iter);
	    }
	  }
	}
      }
    }
  }
  // analyse them
  analyze(particles);
}

void VffGammaAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void VffGammaAnalysis::analyze(tPPtr part) {
  // find the decaying particle
  int imode(-1),id(part->id());
  unsigned int ix(0);
  do {
    if(id==_id[ix]){imode=ix;}
    ++ix;
  }
  while(imode<0&&ix<_id.size());
  if(imode<0) return;
  // find the type of lepton
  imode = 3*imode+(abs(part->children()[0]->id())-11)/2;
  Lorentz5Momentum pphoton;
  unsigned int mult(0);
  Energy emax(ZERO);
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

NoPIOClassDescription<VffGammaAnalysis> VffGammaAnalysis::initVffGammaAnalysis;
// Definition of the static class description member.

void VffGammaAnalysis::Init() {

  static ClassDocumentation<VffGammaAnalysis> documentation
    ("There is no documentation for the VffGammaAnalysis class");

  static ParVector<VffGammaAnalysis,long> interfaceId
    ("Id",
     "PDG codes of the particles to be analysed",
     &VffGammaAnalysis::_id, -1, long(0), 0, 0,
     false, false, Interface::nolimits);

}

inline void VffGammaAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  string titlea,titleb,titlec;
  for(unsigned int ix=0;ix<_id.size();++ix) {
    titlea = " for decay " + getParticleData(_id[ix])->PDGName() + " -> ";
    for(unsigned int iy=0;iy<3;++iy) {
      if(iy==0){titleb = titlea + "e+   e-";}
      else if(iy==1){titleb = titlea + "mu+  mu-";}
      else if(iy==2){titleb = titlea + "tau+ tau-";}
      // number of photons
      titlec = "Photon multiplicity " +titleb;
      using namespace HistogramOptions;
      _nphoton[3*ix+iy]->topdrawOutput(output,Frame|Errorbars|Ylog,
				       "RED",titlec,""
			    "1/SdS/dN0G1",
			    "  G G   XGX",
			    "N0G1",
			    " XGX");
      // fermion masses
      titlec = "Mass of the charged decay products "+titleb;
      _masstotal[3*ix+iy]->topdrawOutput(output,Frame|Errorbars|Ylog,
					  "RED",titlec,""
					  "1/SdS/dm0l2+3l2-31/GeV2-13",
					  "  G G   X X X X XX    X  X",
					  "m0l2+3l2-31",
					  " X X X X XX");
      // total photon energy
      _masstotal[3*ix+iy]->topdrawOutput(output,Frame|Errorbars|Ylog,
					 "RED",titlec,""
					 "1/SdS/dE0G1/GeV2-13",
					 "  G G   XGX    X  X",
					 "E0G1",
					 " XGX");
      titlec="Single photon energy for all events"+titleb;
      _esingle[3*ix+iy]->topdrawOutput(output,Frame|Errorbars|Ylog,
					 "RED",titlec,""
					 "1/SdS/dE0G1/GeV2-13",
					 "  G G   XGX    X  X",
					 "E0G1",
					 " XGX");
	  // all photon
      titlec="All photon energy for all events"+titleb;
      _eall[3*ix+iy]->topdrawOutput(output,Frame|Errorbars|Ylog,
					 "RED",titlec,""
					 "1/SdS/dE0G1/GeV2-13",
					 "  G G   XGX    X  X",
					 "E0G1",
					 " XGX");
    }
  }
}

inline void VffGammaAnalysis::doinitrun() {
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
