// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZPhotonsAnalysis class.
//

#include "ZPhotonsAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void ZPhotonsAnalysis::analyze(tEventPtr event, long, int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  tPVector particles;
  // Rotate to CMS, extract final state particles and call analyze(particles).
  for(unsigned int ix=1, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part(event->primaryCollision()->step(ix)->all());
    ThePEG::ParticleSet::iterator iter(part.begin()),end(part.end());
    for( ;iter!=end;++iter) {
      if((**iter).id()==ParticleID::Z0&&(**iter).children().size()>=2) {
	particles.push_back(*iter);
      }
    }
  }
  // analyse them
  analyze(particles);
}

void ZPhotonsAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void ZPhotonsAnalysis::analyze(tPPtr part) {
  // check we have the right decay
  if(abs(part->children()[0]->id())!=_iferm) return;
  if(abs(part->children()[1]->id())!=_iferm) return;
  // pphoton is the total momentum of all photons, mult is the multiplicity,
  // ix/emax is the index/energy of the most energetic photon emitted:
  Lorentz5Momentum pphoton;
  unsigned int mult=0;
  Energy emax=ZERO;
  int imax=-1;
  // loop over the non-fermionic children i.e. the photons:
  for(unsigned int ix=2;ix<part->children().size();++ix) {
    if(part->children()[ix]->id()!=ParticleID::gamma) return;
    pphoton+=part->children()[ix]->momentum();
    *_etotal[4]+=part->children()[ix]->momentum().e()/MeV;
    if(part->children()[ix]->momentum().e()>emax) {
      emax=part->children()[ix]->momentum().e();
      imax=ix;
    }
    Lorentz5Momentum pf(part->children()[0]->momentum());
    Lorentz5Momentum pfb(part->children()[1]->momentum());
    if(part->children()[0]->id()<0&&part->children()[1]>0) swap(pf,pfb);
    pf.boost(-part->momentum().boostVector());
    pfb.boost(-part->momentum().boostVector());
    // bin the cosine of the angle between each photon and the fermion 
    Lorentz5Momentum pphot(part->children()[ix]->momentum());
    pphot.boost(-part->momentum().boostVector());
    *_cphoton+=pf.vect().cosTheta(pphot.vect());
    ++mult;
  }
  // boost the combined fermion momentum and the total photon momentum
  // to the Z/gamma rest framefermion rest frame:
  Lorentz5Momentum pferm(part->children()[0]->momentum()+
			 part->children()[1]->momentum());
  pferm.boost(-part->momentum().boostVector());
  pphoton.boost(-part->momentum().boostVector());
  pferm.rescaleMass();
  pphoton.rescaleMass();
  // bin the photon multiplicity of this decay:
  *_nphoton+=mult;
  for(unsigned int ix=0;ix<3;++ix) {
    // bin the invariant mass of the fermions and the total photon energy
    // in histograms of varying ranges (ix):  
    *_masstotal[ix]+=pferm.mass()/MeV;
    if(mult>0) *_etotal[ix]+=pphoton.e()/MeV;
  }
  // bin the energy of the most energetic photon:
  if(imax>0) *_etotal[3]+=part->children()[imax]->momentum().e()/MeV;
  if(mult<20) {
    // bin the difermion invariant mass and the total photon energy 
    // according to multiplicity:
    *_mphoton[mult]+=pferm.mass()/MeV;
    *_ephoton[mult]+=pphoton.e()/MeV;
  }
}

ClassDescription<ZPhotonsAnalysis> ZPhotonsAnalysis::initZPhotonsAnalysis;
// Definition of the static class description member.

void ZPhotonsAnalysis::Init() {

  static ClassDocumentation<ZPhotonsAnalysis> documentation
    ("There is no documentation for the ZPhotonsAnalysis class");

  static Parameter<ZPhotonsAnalysis,int> interfaceFermion
    ("Fermion",
     "Id code of fermion",
     &ZPhotonsAnalysis::_iferm, 11, 11, 15,
     false, false, Interface::limited);

}

inline void ZPhotonsAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _nphoton->topdrawOutput(output,Frame|Errorbars|Ylog,
			  "RED",
			  "Photon Multiplicity",
			  "                   ",
			  "1/SdS/dN0G1",
			  "  G G   XGX",
			  "N0G1",
			  " XGX");
  for(unsigned int ix=0;ix<3;++ix) {
    _masstotal[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
			  "RED",
			  "Fermion mass for all events",
			  "                          ",
			  "1/SdS/d/GeV2-13",
			  "  G G      X  X",
			  "m0l2+3l2-31/GeV",
			  " X X X X XX    ");
    _etotal[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
			  "RED",
			  "Photon Energy for all events",
			  "                   ",
			  "1/SdS/dE0G1/GeV2-13",
			  "  G G   XGX    X  X",
			  "E0G1/GeV",
			  " XGX    ");
  }
  _etotal[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "Photon Energy for all events",
			    "                   ",
			    "1/SdS/dE0G1/GeV2-13",
			    "  G G   XGX    X  X",
			    "E0G1/GeV",
			    " XGX    ");
  _etotal[4]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "Photon Energy for all events",
			    "                   ",
			    "1/SdS/dE0G1/GeV2-13",
			    "  G G   XGX    X  X",
			    "E0G1/GeV",
			    " XGX    ");

  _cphoton->topdrawOutput(output,Frame|Errorbars|Ylog,
                          "RED","Photon cosine wrt fermion","",
                          "1/SdS/dc0G1",
                          "  G G   XGX",
                          "c0G1",
                          " XGX    ");
  for(unsigned int ix=0;ix<20;++ix) {
    ostringstream titlea;
    titlea << "Fermion mass for "  << ix << " photons " << flush;
    _mphoton[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",titlea.str(),"",
				"1/SdS/d/GeV2-13",
				"  G G      X  X",
				"m0l2+3l2-31/GeV",
				" X X X X XX    ");
    ostringstream titleb;
    titleb << "photon energy for " << ix << " photons " << flush;
    _ephoton[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",titleb.str(),"",
				"1/SdS/dE0G1/GeV2-13",
				"  G G   XGX    X  X",
				"E0G1/GeV",
				" XGX    ");;
  }
}

inline void ZPhotonsAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  _masstotal.push_back(new_ptr(Histogram(0.,92000.,400)));
  _etotal   .push_back(new_ptr(Histogram(0.,92000.,400)));
  _masstotal.push_back(new_ptr(Histogram(0.,80000.,200)));
  _etotal   .push_back(new_ptr(Histogram(0.,10000.,200)));
  _masstotal.push_back(new_ptr(Histogram(80000.,92000.,200)));
  _etotal   .push_back(new_ptr(Histogram(10000.,92000.,200)));
  for(unsigned int ix=0;ix<20;++ix) {
    _mphoton.push_back(new_ptr(Histogram(0.,92000.,400)));
    _ephoton.push_back(new_ptr(Histogram(0.,92000.,400)));
  }
  _cphoton=new_ptr(Histogram(-1.,1.,100));
  _etotal.push_back(new_ptr(Histogram(0.,92000.,400)));
  _etotal.push_back(new_ptr(Histogram(0.,92000.,400)));
  _nphoton=new_ptr(Histogram(-0.5,100.5,101));
}

void ZPhotonsAnalysis::persistentOutput(PersistentOStream & os) const {
  os <<_iferm;
}

void ZPhotonsAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> _iferm;
}
