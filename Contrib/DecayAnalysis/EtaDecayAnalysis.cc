// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaDecayAnalysis class.
//

#include "EtaDecayAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;
void EtaDecayAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).  
  // find all eta and eta' particles 
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix)
    {
      ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
      ThePEG::ParticleSet::iterator iter=part.begin();
      ThePEG::ParticleSet::iterator end=part.end();
      for( ;iter!=end;++iter)
	{if((**iter).id()==ParticleID::eta||(**iter).id()==ParticleID::etaprime)
	    {particles.push_back(*iter);}}
    }
  // analyse them
  analyze(particles);
}

void EtaDecayAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void EtaDecayAnalysis::analyze(tPPtr part) {
  // ensure 2 or 3 decay products
  if(part->children().size()!=3&&part->children().size()!=2) return;
  int imeson;
  if(part->id()==ParticleID::eta) imeson=0;
  else if(part->id()==ParticleID::etaprime) imeson=1;
  else return;
  // pi0 gamma gamma analysis
  if(part->children().size()==3&&
     part->children()[0]->id()==ParticleID::pi0&&
     part->children()[1]->id()==ParticleID::gamma&&
     part->children()[2]->id()==ParticleID::gamma) {
    Lorentz5Momentum ptemp=
      part->children()[1]->momentum()+part->children()[2]->momentum();
    ptemp.rescaleMass();
    *_mgammagamma[imeson] +=ptemp.mass()/MeV;
    ptemp=part->children()[0]->momentum()+part->children()[1]->momentum();
    ptemp.rescaleMass();
    *_mpi0gamma[imeson]   +=ptemp.mass()/MeV;
    ptemp=part->children()[0]->momentum()+part->children()[2]->momentum();
    ptemp.rescaleMass();
    *_mpi0gamma[imeson]   +=ptemp.mass()/MeV;
  }
  // pi+pi-gamma analysis
  else if((part->children().size()==3&&
	   part->children()[0]->id()==ParticleID::piplus&&
	   part->children()[1]->id()==ParticleID::piminus&&
	   part->children()[2]->id()==ParticleID::gamma)||
	  (part->children().size()==2&&
	   part->children()[0]->id()==ParticleID::rho0&&
	   part->children()[1]->id()==ParticleID::gamma)) {
    Lorentz5Momentum pout[3];
    if(part->children().size()==2) {
      pout[0]=part->children()[0]->children()[0]->momentum();
      pout[1]=part->children()[0]->children()[1]->momentum();
      pout[2]=part->children()[1]->momentum();
    }
    else {
      pout[0]=part->children()[0]->momentum();
      pout[1]=part->children()[1]->momentum();
      pout[2]=part->children()[2]->momentum();
    }
    Lorentz5Momentum ptemp=pout[0]+pout[1];
    ptemp.rescaleMass();
    *_mpippim[imeson]+=ptemp.mass()/MeV;
    Energy egamma = 
      0.5*(part->mass()*part->mass()-ptemp.mass()*ptemp.mass())/part->mass();
    *_photonenergy[imeson]+=egamma/MeV;
    ptemp=pout[imeson]+pout[2];ptemp.rescaleMass();
    *_mpipgamma[imeson]+=ptemp.mass()/MeV;
    ptemp=pout[1]+pout[2];ptemp.rescaleMass();
    *_mpimgamma[imeson]+=ptemp.mass()/MeV;
  }
  else {
    vector<Lorentz5Momentum> ppim,ppip,ppi0,peta; 
    for(unsigned int ix=0;ix<part->children().size();++ix) {
       long id = part->children()[ix]->id();
       if(id==ParticleID::piplus)       
	 ppip  .push_back(part->children()[ix]->momentum());
       else if(id==ParticleID::piminus) 
	 ppim  .push_back(part->children()[ix]->momentum());
       else if(id==ParticleID::pi0)     
	 ppi0  .push_back(part->children()[ix]->momentum());
       else if(id==ParticleID::eta)     
	 peta  .push_back(part->children()[ix]->momentum());
    }
    if(ppi0.size()==3) {
      *_dpi0pi0[imeson]+=(ppi0[0]+ppi0[1]).m()/MeV;
      *_dpi0pi0[imeson]+=(ppi0[0]+ppi0[2]).m()/MeV;
      *_dpi0pi0[imeson]+=(ppi0[1]+ppi0[2]).m()/MeV;
    }
    else if(ppip.size()==1&&ppim.size()==1&&ppi0.size()==1) {
      *_dpi0pip[imeson]+=(ppi0[0]+ppip[0]).m()/MeV;
      *_dpi0pim[imeson]+=(ppi0[0]+ppim[0]).m()/MeV;
      *_dpippim[imeson]+=(ppip[0]+ppim[0]).m()/MeV;
    }
    else if(ppi0.size()==2&&peta.size()==1) {
      *_dpi0pi0[2]+=(ppi0[0]+ppi0[1]).m()/MeV;
      *_dpi0eta[0]+=(ppi0[0]+peta[0]).m()/MeV;
      *_dpi0eta[0]+=(ppi0[1]+peta[0]).m()/MeV;
    }
    else if(ppip.size()==1&&ppim.size()==1&&peta.size()==1) {
      *_dpippim[2]+=(ppip[0]+ppim[0]).m()/MeV;
      *_dpipeta[0]+=(ppip[0]+peta[0]).m()/MeV;
      *_dpimeta[0]+=(ppim[0]+peta[0]).m()/MeV;
    }
  }
}

NoPIOClassDescription<EtaDecayAnalysis> EtaDecayAnalysis::initEtaDecayAnalysis;
// Definition of the static class description member.

void EtaDecayAnalysis::Init() {

  static ClassDocumentation<EtaDecayAnalysis> documentation
    ("The EtaDecayAnalysis class performs the analysis of the decays of eta and eta\'"
     " mesons");

}

void EtaDecayAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _mgammagamma[0]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 "GG mass in HRP203GG",
				 "GG         GWGX XGG",
				 "1/GdG/dm0GG1/MeV2-13",
				 "  F F   XGGX    X  X",
				 "m0GG1/MeV",
				 " XGGX    ");
  _mpi0gamma[0]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 "P203G mass in HRP203GG",
				 "GX XG         GWGX XGG",
				 "1/GdG/dm0P203G1/MeV2-13",
				 "  F F   XGX XGX    X  X",
				 "m0P203G1/MeV",
				 " XGX XGX    ");
  _mgammagamma[1]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 "GG mass in H'RP203GG",
				 "GG         G WGX XGG",
				 "1/GdG/dm0GG1/MeV2-13",
				 "  F F   XGGX    X  X",
				 "m0GG1/MeV",
				 " XGGX    ");
  _mpi0gamma[1]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 "P203G mass in H'RP203GG",
				 "GX XG         G WGX XGG",
				 "1/GdG/dm0P203G1/MeV2-13",
				 "  F F   XGX XGX    X  X",
				 "m0P203G1/MeV",
				 " XGX XGX    ");
  _mpipgamma[0]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 "P2+3G mass in HRP2+3P2-3G",
				 "GX XG         GWGX XGX XG",
				 "1/GdG/dm0P2+3G1/MeV2-13",
				 "  F F   XGX XGX    X  X",
				 "m0P2+3G1/MeV",
				 " XGX XGX    ");
  _mpimgamma[0]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 "P2-3G mass in HRP2+3P2-3G",
				 "GX XG         GWGX XGX XG",
				 "1/GdG/dm0P2-3G1/MeV2-13",
				 "  F F   XGX XGX    X  X",
				 "m0P2-3G1/MeV",
				 " XGX XGX    ");
  _photonenergy[0]->topdrawOutput(output,Frame|Errorbars,
				  "RED",
				  "G Energy in HRP2+3P2-3G",
				  "G           GWGX XGX XG",
				  "1/GdG/dE0G1/MeV2-13",
				  "  F F   XGX    X  X",
				  "E0G1/MeV",
				  " XGX    ");
  _mpippim[0]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P2-3P2+3 mass in HRP2+3P2-3G",
			     "GX XGX X         GWGX XGX XG",
			     "1/GdG/dm0P2-3P2+31/MeV2-13",
			     "  F F   0GX XGX X1    X  X",
			     "m0P2-3P2+31/MeV",
			     " XGX XGX XX    ");
  _mpipgamma[1]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 "P2+3G mass in H'RP2+3P2-3G",
				 "GX XG         G WGX XGX XG",
				 "1/GdG/dm0P2+3G1/MeV2-13",
				 "  F F   XGX XGX    X  X",
				 "m0P2+3G1/MeV",
				 " XGX XGX    ");
  _mpimgamma[1]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 "P2-3G mass in H'RP2+3P2-3G",
				 "GX XG         G WGX XGX XG",
				 "1/GdG/dm0P2-3G1/MeV2-13",
				 "  F F   XGX XGX    X  X",
				 "m0P2-3G1/MeV",
				 " XGX XGX    ");
  _photonenergy[1]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 "G Energy in H'RP2+3P2-3G",
				 "G           G WGX XGX XG",
				 "1/GdG/dE0G1/MeV2-13",
				 "  F F   XGX    X  X",
				 "E0G1/MeV",
				 " XGX    ");
  _mpippim[1]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P2-3P2+3 mass in H'RP2+3P2-3G",
			     "GX XGX X         G WGX XGX XG",
			     "1/GdG/dm0P2-3P2+31/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P2-3P2+31/MeV",
			     " XGX XGX XX    ");
  _dpi0pi0[0]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P203P203 mass in HRP203P203P203",
			     "GX XGX X         GWGX XGX XGX X",
			     "1/GdG/dm0P203P2031/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P203P2031/MeV",
			     " XGX XGX XX    ");
  _dpi0pi0[1]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P203P203 mass in H'RP203P203P203",
			     "GX XGX X         G WGX XGX XGX X",
			     "1/GdG/dm0P203P2031/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P203P2031/MeV",
			     " XGX XGX XX    ");
  _dpippim[0]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P2-3P2+3 mass in HRP2+3P2-3P203",
			     "GX XGX X         GWGX XGX XGX X",
			     "1/GdG/dm0P2-3P2+31/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P2-3P2+31/MeV",
			     " XGX XGX XX    ");
  _dpi0pip[0]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P203P2+3 mass in HRP2+3P2-3P203",
			     "GX XGX X         GWGX XGX XGX X",
			     "1/GdG/dm0P203P2+31/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P203P2+31/MeV",
			     " XGX XGX XX    ");
  _dpi0pim[0]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P203P2-3 mass in HRP2+3P2-3P203",
			     "GX XGX X         GWGX XGX XGX X",
			     "1/GdG/dm0P203P2-31/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P203P2-31/MeV",
			     " XGX XGX XX    ");
  _dpi0pim[1]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P203P2-3 mass in H'RP2+3P2-3P203",
			     "GX XGX X         G WGX XGX XGX X",
			     "1/GdG/dm0P203P2-31/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P203P2-31/MeV",
			     " XGX XGX XX    ");
  _dpi0pip[1]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P203P2+3 mass in H'RP2+3P2-3P203",
			     "GX XGX X         G WGX XGX XGX X",
			     "1/GdG/dm0P203P2+31/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P203P2+31/MeV",
			     " XGX XGX XX    ");
  _dpippim[1]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P2-3P2+3 mass in H'RP2+3P2-3P203",
			     "GX XGX X         G WGX XGX XGX X",
			     "1/GdG/dm0P2-3P2+31/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P2-3P2+31/MeV",
			     " XGX XGX XX    ");
  _dpippim[2]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P2-3P2+3 mass in H'RP2+3P2-3H",
			     "GX XGX X         G WGX XGX XG",
			     "1/GdG/dm0P2-3P2+31/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P2-3P2+31/MeV",
			     " XGX XGX XX    ");
  _dpipeta[0]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "HP2+3 mass in H'RP2+3P2-3H",
			     "GGX X         G WGX XGX XG",
			     "1/GdG/dm0HP2+31/MeV2-13",
			     "  F F   XGGX XX    X  X",
			     "m0HP2+31/MeV",
			     " XGGX XX    ");
  _dpimeta[0]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "HP2-3 mass in H'RP2+3P2-3H",
			     "GGX X         G WGX XGX XG",
			     "1/GdG/dm0HP2-31/MeV2-13",
			     "  F F   XGGX XX    X  X",
			     "m0HP2-31/MeV",
			     " XGGX XX    ");
  _dpi0pi0[2]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "P203P203 mass in H'RP203P203H",
			     "GX XGX X         G WGX XGX XG",
			     "1/GdG/dm0P203P2031/MeV2-13",
			     "  F F   XGX XGX XX    X  X",
			     "m0P203P2031/MeV",
			     " XGX XGX XX    ");
  _dpi0eta[0]->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "HP203 mass in H'RP203P203H",
			     "GGX X         G WGX XGX XG",
			     "1/GdG/dm0HP2031/MeV2-13",
			     "  F F   XGGX XX    X  X",
			     "m0HP2031/MeV",
			     " XGGX XX    ");
}

void EtaDecayAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  double meta[2]={547.45, 957.78};
  for(unsigned int ix=0;ix<2;++ix) {
    _mgammagamma .push_back(new_ptr(Histogram(0.,meta[ix],200)));
    _mpi0gamma   .push_back(new_ptr(Histogram(0.,meta[ix],200)));
    _mpipgamma   .push_back(new_ptr(Histogram(0.,meta[ix],200)));
    _mpimgamma   .push_back(new_ptr(Histogram(0.,meta[ix],200)));
    _photonenergy.push_back(new_ptr(Histogram(0.,meta[ix],200)));
    _mpippim     .push_back(new_ptr(Histogram(0.,meta[ix],200)));

    _dpippim     .push_back(new_ptr(Histogram(200.,meta[ix],200)));
    _dpi0pi0     .push_back(new_ptr(Histogram(200.,meta[ix],200)));
    _dpi0pip     .push_back(new_ptr(Histogram(200.,meta[ix],200)));
    _dpi0pim     .push_back(new_ptr(Histogram(200.,meta[ix],200)));
  }
  _dpi0pi0.push_back(new_ptr(Histogram(200.,500.,200)));
  _dpippim.push_back(new_ptr(Histogram(200.,500.,200)));
  _dpipeta.push_back(new_ptr(Histogram(500.,meta[1],200)));
  _dpimeta.push_back(new_ptr(Histogram(500.,meta[1],200)));
  _dpi0eta.push_back(new_ptr(Histogram(500.,meta[1],200)));
}
