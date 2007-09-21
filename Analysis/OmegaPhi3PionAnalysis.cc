// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaPhi3PionAnalysis class.
//

#include "OmegaPhi3PionAnalysis.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void OmegaPhi3PionAnalysis::analyze(tEventPtr event,long , int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find all omega and phi particles 
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      if((**iter).id()==ParticleID::omega||(**iter).id()==ParticleID::phi) {
	particles.push_back(*iter);
      }
    }
  }
  // analyse them
  analyze(particles);
}

void OmegaPhi3PionAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void OmegaPhi3PionAnalysis::analyze(tPPtr part) {
  Lorentz5Momentum pip,pim,pi0; unsigned int imode;
  bool allowed(false);
  // reorder the particles
  ParticleVector children=part->children();
  if(children.size()==2) {
    // vector first
    if(abs(children[1]->id())%10==3) swap(children[0],children[1]);
    if((children[0]->id()==113||abs(children[0]->id())==213)&&
       (children[1]->id()==111||abs(children[1]->id())==211)) {
      allowed=true;
      vector<tPPtr> temp;
      temp.push_back(children[1]);
      temp.push_back(children[0]->children()[0]);
      temp.push_back(children[0]->children()[1]);
      for(unsigned int ix=0;ix<3;++ix) {
	if(temp[ix]->id()== 111) pi0=temp[ix]->momentum();
	if(temp[ix]->id()== 211) pip=temp[ix]->momentum();
	if(temp[ix]->id()==-211) pim=temp[ix]->momentum();
      }
    }
  }
  else if(children.size()==3) {
    // neutral pion first
    if(children[1]->id()==111) swap(children[0],children[1]);
    if(children[2]->id()==111) swap(children[0],children[2]);
    // postive pion second
    if(children[2]->id()==211) swap(children[1],children[2]);
    if(children[0]->id()== 111&&children[1]->id()==211&&
       children[2]->id()==-211) {
      allowed=true;
      pi0=children[0]->momentum();
      pip=children[1]->momentum();
      pim=children[2]->momentum();
    }
  }
  if(!allowed) return;
  if(part->id()==ParticleID::omega) imode=0;
  else if(part->id()==ParticleID::phi) imode=1;
  else return;
  Boost boostv(-part->momentum().boostVector());
  pi0.boost(boostv);
  pip.boost(boostv);
  pim.boost(boostv);
  Lorentz5Momentum ptemp;
  ptemp = pi0+pip;ptemp.rescaleMass();
  *_mplus[imode]+=ptemp.mass()/MeV;
  ptemp = pi0+pim;ptemp.rescaleMass();
  *_mminus[imode]+=ptemp.mass()/MeV;
  ptemp = pip+pim;ptemp.rescaleMass();
  *_m0[imode]+=ptemp.mass()/MeV;
  Energy x = pip.e()-pim.e();
  Energy y = pi0.e()-pi0.m();
  *_xhist[imode]+=x/MeV;
  *_yhist[imode]+=y/MeV;
  if(_xvalue[imode].size()<_nmax) {
    _xvalue[imode].push_back(x);
    _yvalue[imode].push_back(y);
  }
}

NoPIOClassDescription<OmegaPhi3PionAnalysis> OmegaPhi3PionAnalysis::initOmegaPhi3PionAnalysis;
// Definition of the static class description member.

void OmegaPhi3PionAnalysis::Init() {

  static ClassDocumentation<OmegaPhi3PionAnalysis> documentation
    ("There is no documentation for the OmegaPhi3PionAnalysis class");

  static Parameter<OmegaPhi3PionAnalysis,unsigned int> interfaceNMax
    ("MaxPoints",
     "Maximum number of points for the Dalitz plots",
     &OmegaPhi3PionAnalysis::_nmax, 50000, 100, 1000000,
     false, false, Interface::limited);

}

void OmegaPhi3PionAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + 
    string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _xhist[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "x distribution in WRP2+3P2-3P203",
			   "                  GWGX XGX XGX X",
			   "1/SdS/dx/GeV2-13",
			   "  G G       X  X",
			   "x/MeV",
			   "     ");
  _yhist[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "y distribution in WRP2+3P2-3P203",
			   "                  GWGX XGX XGX X",
			   "1/SdS/dy/GeV2-13",
			   "  G G       X  X",
			   "y/MeV",
			   "     ");
  _mplus[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "R2+3 mass in WRP2+3P2-3P203",
			   "GX X         GWGX XGX XGX X",
			   "1/SdS/dm0R2+31/GeV2-13",
			   "  G G   XGX XX    X  X",
			   "m0R2+31/MeV",
			   " XGX XX    ");
  _mminus[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "R2-3 mass in WRP2+3P2-3P203",
			    "GX X         GWGX XGX XGX X",
			    "1/SdS/dm0R2-31/GeV2-13",
			    "  G G   XGX XX    X  X",
			    "m0R2-31/MeV",
			    " XGX XX    ");
  _m0[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "R203 mass in WRP2+3P2-3P203",
			    "GX X         GWGX XGX XGX X",
			    "1/SdS/dm0R2031/GeV2-13",
			    "  G G   XGX XX    X  X",
			    "m0R2031/MeV",
			    " XGX XX    ");
   output << "new frame\n";
   output << "set font duplex\n";
   output << "set limits x -250 250 y 0 250\n";
   output << "set order x y \n";
   output << "title top \"Dalitz plot for W\"\n";
   output << "case      \"                G\"\n"; 
   for(unsigned int ix=0;ix<_xvalue[0].size();++ix) {
     output << ounit(_xvalue[0][ix],MeV) << "   " << ounit(_yvalue[0][ix],MeV) << "\n";
   }
   output << "plot\n";
  _xhist[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "x distribution in FRP2+3P2-3P203",
			   "                  GWGX XGX XGX X",
			   "1/SdS/dx/GeV2-13",
			   "  G G       X  X",
			   "x/MeV",
			   "     ");
  _yhist[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "y distribution in FRP2+3P2-3P203",
			   "                  GWGX XGX XGX X",
			   "1/SdS/dy/GeV2-13",
			   "  G G       X  X",
			   "y/MeV",
			   "     ");
  _mplus[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "R2+3 mass in FRP2+3P2-3P203",
			   "GX X         GWGX XGX XGX X",
			   "1/SdS/dm0R2+31/GeV2-13",
			   "  G G   XGX XX    X  X",
			   "m0R2+31/MeV",
			   " XGX XX    ");
  _mminus[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "R2-3 mass in FRP2+3P2-3P203",
			    "GX X         GWGX XGX XGX X",
			    "1/SdS/dm0R2-31/GeV2-13",
			    "  G G   XGX XX    X  X",
			    "m0R2-31/MeV",
			    " XGX XX    ");
  _m0[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "R203 mass in FRP2+3P2-3P203",
			    "GX X         GWGX XGX XGX X",
			    "1/SdS/dm0R2031/GeV2-13",
			    "  G G   XGX XX    X  X",
			    "m0R2031/MeV",
			    " XGX XX    ");
   output << "new frame\n";
   output << "set font duplex\n";
   output << "set limits x -400 400 y 0 400\n";
   output << "set order x y \n";
   output << "title top \"Dalitz plot for F\"\n";
   output << "case      \"                G\"\n"; 
   for(unsigned int ix=0;ix<_xvalue[1].size();++ix) {
     output << _xvalue[1][ix]/MeV << "   " << _yvalue[1][ix]/MeV << "\n";
   }
   output << "plot\n";
}

void OmegaPhi3PionAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  for(unsigned int ix=0;ix<2;++ix) {
    _xhist  .push_back(new_ptr(Histogram(-400.,400.  ,200)));
    _yhist  .push_back(new_ptr(Histogram(0.   ,400.  ,200)));
    _mplus  .push_back(new_ptr(Histogram(0.   ,1000.,200)));
    _mminus .push_back(new_ptr(Histogram(0.   ,1000.,200)));
    _m0     .push_back(new_ptr(Histogram(0.   ,1000.,200)));
  }
  _xvalue.resize(2);
  _yvalue.resize(2);
}
