// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauCorrelationAnalysis class.
//

#include "TauCorrelationAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include <ThePEG/PDT/EnumParticles.h>

using namespace Herwig;


void TauCorrelationAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _phi->topdrawOutput(output,Frame|Errorbars,
		      "RED",
		      "F2*3 in h203RT2+3T2-3 with TRPN0T1",
		      "GX X     X XWGX XGX X      GWGGXGX",
		      "1/GdG/dF2*3",
		      "  F F  GX X",
		      "F2*3",
		      "GX X");  
  _delta->topdrawOutput(output,Frame|Errorbars,
			"RED",
			"D2*3 in h203RT2+3T2-3 with TRPN0T1",
			"GX X     X XWGX XGX X      GWGGXGX",
			"1/GdG/dD2*3",
			"  F F  GX X",
			"D2*3",
			"GX X");
  _rhoangle1->topdrawOutput(output,Frame|Errorbars,
		      "RED",
		      "F2*3 in h203RT2+3T2-3 with TRRN0T1 and y011y021>0",
		      "GX X     X XWGX XGX X      GWGGXGX      X X X X  ",
		      "1/GdG/dF2*3",
		      "  F F  GX X",
		      "F2*3",
		      "GX X"); 
  _rhoangle2->topdrawOutput(output,Frame|Errorbars,
		      "RED",
		      "F2*3 in h203RT2+3T2-3 with TRRN0T1 and y011y021<0",
		      "GX X     X XWGX XGX X      GWGGXGX      X X X X  ",
		      "1/GdG/dF2*3",
		      "  F F  GX X",
		      "F2*3",
		      "GX X"); 
}

void TauCorrelationAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  using Constants::pi;
  _phi       = new_ptr(Histogram(0.,pi,100));
  _delta     = new_ptr(Histogram(3.,pi,100));
  _rhoangle1 = new_ptr(Histogram(0.,pi,100));
  _rhoangle2 = new_ptr(Histogram(0.,pi,100));
}

void TauCorrelationAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // find all higgs bosons particles 
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
      ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
      ThePEG::ParticleSet::iterator iter=part.begin();
      ThePEG::ParticleSet::iterator end=part.end();
      for( ;iter!=end;++iter) {
	if((**iter).id()==ParticleID::h0) particles.push_back(*iter);
      }
  }
  // analyse them
  analyze(particles);
}

void TauCorrelationAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

namespace {

tPPtr findTau(tPPtr parent) {
  for(unsigned int ix=0;ix<parent->children().size();++ix) {
    if(parent->children()[ix]->id()==parent->id())
      return findTau(parent->children()[ix]);
  }
  return parent;
}

}

void TauCorrelationAnalysis::analyze(tPPtr part) {
  // check the number of children of the particle
  if(part->children().size()!=2) return;
  // and they are tau's
  if(abs(part->children()[0]->id())!=ParticleID::tauminus||
     abs(part->children()[1]->id())!=ParticleID::tauminus) return;
  ParticleVector children;
  for(unsigned int ix=0;ix<part->children().size();++ix)
    children.push_back(findTau(part->children()[ix]));
  // and number of children
  if(!((children[0]->children().size()==2||
	children[0]->children().size()==3)&&
       (children[1]->children().size()==2||
	children[1]->children().size()==3))) return;
  // call rho and pi specific analysis
  analyzePi(part,children);
  analyzeRho(part,children);
}

NoPIOClassDescription<TauCorrelationAnalysis> TauCorrelationAnalysis::initTauCorrelationAnalysis;
// Definition of the static class description member.

void TauCorrelationAnalysis::Init() {

  static ClassDocumentation<TauCorrelationAnalysis> documentation
    ("The TauCorrelationAnalysis class performs the analysis of correlation effects"
     " in the decays of tau's produced in Higgs decay");

}

void TauCorrelationAnalysis::analyzePi(tPPtr part,
				       ParticleVector children) {
  if(children[0]->children().size()!=2||
     children[1]->children().size()!=2) return;
  // now examine the decay products
  tPPtr taup,taum,pim,pip,nu,nub;
  for(unsigned int ix=0;ix<2;++ix) {
    if(children[ix]->id()==ParticleID::tauminus) {
      taum=children[ix];
      for(unsigned int ix=0;ix<2;++ix) {
	if(taum->children()[ix]->id()==ParticleID::piminus) {
	  pim=taum->children()[ix];
	}
	else if(taum->children()[ix]->id()==ParticleID::nu_tau) {
	  nu=taum->children()[ix];
	}
      }
    }
    else {
      taup=children[ix];
      for(unsigned int ix=0;ix<2;++ix) {
	if(taup->children()[ix]->id()==ParticleID::piplus) {
	  pip=taup->children()[ix];
	}
	else if(taup->children()[ix]->id()==ParticleID::nu_taubar) {
	  nub=taup->children()[ix];
	}
      }
    }
  }
  if(!taup||!taum||!pim||!pip||!nu||!nub){return;}
  Boost bv(-part->momentum().boostVector());
  Lorentz5Momentum ptaup(taup->momentum());ptaup.boost(bv);
  Lorentz5Momentum ptaum(taum->momentum());ptaum.boost(bv);
  Lorentz5Momentum ppim( pim->momentum() );ppim.boost(bv);
  Lorentz5Momentum ppip( pip->momentum() );ppip.boost(bv);
  Lorentz5Momentum pnu(  nu->momentum()  );pnu.boost(bv);
  Lorentz5Momentum pnub( nub->momentum() );pnub.boost(bv);
  ThreeVector<Energy2> norm1(ppip.vect().cross(pnub.vect()));
  ThreeVector<Energy2> norm2(ppim.vect().cross(pnu.vect()));
  double phi=norm1.angle(norm2);
  *_phi   +=phi;
  *_delta +=ppip.vect().angle(ppim.vect());
}

void TauCorrelationAnalysis::analyzeRho(tPPtr part,
					ParticleVector children) {
  // now examine the decay products
  tPPtr taup,taum,pim,pip,nu,nub,pi0a,pi0b,rhop,rhom;
  int idtemp;
  for(unsigned int ix=0;ix<2;++ix) {
    if(children[ix]->id()==ParticleID::tauminus) {
      taum=children[ix];
      for(unsigned int iz=0;iz<taum->children().size();++iz) {
	idtemp=taum->children()[iz]->id();
	if(idtemp==-213||idtemp==-100213||idtemp==-30213) {
	  tPPtr rhom=taum->children()[iz];
	  if(rhom->children().size()!=2) return;
	  for(unsigned int iy=0;iy<2;++iy) {
	    int idb=rhom->children()[iy]->id();
	    if(idb==ParticleID::piminus)  pim=rhom->children()[iy];
	    else if(idb==ParticleID::pi0) pi0a=rhom->children()[iy];
	  }
	}
	else if(idtemp==ParticleID::nu_tau) {
	  nu=taum->children()[iz];
	}
	else if(idtemp==ParticleID::piminus)  pim =taum->children()[iz];
	else if(idtemp==ParticleID::pi0)      pi0a=taum->children()[iz];
      }
    }
    else {
      taup=children[ix];
      for(unsigned int iz=0;iz<taup->children().size();++iz) {
	idtemp=taup->children()[iz]->id();
	if(idtemp==213||idtemp==100213||idtemp==30213) {
	  tPPtr rhop=taup->children()[iz];
	  if(rhop->children().size()!=2) return;
	  for(unsigned int iy=0;iy<2;++iy) {
	    int idb=rhop->children()[iy]->id();
	    if(idb==ParticleID::piplus)   pip=rhop->children()[iy];
	    else if(idb==ParticleID::pi0) pi0b=rhop->children()[iy];
	  }
	}
	else if(idtemp==ParticleID::nu_taubar) {
	  nub=taup->children()[iz];
	}
	else if(idtemp==ParticleID::piplus) pip  = taup->children()[iz];
	else if(idtemp==ParticleID::pi0)    pi0b = taup->children()[iz];
      }
    }
  }
  if(!taup||!taum||!pim||!pip||!pi0a||!pi0b||!nu||!nub) return;
  LorentzMomentum prest(pim ->momentum()+pip ->momentum()+
			pi0a->momentum()+pi0b->momentum());
  Boost bv(-prest.boostVector());
  Lorentz5Momentum ppim( pim->momentum() );ppim.boost(bv);
  Lorentz5Momentum ppip( pip->momentum() );ppip.boost(bv);
  Lorentz5Momentum ppi0a( pi0a->momentum() );ppi0a.boost(bv);
  Lorentz5Momentum ppi0b( pi0b->momentum() );ppi0b.boost(bv);
  ThreeVector<Energy2> norm1(ppip.vect().cross(ppi0b.vect()));
  ThreeVector<Energy2> norm2(ppim.vect().cross(ppi0a.vect()));
  double phi=norm1.angle(norm2);

  Lorentz5Momentum ptaup(taup->momentum());
  Lorentz5Momentum ptaum(taum->momentum());
  bv = -ptaum.boostVector();
  ppim =pim->momentum();ppim.boost(bv);
  ppi0a=pi0a->momentum();ppi0a.boost(bv);
  bv = -ptaup.boostVector();
  ppip =pip->momentum();ppip.boost(bv);
  ppi0b=pi0b->momentum();ppi0b.boost(bv);
  double y1=(ppip.e()-ppi0b.e())/(ppip.e()+ppi0b.e());
  double y2=(ppim.e()-ppi0a.e())/(ppim.e()+ppi0a.e());
  if(y1*y2>0) *_rhoangle1+=phi;
  else        *_rhoangle2+=phi;
}
