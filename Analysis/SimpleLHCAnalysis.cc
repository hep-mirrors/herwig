// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleLHCAnalysis class.
//

#include "SimpleLHCAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void SimpleLHCAnalysis::analyze(tEventPtr event, long, int, int) {
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z
  Lorentz5Momentum pz;
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator stest =event->primaryCollision()->steps().end();
  StepVector::const_iterator send=sit;
  ++send;
  if(send==stest) --send;
  ++send;
  if(send==stest) --send;
  ++send;
  for(;sit!=send;++sit) {
    ParticleSet part=(**sit).all();
    ParticleSet::const_iterator iter = part.begin(), end = part.end();
    for( ;iter!=end;++iter) {
      if((**iter).children().size()!=2) continue;
      if((**iter).id()==ParticleID::Z0||(**iter).id()==ParticleID::gamma) {
	pz=(*iter)->momentum();
	double pt = pz.perp()/GeV;
	double mz = pz.mass()/GeV;
	if(mz>20.&&mz<80.)        (_ptZ[1])+=(pt);
	else if (mz>80.&&mz<100.) (_ptZ[2])+=(pt);
	else if (mz>100.)         (_ptZ[3])+=(pt);
	(_ptZ[0])+=(pt);
	(_mZ)+=(mz);
	(_rapZ)+=pz.rapidity();
	(_phiZ)+=pz.phi();
      } 
      else if ((**iter).id()==ParticleID::Wplus) {
	pz=(*iter)->momentum();
	double pt = pz.perp()/GeV;
	double mz = pz.mass()/GeV;
	if(mz>20.&&mz<80.)        (_ptWp[1])+=(pt);
	else if (mz>80.&&mz<100.) (_ptWp[2])+=(pt);
	else if (mz>100.)         (_ptWp[3])+=(pt);
	(_ptWp[0])+=(pt);
	(_mWp)+=(mz);
	(_rapWp)+=pz.rapidity();
	(_phiWp)+=pz.phi();
      } 
      else if ((**iter).id()==ParticleID::Wminus) {
	pz=(*iter)->momentum();
	double pt = pz.perp()/GeV;
	double mz = pz.mass()/GeV;
	if(mz>20.&&mz<80.)        (_ptWm[1])+=(pt);
	else if (mz>80.&&mz<100.) (_ptWm[2])+=(pt);
	else if (mz>100.)         (_ptWm[3])+=(pt);
	(_ptWm[0])+=(pt);
	(_mWm)+=(mz);
	(_rapWm)+=pz.rapidity();
	(_phiWm)+=pz.phi();
      }
    }
  }
}

NoPIOClassDescription<SimpleLHCAnalysis> SimpleLHCAnalysis::initSimpleLHCAnalysis;
// Definition of the static class description member.

void SimpleLHCAnalysis::Init() {

  static ClassDocumentation<SimpleLHCAnalysis> documentation
    ("The SimpleLHCAnalysis class performs a simple analysis of W and"
     " Z production in hadron-hadron collisions");

}

void SimpleLHCAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  string title;
  using namespace HistogramOptions;
  for(unsigned int ix=0;ix<4;++ix) {
    if(ix==0){title="pt of Z for all masses ";}
    else if(ix==1){title="pt of Z for mass 40-80 GeV";}
    else if(ix==2){title="pt of Z for mass 80-100 GeV";}
    else if(ix==3){title="pt of Z for mass 100- GeV";}
    _ptZ[ix].topdrawOutput(outfile,Frame,"BLACK",title);
    _ptZ[ix].topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
    if(ix==0){title="pt of Wp for all masses ";}
    else if(ix==1){title="pt of Wp for mass 40-80 GeV";}
    else if(ix==2){title="pt of Wp for mass 80-100 GeV";}
    else if(ix==3){title="pt of Wp for mass 100- GeV";}
    _ptWp[ix].topdrawOutput(outfile,Frame,"BLACK",title);
    _ptWp[ix].topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
    if(ix==0){title="pt of Wm for all masses ";}
    else if(ix==1){title="pt of Wm for mass 40-80 GeV";}
    else if(ix==2){title="pt of Wm for mass 80-100 GeV";}
    else if(ix==3){title="pt of Wm for mass 100- GeV";}
    _ptWm[ix].topdrawOutput(outfile,Frame,"BLACK",title);
    _ptWm[ix].topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  }
  _mZ.topdrawOutput(outfile,Frame,"BLACK","Mass of Z");
  _mZ.topdrawOutput(outfile,Frame|Ylog,"BLACK", "Mass of Z");
  _mWp.topdrawOutput(outfile,Frame,"BLACK","Mass of Wp");
  _mWp.topdrawOutput(outfile,Frame|Ylog,"BLACK", "Mass of Wp");
  _mWm.topdrawOutput(outfile,Frame,"BLACK","Mass of Wm");
  _mWm.topdrawOutput(outfile,Frame|Ylog,"BLACK", "Mass of Wm");
  _rapZ.topdrawOutput(outfile,Frame,"BLACK","Rapidity of Z");
  _rapZ.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of Z");
  _rapWp.topdrawOutput(outfile,Frame,"BLACK","Rapidity of Wp");
  _rapWp.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of Wp");
  _rapWm.topdrawOutput(outfile,Frame,"BLACK","Rapidity of Wm");
  _rapWm.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of Wm");

  _phiZ.topdrawOutput(outfile,Frame,"BLACK","Azimuth of Z");
  _phiWp.topdrawOutput(outfile,Frame,"BLACK","Azimuth of Wp");
  _phiWm.topdrawOutput(outfile,Frame,"BLACK","Azimuth of Wm");
}

