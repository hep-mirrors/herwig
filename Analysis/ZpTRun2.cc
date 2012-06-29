// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZpTRun2 class.
//

#include "ZpTRun2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

inline IBPtr ZpTRun2::clone() const {
  return new_ptr(*this);
}

inline IBPtr ZpTRun2::fullclone() const {
  return new_ptr(*this);
}

void ZpTRun2::analyze(tEventPtr event, long , int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  tParticleVector outgoing = event->primaryCollision()->step(1)->getFinalState();
  ParticleVector Zdecay;
  // search for electrons and positions from Z/gamma
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    int id=outgoing[ix]->id();
    if(abs(id)==ParticleID::eminus) {
      PPtr part=outgoing[ix];
      do {
	part=part->parents()[0];
      }
      while (part->id()==id);
      if(part->id()==ParticleID::gamma||part->id()==ParticleID::Z0) 
	Zdecay.push_back(outgoing[ix]);
    }
  }
  if(Zdecay.size()!=2) return;
  Lorentz5Momentum pZ=Zdecay[0]->momentum()+Zdecay[1]->momentum();
  pZ.rescaleMass();
  if(pZ.mass()<40.*GeV||pZ.mass()>200.*GeV) return;
  Energy pT = pZ.perp();
  if(pT>260.*GeV) return;
  *_pt += pT/GeV;
}

NoPIOClassDescription<ZpTRun2> ZpTRun2::initZpTRun2;
// Definition of the static class description member.

void ZpTRun2::Init() {

  static ClassDocumentation<ZpTRun2> documentation
    ("There is no documentation for the ZpTRun2 class");

}

void ZpTRun2::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _pt->topdrawOutput(output,Frame|Errorbars|Ylog,"RED",
		     "pt of Z (mass 40 GeV to 200 GeV) compared to D0 run II data", 
		     "                           ",
		     "1/SdS/dpT",
		     "  G G       ",
		     "pT",
		     "   ");
}

void ZpTRun2::doinitrun() {
  AnalysisHandler::doinitrun();
  //data from hep-ex 07120803
  //normalised to one
  double abins[] ={ 0.0 , 2.5, 5.0, 7.5,10.0,
		    12.5,15.0,17.5,20.0,22.5,
		    25.0,27.5,30.0,40.0,50.0,
		    60.0,70.0,80.0,90.0,100.0,
		    140.0,180.0,220.0,260.0};
  
  double adata[] = { 0.0533419,	0.0810158,	0.0634688,	0.0444179,
		     0.0314839,	0.0246659,	0.01865,	0.014238,
		     0.010929,	0.00942508,	0.00691838,	0.00551469,
		     0.00391039,	0.00210559,	0.0011029,	
		     0.000731948,	0.000421119,	0.000250669,	
		     0.00016043,	6.01598e-05,	1.1029e-05,	
		     3.00799e-06,	7.11888e-07 };

  double aerror[] = {0.00273679,	0.00225319,	0.0017852,	
		     0.0014251, 0.0011344,	0.000924408,	0.000783108,	0.000708988,
		     0.000501329,	0.000448409,	0.000361519,	0.000317069,
		     0.0001418,	9.24408e-05,	5.84649e-05,	4.48409e-05,
		     3.61519e-05,	2.24199e-05,	1.8838e-05,	5.84649e-06,
		     2.12459e-06,	1.0468e-06,	6.14578e-07 };

  vector<double>  bins(abins ,abins +24);
  vector<double>  data(adata ,adata +23);
  vector<double> error(aerror,aerror+23);
  _pt = new_ptr(Histogram(bins,data,error));
}
