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

void ZpTRun2::analyze(tEventPtr event, long ieve, int loop, int state) {
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
  if(pT>280.*GeV) return;
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
  _pt->topdrawOutput(output,Frame|Ylog,"RED",
		     "pt of Z (mass 40 GeV to 200 GeV) compared to D0 run II data", 
		     "                           ",
		     "1/SdS/dpT",
		     "  G G       ",
		     "pT",
		     "   ");
}

void ZpTRun2::doinitrun() {
  AnalysisHandler::doinitrun();
  double abins[] ={ 0.0 , 2.5, 5.0, 7.5,10.0,
		    12.5,15.0,17.5,20.0,22.5,
		    25.0,27.5,30.0,40.0,50.0,
		    60.0,70.0,80.0,90.0,100.0,
		    140.0,180.0,220.0,260.0};
  double adata[] ={0.53342E-01,0.81016E-01,0.63469E-01,0.44418E-01,0.31484E-01,
		   0.24666E-01,0.18650E-01,0.14238E-01,0.10929E-01,0.94251E-02,
		   0.69184E-02,0.55147E-02,0.39104E-02,0.21056E-02,0.11029E-02,
		   0.73195E-03,0.42112E-03,0.25067E-03,0.16043E-03,0.60160E-04,
		   0.11029E-04,0.30080E-05,0.71189E-06};
  double aerror[]={0.27368E-02,0.22532E-02,0.17852E-02,0.14251E-02,0.11344E-02,
		   0.92441E-03,0.78311E-03,0.70899E-03,0.50133E-03,0.44841E-03,
		   0.36152E-03,0.31707E-03,0.14180E-03,0.92441E-04,0.58465E-04,
		   0.44841E-04,0.36152E-04,0.22420E-04,0.18838E-04,0.58465E-05,
		   0.21246E-05,0.10468E-05,0.61458E-06};
  vector<double>  bins(abins ,abins +24);
  vector<double>  data(adata ,adata +23);
  vector<double> error(aerror,aerror+23);
  _pt = new_ptr(Histogram(bins,data,error));
}
