// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauTo2MesonAnalysis class.
//

#include "TauTo2MesonAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

inline void TauTo2MesonAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  // BELLE data for the pipi distribution
  double valsBELLE[]={0.076,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,
		      0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,
		      1.050,1.100,1.150,1.200,1.250,1.300,1.350,1.400,1.450,1.500,
		      1.550,1.600,1.650,1.700,1.750,1.800,1.850,1.900,1.950,2.000,
		      2.050,2.100,2.150,2.200,2.250,2.300,2.350,2.400,2.450,2.500,
		      2.550,2.600,2.650,2.700,2.750,2.800,2.850,2.900,2.950,3.000,
		      3.050,3.100,3.150};
  double dvalBELLE[]={  10.3,  64.5, 123.2, 191.2, 284.3,
			443.5, 690.8,1113.9,1781.8,2716.5,
			3260.2,2765.4,1879.6,1222.8, 816.7,
			590.1, 424.4, 324.9, 259.1, 202.0,
			175.7, 147.9, 122.6, 110.6,  97.3,
			82.1,  74.0,  63.2,  55.7,  49.3,
			43.2,331.93,298.43,263.91,198.85,
			176.86,145.50,110.07, 81.40, 64.57,
			43.10, 27.46, 20.27, 14.58,  8.89,
			3.77,  0.31, -0.69,  0.07,  1.52,
			2.62,  3.64,  4.17,  3.85,  3.16,
			3.44,  1.96,  1.86,  1.55,  0.89,
			0.32,  0.04};
  double errorstatBELLE[]={ 5.0 , 3.0 , 3.0 , 3.4 , 4.3 , 5.9 , 8.6 ,13.3 ,20.8 ,31.1 ,
			    37.2 ,31.7 ,21.9 ,14.6 ,10.1 , 7.7 , 5.9 , 4.9 , 4.2 , 3.6 ,
			    3.3 , 3.0 , 2.7 , 2.5 , 2.3 , 2.1 , 1.9 , 1.7 , 1.5 , 1.4 ,
			    1.3 ,11.43, 9.26, 8.89, 8.09, 6.14, 5.35, 4.61, 3.76, 3.45,
			    2.76, 2.53, 2.17, 2.03, 1.82, 1.58, 1.42, 1.34, 1.26, 1.05,
			    0.97, 1.01, 2.95, 1.12, 1.32, 0.99, 0.60, 0.61, 0.43, 0.42,
			    0.11, 0.01};
  double errorsystBELLE[]={ 4.1 , 4.6 , 3.8 , 3.8 , 3.6 , 3.8 , 3.8 , 4.0 , 4.6 , 4.1 ,
		       3.8 , 3.1 , 2.3 , 1.8 , 1.3 , 1.1 , 1.0 , 0.8 , 0.9 , 0.6 ,
		       0.8 , 0.5 , 0.7 , 0.4 , 0.5 , 0.4 , 0.4 , 0.4 , 0.3 , 0.4 ,
		       0.3 , 2.33, 2.12, 2.54, 3.16, 2.85, 2.74, 2.11, 1.76, 1.78,
		       1.87, 1.93, 2.11, 2.17, 2.21, 1.98, 1.45, 1.18, 1.27, 1.41,
		       1.55, 1.81, 1.89, 1.66, 1.30, 1.37, 0.76, 0.71, 0.57, 0.32,
		       0.12, 0.01};
  double evalBELLE[62];
  for(unsigned int ix=0;ix<62;++ix) {
    evalBELLE[ix]=sqrt(sqr(errorstatBELLE[ix])+sqr(errorsystBELLE[ix]));
    if(ix<31) {
      dvalBELLE[ix]*=1e-3;
      evalBELLE[ix]*=1e-3;
    }
    else {
      dvalBELLE[ix]*=1e-4;
      evalBELLE[ix]*=1e-4;
    }
  }
  vector<double> bins(valsBELLE ,valsBELLE +63);
  vector<double> data  = vector<double>(dvalBELLE ,dvalBELLE +62);
  vector<double> error = vector<double>(evalBELLE ,evalBELLE +62);
  _m2pipiBELLE= new_ptr(Histogram(bins,data,error));
  double valsCLEO[]={0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,
		     0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,0.750,
		     0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,1.000,
		     1.050,1.100,1.150,1.200,1.250,1.300,1.350,1.400,1.450,1.500,
		     1.550,1.600,1.650,1.700};
  double dvalCLEO[]={   1.5 ,   8.0 ,   6.0 ,   8.5 ,   15.6 ,
		       16.2 ,  24.9 ,  41.4 ,  50.6 ,   60.9 ,
		       79.8 , 107.4 , 144.3 , 204.5 ,  269.1 ,
		      385.8 , 571.5 , 826.8 ,1078.4 , 1228.1 ,
		     1114.7 , 878.1 , 629.3 , 446.7 ,  326.2 ,
		      262.1 , 207.3 , 158.8 , 129.6 ,  202.8 ,
		      151.0 , 111.0 ,  87.0 ,  63.9 ,   42.7 ,
		       29.2 ,  18.1 ,   6.98,   2.91,    0.71,
		       0.59,    0.68,   0.28};
  double evalCLEO[]={ 1.4 , 2.5 , 2.6 , 2.3 , 2.6 , 2.9 , 3.1 , 3.4 , 3.7 , 4.0 ,
		      4.4 , 4.7 , 5.2 , 5.9 , 6.5 , 7.5 , 8.7 ,10.1 ,11.3 ,11.8 ,
		     11.0 , 9.6 , 7.9 , 6.6 , 5.5 , 4.9 , 4.3 , 3.7 , 3.4 , 4.8 ,
		      4.1 , 3.4 , 3.0 , 2.5 , 2.0 , 1.7 , 1.3 , 0.84, 0.59, 0.32,
		      0.25, 0.26, 0.21};
  for(unsigned int ix=0;ix<43;++ix) {
    if(ix<29) {
      dvalCLEO[ix]*=1e-4/0.025;
      evalCLEO[ix]*=1e-4/0.025;
    }
    else {
      dvalCLEO[ix]*=1e-4/0.05;
      evalCLEO[ix]*=1e-4/0.05;
    }
  }
  bins  = vector<double>(valsCLEO ,valsCLEO +44);
  data  = vector<double>(dvalCLEO ,dvalCLEO +43);
  error = vector<double>(evalCLEO ,evalCLEO +43);
  _mpipiCLEO= new_ptr(Histogram(bins,data,error));
  _m2KpiA=new_ptr(Histogram(0.,3.15,200));
  _mKpiA =new_ptr(Histogram(0.,1.8,200));
  _m2KpiB=new_ptr(Histogram(0.,3.15,200));
  _mKpiB =new_ptr(Histogram(0.,1.8,200));
  _m2KpiC=new_ptr(Histogram(0.,3.15,200));
  _m2KpiD=new_ptr(Histogram(0.,3.15,200));
  _m2Keta=new_ptr(Histogram(0.,3.15,200));
  _mKeta =new_ptr(Histogram(0.,1.8,200));
  _m2KK  =new_ptr(Histogram(0.,3.15,200));
  _mKK   =new_ptr(Histogram(0.,1.8,200));
}

void TauTo2MesonAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  tPVector hadrons=event->getFinalState();
  map<tPPtr,ParticleVector> taus;
  for(unsigned int ix=0;ix<hadrons.size();++ix) {
    PPtr mother=hadrons[ix];
    do {
      if(!mother->parents().empty()) mother=mother->parents()[0];
      else                           mother=tPPtr();
    }
    while(mother&&abs(mother->id())!=ParticleID::tauminus);
    if(mother&&abs(mother->id())==ParticleID::tauminus) {
      if(taus.find(mother)==taus.end()) {
	taus.insert(make_pair(mother,ParticleVector()));
      }
      taus[mother].push_back(hadrons[ix]);
    }
  }
  map<tPPtr,ParticleVector>::const_iterator tit;
  for(tit=taus.begin();tit!=taus.end();++tit) {
    if(tit->second.size()!=3) continue;
    ParticleVector decay=tit->second;
    using Constants::pi;
    unsigned int npi0(0),npip(0),nkp(0),nk0(0),neta(0);
    Lorentz5Momentum pdecay;
    for(unsigned int ix=0;ix<decay.size();++ix) {
      long id = abs(decay[ix]->id());
      if(id!=ParticleID::nu_tau) pdecay+=decay[ix]->momentum();
      if(id==ParticleID::pi0)         ++npi0;
      else if(id==ParticleID::piplus) ++npip;
      else if(id==ParticleID::Kplus)  ++nkp;
      else if(id==ParticleID::K0)     ++nk0;
      else if(id==ParticleID::eta)    ++neta;
    }
    double mass2 = pdecay.m2()/GeV2;
    double mass = sqrt(mass2);
    if(npi0==1&&npip==1) {
      *_m2pipiBELLE+=mass2;
      *_mpipiCLEO+=mass;
    }
    else if(nkp==1&&npi0==1) {
      *_m2KpiA+=mass2;
      *_mKpiA+=mass;
      double mtau=getParticleData(ParticleID::tauminus)->mass()/GeV;
      double m1=getParticleData(ParticleID::Kplus)->mass()/GeV;
      double m2=getParticleData(ParticleID::pi0)->mass()/GeV;
      double q=0.5*sqrt((sqr(mass2-sqr(m1)-sqr(m2))-4*sqr(m1*m2))/mass2);
      double fact = sqr(mtau*mtau-mass2)*q/mass2/mass/mtau*
	(2.*mass2+sqr(mtau))/3./sqr(mtau)
	*generator()->standardModel()->CKM(0,1)/pow(2.,8)/pow(pi,3)*sqr(1.16637e-5/GeV2)
	/(2.26501e-11*MeV)*GeV2*GeV2*GeV;
      _m2KpiC->addWeighted(mass2,1./fact);
    }
    else if(nk0==1&&npip==1) {
      *_m2KpiB+=mass2;
      *_mKpiB+=mass;
      double mtau=getParticleData(ParticleID::tauminus)->mass()/GeV;
      double m1=getParticleData(ParticleID::K0)->mass()/GeV;
      double m2=getParticleData(ParticleID::piplus)->mass()/GeV;
      double q=0.5*sqrt((sqr(mass2-sqr(m1)-sqr(m2))-4*sqr(m1*m2))/mass2);
      double fact = sqr(mtau*mtau-mass2)*q/mass2/mass/mtau*
	(2.*mass2+sqr(mtau))/3./sqr(mtau)
	*generator()->standardModel()->CKM(0,1)/pow(2.,8)/pow(pi,3)*sqr(1.16637e-5/GeV2)
	/(2.26501e-11*MeV)*GeV2*GeV2*GeV;
      _m2KpiD->addWeighted(mass2,1./fact);
    }
    else if(neta==1&&nkp==1) {
      *_m2Keta+=mass2;
      *_mKeta+=mass;
    }
    else if(nkp==1&&nk0==1) {
      *_m2KK  +=mass2;
      *_mKK  +=mass;
    }
  }
}

NoPIOClassDescription<TauTo2MesonAnalysis> TauTo2MesonAnalysis::initTauTo2MesonAnalysis;
// Definition of the static class description member.

void TauTo2MesonAnalysis::Init() {

  static ClassDocumentation<TauTo2MesonAnalysis> documentation
    ("The TauTo2MesonAnalysis class plots the mass distributions of tau decays to"
     "two mesons.");

}

void TauTo2MesonAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _m2pipiBELLE->topdrawOutput(output,Frame|Errorbars|Ylog,
			 "RED",
			 "P2+3P203 mass",
			 "GX XGX X     ",
			 "1/SdS/dm2230P2+3P2031/GeV2-23",
			 "  G G   X XXGX XGX XX    X  X",
			 "m2230P2+3P2031/GeV223",
			 " X XXGX XGX XX    X X");
  _mpipiCLEO->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2+3P203 mass",
			    "GX XGX X     ",
			    "1/SdS/dm0P2+3P2031/GeV2-13",
			    "  G G   XGX XGX XX    X  X",
			    "m0P2+3P2031/GeV",
			    " XGX XGX XX    ");
  _m2KpiA->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"K2+3P203 mass",
			" X XGX X     ",
			"1/SdS/dm2230K2+3P2031/GeV2-23",
			"  G G   X XX X XGX XX    X  X",
			"m2230K2+3P2031/GeV223",
			" X XX X XGX XX    X X");
  _mKpiA->topdrawOutput(output,Frame|Errorbars|Ylog,
			 "RED",
			 "K2+3P203 mass",
			 " X XGX X     ",
			 "1/SdS/dm0K2+3P2031/GeV2-13",
			 "  G G   X X XGX XX    X  X",
			 "m0K2+3P2031/GeV",
			 " X X XGX XX    ");
  _m2KpiC->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"K2+3P203 mass",
			" X XGX X     ",
			"1/SdS/dm2230K2+3P2031/GeV2-23",
			"  G G   X XX X XGX XX    X  X",
			"m2230K2+3P2031/GeV223",
			" X XX X XGX XX    X X");
  _m2KpiB->topdrawOutput(output,Frame|Errorbars|Ylog,
			 "RED",
			 "P2+3K203 mass",
			 "GX X X X     ",
			 "1/SdS/dm2230P2+3K2031/GeV2-23",
			 "  G G   X XXGX X X XX    X  X",
			 "m2230P2+3K2031/GeV223",
			 " X XXGX X X XX    X X");
  _mKpiB->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"P2+3K203 mass",
			"GX X X X     ",
			"1/SdS/dm0P2+3K2031/GeV2-13",
			"  G G   XGX X X XX    X  X",
			"m0P2+3K2031/GeV",
			" XGX X X XX    ");
  _m2KpiD->topdrawOutput(output,Frame|Errorbars|Ylog,
			 "RED",
			 "P2+3K203 mass",
			 "GX X X X     ",
			 "1/SdS/dm2230P2+3K2031/GeV2-23",
			 "  G G   X XXGX X X XX    X  X",
			 "m2230P2+3K2031/GeV223",
			 " X XXGX X X XX    X X");
  _m2Keta->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"K2+3H mass",
			" X XG     ",
			"1/SdS/dm2230K2+3H1/GeV2-23",
			"  G G   X XX X XGX    X  X",
			"m2230K2+3H1/GeV223",
			" X XX X XGX    X X");
  _mKeta->topdrawOutput(output,Frame|Errorbars|Ylog,
			 "RED",
			 "K2+3H mass",
			 " X XG     ",
			 "1/SdS/dm0K2+3H1/GeV2-13",
			 "  G G   X X XGX    X  X",
			 "m0K2+3H1/GeV",
			 " X X XGX    ");
  _m2KK->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"K2+3K203 mass",
			" X X X X     ",
			"1/SdS/dm2230K2+3K2031/GeV2-23",
			"  G G   X XX X X X XX    X  X",
			"m2230K2+3K2031/GeV223",
			" X XX X X X XX    X X");
  _mKK->topdrawOutput(output,Frame|Errorbars|Ylog,
			 "RED",
			 "K2+3K203 mass",
			 " X X X X     ",
			 "1/SdS/dm0K2+3K2031/GeV2-13",
			 "  G G   X X X X XX    X  X",
			 "m0K2+3K2031/GeV",
			 " X X X X XX    ");
}


