// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QuickVBFHadron class.
//

#include "QuickVBFHadron.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"

#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/PDF.h"

using namespace Herwig;

void QuickVBFHadron::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tPVector part = event->getFinalState();
  Lorentz5Momentum pjj;
  Energy2 mHsq(ZERO);
  for(tPVector::const_iterator iter = part.begin(), end = part.end();
      iter!=end;++iter) {
    if((**iter).id()==ParticleID::h0) {
      *_mH     += (**iter).momentum().m()/GeV;
      mHsq=(**iter).momentum().m2();
      *_yH     += (**iter).momentum().rapidity();
      *_phiH   += (**iter).momentum().phi()+Constants::pi;
      *_pTH[0] += (**iter).momentum().perp()/GeV;
      *_pTH[1] += (**iter).momentum().perp()/GeV;
    }
    else if((**iter).id()!=82) {
      *_yjet     += (**iter).momentum().rapidity();
      *_phijet   += (**iter).momentum().phi()+Constants::pi;
      *_pTjet[0] += (**iter).momentum().perp()/GeV;
      *_pTjet[1] += (**iter).momentum().perp()/GeV;
      pjj+=(**iter).momentum();
    }
  }
  *_mjj += pjj.m()/GeV;

  if (!_doOnce) {
    // Find the pdf information
    tSubProPtr sub = event->primarySubProcess();
    // get the event handler
    tcEHPtr eh = dynamic_ptr_cast<tcEHPtr>(event->handler());
    // get the pdfs
    pair<PDF,PDF> pdfs;
    pdfs.first  =  eh->pdf<PDF>(sub->incoming().first );
    pdfs.second =  eh->pdf<PDF>(sub->incoming().second);
    // Make a semi-analytical prediction of the convolution over pdfs
    const double num=1e7;
    double allsum[14][14];
    double allsqr[14][14];
    for (int j=0; j<=13 ; ++j) {
      for (int k=0; k<=13 ; ++k) {
	allsum[j][k]=0;
	allsqr[j][k]=0;
      }
    }
    double zzsum=0;
    double zzsqr=0;
    double wwsum=0;
    double wwsqr=0;
    Energy2 scale=mHsq;
    Energy2 s=event->primaryCollision()->m2();
    cerr << "Starting semi-analytical calculation.  Higgs mass/GeV=" <<
      sqrt(mHsq)/GeV << ", sqrt(s)/GeV=" << sqrt(s)/GeV << "\n";
    double tmin=mHsq/s;
    double tlmin=log(tmin);
    for (int i=0 ; i<int(num+0.5) ; ++i) {
      double xxmin=exp(tlmin*UseRandom::rnd());
      double xlmin=log(xxmin);
      double x1=exp(xlmin*UseRandom::rnd());
      double x2=xxmin/x1;
      // out of laziness, just don't use the 0'th elements, to be like fortran
      double disf[14][3];
      for (int j=1 ; j<=13; ++j) {
	for (int k=1 ; k<=2 ; ++k) {
	  disf[j][k]=0;
	}
      }
      for (int j=1 ; j<=6 ; ++j) {
	disf[j][1]=pdfs.first .xfx((tcPDPtr)getParticleData(j),scale,x1);
	disf[j+6][1]=pdfs.first .xfx((tcPDPtr)getParticleData(-j),scale,x1);
      }
      disf[13][1]=pdfs.first .xfx((tcPDPtr)getParticleData(21),scale,x1);
      for (int j=1 ; j<=6 ; ++j) {
	disf[j][2]=pdfs.second.xfx((tcPDPtr)getParticleData(j),scale,x2);
	disf[j+6][2]=pdfs.second.xfx((tcPDPtr)getParticleData(-j),scale,x2);
      }
      disf[13][2]=pdfs.second.xfx((tcPDPtr)getParticleData(21),scale,x2);
      double w=tlmin*xlmin*
	(disf[1][1]+disf[2][1]+disf[3][1]+disf[4][1]+
	 disf[7][1]+disf[8][1]+disf[9][1]+disf[10][1])*
	(disf[1][2]+disf[2][2]+disf[3][2]+disf[4][2]+
	 disf[7][2]+disf[8][2]+disf[9][2]+disf[10][2]);
      zzsum+=w;
      zzsqr+=sqr(w);
      w=tlmin*xlmin*
	((disf[1][1]+disf[3][1]+disf[8][1]+disf[10][1])*
	 (disf[2][2]+disf[4][2]+disf[7][2]+disf[9][2])+
	 (disf[2][1]+disf[4][1]+disf[7][1]+disf[9][1])*
	 (disf[1][2]+disf[3][2]+disf[8][2]+disf[10][2]));
      // Extra factor because Hw++ test program sums over 4 Cabibbo possibilities
      w*=4;
      wwsum+=w;
      wwsqr+=sqr(w);
      for (int j=1; j<=13; ++j) {
	for (int k=1; k<=13; ++k) {
	  w=tlmin*xlmin*disf[j][1]*disf[k][2];
	  allsum[j][k]+=w;
	  allsqr[j][k]+=sqr(w);
	}
      }
    }
    cerr << "Analytical result for ZZ=" << zzsum/num << "+-" << sqrt(zzsqr-sqr(zzsum)/num)/num << "\n";
    cerr << "Analytical result for WW=" << wwsum/num << "+-" << sqrt(wwsqr-sqr(wwsum)/num)/num << "\n";
    for (int j=1; j<=13; ++j) {
      int idj=0;
      if (j>=1&&j<=5) idj=j;
      if (j>=7&&j<=11)idj=6-j;
      if (idj!=0) {
	for (int k=1; k<=13; ++k) {
	  int idk=0;
	  if (k>=1&&k<=5) idk=k;
	  if (k>=7&&k<=11)idk=6-k;
	  if (idk!=0 && allsum[j][k]>0) {
	    cerr << getParticleData(idj)->PDGName() << " " <<
	      getParticleData(idk)->PDGName() << "\t" <<
	      allsum[j][k]/num << "+-" << sqrt(allsqr[j][k]-sqr(allsum[j][k])/num)/num << "\n";
	  }
	}
      }
    }
    _doOnce=true;
  }
  /*
  // ids of the partons going into the primary sub process
  tSubProPtr sub = event->primarySubProcess();
  int id1 = sub->incoming().first ->id();
  int id2 = sub->incoming().second->id();
  // get the event handler
  tcEHPtr eh = dynamic_ptr_cast<tcEHPtr>(event->handler());
  // get the values of x
  double x1 = eh->lastX1();
  double x2 = eh->lastX2();
  // get the pdfs
  pair<PDF,PDF> pdfs;
  pdfs.first  =  eh->pdf<PDF>(sub->incoming().first );
  pdfs.second =  eh->pdf<PDF>(sub->incoming().second);
  // get the scale
  Energy2 scale = eh->lastScale();
  // get the values of the pdfs
  double pdf1 = pdfs.first .xfx(sub->incoming().first ->dataPtr(),scale,x1);
  double pdf2 = pdfs.second.xfx(sub->incoming().second->dataPtr(),scale,x2);
  // print them to cerr
  if (x1<3e-4 || x2<3e-4) {
  cerr << "The x values are:" << x1 << " " << x2 << "\n";
  cerr << "The scale is:" << sqrt(scale)/GeV << " GeV \n";
  cerr << "The parton ids are:" << id1 << " " << id2 << "\n";
  cerr << "The pdfs are:" << pdf1 << " " << pdf2 << "\n\n";
  }
  for (int i=0;i<=100;++i) {
    x1 = pow(10.0,4*(double(i)/100-1));
    pdf1 = pdfs.first .xfx(sub->incoming().first ->dataPtr(),scale,x1);
    cerr << x1 << " " << pdf1 << "\n";
  }
  */
}

bool QuickVBFHadron::_doOnce = false;

IBPtr QuickVBFHadron::clone() const {
  return new_ptr(*this);
}

IBPtr QuickVBFHadron::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<QuickVBFHadron> QuickVBFHadron::initQuickVBFHadron;
// Definition of the static class description member.

void QuickVBFHadron::Init() {

  static ClassDocumentation<QuickVBFHadron> documentation
    ("There is no documentation for the QuickVBFHadron class");

}

void QuickVBFHadron::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  string title,species;
  title = "mass of H";
  _mH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "rapidity of H";
  _yH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "pT of H";
  _pTH[0]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  _pTH[1]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  title = "azimuth of H";
  _phiH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "rapidity of jet";
  _yjet->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "pT of jet";
  _pTjet[0]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  _pTjet[1]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  title = "azimuth of jet";
  _phijet->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "mjj";
  _mjj->topdrawOutput(outfile,Frame,"BLACK",title);
}

void QuickVBFHadron::doinitrun() {
  AnalysisHandler::doinitrun();
  if(getParticleData(ParticleID::h0)->mass()>200.*GeV) 
    _mH     = new_ptr(Histogram(200.,            400.,200));
  else
    _mH     = new_ptr(Histogram(114.,            116.0,200));
  _yH       = new_ptr(Histogram( -10.0,              10.0,200));
  _phiH     = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _pTH[0]      = new_ptr(Histogram(  0.0,1000.,1000));
  _pTH[1]      = new_ptr(Histogram(  0.0,1000.,100));
  _yjet     = new_ptr(Histogram( -10.0,              10.0,200));
  _phijet   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _pTjet[0] = new_ptr(Histogram(  0.0,1000.,1000));
  _pTjet[1] = new_ptr(Histogram(  0.0,1000.,100));
  _mjj = new_ptr(Histogram(0.0,200.,100));
}
