// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DDalitzAnalysis class.
//

#include "DDalitzAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

  
void DDalitzAnalysis::
findChildren(tPPtr part,ParticleVector & prod) {
  if(part->children().empty()) {
    prod.push_back(part);
  }
  else {
    for(unsigned ix=0;ix<part->children().size();++ix) {
      findChildren(part->children()[ix],prod);
    }
  }
}

void DDalitzAnalysis::analyze(tEventPtr event, long , int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find all D0/Dbar0 and D+/-
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      if(abs((**iter).id())==ParticleID::D0||
	 abs((**iter).id())==ParticleID::Dplus||
	 abs((**iter).id())==ParticleID::D_splus) {
	particles.push_back(*iter);
      }
    }
  }
  // analyse them
  analyze(particles);
}

LorentzRotation DDalitzAnalysis::transform(tEventPtr ) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void DDalitzAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void DDalitzAnalysis::analyze(tPPtr part) {
  // find the stable decay products
  ParticleVector products;
  findChildren(part,products);
  vector<Lorentz5Momentum> ppim,ppip,pk0,pkp,pkm,ppi0;
  if(products.size()!=3) return;
  for(unsigned int ix=0;ix<products.size();++ix) {
    int id0=products[ix]->id();
    if(id0==ParticleID::piplus)           ppip.push_back(products[ix]->momentum());
    else if(id0==ParticleID::piminus)     ppim.push_back(products[ix]->momentum());
    else if(id0==ParticleID::pi0)         ppi0.push_back(products[ix]->momentum());
    else if(abs(id0)==ParticleID::K0)     pk0 .push_back(products[ix]->momentum());
    else if(id0==ParticleID::K_L0)        pk0 .push_back(products[ix]->momentum());
    else if(id0==ParticleID::K_S0)        pk0 .push_back(products[ix]->momentum());
    else if(id0==ParticleID::Kplus)  pkp.push_back(products[ix]->momentum());
    else if(id0==ParticleID::Kminus) pkm.push_back(products[ix]->momentum());
  }
  Energy2 mplus,mminus,mpipi,mneut;
  if(ppim.size()==1&&ppip.size()==1&&pk0.size()==1) {
    if(part->id()==ParticleID::D0) {
      mminus = (ppim.back()+pk0.back() ).m2();
      mplus  = (ppip.back()+pk0.back() ).m2();
      mpipi  = (ppip.back()+ppim.back()).m2();
    }
    else {
      mminus = (ppip.back()+pk0.back()).m2();
      mplus  = (ppim.back()+pk0.back()).m2();
      mpipi  = (ppip.back()+ppim.back()).m2();
    }
    *_m2plus1  +=mplus/GeV2;
    *_m2minus1 +=mminus/GeV2;
    *_m2pipi1  +=mpipi/GeV2;
    if(_points1.size()<200000) {
      _points1.push_back(make_pair(mplus,mminus));
    }
  }
  else if ((ppim.size()==1&&pkp.size()==1&&ppi0.size()==1)||
	   (ppip.size()==1&&pkm.size()==1&&ppi0.size()==1)) {
    if(part->id()==ParticleID::D0) {
      mneut =(pkm .back()+ppip.back()).m2();
      mminus=(pkm .back()+ppi0.back()).m2();
      mpipi =(ppip.back()+ppi0.back()).m2();
    }
    else {
      mneut =(pkp .back()+ppim.back()).m2();
      mminus=(pkp .back()+ppi0.back()).m2();
      mpipi =(ppim.back()+ppi0.back()).m2();
    }
    *_m2neutral2  +=mneut/GeV2;
    *_m2minus2 +=mminus/GeV2;
    *_m2pipi2  +=mpipi/GeV2;
    if(_points2.size()<200000) {
      _points2.push_back(make_pair(mminus,mneut));
    }
  }
  else if ((ppip.size()==2&&pkm.size()==1)||
	   (ppim.size()==2&&pkp.size()==1)) {
    if(part->id()==ParticleID::Dplus) {
      mplus  = (pkm[0] +ppip[0]).m2();
      mminus = (pkm[0] +ppip[1]).m2();
      mpipi  = (ppip[0]+ppip[1]).m2();
    }
    else {
      mplus  = (pkp[0] +ppim[0]).m2();
      mminus = (pkp[0] +ppim[1]).m2();
      mpipi  = (ppim[0]+ppim[1]).m2();
    }
    if(mplus<mminus) swap(mplus,mminus);
    *_mKpilow3 += mminus/GeV2;
    *_mKpihigh3+= mplus/GeV2;
    *_mKpiall3 += mminus/GeV2;
    *_mKpiall3 += mplus/GeV2;
    *_mpipi3   += mpipi/GeV2;
    if(_points3.size()<200000) {
      _points3.push_back(make_pair(mminus,mpipi));
    }
  }
  else if ((ppip.size()==1&&ppi0.size()==1&&pk0.size()==1)||
	   (ppim.size()==1&&ppi0.size()==1&&pk0.size()==1)) {
    if(part->id()==ParticleID::Dplus) {
      mminus = (pk0[0]+ppip[0]).m2();
      mplus  = (pk0[0]+ppi0[0]).m2();
      mpipi  = (ppip[0]+ppi0[0]).m2();
    }
    else {
      mminus = (pk0[0]+ppim[0]).m2();
      mplus  = (pk0[0]+ppi0[0]).m2();
      mpipi  = (ppim[0]+ppi0[0]).m2();
    }
    *_m2Kpip4 += mminus/GeV2;
    *_m2pipi4 += mpipi /GeV2;
    *_m2Kpi04 += mplus /GeV2;
    if(_points4.size()<200000) {
      _points4.push_back(make_pair(mplus,mpipi));
    }
  }
  else if ((ppim.size()==1&&pkp.size()==1&&ppip.size()==1)||
	   (ppim.size()==1&&pkm.size()==1&&ppip.size()==1)) {
    unsigned int itype(0);
    if(part->id()==ParticleID::Dplus) {
      itype=1;
      mplus  = (pkp [0]+ppip[0]).m2();
      mminus = (pkp [0]+ppim[0]).m2();
      mpipi  = (ppip[0]+ppim[0]).m2();
    }
    else if(part->id()==ParticleID::Dminus) {
      itype=1;
      mplus  = (pkm [0]+ppim[0]).m2();
      mminus = (pkm [0]+ppip[0]).m2();
      mpipi  = (ppip[0]+ppim[0]).m2();
    }
    else if(part->id()==ParticleID::D_splus) {
      itype=2;
      mplus  = (pkp [0]+ppip[0]).m2();
      mminus = (pkp [0]+ppim[0]).m2();
      mpipi  = (ppip[0]+ppim[0]).m2();
    }
    else if(part->id()==ParticleID::D_sminus) {
      itype=2;
      mplus  = (pkm [0]+ppim[0]).m2();
      mminus = (pkm [0]+ppip[0]).m2();
      mpipi  = (ppip[0]+ppim[0]).m2();
    }
    else return;
    if(itype==1) {
      *_mkppim5 +=mminus/GeV2;
      *_mkppip5 +=mplus /GeV2;
      *_mpippim5+=mpipi /GeV2;
      _points5.push_back(make_pair(mminus,mpipi));
    }
    else {
      *_mkppim6 +=mminus/GeV2;
      *_mkppip6 +=mplus /GeV2;
      *_mpippim6+=mpipi /GeV2;
      _points6.push_back(make_pair(mminus,mpipi));
    }
  }


}

void DDalitzAnalysis::persistentOutput(PersistentOStream & ) const {
}

void DDalitzAnalysis::persistentInput(PersistentIStream & , int) {
}

ClassDescription<DDalitzAnalysis> DDalitzAnalysis::initDDalitzAnalysis;
// Definition of the static class description member.

void DDalitzAnalysis::Init() {

  static ClassDocumentation<DDalitzAnalysis> documentation
    ("There is no documentation for the DDalitzAnalysis class");

}

void DDalitzAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _m2plus1->topdrawOutput(output,Frame|Errorbars,"RED",
			 "m0+1223 for D203RK203P2+3P2-3",
			 " X XX X      X XW X XGX XGX X",
			 "1/GdG/m0+1223/GeV2-23",
			 "  F F  X XX X    X  X",
			 "m0+1223/GeV223",
			 " X XX X    X X");
  _m2minus1->topdrawOutput(output,Frame|Errorbars,"RED",
			  "m0-1223 for D203RK203P2+3P2-3",
			  " X XX X      X XW X XGX XGX X",
			  "1/GdG/m0-1223/GeV2-23",
			  "  F F  X XX X    X  X",
			  "m0-1223/GeV223",
			  " X XX X    X X");
  _m2pipi1->topdrawOutput(output,Frame|Errorbars,"RED",
			 "m0PP1223 for D203RK203P2+3P2-3",
			 " XGGXX X      X XW X XGX XGX X",
			 "1/GdG/m0PP1223/GeV2-23",
			 "  F F  XGGXX X    X  X",
			 "m0PP1223/GeV223",
			 " XGGXX X    X X");
  _m2minus2->topdrawOutput(output,Frame|Errorbars,"RED",
			   "m0K2-3P2031223 for D203RK2-3P2+3P203",
			   " X X XGX XXX X      X XW X XGX XGX X",
			   "1/GdG/m0K2-3P2031223/GeV2-23",
			   "  F F  X X XGX XXX X    X  X",
			   "m0K2-3P2031223/GeV223",
			   " X X XGX XXX X    X X");
  _m2neutral2->topdrawOutput(output,Frame|Errorbars,"RED",
			     "m0K2-3P2+31223 for D203RK2-3P2+3P203",
			     " X X XGX XXX X      X XW X XGX XGX X",
			     "1/GdG/m0K2-3P2+31223/GeV2-23",
			     "  F F  X X XGX XXX X    X  X",
			     "m0K2-3P2+31223/GeV223",
			     " X X XGX XXX X    X X");
  _m2pipi2->topdrawOutput(output,Frame|Errorbars,"RED",
			 "m0PP1223 for D203RK2-3P2+3P203",
			 " XGGXX X      X XW X XGX XGX X",
			 "1/GdG/m0PP1223/GeV2-23",
			 "  F F  XGGXX X    X  X",
			 "m0PP1223/GeV223",
			 " XGGXX X    X X");
  _mKpilow3->topdrawOutput(output,Frame|Errorbars,"RED",
			   "m0K2-3P2+31223 for D2+3RK2-3P2+3P2+3 low",
			   " X X XGX XXX X      X XW X XGX XGX X    ",
			   "1/GdG/m0K2-3P2+31223/GeV2-23",
			   "  F F  X X XGX XXX X    X  X",
			   "m0K2-3P2+31223/GeV223",
			   " X X XGX XXX X    X X");
  _mKpihigh3->topdrawOutput(output,Frame|Errorbars,"RED",
			   "m0K2-3P2+31223 for D2+3RK2-3P2+3P2+3 high",
			   " X X XGX XXX X      X XW X XGX XGX X     ",
			   "1/GdG/m0K2-3P2+31223/GeV2-23",
			   "  F F  X X XGX XXX X    X  X",
			   "m0K2-3P2+31223/GeV223",
			   " X X XGX XXX X    X X");
  _mKpiall3->topdrawOutput(output,Frame|Errorbars,"RED",
			   "m0K2-3P2+31223 for D2+3RK2-3P2+3P2+3 all",
			   " X X XGX XXX X      X XW X XGX XGX X    ",
			   "1/GdG/m0K2-3P2+31223/GeV2-23",
			   "  F F  X X XGX XXX X    X  X",
			   "m0K2-3P2+31223/GeV223",
			   " X X XGX XXX X    X X");
  _mpipi3->topdrawOutput(output,Frame|Errorbars,"RED",
			 "m0P2+3P2+31223 for D2+3RK2-3P2+3P2+3 all",
			 " XGX XGX XXX X      X XW X XGX XGX X    ",
			 "1/GdG/m0P2+3P2+31223/GeV2-23",
			 "  F F  XGX XGX XXX X    X  X",
			 "m0P2+3P2+31223/GeV223",
			 " XGX XGX XXX X    X X");
  _m2Kpip4->topdrawOutput(output,Frame|Errorbars,"RED",
			 "m0K203P2+31223 for D2+3RK203P2+3P203",
			 " X X XGX XXX X      X XW X XGX XGX X",
			 "1/GdG/m0K203P2+31223/GeV2-23",
			 "  F F  X X XGX XXX X    X  X",
			 "m0K203P2+31223/GeV223",
			 " X X XGX XXX X    X X");
  _m2pipi4->topdrawOutput(output,Frame|Errorbars,"RED",
			  "m0P203P2+31223 for D2+3RK203P2+3P203",
			  " XGX XGX XXX X      X XW X XGX XGX X",
			  "1/GdG/m0P203P2+31223/GeV2-23",
			  "  F F  XGX XGX XXX X    X  X",
			  "m0P203P2+31223/GeV223",
			  " XGX XGX XXX X    X X");
  _m2Kpi04->topdrawOutput(output,Frame|Errorbars,"RED",
			 "m0K203P2-31223 for D2+3RK203P2+3P203",
			 " X X XGX XXX X      X XW X XGX XGX X",
			 "1/GdG/m0K203P2-31223/GeV2-23",
			 "  F F  X X XGX XXX X    X  X",
			 "m0K203P2-31223/GeV223",
			 " X X XGX XXX X    X X");
  _mkppim5->topdrawOutput(output,Frame|Errorbars,"RED",
			  "m0K2+3P2-31223 for D2+3RK2+3P2-3P2+3",
			  " X X XGX XXX X      X XW X XGX XGX X",
			  "1/GdG/m0K2+3P2-31223/GeV2-23",
			  "  F F  X X XGX XXX X    X  X",
			  "m0K2+3P2-31223/GeV223",
			  " X X XGX XXX X    X X");
  _mkppip5->topdrawOutput(output,Frame|Errorbars,"RED",
			  "m0K2+3P2+31223 for D2+3RK2+3P2-3P2+3",
			  " X X XGX XXX X      X XW X XGX XGX X",
			  "1/GdG/m0K2+3P2+31223/GeV2-23",
			  "  F F  X X XGX XXX X    X  X",
			  "m0K2+3P2+31223/GeV223",
			  " X X XGX XXX X    X X");
  _mpippim5->topdrawOutput(output,Frame|Errorbars,"RED",
			   "m0P2+3P2-31223 for D2+3RK2+3P2-3P2+3",
			   " XGX XGX XXX X      X XW X XGX XGX X",
			   "1/GdG/m0P2+3P2-31223/GeV2-23",
			   "  F F  XGX XGX XXX X    X  X",
			   "m0P2+3P2-31223/GeV223",
			   " XGX XGX XXX X    X X");
  _mkppim6->topdrawOutput(output,Frame|Errorbars,"RED",
			  "m0K2+3P2-31223 for D0s12+3RK2+3P2-3P2+3",
			  " X X XGX XXX X      X XX XW X XGX XGX X",
			  "1/GdG/m0K2+3P2-31223/GeV2-23",
			  "  F F  X X XGX XXX X    X  X",
			  "m0K2+3P2-31223/GeV223",
			  " X X XGX XXX X    X X");
  _mkppip6->topdrawOutput(output,Frame|Errorbars,"RED",
			  "m0K2+3P2+31223 for D0s12+3RK2+3P2-3P2+3",
			  " X X XGX XXX X      X XX XW X XGX XGX X",
			  "1/GdG/m0K2+3P2+31223/GeV2-23",
			  "  F F  X X XGX XXX X    X  X",
			  "m0K2+3P2+31223/GeV223",
			  " X X XGX XXX X    X X");
  _mpippim6->topdrawOutput(output,Frame|Errorbars,"RED",
			   "m0P2+3P2-31223 for D0s12+3RK2+3P2-3P2+3",
			   " XGX XGX XXX X      X XX XW X XGX XGX X",
			   "1/GdG/m0P2+3P2-31223/GeV2-23",
			   "  F F  XGX XGX XXX X    X  X",
			   "m0P2+3P2-31223/GeV223",
			   " XGX XGX XXX X    X X");
  if(!_points1.empty()) {
    output << "new frame\n";
    output << "set font duplex\n";
    output << "set limits x 0 3 y 0 3\n";
    output << "set order x y \n";
    output << "title top \"Dalitz plot for D203RK203P2+3P2-3\"\n";
    output << "case      \"                 X XW X XGX XGX X\"\n";
    output << "title left   \"m0-1223/GeV223\"\n";
    output << "case         \" X XX X    X X\"\n";
    output << "title bottom \"m0+1223/GeV223\"\n";
    output << "case         \" X XX X    X X\"\n";
    for(unsigned int ix=0;ix<_points1.size();++ix) {
      output << _points1[ix].first /GeV2 << " "
	     << _points1[ix].second/GeV2 << "\n"; 
      if(ix%50000==0) output << "plot\n";
    }
    output << "plot\n";
  }
  if(!_points2.empty()) {
    output << "new frame\n";
    output << "set font duplex\n";
    output << "set limits x 0 3 y 0 3\n";
    output << "set order x y \n";
    output << "title top \"Dalitz plot for D203RK2-3P2+3P203\"\n";
    output << "case      \"                 X XW X XGX XGX X\"\n";
    output << "title left   \"m0-1223/GeV223\"\n";
    output << "case         \" X XX X    X X\"\n";
    output << "title bottom \"m0+1223/GeV223\"\n";
    output << "case         \" X XX X    X X\"\n";
    for(unsigned int ix=0;ix<_points2.size();++ix) {
      output << _points2[ix].first /GeV2 << " "
	     << _points2[ix].second/GeV2 << "\n"; 
      if(ix%50000==0) output << "plot\n";
    }
    output << "plot\n";
  }
  if(!_points3.empty()) {
    output << "new frame\n";
    output << "set font duplex\n";
    output << "set limits x 0 3 y 0 3\n";
    output << "set order x y \n";
    output << "title top \"Dalitz plot for D2+3RK2-3P2+3P2+3\"\n";
    output << "case      \"                 X XW X XGX XGX X\"\n";
    output << "title left   \"m0K2-3P2+31223\"\n";
    output << "case         \" X X XGX XXX X\"\n";
    output << "title bottom \"m0P2+3P2+31223\"\n";
    output << "case         \" XGX XGX XXX X\"\n";
    for(unsigned int ix=0;ix<_points3.size();++ix) {
      output << _points3[ix].first /GeV2 << " "
	     << _points3[ix].second/GeV2 << "\n"; 
      if(ix%50000==0) output << "plot\n";
    }
    output << "plot\n";
  }
  if(!_points4.empty()) {
    output << "new frame\n";
    output << "set font duplex\n";
    output << "set limits x 0 3 y 0 3\n";
    output << "set order x y \n";
    output << "title top \"Dalitz plot for D2+3RK203P2+3P203\"\n";
    output << "case      \"                 X XW X XGX XGX X\"\n";
    output << "title left   \"m0P2+3P2031223\"\n";
    output << "case         \" XGX XGX XXX X\"\n";
    output << "title bottom \"m0K203P2031223\"\n";
    output << "case         \" X X XGX XXX X\"\n";
    for(unsigned int ix=0;ix<_points4.size();++ix) {
      output << _points4[ix].first /GeV2 << " "
	     << _points4[ix].second/GeV2 << "\n"; 
      if(ix%50000==0) output << "plot\n";
    }
    output << "plot\n";
  }
  if(!_points5.empty()) {
    output << "new frame\n";
    output << "set font duplex\n";
    output << "set limits x 0 3 y 0 3\n";
    output << "set order x y \n";
    output << "title top \"Dalitz plot for D2+3RK2+3P2-3P2+3\"\n";
    output << "case      \"                 X XW X XGX XGX X\"\n";
    output << "title left   \"m0P2+3P2-31223\"\n";
    output << "case         \" XGX XGX XXX X\"\n";
    output << "title bottom \"m0K2+3P2-31223\"\n";
    output << "case         \" X X XGX XXX X\"\n";
    for(unsigned int ix=0;ix<_points5.size();++ix) {
      output << _points5[ix].first /GeV2 << " "
	     << _points5[ix].second/GeV2 << "\n"; 
      if(ix%50000==0) output << "plot\n";
    }
    output << "plot\n";
  }
  if(!_points6.empty()) {
    output << "new frame\n";
    output << "set font duplex\n";
    output << "set limits x 0 3.5 y 0 2.5\n";
    output << "set order x y \n";
    output << "title top \"Dalitz plot for D0s12+3RK2+3P2-3P2+3\"\n";
    output << "case      \"                 X XX XW X XGX XGX X\"\n";
    output << "title left   \"m0P2+3P2-31223\"\n";
    output << "case         \" XGX XGX XXX X\"\n";
    output << "title bottom \"m0K2+3P2-31223\"\n";
    output << "case         \" X X XGX XXX X\"\n";
    for(unsigned int ix=0;ix<_points6.size();++ix) {
      output << _points6[ix].first /GeV2 << " "
	     << _points6[ix].second/GeV2 << "\n"; 
      if(ix%50000==0) output << "plot\n";
    }
    output << "plot\n";
  }
}
 
void DDalitzAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  _m2plus1    = new_ptr(Histogram(0.,3.,200));
  _m2minus1   = new_ptr(Histogram(0.,3.,200));
  _m2pipi1    = new_ptr(Histogram(0.,2.,200));
  _m2minus2   = new_ptr(Histogram(0.,3.,200));
  _m2neutral2 = new_ptr(Histogram(0.,3.,200));
  _m2pipi2    = new_ptr(Histogram(0.,2.,200));
  _mKpilow3   = new_ptr(Histogram(0.,3.,200));
  _mKpihigh3  = new_ptr(Histogram(0.,3.,200));
  _mKpiall3   = new_ptr(Histogram(0.,3.,200));
  _mpipi3     = new_ptr(Histogram(0.,2.,200));
  _m2Kpip4    = new_ptr(Histogram(0.,3.,200));
  _m2pipi4    = new_ptr(Histogram(0.,2.,200));
  _m2Kpi04    = new_ptr(Histogram(0.,3.,200));
  _mkppim5    = new_ptr(Histogram(0.,3.,200));
  _mkppip5    = new_ptr(Histogram(0.,3.,200));
  _mpippim5   = new_ptr(Histogram(0.,3.,200));
  _mkppim6    = new_ptr(Histogram(0.,3.5,200));
  _mkppip6    = new_ptr(Histogram(0.,3.5,200));
  _mpippim6   = new_ptr(Histogram(0.,3.,200));
}
