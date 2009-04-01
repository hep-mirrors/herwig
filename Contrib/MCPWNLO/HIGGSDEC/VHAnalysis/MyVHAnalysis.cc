// -*- C++ -*-
//
// MyVHAnalysis.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MyVHAnalysis class.
//

#include "MyVHAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MyVHAnalysis::MyVHAnalysis() :
   _Mbb(0.,150.,150),
   _Lbb(-350.,350.,175), _rapbb(-5.,5.,100), 
   _PTbmH(0.,300.,125),
   _PTVH(0.,600.,150), _Ebb(0.,800., 200) {}

  void MyVHAnalysis::analyze(tEventPtr event, long, int, int) {

  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z
    string infile = "PPVH.dat";
    ifstream normin;
    string stringin = "";
    int IPROC = 24;
    int stringdoub = 0;
    bool next = 0, filerror = 0;
    normin.open(infile.c_str());
    if(!normin) { cerr << "Error: Failed to open file " << infile << endl; filerror = 1; IPROC = 24; }
    if(!filerror) {
      while(normin) {

	normin >> stringin;

	if(next) {   normin >> stringdoub; IPROC = stringdoub; break; }
	if(stringin == "IPROC") { next = 1; }

      }
    }

  Lorentz5Momentum pz, pVd;
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator stest =event->primaryCollision()->steps().end();
  StepVector::const_iterator send=sit;
  ++send;
  if(send==stest) --send;
  ++send;
  if(send==stest) --send;
  ++send;
  double ptbb, mzbb, pzz, rapbb, ptbb2, pggg, pzemax, pzemax2, ptbb3;
  pzemax = 0.; pzemax2 = 0.;
  Lorentz5Momentum pzbb, pzV, pg, pgg, pzb, pzb2;
  for(sit=send-1;sit!=send;++sit) {
    ParticleSet part=(**sit).all();
    ParticleSet::const_iterator iter = part.begin(), end = part.end();
    for( ;iter!=end;++iter) {
      // if((**iter).children().size()==2) continue;
      //   if(abs((**iter).id())!= 5 || abs((**iter).id())!=11 || abs((**iter).id())!=12) continue;
      if((**iter).id() == 5) {
	pz =(*iter)->momentum();
	if (pz.e() > pzemax) {pzb = pz;}}
      if((**iter).id() == -5) {
        pz =(*iter)->momentum();
        if (pz.e() > pzemax2) {pzb2 = pz;}}
      if (IPROC == 24) {
      if (abs((**iter).id())== 11 || abs((**iter).id()) == 12) {
	pVd = (*iter)->momentum();
	pzV += pVd;
      }} else if (IPROC == 23) {
      if (abs((**iter).id())== 11) {
        pVd = (*iter)->momentum();
        pzV += pVd;
      }}

     
    }
     pzbb = pzb + pzb2 ;
     ptbb = sqrt(pow(pzbb.x(),2)+pow(pzbb.y(),2)+pow(pzbb.z(),2))*sin(acos((pzV.x()*pzbb.x()+pzV.y()*pzbb.y()+pzV.z()*pzbb.z())/(sqrt(pow(pzV.x(),2)+pow(pzV.y(),2)+pow(pzV.z(),2))*sqrt(pow(pzbb.x(),2)+pow(pzbb.y(),2)+pow(pzbb.z(),2)))));
      mzbb = sqrt(pow(pzbb.e(),2)- pow(pzbb.x(),2)-pow(pzbb.y(),2)-pow(pzbb.z(),2));
      ptbb = abs(ptbb)/GeV;
      ptbb2 = pzbb.perp()/GeV;
      mzbb = mzbb/GeV;
      pzz = pzbb.z()/GeV;
      rapbb = pzbb.rapidity();
       if (mzbb < 10.0/GeV) continue;
      _PTVH.addWeighted(ptbb           ,event->weight());
      _Mbb    .addWeighted(mzbb           ,event->weight());
      _Lbb   .addWeighted(pzz           ,event->weight());
      _rapbb   .addWeighted(rapbb           ,event->weight());
      _PTbmH .addWeighted(ptbb2           ,event->weight());
      _Ebb.addWeighted(pzbb.e()/GeV           ,event->weight());

      //cout << ptbb << "\t" << mzbb << "\t" << ptbb2 << endl;
	//if(mz>20.&&mz<80.)        _ptZ[1].addWeighted(pt,event->weight());
	//	else if (mz>80.&&mz<100.) _ptZ[2].addWeighted(pt,event->weight());
	//else if (mz>100.)         _ptZ[3].addWeighted(pt,event->weight());
	//_ptZ[0].addWeighted(pt           ,event->weight());
	//_mZ    .addWeighted(mz           ,event->weight());
	//	_rapZ  .addWeighted(pz.rapidity(),event->weight());
	//_phiZ  .addWeighted(pz.phi()     ,event->weight());
	//    } 
      //else if ((**iter).id()==ParticleID::Wplus) {
      //	pz=(*iter)->momentum();
      //	double pt = pz.perp()/GeV;
      //	double mz = pz.mass()/GeV;
      //	if(mz>20.&&mz<80.)        _ptWp[1].addWeighted(pt,event->weight());
      //	else if (mz>80.&&mz<100.) _ptWp[2].addWeighted(pt,event->weight());
      //	else if (mz>100.)         _ptWp[3].addWeighted(pt,event->weight());
      //	_ptWp[0].addWeighted(pt           ,event->weight());
      //	_mWp    .addWeighted(mz           ,event->weight());
      //	_rapWp  .addWeighted(pz.rapidity(),event->weight());
      //	_phiWp  .addWeighted(pz.phi()     ,event->weight());
      // } 
      // else if ((**iter).id()==ParticleID::Wminus) {
      //	pz=(*iter)->momentum();
      //double pt = pz.perp()/GeV;
      //	double mz = pz.mass()/GeV;
      //	if(mz>20.&&mz<80.)        (_ptWm[1]).addWeighted(pt,event->weight());
      //	else if (mz>80.&&mz<100.) (_ptWm[2]).addWeighted(pt,event->weight());
      //	else if (mz>100.)         (_ptWm[3]).addWeighted(pt,event->weight());
      //	_ptWm[0].addWeighted(pt           ,event->weight());
      //	_mWm    .addWeighted(mz           ,event->weight());
      //	_rapWm  .addWeighted(pz.rapidity(),event->weight());
      //	_phiWm  .addWeighted(pz.phi()     ,event->weight());
	// }
      // }
    
  }

}
//  Lorentz5Momentum ps;
// double pss;
//  for (sit = send; sit != send+1; sit++) {
//    ParticleSet part2(**sit).all();
//    ParticleSet::const_iterator iter2= part2.begin(), end2 = part2.end();
//    for( ;iter2!=end2;++iter2) {
//     ps = (*iter2)->momentum();
//     pss += abs(ps.perp()/GeV);
// }
//   _rapWm.addWeighted(pss,event->weight());
// }

NoPIOClassDescription<MyVHAnalysis> MyVHAnalysis::initMyVHAnalysis;
// Definition of the static class description member.

void MyVHAnalysis::Init() {

  static ClassDocumentation<MyVHAnalysis> documentation
    ("The MyVHAnalysis class performs a simple analysis of W and"
     " Z production in hadron-hadron collisions");

}

void MyVHAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  string title;
  using namespace HistogramOptions;
  //  for(unsigned int ix=0;ix<4;++ix) {
  // if(ix==0){title="pt of Z for all masses ";}
  // else if(ix==1){title="pt of Z for mass 40-80 GeV";}
  // else if(ix==2){title="pt of Z for mass 80-100 GeV";}
  // else if(ix==3){title="pt of Z for mass 100- GeV";}
  // _ptZ[ix].topdrawOutput(outfile,Frame,"BLACK",title);
  // _ptZ[ix].topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  // if(ix==0){title="pt of Wp for all masses ";}
  // else if(ix==1){title="pt of Wp for mass 40-80 GeV";}
  // else if(ix==2){title="pt of Wp for mass 80-100 GeV";}
  // else if(ix==3){title="pt of Wp for mass 100- GeV";}
  // _ptWp[ix].topdrawOutput(outfile,Frame,"BLACK",title);
  // _ptWp[ix].topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  // if(ix==0){title="pt of Wm for all masses ";}
  // else if(ix==1){title="pt of Wm for mass 40-80 GeV";}
  // else if(ix==2){title="pt of Wm for mass 80-100 GeV";}
  // else if(ix==3){title="pt of Wm for mass 100- GeV";}
  // _ptWm[ix].topdrawOutput(outfile,Frame,"BLACK",title);
  // _ptWm[ix].topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  // }
   _Mbb.topdrawOutput(outfile,Frame,"BLACK","Mass of b anti-b pair");
   _Ebb.topdrawOutput(outfile,Frame|Ylog,"BLACK","Energy of b anti-b pair");
   _PTbmH.topdrawOutput(outfile,Frame|Ylog,"BLACK","pT of b anti-b pair wrt to beam axis");
   _PTVH.topdrawOutput(outfile,Frame|Ylog,"BLACK","pT of b anti-b pair wrt to the vector boson");
   _Lbb.topdrawOutput(outfile,Frame,"BLACK", "Longtitudinal mmt of b anti-b pair");
   _rapbb.topdrawOutput(outfile,Frame|Ylog,"BLACK", "Rapidity of b anti-b pair");
}

