// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MySimpleAnalysis class.
//

#include "MySimpleAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MySimpleAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MySimpleAnalysis::MySimpleAnalysis() :
  _ptZ(4,Histogram(0.,250.,250)), _ptWp(4,Histogram(0.,250.,250)),_ptWm(4,Histogram(0.,250.,250)), 
  _mWp(-250.,250.,100),_mZ(-250.,250.,100),_mWm(-250.,250.,100),
  _rapZ(-10.,10., 100), _rapWp(-10.,10., 100),_rapWm(-10.,10., 100),
  _rapem(-10.,10., 100),_rapep(-10.,10., 100),_rapnu(-10.,10., 100),_rapanu(-10.,10., 100),
 
     _ptW4(0.,20.,10),
     _ptW5(20.,40.,4),
     _ptW6(40.,80.,4),
     _ptW7(80.,100.,1),
     _ptW8(100.,120.,1),
     _ptW9(120.,160.,1),
     _ptW10(160.,200.,1),
     _ptZ4(0.,12.,24),
     _ptZ5(12.,20.,8),
     _ptZ6(20.,30.,5),
     _ptZ7(30.,50.,5),
     _ptZ8(50.,100.,5),
     _ptZ9(100.,150.,2),
     _ptZ10(150.,200.,1)


{}

//MySimpleAnalysis::~MySimpleAnalysis() {}

// namespace {
//   inline Lorentz5Momentum getMomentum(tcPPtr particle) {
//     return particle->momentum();
//     //Lorentz5Momentum tmp = particle->children()[0]->next()->momentum();
//     //tmp += particle->children()[1]->next()->momentum();
//     //tmp.rescaleMass();
//     //return tmp;

//   }
// }
   
void MySimpleAnalysis::analyze(tEventPtr event, long, int, int) {
string infile = "DYPP.dat";

  ifstream normin;
  string stringin = "";
  double normfac = 1;
  double stringdoub = 0;
  bool next = 0, filerror = 0;
  normin.open(infile.c_str());
  if(!normin) { cerr << "Error: Failed to open file " << infile << endl; filerror = 1; normfac = 1.; }
  if(!filerror) {
    while(normin) { 
      
      normin >> stringin;
     
      if(next) {   normin >> stringdoub; normfac = stringdoub; break; }
      if(stringin == "NORMFACTOR") { next = 1; }

    }
   
  }
 
    eventweight_=event->weight()*normfac;
   
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z or W
  Lorentz5Momentum pz;
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator send=event->primaryCollision()->steps().end();
  for(;sit!=send;++sit)
    {
      ParticleSet part=(**sit).all();
      ParticleSet::const_iterator iter=part.begin();
      ParticleSet::const_iterator end =part.end();
      for( ;iter!=end;++iter)
	{
	  if(((**iter).id()==ParticleID::Z0||(**iter).id()==ParticleID::gamma)&&
	     (**iter).children().size()==2)
	    {
	      pz=(*iter)->momentum();
	      Energy pt = pz.perp()/GeV;
	      Energy mz = pz.mass()/GeV;
	      if(mz>20.&&mz<80.) _ptZ[1].addWeighted(pt,eventweight_);
	      else if (mz>80.&&mz<100.) _ptZ[2].addWeighted(pt,eventweight_);
	      else if (mz>100.) _ptZ[3].addWeighted(pt,eventweight_);
	      if(mz>66.&&mz<116.) {
               _ptZ4.addWeighted(pt,eventweight_);
	       _ptZ5.addWeighted(pt,eventweight_);
	       _ptZ6.addWeighted(pt,eventweight_);
	       _ptZ7.addWeighted(pt,eventweight_);
	       _ptZ8.addWeighted(pt,eventweight_);
	       _ptZ9.addWeighted(pt,eventweight_);
	       _ptZ10.addWeighted(pt,eventweight_);}
              _ptZ[0].addWeighted(pt,eventweight_);
	      _mZ.addWeighted(mz,eventweight_);
	      double rap = 0.5*log((pz.e()+pz.z())/(pz.e()-pz.z()));
	      _rapZ.addWeighted(rap,eventweight_);
	      _phiZ.addWeighted(pz.phi(),eventweight_);
	    } else if ((**iter).id()==ParticleID::Wplus) {
	      pz=(*iter)->momentum();
	      Energy pt = pz.perp()/GeV;
	      Energy mz = pz.mass()/GeV;
	      if(mz>20.&&mz<80.) _ptWp[1].addWeighted(pt,eventweight_);
	      else if (mz>80.&&mz<100.) _ptWp[2].addWeighted(pt,eventweight_);
	      else if (mz>100.) _ptWp[3].addWeighted(pt,eventweight_);
	       _ptWp[0].addWeighted(pt,eventweight_);
	       _mWp.addWeighted(mz,eventweight_);
	      double rap = 0.5*log((pz.e()+pz.z())/(pz.e()-pz.z()));
	      _rapWp.addWeighted(rap,eventweight_);
	      _phiWp.addWeighted(pz.phi(),eventweight_);
	    } else if ((**iter).id()==ParticleID::Wminus) {
	      pz=(*iter)->momentum();
	      Energy pt = pz.perp()/GeV;
	      Energy mz = pz.mass()/GeV;
	      if(mz>20.&&mz<80.) _ptWm[1].addWeighted(pt,eventweight_);
	      else if (mz>80.&&mz<100.) _ptWm[2].addWeighted(pt,eventweight_);
	      else if (mz>100.) _ptWm[3].addWeighted(pt,eventweight_);
	      _ptWm[0].addWeighted(pt,eventweight_);
	      _mWm.addWeighted(mz,eventweight_);
	      double rap = 0.5*log((pz.e()+pz.z())/(pz.e()-pz.z()));
	      _rapWm.addWeighted(rap,eventweight_);
	      _phiWm.addWeighted(pz.phi(),eventweight_);
	    }
	
	  if ((**iter).id()==ParticleID::Wplus || (**iter).id()==ParticleID::Wminus) {
	      pz=(*iter)->momentum();
	      Energy pt = pz.perp()/GeV;
	      Energy mz = pz.mass()/GeV;
	        if(mz>70.&&mz<100.){ 
	       _ptW4.addWeighted(pt,eventweight_);
	       _ptW5.addWeighted(pt,eventweight_);
	       _ptW6.addWeighted(pt,eventweight_);
	       _ptW7.addWeighted(pt,eventweight_);
	       _ptW8.addWeighted(pt,eventweight_);
	       _ptW9.addWeighted(pt,eventweight_);
	       _ptW10.addWeighted(pt,eventweight_);
	      }}
	}
    }

  
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the leptons

 
  ParticleVector lept;
  lept.clear();
 

  for(sit=send-1;sit!=send;++sit)
  {
  ParticleSet part=(**sit).all();
  ParticleSet::const_iterator iter1=part.begin();
  ParticleSet::const_iterator end1 =part.end();

  for( ;iter1!=end1;++iter1) {
     if ((fabs((**iter1).id())==11 || fabs((**iter1).id())==12) && (**iter1).children().size()==0) {
  
      lept.push_back((*iter1));
     } }
	  


  Lorentz5Momentum pzep(0.0,0.0,0.0,0.0,0.0);
  Lorentz5Momentum pzem(0.0,0.0,0.0,0.0,0.0);
  Lorentz5Momentum pznu(0.0,0.0,0.0,0.0,0.0);
  Lorentz5Momentum pzanu(0.0,0.0,0.0,0.0,0.0);
  Lorentz5Momentum fullmom(0.0,0.0,0.0,0.0,0.0);
  ParticleVector::const_iterator iter2=lept.begin();
  ParticleVector::const_iterator end2 =lept.end();
  for( ;iter2!=end2;++iter2) {
    if((**iter2).id()==-11) {
      pzep=(*iter2)->momentum();}
    else if ((**iter2).id()==11) {
      pzem=(*iter2)->momentum();} 
    else if ((**iter2).id()==12) {
      pznu=(*iter2)->momentum();} 
    else if ((**iter2).id()==-12) {
      pzanu=(*iter2)->momentum();} 

  }
  
    double rapem = 0.5*log((pzem.e()+pzem.z())/(pzem.e()-pzem.z()));
    double rapep = 0.5*log((pzep.e()+pzep.z())/(pzep.e()-pzep.z()));
    double rapnu = 0.5*log((pznu.e()+pznu.z())/(pznu.e()-pznu.z()));
    double rapanu = 0.5*log((pzanu.e()+pzanu.z())/(pzanu.e()-pzanu.z()));
    

    _rapem.addWeighted(rapem,eventweight_); 
    _rapep.addWeighted(rapep,eventweight_);
    _rapnu.addWeighted(rapnu,eventweight_);
    _rapanu.addWeighted(rapanu,eventweight_);



  }
}




//LorentzRotation MySimpleAnalysis::transform(tEventPtr) const {
//  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
//}

//void MySimpleAnalysis::analyze(const tPVector & particles) {
//  AnalysisHandler::analyze(particles);
//}

//void MySimpleAnalysis::analyze(tPPtr) {}

//void MySimpleAnalysis::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
//}

//void MySimpleAnalysis::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
//}

NoPIOClassDescription<MySimpleAnalysis> MySimpleAnalysis::initMySimpleAnalysis;
// Definition of the static class description member.

void MySimpleAnalysis::Init() {

  static ClassDocumentation<MySimpleAnalysis> documentation
    ("There is no documentation for the MySimpleAnalysis class");

}

void MySimpleAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  string title;
  using namespace HistogramOptions;
     title="pT of Z boson (bin-sized for data comparison with CDF Run I at 1.8TeV)";
     _ptZ4.topdrawOutput(outfile,Frame,"",title);
     _ptZ5.topdrawOutput(outfile,None,"", title="");
    _ptZ6.topdrawOutput(outfile,None,"",title="");
     _ptZ7.topdrawOutput(outfile,None,"",title="");
    _ptZ8.topdrawOutput(outfile,None,"",title="");
    _ptZ9.topdrawOutput(outfile,None,"",title="");
    _ptZ10.topdrawOutput(outfile,None,"BLACK",title="");
     title="pT of W boson (bin-sized for data comparison with D0 Run I at 1.8TeV)";
     _ptW4.topdrawOutput(outfile,Frame,"",title);
     _ptW5.topdrawOutput(outfile,None,"", title="");
    _ptW6.topdrawOutput(outfile,None,"",title="");
     _ptW7.topdrawOutput(outfile,None,"",title="");
    _ptW8.topdrawOutput(outfile,None,"",title="");
    _ptW9.topdrawOutput(outfile,None,"",title="");
    _ptW10.topdrawOutput(outfile,None,"BLACK",title="");

    //   _phianu.topdrawOutput(outfile,Frame,"BLACK","Azimuth of anitneutrino");
    // _phinu.topdrawOutput(outfile,Frame,"BLACK","Azimuth of neutrino");
    // _phiep.topdrawOutput(outfile,Frame,"BLACK","Azimuth of positron");
    // _phiem.topdrawOutput(outfile,Frame,"BLACK","Azimuth of electron");
   
   _rapanu.topdrawOutput(outfile,Frame,"BLACK","Rapidity of anitneutrino");
   _rapnu.topdrawOutput(outfile,Frame,"BLACK","Rapidity of neutrino");
   _rapep.topdrawOutput(outfile,Frame,"BLACK","Rapidity of positron");
   _rapem.topdrawOutput(outfile,Frame,"BLACK","Rapidity of electron");
    for(unsigned int ix=0;ix<1;++ix)
       {
	       if(ix==0){title="pt of Z for all masses ";}
      else if(ix==1){title="pt of Z for mass 40-80 GeV";}
      else if(ix==2){title="pt of Z for mass 80-100 GeV";}
      else if(ix==3){title="pt of Z for mass 100- GeV";}
	       // _ptZ[ix].topdrawOutput(outfile,Frame,"BLACK",title);
       _ptZ[ix].topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
         if(ix==0){title="pt of Wp for all masses ";}
      else if(ix==1){title="pt of Wp for mass 40-80 GeV";}
      else if(ix==2){title="pt of Wp for mass 80-100 GeV";}
      else if(ix==3){title="pt of Wp for mass 100- GeV";}
	 //   _ptWp[ix].topdrawOutput(outfile,Frame,"BLACK",title);
       _ptWp[ix].topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
          if(ix==0){title="pt of Wm for all masses ";}
      else if(ix==1){title="pt of Wm for mass 40-80 GeV";}
      else if(ix==2){title="pt of Wm for mass 80-100 GeV";}
      else if(ix==3){title="pt of Wm for mass 100- GeV";}
	  //   _ptWm[ix].topdrawOutput(outfile,Frame,"BLACK",title);
       _ptWm[ix].topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
       }
      _mWp.topdrawOutput(outfile,Frame,"BLACK","Mass of Wp");
      _mZ.topdrawOutput(outfile,Frame,"BLACK","Mass of Z");
      //_mZ.topdrawOutput(outfile,Frame|Ylog,"BLACK", "Mass of Z");
    
 //  _mWp.topdrawOutput(outfile,Frame|Ylog,"BLACK", "Mass of Wp");
  _mWm.topdrawOutput(outfile,Frame,"BLACK","Mass of Wm");
//   _mWm.topdrawOutput(outfile,Frame|Ylog,"BLACK", "Mass of Wm");
    _rapZ.topdrawOutput(outfile,Frame,"BLACK","Rapidity of Z");
//   _rapZ.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of Z");
  _rapWp.topdrawOutput(outfile,Frame,"BLACK","Rapidity of Wp");
//   _rapWp.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of Wp");
  _rapWm.topdrawOutput(outfile,Frame,"BLACK","Rapidity of Wm");
//   _rapWm.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of Wm");

//   _phiZ.topdrawOutput(outfile,Frame,"BLACK","Azimuth of Z");
//   _phiWp.topdrawOutput(outfile,Frame,"BLACK","Azimuth of Wp");
//   _phiWm.topdrawOutput(outfile,Frame,"BLACK","Azimuth of Wm");
}



// LorentzRotation MySimpleAnalysis::transform(tEventPtr event) const {
//   return LorentzRotation();
//   // Return the Rotation to the frame in which you want to perform the analysis.
// }

// void MySimpleAnalysis::analyze(const tPVector & particles) {
//   AnalysisHandler::analyze(particles);
// }

// void MySimpleAnalysis::analyze(tPPtr) {}

// void MySimpleAnalysis::persistentOutput(PersistentOStream & os) const {
//   // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
// }

// void MySimpleAnalysis::persistentInput(PersistentIStream & is, int) {
//   // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
// }

// ClassDescription<MySimpleAnalysis> MySimpleAnalysis::initMySimpleAnalysis;
// // Definition of the static class description member.

// void MySimpleAnalysis::Init() {

//   static ClassDocumentation<MySimpleAnalysis> documentation
//     ("There is no documentation for the MySimpleAnalysis class");

// }

