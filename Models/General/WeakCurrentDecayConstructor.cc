// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakCurrentDecayConstructor class.
//

#include "WeakCurrentDecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Decay/General/GeneralCurrentDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;
using ThePEG::Helicity::VertexBasePtr;

void WeakCurrentDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << _theExistingDecayers << _init << _iteration << _points << ounit(_masscut,GeV)
     << _part1 << _part2 << _part3 << _part4 << _part5 << _norm << _current;
}

void WeakCurrentDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >>_theExistingDecayers >> _init >> _iteration >> _points >> iunit(_masscut,GeV)
     >> _part1 >> _part2 >> _part3 >> _part4 >> _part5 >> _norm >> _current;
}

ClassDescription<WeakCurrentDecayConstructor> WeakCurrentDecayConstructor::initWeakCurrentDecayConstructor;
// Definition of the static class description member.

void WeakCurrentDecayConstructor::Init() {

  static ClassDocumentation<WeakCurrentDecayConstructor> documentation
    ("There is no documentation for the WeakCurrentDecayConstructor class");
  
  static Switch<WeakCurrentDecayConstructor,bool> interfaceInitializeDecayers
    ("InitializeDecayers",
     "Initialize new decayers",
     &WeakCurrentDecayConstructor::_init, true, false, false);
  static SwitchOption interfaceInitializeDecayersInitializeDecayersOn
    (interfaceInitializeDecayers,
     "On",
     "Initialize new decayers to find max weights",
     1);
  static SwitchOption interfaceInitializeDecayersoff
    (interfaceInitializeDecayers,
     "Off",
     "Use supplied weights for integration",
     0);
  
  static Parameter<WeakCurrentDecayConstructor,int> interfaceInitIteration
    ("InitIteration",
     "Number of iterations to optimise integration weights",
     &WeakCurrentDecayConstructor::_iteration, 5, 0, 5,
     false, false, true);

  static Parameter<WeakCurrentDecayConstructor,int> interfaceInitPoints
    ("InitPoints",
     "Number of points to generate when optimising integration",
     &WeakCurrentDecayConstructor::_points, 10000, 100, 100000000,
     false, false, true);

  static ParVector<WeakCurrentDecayConstructor,long> interfaceParticle1
    ("Particle1",
     "The first decay product",
     &WeakCurrentDecayConstructor::_part1, -1, long(0), long(-1000000), long(1000000),
     false, false, Interface::limited);

  static ParVector<WeakCurrentDecayConstructor,long> interfaceParticle2
    ("Particle2",
     "The second decay product",
     &WeakCurrentDecayConstructor::_part2, -1, long(0), long(-1000000), long(1000000),
     false, false, Interface::limited);

  static ParVector<WeakCurrentDecayConstructor,long> interfaceParticle3
    ("Particle3",
     "The third   decay product",
     &WeakCurrentDecayConstructor::_part3, -1, long(0), long(-1000000), long(1000000),
     false, false, Interface::limited);

  static ParVector<WeakCurrentDecayConstructor,long> interfaceParticle4
    ("Particle4",
     "The fourth decay product",
     &WeakCurrentDecayConstructor::_part4, -1, long(0), long(-1000000), long(1000000),
     false, false, Interface::limited);

  static ParVector<WeakCurrentDecayConstructor,long> interfaceParticle5
    ("Particle5",
     "The fifth decay product",
     &WeakCurrentDecayConstructor::_part5, -1, long(0), long(-1000000), long(1000000),
     false, false, Interface::limited);

  static ParVector<WeakCurrentDecayConstructor,double> interfaceNormalisation
    ("Normalisation",
     "The normalisation of the different modes",
     &WeakCurrentDecayConstructor::_norm, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static RefVector<WeakCurrentDecayConstructor,WeakDecayCurrent> interfaceCurrent
    ("Current",
     "The current for the decay mode",
     &WeakCurrentDecayConstructor::_current, -1, false, false, true, false, false);

  static Parameter<WeakCurrentDecayConstructor,Energy> interfaceMassCut
    ("MassCut",
     "The maximum mass difference for the decay",
     &WeakCurrentDecayConstructor::_masscut, GeV, 5.0*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

void WeakCurrentDecayConstructor::DecayList(const vector<PDPtr> & part) {
  _theModel->init();
  unsigned int nv(_theModel->numberOfVertices());
  // make sure vertices are initialized
  for(unsigned int i = 0; i < nv; ++i) _theModel->vertex(i)->init();
  // resize the vectors
  _theExistingDecayers.
    resize(nv,vector<map<WeakDecayCurrentPtr,GeneralCurrentDecayerPtr> >
	   (3,map<WeakDecayCurrentPtr,GeneralCurrentDecayerPtr>()));
  PDVector decays;
  unsigned int np = part.size();
  for(unsigned int ipart = 0; ipart < np; ++ipart) {
    for(unsigned int iv = 0; iv < nv; ++iv) {
      for(unsigned int ilist = 0; ilist < 3; ++ilist) { 
	decays = createModes(part[ipart], _theModel->vertex(iv),
			     ilist, iv);
	if(decays.size() > 0){
	  tPDPtr incpart = part[ipart]->CC() ? part[ipart]->CC() : tPDPtr(part[ipart]);
	  createDecayMode(incpart, decays, _theExistingDecayers[iv][ilist]);
	}
      }
    }
  }
  //ParticleData objects need updating 
  for(unsigned int ip= 0; ip < np; ++ip) {
    PDPtr pdp = part[ip];
    pdp->touch();
    pdp->update();
    if(pdp->CC()) pdp->CC()->synchronize();
  }
}
  
vector<PDPtr> WeakCurrentDecayConstructor::createModes(const PDPtr inpart,
						       const VertexBasePtr vert,
						       unsigned int ilist,
						       unsigned int iv) {
  int id = inpart->id();
  if(id < 0) return vector<PDPtr>(0);
  Energy m1(inpart->mass());
  vector<PDPtr> decaylist;
  if(vert->getNpoint()==3 && vert->incoming(id)) {
    decaylist = vert->search(ilist,id);
    for(PDVector::iterator iter=decaylist.begin();iter!=decaylist.end();) {
      Energy m2,m3;
      bool cc1((*iter)->CC()),cc2((*(iter+1))->CC()),cc3((*(iter+2))->CC());
      int id2,id3;
      if((*iter)->id()==id) {
 	m2  = (*(iter+1))->mass();
 	m3  = (*(iter+2))->mass();
	id2 = abs((*(iter+1))->id());
	id3 = abs((*(iter+2))->id());
 	if(cc1) {
 	  if(cc2) *(iter+1) = (*(iter+1))->CC();	  	  
 	  if(cc3) *(iter+2) = (*(iter+2))->CC();	  
 	}
      }
      else if((*(iter+1))->id()==id) {
 	m2 = (*iter)    ->mass();
 	m3 = (*(iter+2))->mass();
	id2 = abs((*iter)    ->id());
	id3 = abs((*(iter+2))->id());
  	if(cc2) {
 	  if(cc1) *iter     = (*iter)->CC();
 	  if(cc3) *(iter+2) = (*(iter+2))->CC();
 	}
      }
      else {
 	m2 = (*iter)    ->mass();
 	m3 = (*(iter+1))->mass();
	id2 = abs((*iter)    ->id());
	id3 = abs((*(iter+1))->id());
 	if(cc3) {
	  if(cc1) *iter     = (*iter)->CC();
	  if(cc2) *(iter+1) = (*(iter+1))->CC();
	}
      }
      if(id2==ParticleID::Wplus&&(m1-m3>=0.*MeV && m1-m3<=_masscut))      iter+=3;
      else if(id3==ParticleID::Wplus&&(m1-m2>=0.*MeV && m1-m2<=_masscut)) iter+=3;
      else decaylist.erase(iter,iter+3);
    }
    if(decaylist.size() > 0) createDecayer(vert,ilist,iv);
  }
  return decaylist;
}

void WeakCurrentDecayConstructor::createDecayer(const VertexBasePtr vert,
						unsigned int icol,
						unsigned int ivert) {
  if(_theExistingDecayers[ivert][icol].empty()) {
    string name;
    switch(vert->getName()) {
    case FFV : 
      name = "FFVCurrentDecayer";
      break;
    default : throw NBodyDecayConstructorError() << "Cannot find appropriate "
 						 << "vertex to create "
 						 << "decayer\n";
    }
//     case VVS : {
//       if( icol == 0 || icol == 1) 
// 	name = "VVSDecayer";
//       else 
// 	name = "SVVDecayer";
//     } 
//       break;
//     case GeneralSVV : {
//       if(icol == 0) 
// 	name = "SVVLoopDecayer";
//     }
//       break;
//     case VSS : {
//       if(icol == 0)
// 	name = "VSSDecayer";
//       else 
// 	name = "SSVDecayer";
//     }
//       break;
//     case VVT : {
//       if(icol == 2)
// 	name = "TVVDecayer";
//     }
    ostringstream fullname;
    fullname << "/Defaults/Decays/" << name << "_" 
 	     << ivert << "_" << icol;
    string classname = "Herwig++::" + name;
    ostringstream cut;
    cut << _masscut/GeV;
    for(unsigned int ix=0;ix<_part1.size();++ix) {
      GeneralCurrentDecayerPtr decayer;
      ostringstream fullname2;
      fullname2 << fullname.str() << "_" << ix;
      if(_theExistingDecayers[ivert][icol].find(_current[ix])==
	 _theExistingDecayers[ivert][icol].end()) {
	decayer = dynamic_ptr_cast<GeneralCurrentDecayerPtr>
	  (generator()->preinitCreate(classname,fullname2.str()));
	string msg = generator()->preinitInterface(decayer, "DecayVertex", 
						   "set", vert->fullName());
	if(msg.find("Error:") != string::npos)
	  throw NBodyDecayConstructorError() 
	    << "WeakCurrentDecayConstructor::createDecayer - An error occurred while "
	    << "setting the vertex for " << decayer->fullName()
	    << " - " << msg
	    << Exception::abortnow;
	msg = generator()->preinitInterface(decayer, "Current","set",
					    _current[ix]->fullName());
	if(msg.find("Error:") != string::npos)
	  throw NBodyDecayConstructorError() 
	    << "WeakCurrentDecayConstructor::createDecayer - An error occurred while "
	    << "setting the current for " << decayer->fullName()
	    << " - " << msg
	    << Exception::abortnow;
	msg = generator()->preinitInterface(decayer, "MaximumMass","set",cut.str());
	if(msg.find("Error:") != string::npos)
	  throw NBodyDecayConstructorError() 
	    << "WeakCurrentDecayConstructor::createDecayer - An error occurred while "
	    << "setting the cut-off for " << decayer->fullName()
	    << " - " << msg
	    << Exception::abortnow;
	if(_init) initializeDecayers(fullname2.str());
	decayer->init();
	_theExistingDecayers[ivert][icol][_current[ix]]=decayer;
      }
    }
  }
  else 
    return;
}

void WeakCurrentDecayConstructor::
createDecayMode(PDPtr inpart, const PDVector & decays,
		map<WeakDecayCurrentPtr,GeneralCurrentDecayerPtr> decayers) {
  if(decays.empty()) {
    throw NBodyDecayConstructorError() 
      << "WeakCurrentDecayConstructor::createDecayMode - No decayers\n"
      << Exception::runerror;
  }
  // the partial widths
  Energy totalWidth(0.*MeV);
  vector<vector<Energy> > pWidths(decays.size()/3);
  vector<vector<string> > tags(decays.size()/3);
  vector<vector<WeakDecayCurrentPtr> > currents(decays.size()/3);
  PDVector particles(3);
  if(inpart->CC()) inpart = inpart->CC();
  inpart->stable(false);
  particles[0] = inpart;
  string dmtag,dmtagb;
  bool Wplus;
  for(unsigned int ix = 0; ix < decays.size(); ix += 3) {
    if(decays[ix]->id() == inpart->id()) {
      particles[1] = decays[ix+1];
      particles[2] = decays[ix+2];
    }
    else if(decays[ix+1]->id() == inpart->id()) {
      particles[1] = decays[ix];
      particles[2] = decays[ix+2];
    }
    else {
      particles[1] = decays[ix];
      particles[2] = decays[ix+1];
    }
    if(abs(particles[1]->id())==ParticleID::Wplus) swap(particles[1],particles[2]);
    Wplus=particles[2]->id()==ParticleID::Wplus;
    particles.resize(2);
    dmtag = particles[0]->PDGName() + "->" + particles[1]->PDGName();
    for(unsigned int iy=0;iy<_current.size();++iy) {
       particles.resize(2);
      dmtagb=dmtag;
      vector<tPDPtr> wprod;
      if(_part1[iy]!=0) wprod.push_back(getParticleData(_part1[iy]));
      if(_part2[iy]!=0) wprod.push_back(getParticleData(_part2[iy]));
      if(_part3[iy]!=0) wprod.push_back(getParticleData(_part3[iy]));
      if(_part4[iy]!=0) wprod.push_back(getParticleData(_part4[iy]));
      if(_part5[iy]!=0) wprod.push_back(getParticleData(_part5[iy]));
      int icharge=0;
      for(unsigned int iz=0;iz<wprod.size();++iz) icharge+=wprod[iz]->iCharge();
      bool cc = (Wplus&&icharge==-3)||(!Wplus&&icharge==3);
      for(unsigned int iz=0;iz<wprod.size();++iz) {
 	if(cc&&wprod[iz]->CC())  wprod[iz]=wprod[iz]->CC();
 	dmtagb+="," + wprod[iz]->      PDGName();
      }
      dmtagb += ";";
      pWidths[ix/3].push_back(decayers.find(_current[iy])
			      ->second->partialWidth(inpart,particles[1],wprod));
      tags[ix/3].push_back(dmtagb);
      currents[ix/3].push_back(_current[iy]);
      totalWidth += pWidths[ix/3][iy];
    }
  }
  for(unsigned int ix=0;ix<tags.size();++ix) {
    for(unsigned int iy=0;iy<tags[ix].size();++iy) {
      double tbr = pWidths[ix][iy]/totalWidth;
      tDMPtr thedm= generator()->findDecayMode(tags[ix][iy]);
      if( !thedm ) {
	tDMPtr ndm = generator()->preinitCreateDecayMode(tags[ix][iy]);
	if(ndm) {
	  generator()->preinitInterface(ndm, "Decayer", "set",
					decayers.find(currents[ix][iy])
					->second->fullName());
	  generator()->preinitInterface(ndm, "OnOff", "set", "1");
	  ostringstream br;
	  br << tbr;
	  if(!br)
	    throw NBodyDecayConstructorError()
	      << "Error with branching ratio stream. "
	      << "Branching ratio set to zero for decay mode "
	      << tags[ix][iy]
	      << Exception::warning;
	  else {
	    generator()->preinitInterface(ndm, "BranchingRatio",
					  "set", br.str());
	  }
	}
	else 
	  throw NBodyDecayConstructorError() 
	    << "WeakCurrentDecayConstructor::createDecayMode - Needed to create "
	    << "new decaymode but one could not be created for the tag " 
	    << tags[ix][iy] 
	    << Exception::warning;
      }
      else {
	string::size_type idx = (thedm->decayer()->fullName()).find("Mambo");
	if(idx != string::npos)
	  generator()->preinitInterface(thedm, "Decayer", "set", 
					decayers.find(currents[ix][iy])
					->second->fullName());
      }
    }
  }
}

void WeakCurrentDecayConstructor::initializeDecayers(string fullname) const {
  ostringstream value;
  value << _init;
  string msg = generator()->preinitInterface(fullname, "Initialize", "set",
					     value.str());
  value.str("");
  value << _iteration;
  msg=generator()->preinitInterface(fullname, "Iteration", "set",
				value.str());
  value.str("");
  value << _points;
  msg=generator()->preinitInterface(fullname, "Points", "set",
				value.str());
}
