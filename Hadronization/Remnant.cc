// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Remnant class.
//

#include "Remnant.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/VSelector.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Remnant.tcc"
#endif


using namespace Herwig;
using namespace ThePEG;

ClassDescription<Remnant> Remnant::initRemnant;
// Definition of the static class description member.

void Remnant::Init() {

  static ClassDocumentation<Remnant> documentation
    ("There is no documentation for the Remnant class");

}

Remnant::Remnant(tcEventPDPtr x) :Particle(x)
{}

Remnant::Remnant(PartonBinInstance & pb,const LorentzMomentum & p) : 
Particle(CurrentGenerator::current().getParticleData(ExtraParticleID::Remnant))
{ 
  // set the momentum of the remnant
  set5Momentum(p);
  // id of the particle
  _parent=pb.particleData();
  int pid(_parent->id());
  // get the valence flavours
  _sign = pid < 0 ? -1 : 1;
  // beam particle is a baryon
  if(BaryonMatcher::Check(*(pb.particleData())))
    {
      // get the valence flavours
      _valence.resize(3);
      _valence[0]=(abs(pid)/1000)%10;
      _valence[1]=(abs(pid)/100)%10;
      _valence[2]=(abs(pid)/10)%10;
    }
  // beam particle is a meson
  else if(MesonMatcher::Check(*(pb.particleData())))
    {throw Exception() << "Meson requested in Remant::Remnant() but not implemented "
		       << Exception::runerror;}
  // unknown type of beam particle
  else
    {throw Exception() << " requested in Remant::Remnant() but not implemented "
			<< Exception::runerror;}
  // work out the flavours of the remnants
  obtainConstituents(pb.partonData()->id());
}

PPtr Remnant::clone() const {
  return ptr_new<PPtr>(*this);
}

PPtr Remnant::fullclone() const {
  return clone();
}

void Remnant::obtainConstituents(int extracted)
{
  // set the code of the extracted particle
  _extracted=extracted;
  // construct the remnant
  _constituents.resize(0);
  // copy of the valence partons to construct the remnant
  vector<int> vtemp(_valence);
  // see if the parton is one of the valence ones
  vector<int>::iterator v=find(vtemp.begin(),vtemp.end(),_sign*_extracted);
  // if it is
  bool isvalence(false);
  if(v!=vtemp.end())
    {
      vtemp.erase(v);
      isvalence=true;
    }
  // if valence then the remnant is a diquark
  if(isvalence)
    {
      // this is the spin 1 diquark
      int idqr = 1000*max(vtemp[0],vtemp[1])+100*min(vtemp[0],vtemp[1])+3;
      // if flavours the same could be spin-0 (makes no difference in Hw++)
      if(vtemp[0]!=vtemp[1] && UseRandom::rnd() < 0.25) idqr-=2;
      _constituents.push_back(CurrentGenerator::current().getParticleData(_sign*idqr)->produceParticle());
    }
  // otherwise all constituents
  else
    {
      // obtain the possible quarks and diquarks for the valence bit
      VSelector< pair< int, int > > valenceselector;
      int iq1,iq2,iq3;
      for(iq1 = 0; iq1 < 3; iq1++)
	{
	  iq2 = (iq1+1)%3;
	  iq3 = (iq2+1)%3;
	  // This is the id of the diquark (spin 1) that accompanies iq1
	  int idq = 1000*max(vtemp[iq2], vtemp[iq3]) +
	    100*min(vtemp[iq2], vtemp[iq3]) + 3;
	  valenceselector.insert(3.0, make_pair(vtemp[iq1], idq));
	  if(vtemp[iq2] == vtemp[iq3]) continue;
	  // If they are different, we have spin 0 combination too
	  valenceselector.insert(1.0, make_pair(vtemp[iq1], idq-2));
	}
      // select a quark-diquark pair and add to remnant
      pair<int,int> rr = valenceselector.select(UseRandom::current());
      _constituents.push_back(CurrentGenerator::current().getParticleData(rr.first *_sign)->produceParticle());
      _constituents.push_back(CurrentGenerator::current().getParticleData(rr.second*_sign)->produceParticle());
      // if we extracted a sea quark/antiquark then we to add the antiparticle
      // as well
      if(_extracted!=ParticleID::g)
	{_constituents.push_back(CurrentGenerator::current().getParticleData(-_extracted)->produceParticle());}
    }
}

void Remnant::regenerate(tPPtr extracted,Lorentz5Momentum ptotal)
{
  // change the momentum
  set5Momentum(ptotal);
  // change the constituents
  obtainConstituents(extracted->id());
  // remake the colour connections
  // remove old colour connection
  if(this->colourLine())
    this->colourLine()->removeColoured(this,false);
  if(this->antiColourLine())
    this->antiColourLine()->removeColoured(this,true);
  // make the new colour lines
  if(extracted->colourLine())
    extracted->colourLine()->addColoured(this,true);
  if(extracted->antiColourLine())
    extracted->antiColourLine()->addColoured(this,false);
}

void Remnant::createRemnant(tStepPtr pstep)
{
  cout << "testing in create remnant\n";
  for(unsigned int ix=0;ix<_valence.size();++ix)
    {cout << "testing valence " << ix << " " << _valence[ix] << '\n';}
  for(unsigned int ix=0;ix<_constituents.size();++ix)
    {cout << "testing consitituents " << _constituents[ix]->PDGName() << '\n';}
  // if only one constituent just add it
  if(_constituents.size()==1)
    {
      _constituents[0]->set5Momentum(momentum());
      // set up the colours
      if(this->colourLine()    )
	this->colourLine()->addColoured(_constituents[0],true);
      if(this->antiColourLine())
	this->antiColourLine()->addColoured(_constituents[0],true);
      this->addChild(_constituents[0]);
      pstep->addDecayProduct(_constituents[0]);
      return;
    }
  // if two constituents
  else if(_constituents.size()==2)
    {
      cerr << "Remnant::createRemnant() testing two constituents " << _extracted << '\n';
      exit(1);
    }
  else if(_constituents.size()==3)
    {
      cerr << "Remnant::createRemnant() testing three constituents " << _extracted << '\n';
      // first we need a forced splitting of the 

      exit(1);
    }
  else
    {
      cerr << "Remnant::createRemnant() testing #constituents != 1,2 or 3 " << _extracted << '\n';
      exit(1);
    }
}
