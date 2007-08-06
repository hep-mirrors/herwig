// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericWidthGenerator class.
//

#include "GenericWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "TwoBodyAllOnCalculator.h"
#include "OneOffShellCalculator.h"
#include "TwoOffShellCalculator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/Repository.h"

using namespace Herwig;

void GenericWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << _theParticle << ounit(_mass,GeV) << _prefactor << _MEtype << _MEcode
     << ounit(_MEmass1,GeV) << ounit(_MEmass2,GeV) << _MEcoupling << _modeon
     << ounit(_intermasses,GeV) << ounit(_interwidths,GeV) 
     << _noofentries << _initialize << _BRnorm
     << _npoints << _decaymodes << _decaytags << ounit(_minmass,GeV) 
     << _BRminimum << _intorder;
}

void GenericWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _theParticle >> iunit(_mass,GeV) >> _prefactor >> _MEtype >> _MEcode 
     >> iunit(_MEmass1,GeV) >> iunit(_MEmass2,GeV) >> _MEcoupling >>_modeon
     >> iunit(_intermasses,GeV) >> iunit(_interwidths,GeV)
     >> _noofentries >> _initialize >> _BRnorm
     >> _npoints >> _decaymodes >> _decaytags >> iunit(_minmass,GeV)
     >> _BRminimum >> _intorder;
}

ClassDescription<GenericWidthGenerator> GenericWidthGenerator::initGenericWidthGenerator;
// Definition of the static class description member.

void GenericWidthGenerator::Init() {

  static ClassDocumentation<GenericWidthGenerator> documentation
    ("The GenericWidthGenerator class is the base class for running widths");

  static Reference<GenericWidthGenerator,ParticleData> interfaceParticle
    ("Particle",
     "The particle for which this is the width generator",
     &GenericWidthGenerator::_theParticle, false, false, true, false, false);

  static Switch<GenericWidthGenerator,bool> interfaceInitialize
    ("Initialize",
     "Initialize the width using the particle data object",
     &GenericWidthGenerator::_initialize, false, false, false);
  static SwitchOption interfaceInitializeInitialization
    (interfaceInitialize,
     "Initialization",
     "Do the initialization",
     true);
  static SwitchOption interfaceInitializeNoInitialization
    (interfaceInitialize,
     "NoInitialization",
     "Don't do the initalization",
     false);

  static ParVector<GenericWidthGenerator,int> interfacemetype
    ("MEtype",
     "The type of matrix element either 2-body from this class or higher from"
     " class inheriting from this",
     &GenericWidthGenerator::_MEtype,
     0, 0, 0, 0, 3, false, false, true);

  static ParVector<GenericWidthGenerator,int> interfacemecode
    ("MEcode",
     "The code of matrix element either 2-body from this class or higher from"
     " class inheriting from this",
     &GenericWidthGenerator::_MEcode,
     0, 0, 0, -1, 200, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceMinimumMasses
    ("MinimumMasses",
     "The minimum mass of the decay products",
     &GenericWidthGenerator::_minmass,
     GeV, 0, 0*GeV, 0*GeV,  1.E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,double> interfaceMEcoupling
    ("MEcoupling",
     "The coupling for a given ME",
     &GenericWidthGenerator::_MEcoupling,
     0, 0, 0, 0, 1.E12, false, false, true);

  static ParVector<GenericWidthGenerator,bool> interfaceModeOn
    ("ModeOn",
     "Is this mode included in the total width calculation",
     &GenericWidthGenerator::_modeon,
     0, 0, 0, 0, 1, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceMEmass1
    ("MEmass1",
     "The mass for first particle in a two body mode",
     &GenericWidthGenerator::_MEmass1,
     GeV, 0, 0*GeV, 0*GeV,  1.E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceMEmass2
    ("MEmass2",
     "The mass for second particle in a two body mode",
     &GenericWidthGenerator::_MEmass2,
     GeV, 0, 0*GeV, 0*GeV,  1.E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceInterpolationMasses
    ("InterpolationMasses",
     "The masses for interpolation table",
     &GenericWidthGenerator::_intermasses,
     GeV, 0, 0*GeV, 0*GeV,  1E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceInterpolationWidths
    ("InterpolationWidths",
     "The widths for interpolation table",
     &GenericWidthGenerator::_interwidths,
     GeV, 0, 0*GeV, 0*GeV,  1E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,int> interfacenoofenteries
    ("NumberofEntries",
     "The number of entries in the table after this mode",
     &GenericWidthGenerator::_noofentries,
     0, 0, 0, 0, 100000000, false, false, true);

  static Switch<GenericWidthGenerator,bool> interfaceBRNormalize
    ("BRNormalize",
     "Normalize the partial widths so that they have the value BR*Total Width"
     " for an on-shell particle",
     &GenericWidthGenerator::_BRnorm, false, false, false);
  static SwitchOption interfaceBRNormalizeNormalize
    (interfaceBRNormalize,
     "Normalize",
     "Perform the normalization",
     true);
  static SwitchOption interfaceBRNormalizeNoNormalisation
    (interfaceBRNormalize,
     "NoNormalisation",
     "Do not perform the normalization",
     false);

  static Parameter<GenericWidthGenerator,double> interfaceBRMinimum
    ("BRMinimum",
     "Minimum branching ratio for inclusion in the running width calculation.",
     &GenericWidthGenerator::_BRminimum, 0.01, 0.0, 1.0,
     false, false, true);

  static Parameter<GenericWidthGenerator,double> interfacePrefactor
    ("Prefactor",
     "The prefactor to get the correct on-shell width",
     &GenericWidthGenerator::_prefactor, 1.0, 0., 1000.,
     false, false, false);

  static Parameter<GenericWidthGenerator,int> interfacePoints
    ("Points",
     "Number of points to use for interpolation tables when needed",
     &GenericWidthGenerator::_npoints, 50, 5, 1000,
     false, false, true);

  static ParVector<GenericWidthGenerator,string> interfaceDecayModes
    ("DecayModes",
     "The tags for the decay modes used in the width generator",
     &GenericWidthGenerator::_decaytags, -1, "", "", "",
     false, false, Interface::nolimits);

  static Parameter<GenericWidthGenerator,unsigned int> interfaceInterpolationOrder
    ("InterpolationOrder",
     "The interpolation order for the tables",
     &GenericWidthGenerator::_intorder, 1, 1, 5,
     false, false, Interface::limited);

}

Energy GenericWidthGenerator::width(const ParticleData &, Energy m) const {
  Energy gamma= Energy();
  for(unsigned int ix =0;ix<_MEcoupling.size();++ix) {
    if(_modeon[ix]) gamma +=partialWidth(ix,m);
  }
  return gamma*_prefactor;
}

void GenericWidthGenerator::doinit() throw(InitException) {
  WidthGenerator::doinit();
  // make sure the particle data object was initialized
  _theParticle->init();
  tDecayIntegratorPtr decayer;
  // mass of the decaying particle
  _mass = _theParticle->mass();
  if(_initialize) {
    // the initial prefactor
    _prefactor=1.;
    // resize all the storage vectors
    _MEtype.clear();
    _MEcode.clear();
    _MEmass1.clear();
    _MEmass2.clear();
    _MEcoupling.clear(); 
    _modeon.clear();
    _minmass.clear();
    _intermasses.clear();
    _interwidths.clear();
    _noofentries.clear();
    _decaytags.clear();
    // integrators that we may need
    WidthCalculatorBasePtr widthptr;
    // get the list of decay modes as a decay selector
    DecayMap modes=_theParticle->decaySelector();
    DecayMap::const_iterator start=modes.begin();
    DecayMap::const_iterator end=modes.end();
    tPDPtr part1,part2;
    tGenericMassGeneratorPtr massgen1,massgen2;
    // loop over the decay modes to get the partial widths
    for(;start!=end;++start) {
      // the decay mode
      tcDMPtr mode=(*start).second;      
      _decaymodes.push_back(const_ptr_cast<DMPtr>(mode));
      _decaytags.push_back(_decaymodes.back()->tag());
      ParticleMSet::const_iterator pit(mode->products().begin());
      // minimum mass for the decaymode
      Energy minmass = Energy();
      for(;pit!=mode->products().end();++pit) {
	(**pit).init();
	minmass+=(**pit).massMin();
      }
      _minmass.push_back(minmass);
      pit=mode->products().begin();
      // its decayer
      decayer=dynamic_ptr_cast<tDecayIntegratorPtr>(mode->decayer());
      if(decayer) decayer->init();
      // if there's no decayer then set the partial width to the br times the
      // on-shell value
      if(!decayer) {
	_MEtype.push_back(0);
	_MEcode.push_back(0);
	_MEcoupling.push_back(mode->brat());
	_MEmass1.push_back(0.*GeV);
	_MEmass2.push_back(0.*GeV);
	_noofentries.push_back(_intermasses.size());
	_modeon.push_back(mode->brat()>_BRminimum);
	setupMode(mode,decayer,_MEtype.size()-1);
      }
      else if(mode->products().size()==2) {
	// the outgoing particles
	ParticleMSet::const_iterator pit = mode->products().begin();
	part1=*pit;++pit;
	part2=*pit;
	// mass generators
	if(part1->massGenerator())
	  massgen1=dynamic_ptr_cast<tGenericMassGeneratorPtr>(part1->massGenerator());
	else
	  massgen1=tGenericMassGeneratorPtr();
	if(part2->massGenerator())
	  massgen2=dynamic_ptr_cast<tGenericMassGeneratorPtr>(part2->massGenerator());
	else
	  massgen2=tGenericMassGeneratorPtr();
	if(massgen1) massgen1->init();
	if(massgen2) massgen2->init();
	double coupling(0.);
	int mecode(-1);
	bool order(decayer->twoBodyMEcode(*mode,mecode,coupling));
	_MEcode.push_back(mecode);
	_MEcoupling.push_back(coupling);
	_modeon.push_back(mode->brat()>_BRminimum);
	if(order) {
	  _MEmass1.push_back(part1->mass());
	  _MEmass2.push_back(part2->mass());
	}
	else {
	  _MEmass1.push_back(part2->mass());
	  _MEmass2.push_back(part1->mass());
	}
	// perform setup in the inheriting class
	setupMode(mode,decayer,_MEcode.size()-1);
	// both particles on shell
	if(!massgen1&&!massgen2) {
	  _MEtype.push_back(1);
	  _noofentries.push_back(_intermasses.size());
	  if(_BRnorm) {
	    if(_mass>_MEmass1[_MEtype.size()-1]+_MEmass2[_MEtype.size()-1]) {
	      Energy gamma(partial2BodyWidth(_MEtype.size()-1,_mass));
	      double ratio(mode->brat()*mode->parent()->width()/gamma);
	      ratio=sqrt(ratio);
	      _MEcoupling.back() *=ratio;
	    }
	  }
	}
	else {
	  // one off-shell particle
	  if(!massgen1||!massgen2) {
	    // create the width calculator
	    tGenericWidthGeneratorPtr 
	      ttthis(const_ptr_cast<tGenericWidthGeneratorPtr>(this));
	    WidthCalculatorBasePtr twobody
	      (new_ptr(TwoBodyAllOnCalculator(ttthis,_MEcode.size()-1,
					      _MEmass1[_MEcode.size()-1],
					      _MEmass2[_MEcode.size()-1])));
	    int ioff = ((part1->massGenerator()&&!order)||
			(part2->massGenerator()&&order)) ? 2 : 1;
	    if(massgen1)
	      widthptr=new_ptr(OneOffShellCalculator(ioff,twobody,massgen1,0.*GeV));
	    else
	      widthptr=new_ptr(OneOffShellCalculator(ioff,twobody,massgen2,0.*GeV));
	  }
	  else {
	    int ioff   = order ? 1 : 2;
	    int iother = order ? 2 : 1;
	    // create the width calculator
	    tGenericWidthGeneratorPtr 
	      ttthis(const_ptr_cast<tGenericWidthGeneratorPtr>(this));
	    // this is the both on-shell case
	    WidthCalculatorBasePtr twobody
	      (new_ptr(TwoBodyAllOnCalculator(ttthis,_MEcode.size()-1,
					      _MEmass1[_MEcode.size()-1],
					      _MEmass2[_MEcode.size()-1])));
	    // this is the first off-shell
	    WidthCalculatorBasePtr widthptr2=
	      new_ptr(OneOffShellCalculator(ioff,twobody,massgen1,0*GeV));
	    widthptr=new_ptr(TwoOffShellCalculator(iother,widthptr2,massgen2,
						   0*GeV,massgen1->lowerLimit()));
	  }
	  // set up the interpolation table
	  Energy test(part1->massMin()+part2->massMin());
	  Energy min(max(_theParticle->massMin(),test)),upp(_theParticle->massMax());
	  Energy step((upp-min)/(_npoints-1));
	  Energy moff(min);
	  Energy2 moff2;
	  // additional points to improve the interpolation
	  if(min==test) {
	    _intermasses.push_back(moff-2.*step);_interwidths.push_back(0*GeV);
	    _intermasses.push_back(moff-   step);_interwidths.push_back(0*GeV);
	    _intermasses.push_back(moff        );_interwidths.push_back(0*GeV);
	    double fact(exp(0.1*log(1.+step/moff)));
	    for(unsigned int ix=0;ix<10;++ix) {
	      moff*=fact;
	      moff2=sqr(moff);
	      _intermasses.push_back(moff);
	      _interwidths.push_back(widthptr->partialWidth(moff2));
	    }
	    moff+=step;
	  }
	  else if(test>min-2.*step) {
	    _intermasses.push_back(moff-2.*step);_interwidths.push_back(0*GeV);
	    _intermasses.push_back(test        );_interwidths.push_back(0*GeV);
	  }
	  else {
	    _intermasses.push_back(moff-2.*step);
	    _interwidths.push_back(widthptr->partialWidth((moff-2.*step)*
							  (moff-2.*step)));
	    _intermasses.push_back(moff-   step);
	    _interwidths.push_back(widthptr->partialWidth((moff-   step)*
							  (moff-   step)));
	  }
	  for(; moff<upp+2.5*step;moff+=step) {
	    moff2=moff*moff;
	    _intermasses.push_back(moff);
	    _interwidths.push_back(widthptr->partialWidth(moff2));
	  }
	  coupling=1.;
	  if(_BRnorm) {
	    double ratio(1.);
	    if((massgen1&&massgen2&&
		_mass>massgen1->lowerLimit()+massgen2->lowerLimit())||
	       (massgen1&&!massgen2&&
		_mass>massgen1->lowerLimit()+part2->mass())||
	       (massgen2&&!massgen1&&
		_mass>massgen2->lowerLimit()+part1->mass())||
	       (!massgen1&&!massgen2&&
		_mass>part1->mass()+part2->mass())) {
	      Energy gamma(widthptr->partialWidth(_mass*_mass));
	      ratio=mode->brat()*mode->parent()->width()/gamma;
	    }
	    _MEcoupling.back()=ratio;
	  }
	  _MEtype.push_back(2);
	  _MEcode.back()=0;
	  unsigned int ix=0;
	  if(_MEtype.size()>1){ix=_noofentries[_MEtype.size()-2];}
	  _noofentries.push_back(_intermasses.size());
	  _interpolators.resize(_MEtype.size());
	  // get the vectors we will need
	  vector<Energy>::iterator istart= _intermasses.begin();
	  if(_MEtype.size()>1){istart+=_noofentries[_MEtype.size()-2];}
	  vector<Energy>::iterator iend=_intermasses.end();
	  vector<Energy> masses(istart,iend);

	  istart= _interwidths.begin();
	  if(_MEtype.size()>1){istart+=_noofentries[_MEtype.size()-2];}
	  iend=_interwidths.end();
	  vector<Energy> widths(istart,iend);
	  _interpolators.back() = make_InterpolatorPtr(widths,masses,_intorder);
	}
      }
      // higher multiplicities
      else {
	setupMode(mode,decayer,_MEcode.size());
	widthptr=decayer->threeBodyMEIntegrator(*mode);
	if(!widthptr) {
	  _MEtype.push_back(0);
	  _MEcode.push_back(0);
	  _MEcoupling.push_back(mode->brat());
	  _MEmass1.push_back(0.*GeV);
	  _MEmass2.push_back(0.*GeV);
	  _noofentries.push_back(_intermasses.size());
	  _modeon.push_back(mode->brat()>_BRminimum);
	}
	else {
	  Energy step((_theParticle->widthUpCut()+_theParticle->widthLoCut())/
		      (_npoints-1));
	  Energy moff(_theParticle->massMin()),upp(_theParticle->massMax());
	  for( ; moff<upp+0.5*step;moff+=step) {
	    Energy2 moff2=sqr(moff);
	    Energy wtemp=widthptr->partialWidth(moff2);
	    _intermasses.push_back(moff);
	    _interwidths.push_back(wtemp);
	  }
	  double coupling(1.);
	  if(_BRnorm) {
	    Energy gamma = Energy();
	    gamma=widthptr->partialWidth(_mass*_mass);
	    double ratio(mode->brat()*mode->parent()->width()/gamma);
	    coupling *=ratio;
	  }
	  _MEtype.push_back(2);
	  _MEcode.push_back(0);
	  _MEcoupling.push_back(coupling);
	  _MEmass1.push_back(0*GeV);
	  _MEmass2.push_back(0*GeV);
	  _modeon.push_back(mode->brat()>_BRminimum);
	  unsigned int ix=0;
	  if(_MEtype.size()>1){ix=_noofentries[_MEtype.size()-2];}
	  _noofentries.push_back(_intermasses.size());
	  _interpolators.resize(_MEtype.size());
	  // get the vectors we will need
	  vector<Energy>::iterator istart( _intermasses.begin()),
	    iend(_intermasses.end());
	  if(_MEtype.size()>1){istart+=_noofentries[_MEtype.size()-2];}
	  vector<Energy> masses(istart,iend);
	  
	  istart= _interwidths.begin();
	  if(_MEtype.size()>1){istart+=_noofentries[_MEtype.size()-2];}
	  iend=_interwidths.end();
	  vector<Energy> widths(istart,iend);
	  _interpolators.back() = make_InterpolatorPtr(widths,masses,_intorder);
	}
      }
    }
    // now check the overall normalisation of the running width
    Energy gamma = width(*_theParticle,_mass);
    if(gamma>Energy()) _prefactor = _theParticle->width()/gamma;
    // output the info so it can be read back in
  }
  else {
    // get the decay modes from the tags
    if(_decaytags.size()!=0) {
      _decaymodes.clear();
      for(unsigned int ix=0;ix<_decaytags.size();++ix) {
	_decaymodes.push_back(Repository::findDecayMode(_decaytags[ix]));
      }
    }
    // otherwise just use the modes from the selector
    else {
      DecayMap modes(_theParticle->decaySelector());
      DecayMap::const_iterator start(modes.begin()),end(modes.end());
      tcDMPtr mode;
      for(;start!=end;++start) {
	mode=(*start).second;      
	_decaymodes.push_back(const_ptr_cast<DMPtr>(mode));
      }
    }
    // set up the interpolators
    setInterpolators();
  }
  // setup the partial widths in the decayers for normalization
  tDecayIntegratorPtr temp;
  for(unsigned int ix=0;ix<_decaymodes.size();++ix) {
    decayer=dynamic_ptr_cast<tDecayIntegratorPtr>(_decaymodes[ix]->decayer());
    if(decayer) {
      decayer->init();
      decayer->setPartialWidth(*_decaymodes[ix],ix);
    }
  }
  // code to output plots
  string fname = CurrentGenerator::current().filename() + 
    string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  Energy step = step=(_theParticle->massMax()-_theParticle->massMin())/100.;
  output << "SET FONT DUPLEX\n";
  output << "TITLE TOP \"Width for " << _theParticle->name() << "\"\n";
  output << "TITLE BOTTOM \"m/GeV\"\n";
  output << "TITLE LEFT \"G/GeV\"\n";
  output << "CASE       \"F    \"\n";
  output << "SET LIMITS X " 
	 << (_theParticle->massMin()-10.*step)/GeV << " " 
	 << _theParticle->massMax()/GeV << "\n";
  Energy upper(0.*GeV);
  for(Energy etest=_theParticle->massMin();etest<_theParticle->massMax();etest+=step) {
    Energy gamma=width(*_theParticle,etest);
    upper = max(gamma,upper);
    output << etest/GeV << "\t" << gamma/GeV << "\n";
  }
  output << "SET LIMITS Y 0. " << upper/GeV << "\n";
  output << "JOIN\n";
  output << (_theParticle->massMin()-9.*step)/GeV << "\t" 
	 <<  upper*(_MEcode.size()+1)/(_MEcode.size()+2)/GeV << "\n";
  output << (_theParticle->massMin()-7.*step)/GeV << "\t" 
	 <<  upper*(_MEcode.size()+1)/(_MEcode.size()+2)/GeV << "\n";
  output << "JOIN\n";
  output << "TITLE DATA " 
	 << (_theParticle->massMin()-6.*step)/GeV << "\t" 
	 <<  upper*(_MEcode.size()+1)/(_MEcode.size()+2)/GeV 
	 << " \"total\"\n";
  for(unsigned int ix=0;ix<_MEcode.size();++ix) {
    for(Energy etest=_theParticle->massMin();etest<_theParticle->massMax();etest+=step) {
      output << etest/GeV << "\t" << partialWidth(ix,etest)*_prefactor/GeV << "\n";
    }
    switch(ix) {
    case 0:  output << "join red\n"    ; break;
    case 1:  output << "join blue\n"   ; break;
    case 2:  output << "join green\n"  ; break;
    case 3:  output << "join yellow\n" ; break;
    case 4:  output << "join magenta\n"; break;
    case 5:  output << "join cyan\n"   ; break;
    case 6:  output << "join dashes\n" ; break;
    case 7:  output << "join dotted\n" ; break;
    case 8:  output << "join dotdash\n"; break;
    default: output << "join daashes space\n";  break;
    }
    output << (_theParticle->massMin()-9.*step)/GeV << "\t" 
	   <<  upper*(_MEcode.size()-ix)/(_MEcode.size()+2)/GeV << "\n";
    output << (_theParticle->massMin()-7.*step)/GeV << "\t" 
	   <<  upper*(_MEcode.size()-ix)/(_MEcode.size()+2)/GeV << "\n"; 
    switch(ix) {
    case 0:  output << "join red\n"    ; break;
    case 1:  output << "join blue\n"   ; break;
    case 2:  output << "join green\n"  ; break;
    case 3:  output << "join yellow\n" ; break;
    case 4:  output << "join magenta\n"; break;
    case 5:  output << "join cyan\n"   ; break;
    case 6:  output << "join dashes\n" ; break;
    case 7:  output << "join dotted\n" ; break;
    case 8:  output << "join dotdash\n"; break;
    default: output << "join daashes space\n";  break;
    }
    output << "TITLE DATA " 
	   << (_theParticle->massMin()-6.*step)/GeV << "\t" 
	   <<  upper*(_MEcode.size()-ix)/(_MEcode.size()+2)/GeV 
	   << " \"" << _decaytags[ix] << "\"\n";
  }
}
 
void GenericWidthGenerator::setInterpolators() {
  // create the interpolators
  _interpolators.resize(_MEtype.size());
  vector<Energy>::iterator estart(_intermasses.begin()),eend;
  vector<Energy>::iterator wstart(_interwidths.begin()),wend;
  vector<Energy> masses,widths;
  for(unsigned int ix=0;ix<_MEtype.size();++ix) {
    eend=_intermasses.begin()+_noofentries[ix];
    wend=_interwidths.begin()+_noofentries[ix];
    if(_MEtype[ix]==2) {
      masses.assign(estart,eend);
      widths.assign(wstart,wend);
      _interpolators[ix]= make_InterpolatorPtr(widths,masses,_intorder);
    }
    estart=eend;
    wstart=wend;
  }
}

void GenericWidthGenerator::dataBaseOutput(ofstream & output, bool header) {
  if(header) output << "update Width_Generators set parameters=\"";
  // prefactor and general switiches
  output << "set " << fullName() << ":Prefactor "   << _prefactor << "\n";
  output << "set " << fullName() << ":BRNormalize " << _BRnorm    << "\n";
  output << "set " << fullName() << ":BRMinimum "   << _BRminimum << "\n";
  output << "set " << fullName() << ":Points "      << _npoints   << "\n";
  output << "set " << fullName() << ":InterpolationOrder " << _intorder << "\n";
  // the type of the matrix element
  for(unsigned int ix=0;ix<_MEtype.size();++ix) {
    output << "insert " << fullName() << ":MEtype " << ix << " " 
	   << _MEtype[ix] << "\n";
  }
  // the code for thew two body matrix elements
  for(unsigned int ix=0;ix<_MEcode.size();++ix) {
    output << "insert " << fullName() << ":MEcode " 
	   << ix << " " << _MEcode[ix] << "\n";
  }
  // the coupling for trhe two body matrix elements
  for(unsigned int ix=0;ix<_MEcoupling.size();++ix) {
    output << "insert " << fullName() << ":MEcoupling " 
	   << ix << " " << _MEcoupling[ix] << "\n";
  }
  // use this mode for the running width
  for(unsigned int ix=0;ix<_modeon.size();++ix) {
    output << "insert " << fullName() << ":ModeOn " 
	   << ix << " " << _modeon[ix] << "\n";
  }
  // first outgoing mass
  for(unsigned int ix=0;ix<_minmass.size();++ix) {
    output << "insert " << fullName() << ":MinimumMasses " 
	   << ix << " " << _minmass[ix]/GeV << "\n";
  }
  // first outgoing mass
  for(unsigned int ix=0;ix<_MEmass1.size();++ix) {
    output << "insert " << fullName() << ":MEmass1 " 
	   << ix << " " << _MEmass1[ix]/GeV << "\n";
  }
  // second outgoing mass
  for(unsigned int ix=0;ix<_MEmass2.size();++ix) {
    output << "insert " << fullName() << ":MEmass2 " 
	   << ix << " " << _MEmass2[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_decaymodes.size();++ix) {
    output << "insert " << fullName() << ":DecayModes "
	   << ix << " " << _decaytags[ix] << " \n";
  }
  // data for the interpolation tables
  std::streamsize curpre=output.precision();
  output.precision(curpre+2);
  for(unsigned int ix=0;ix<_intermasses.size();++ix) {
    output << "insert " << fullName() 
	   << ":InterpolationMasses " 
	   << ix << " " << _intermasses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_interwidths.size();++ix) {
    output << "insert " << fullName() 
	   << ":InterpolationWidths " 
	   << ix << " " << _interwidths[ix]/GeV << "\n";
  }
  output.precision(curpre);
  for(unsigned int ix=0;ix<_noofentries.size();++ix) {
    output << "insert " << fullName() 
	   << ":NumberofEntries " 
	   << ix << " " << _noofentries[ix] << "\n";
  }  
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() 
		    << "\";" << endl;
}

DecayMap GenericWidthGenerator::rate(const Particle & p) {
  Energy scale(p.mass());
  DecayMap dm;
  // use the running widths to generate the branching ratios
  if(_theParticle->variableRatio()) {
    DecayMap newmap;
    Energy width = _theParticle->width();
    for(unsigned int ix=0;ix<_decaymodes.size();++ix) {
      // DGRELL units?
      if(p.id()==_theParticle->id())
	newmap.insert(partialWidth(ix,scale)/width,_decaymodes[ix]);
      else
	newmap.insert(partialWidth(ix,scale)/width,_decaymodes[ix]->CC());
    }
    dm=newmap;
  }
  // if we are not varying the width return the default
  else dm=p.data().decaySelector();
  return dm;
}

void GenericWidthGenerator::setupMode(tcDMPtr, tDecayIntegratorPtr,
				      unsigned int)
{}

Energy GenericWidthGenerator::partialWidth(int imode,Energy q) const {
  if(q<_minmass[imode]) return Energy();
  Energy gamma;
  if(_MEtype[imode]==0) {
    gamma=_MEcoupling[imode]*_theParticle->width();
  }
  else if(_MEtype[imode]==1) {
    gamma=partial2BodyWidth(imode,q);
  }
  else if(_MEtype[imode]==2) {
    gamma=_MEcoupling[imode]*(*_interpolators[imode])(q);
  }
  else {
    throw Exception() << "Unknown type of mode " << _MEtype[imode] 
		      << "in GenericWidthGenerator::partialWidth()"
		      << Exception::runerror;
  }
  return max(gamma,Energy());
}

void GenericWidthGenerator::dofinish() {
  if(_initialize) {
    string fname = CurrentGenerator::current().filename() + 
      string("-") + name() + string(".output");
    ofstream output(fname.c_str());
    dataBaseOutput(output,true);
  }
  WidthGenerator::dofinish();
}

void GenericWidthGenerator::doinitrun() {
  WidthGenerator::doinitrun();
  // set up the interpolators
  setInterpolators();
}
