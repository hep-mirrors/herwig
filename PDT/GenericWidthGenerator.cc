// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericWidthGenerator class.
//

#include "GenericWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "TwoBodyAllOnCalculator.h"
#include "OneOffShellCalculator.h"
#include "TwoOffShellCalculator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GenericWidthGenerator.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
typedef Selector<tDMPtr> DecayMap;

GenericWidthGenerator::~GenericWidthGenerator() {}

void GenericWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << _theParticle << _mass << _prefactor << _MEtype << _MEcode
     << _MEmass1 << _MEmass2 << _MEcoupling << _modeon
     << _intermasses << _interwidths << _noofentries << _initialize << _BRnorm
     << _npoints << _decaymodes << _minmass << _BRminimum;
}

void GenericWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _theParticle >> _mass >> _prefactor >> _MEtype >> _MEcode 
     >> _MEmass1 >> _MEmass2 >> _MEcoupling >>_modeon
     >> _intermasses >> _interwidths >> _noofentries >> _initialize >> _BRnorm
     >> _npoints >> _decaymodes >> _minmass >> _BRminimum;
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
     0, 0, 0, 0,  1.E12, false, false, true);

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
     0, 0, 0, 0,  1.E12, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceMEmass2
    ("MEmass2",
     "The mass for second particle in a two body mode",
     &GenericWidthGenerator::_MEmass2,
     0, 0, 0, 0,  1.E12, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceInterpolationMasses
    ("InterpolationMasses",
     "The masses for interpolation table",
     &GenericWidthGenerator::_intermasses,
     0, 0, 0, 0,  1E12, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceInterpolationWidths
    ("InterpolationWidths",
     "The widths for interpolation table",
     &GenericWidthGenerator::_interwidths,
     0, 0, 0, 0,  1E12, false, false, true);

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

  static RefVector<GenericWidthGenerator,DecayMode> interfaceDecayModes
    ("DecayModes",
     "The decay modes used in the width generator",
     &GenericWidthGenerator::_decaymodes, -1, false, false, true, false, false);

}

bool GenericWidthGenerator::accept(const ParticleData & in) const
{
  if(!_theParticle){return false;}
  return &in==_theParticle||(in.CC()&&in.CC()==_theParticle);
 }

Energy GenericWidthGenerator::width(const ParticleData & in, Energy m) const
{
  Energy gamma(0.);
  for(unsigned int ix =0;ix<_MEcoupling.size();++ix)
    {if(_modeon[ix]){gamma +=partialWidth(ix,m);}}
  return gamma*_prefactor;
}

DecayMap GenericWidthGenerator::rate(const ParticleData & in) const
{return in.decaySelector();}

void GenericWidthGenerator::doinit() throw(InitException) {
  WidthGenerator::doinit();
  // make sure the particle data object was initialized
  _theParticle->init();
  tDecayIntegratorPtr decayer;
  // mass of the decaying particle
  _mass = _theParticle->mass();
  if(_initialize)
    {
      // the initial prefactor
      _prefactor=1.;
      // resize all the storage vectors
      _MEtype.resize(0);_MEcode.resize(0);
      _MEmass1.resize(0);_MEmass2.resize(0);
      _MEcoupling.resize(0); 
      _modeon.resize(0);
      _minmass.resize(0);
      _intermasses.resize(0);_interwidths.resize(0);
      _noofentries.resize(0);_decaymodes.resize(0);
      // integrators that we may need
      WidthCalculatorBasePtr widthptr;
      // get the list of decay modes as a decay selector
      DecayMap modes=_theParticle->decaySelector();
      DecayMap::const_iterator start=modes.begin();
      DecayMap::const_iterator end=modes.end();
      tPDPtr part1,part2;
      tGenericMassGeneratorPtr massgen1,massgen2;
      // loop over the decay modes to get the partial widths
      for(;start!=end;++start)
	{
	  // the decay mode
	  tcDMPtr mode=(*start).second;      
	  _decaymodes.push_back(const_ptr_cast<DMPtr>(mode));
	  ParticleMSet::const_iterator pit(mode->products().begin());
	  // minimum mass for the decaymode
	  Energy minmass(0.);
	  for(;pit!=mode->products().end();++pit)
	    {
	      (**pit).init();
	      minmass+=(**pit).massMin();
	    }
	  _minmass.push_back(minmass);
	  pit=mode->products().begin();
	  // its decayer
	  decayer=dynamic_ptr_cast<tDecayIntegratorPtr>(mode->decayer());
	  if(decayer){decayer->init();}
	  // if there's no decayer then set the partial width to the br times the
	  // on-shell value
	  if(!decayer)
	    {
	      _MEtype.push_back(0);
	      _MEcode.push_back(0);
	      _MEcoupling.push_back(mode->brat());
	      _MEmass1.push_back(0.);
	      _MEmass2.push_back(0.);
	      _noofentries.push_back(_intermasses.size());
	      _modeon.push_back(mode->brat()>_BRminimum);
	      setupMode(mode,decayer,_MEtype.size()-1);
	    }
	  else if(mode->products().size()==2)
	    {
	      // the outgoing particles
	      ParticleMSet::const_iterator pit = mode->products().begin();
	      part1=*pit;++pit;
	      part2=*pit;
	      // mass generators
	      if(part1->massGenerator()){massgen1=
		  dynamic_ptr_cast<tGenericMassGeneratorPtr>(part1->massGenerator());}
	      else{massgen1=tGenericMassGeneratorPtr();}
	      if(part2->massGenerator()){massgen2=
		  dynamic_ptr_cast<tGenericMassGeneratorPtr>(part2->massGenerator());}
	      else{massgen2=tGenericMassGeneratorPtr();}
	      if(massgen1){massgen1->init();}
	      if(massgen2){massgen2->init();}
	      double coupling(0.);
	      int mecode(-1);
	      bool order(decayer->twoBodyMEcode(*mode,mecode,coupling));
	      _MEcode.push_back(mecode);
	      _MEcoupling.push_back(coupling);
	      _modeon.push_back(mode->brat()>_BRminimum);
	      if(order)
		{_MEmass1.push_back(part1->mass());_MEmass2.push_back(part2->mass());}
	      else
		{_MEmass1.push_back(part2->mass());_MEmass2.push_back(part1->mass());}
	      // perform setup in the inheriting class
	      setupMode(mode,decayer,_MEcode.size()-1);
	      // both particles on shell
	      if(!massgen1&&!massgen2)
		{
		  _MEtype.push_back(1);
		  _noofentries.push_back(_intermasses.size());
		  if(_BRnorm)
		    {
		      if(_mass>_MEmass1[_MEtype.size()-1]+_MEmass2[_MEtype.size()-1])
			{
			  Energy gamma(partial2BodyWidth(_MEtype.size()-1,_mass));
			  double ratio(mode->brat()*mode->parent()->width()/gamma);
			  ratio=sqrt(ratio);
			  _MEcoupling.back() *=ratio;
			}
		    }
		}
	      else
		{
		  // one off-shell particle
		  if(!massgen1||!massgen2)
		    {
		      int ioff(1);
		      // create the width calculator
		      tGenericWidthGeneratorPtr 
			ttthis(const_ptr_cast<tGenericWidthGeneratorPtr>(this));
		      WidthCalculatorBasePtr twobody
			(new_ptr(TwoBodyAllOnCalculator(ttthis,_MEcode.size()-1,
						       _MEmass1[_MEcode.size()-1],
							_MEmass2[_MEcode.size()-1])));
		      if((part1->massGenerator()&&!order)||
			 (part2->massGenerator()&&order)){ioff=2;}
		      if(massgen1)
			{widthptr=
			    new_ptr(OneOffShellCalculator(ioff,twobody,massgen1,0.));}
		      else
			{widthptr=
			    new_ptr(OneOffShellCalculator(ioff,twobody,massgen2,0.));}
		    }
		  else 
		    {
		      int ioff(1),iother(2);if(!order){ioff=2;iother=1;}
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
			new_ptr(OneOffShellCalculator(ioff,twobody,massgen1,0.));
		      widthptr=new_ptr(TwoOffShellCalculator(iother,widthptr2,massgen2,
							     0.,massgen1->lowerLimit()));
		    }
		  // set up the interpolation table
		  Energy min(_theParticle->massMin()),upp(_theParticle->massMax());
		  Energy test(part1->massMin()+part2->massMin());
		  if(min<test){min=test;}
		  Energy step((upp-min)/(_npoints-1));
		  Energy moff(min);
		  Energy2 moff2;
		  // additional points to improve the interpolation
		  if(min==test)
		    {
		      _intermasses.push_back(moff-2.*step);_interwidths.push_back(0.);
		      _intermasses.push_back(moff-   step);_interwidths.push_back(0.);
		      _intermasses.push_back(moff        );_interwidths.push_back(0.);
		      double fact(exp(0.1*log(1.+step/moff)));
		      for(unsigned int ix=0;ix<10;++ix)
			{
			  moff*=fact;
			  moff2=moff*moff;
			  _intermasses.push_back(moff);
			  _interwidths.push_back(widthptr->partialWidth(moff2));
			}
		    }
		  else if(test>min-2.*step)
		    {
		      _intermasses.push_back(moff-2.*step);_interwidths.push_back(0.);
		      _intermasses.push_back(test        );_interwidths.push_back(0.);
		    }
		  else
		    {	      
		      _intermasses.push_back(moff-2.*step);
		      _interwidths.push_back(widthptr->partialWidth((moff-2.*step)*
								    (moff-2.*step)));
		      _intermasses.push_back(moff-   step);
		      _interwidths.push_back(widthptr->partialWidth((moff-   step)*
								    (moff-   step)));
		    }
		  for(; moff<upp+2.5*step;moff+=step)
		    {
		      moff2=moff*moff;
		      _intermasses.push_back(moff);
		      _interwidths.push_back(widthptr->partialWidth(moff2));
		    }
		  coupling=1.;
		  if(_BRnorm)
		    {
		      double ratio(1.);
		      if((massgen1&&massgen2&&
			  _mass>massgen1->lowerLimit()+massgen2->lowerLimit())||
			 (massgen1&&!massgen2&&
			  _mass>massgen1->lowerLimit()+part2->mass())||
			 (massgen2&&!massgen1&&
			  _mass>massgen2->lowerLimit()+part1->mass())||
			 (!massgen1&&!massgen2&&
			  _mass>part1->mass()+part2->mass()))
			{
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
		  vector<Energy> masses=vector<Energy>(istart,iend);
		  istart= _interwidths.begin();
		  if(_MEtype.size()>1){istart+=_noofentries[_MEtype.size()-2];}
		  iend=_interwidths.end();
		  vector<Energy> widths=vector<Energy>(istart,iend);
		  _interpolators.back()=new Interpolator(widths,masses,3);
		}
	    }
	  // higher multiplicities
	  else
	    {
	      // perform setup in the inheriting class
	      setupMode(mode,decayer,_MEcode.size());
	      // see how many off-shell particles there are
	      ParticleMSet::const_iterator pit = mode->products().begin();
	      widthptr=decayer->threeBodyMEIntegrator(*mode);
	      Energy step((_theParticle->widthUpCut()+_theParticle->widthLoCut())/
			  (_npoints-1));
	      Energy moff(_theParticle->massMin()),upp(_theParticle->massMax());
	      Energy2 moff2,wtemp;
	      for( ; moff<upp+0.5*step;moff+=step)
		{
		  moff2=moff*moff;
		  wtemp=widthptr->partialWidth(moff2);
		  _intermasses.push_back(moff);
		  _interwidths.push_back(wtemp);
		}
	      double coupling(1.);
	      if(_BRnorm)
		{
		  Energy gamma(0.);
		  gamma=widthptr->partialWidth(_mass*_mass);
		  double ratio(mode->brat()*mode->parent()->width()/gamma);
		  coupling *=ratio;
		}
	      _MEtype.push_back(2);
	      _MEcode.push_back(0);
	      _MEcoupling.push_back(coupling);
	      _MEmass1.push_back(0.);
	      _MEmass2.push_back(0.);
	      _modeon.push_back(mode->brat()>_BRminimum);
	      unsigned int ix=0;
	      if(_MEtype.size()>1){ix=_noofentries[_MEtype.size()-2];}
	      _noofentries.push_back(_intermasses.size());
	      _interpolators.resize(_MEtype.size());
	      // get the vectors we will need
	      vector<Energy>::iterator istart( _intermasses.begin()),
		iend(_intermasses.end());
	      if(_MEtype.size()>1){istart+=_noofentries[_MEtype.size()-2];}
	      vector<Energy> masses=vector<Energy>(istart,iend);
	      istart= _interwidths.begin();
	      if(_MEtype.size()>1){istart+=_noofentries[_MEtype.size()-2];}
	      iend=_interwidths.end();
	      vector<Energy> widths=vector<Energy>(istart,iend);
	      _interpolators.back()=new Interpolator(widths,masses,3);
	    }
	}
      // now check the overall normalisation of the running width
      Energy gamma = width(*_theParticle,_mass);
      if(gamma>0){_prefactor = _theParticle->width()/gamma;}
      // output the info so it can be read back in
    }
  else
    {
      setInterpolators();
      if(_decaymodes.size()==0)
	{
	  DecayMap modes(_theParticle->decaySelector());
	  DecayMap::const_iterator start(modes.begin()),end(modes.end());
	  tcDMPtr mode;
	  for(;start!=end;++start)
	    {
	      mode=(*start).second;      
	      _decaymodes.push_back(const_ptr_cast<DMPtr>(mode));
	      ParticleMSet::const_iterator pit = mode->products().begin();
	    }
	}
    }
  // setup the partial widths in the decayers for normalization
  tDecayIntegratorPtr temp;
  for(unsigned int ix=0;ix<_decaymodes.size();++ix)
    {
      decayer=dynamic_ptr_cast<tDecayIntegratorPtr>(_decaymodes[ix]->decayer());
      if(decayer)
	{
	  decayer->init();
	  decayer->setPartialWidth(*_decaymodes[ix],ix);
	}
    }
}
  
void GenericWidthGenerator::setInterpolators()
{
  // create the interpolators
  _interpolators.resize(_MEtype.size());
  vector<Energy>::iterator estart(_intermasses.begin()),eend;
  vector<Energy>::iterator wstart(_interwidths.begin()),wend;
  vector<Energy> masses,widths;
  for(unsigned int ix=0;ix<_MEtype.size();++ix)
    {
      eend=_intermasses.begin()+_noofentries[ix];
      wend=_interwidths.begin()+_noofentries[ix];
      if(_MEtype[ix]==2)
	{
	  masses=vector<Energy>(estart,eend);
	  widths=vector<Energy>(wstart,wend);
	  _interpolators[ix]=new Interpolator(widths,masses,3);
	}
      estart=eend;
      wstart=wend;
    }
}

void GenericWidthGenerator::dataBaseOutput(ofstream & output, bool header)
{
  if(header){output << "update Width_Generators set parameters=\"";}
  // prefactor and general switiches
  output << "set " << fullName() << ":Prefactor "   << _prefactor << "\n";
  output << "set " << fullName() << ":BRNormalize " << _BRnorm    << "\n";
  output << "set " << fullName() << ":BRMinimum "   << _BRminimum << "\n";
  output << "set " << fullName() << ":Points "      << _npoints   << "\n";
  // the type of the matrix element
  for(unsigned int ix=0;ix<_MEtype.size();++ix)
    {output << "insert " << fullName() << ":MEtype " << ix << " " 
	    << _MEtype[ix] << "\n";}
  // the code for thew two body matrix elements
  for(unsigned int ix=0;ix<_MEcode.size();++ix)
    {output << "insert " << fullName() << ":MEcode " 
	    << ix << " " << _MEcode[ix] << "\n";}
  // the coupling for trhe two body matrix elements
  for(unsigned int ix=0;ix<_MEcoupling.size();++ix)
    {output << "insert " << fullName() << ":MEcoupling " 
	    << ix << " " << _MEcoupling[ix] << "\n";}
  // use this mode for the running width
  for(unsigned int ix=0;ix<_modeon.size();++ix)
    {output << "insert " << fullName() << ":ModeOn " 
	    << ix << " " << _modeon[ix] << "\n";}
  // first outgoing mass
  for(unsigned int ix=0;ix<_minmass.size();++ix)
    {output << "insert " << fullName() << ":MinimumMasses " 
	    << ix << " " << _minmass[ix] << "\n";}
  // first outgoing mass
  for(unsigned int ix=0;ix<_MEmass1.size();++ix)
    {output << "insert " << fullName() << ":MEmass1 " 
	    << ix << " " << _MEmass1[ix] << "\n";}
  // second outgoing mass
  for(unsigned int ix=0;ix<_MEmass2.size();++ix)
    {output << "insert " << fullName() << ":MEmass2 " 
	    << ix << " " << _MEmass2[ix] << "\n";}
  for(unsigned int ix=0;ix<_decaymodes.size();++ix)
    {output << "insert " << fullName() << ":DecayModes "
	    << ix << " " << _decaymodes[ix]->fullName() << " \n";}
  // data for the interpolation tables
  std::streamsize curpre=output.precision();
  output.precision(curpre+2);
  for(unsigned int ix=0;ix<_intermasses.size();++ix)
    {output << "insert " << fullName() 
	    << ":InterpolationMasses " 
	    << ix << " " << _intermasses[ix] << "\n";}
  for(unsigned int ix=0;ix<_interwidths.size();++ix)
    {output << "insert " << fullName() 
	    << ":InterpolationWidths " 
	    << ix << " " << _interwidths[ix] << "\n";}
  output.precision(curpre);
  for(unsigned int ix=0;ix<_noofentries.size();++ix)
    {output << "insert " << fullName() 
	    << ":NumberofEntries " 
	    << ix << " " << _noofentries[ix] << "\n";}  
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

DecayMap GenericWidthGenerator::rate(const Particle & p) {
  Energy scale(p.mass());
  DecayMap dm;
  // use the running widths to generate the branching ratios
  if(_theParticle->variableRatio())
    {
      DecayMap newmap;
      for(unsigned int ix=0;ix<_decaymodes.size();++ix)
	{
	  if(p.id()==_theParticle->id())
	    {newmap.insert(partialWidth(ix,scale),_decaymodes[ix]);}
	  else
	    {newmap.insert(partialWidth(ix,scale),_decaymodes[ix]->CC());}
	}
      dm=newmap;
    }
  // if we are not varying the width return the default
  else
    {dm=p.data().decaySelector();}
  return dm;
} 
  
void GenericWidthGenerator::setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer,
				      unsigned int imode)
{}

Energy GenericWidthGenerator::partialWidth(int imode,Energy q) const
{
  if(q<_minmass[imode]){return 0.;}
  Energy gamma;
  if(_MEtype[imode]==0){gamma=_MEcoupling[imode]*_theParticle->width();}
  else if(_MEtype[imode]==1){gamma=partial2BodyWidth(imode,q);}
  else if(_MEtype[imode]==2){gamma=_MEcoupling[imode]*(*_interpolators[imode])(q);}
  else{cout << "Unknown type of mode " << _MEtype[imode] << endl;gamma=0.;}
  return max(gamma,0.);
}

}
