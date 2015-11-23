// -*- C++ -*-
//
// GenericWidthGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericWidthGenerator class.
//

#include "GenericWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "TwoBodyAllOnCalculator.h"
#include "OneOffShellCalculator.h"
#include "TwoOffShellCalculator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <ctime>

using namespace Herwig;

DescribeClass<GenericWidthGenerator,WidthGenerator>
describeHerwigGenericWidthGenerator("Herwig::GenericWidthGenerator","");
HERWIG_INTERPOLATOR_CLASSDESC(GenericWidthGenerator,Energy,Energy)


void GenericWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << particle_ << ounit(mass_,GeV) << prefactor_ << MEtype_ << MEcode_
     << ounit(MEmass1_,GeV) << ounit(MEmass2_,GeV) << MEcoupling_ << modeOn_
     << ounit(interMasses_,GeV) << ounit(interWidths_,GeV) 
     << noOfEntries_ << initialize_ << output_ << BRnorm_ << twoBodyOnly_
     << npoints_ << decayModes_ << decayTags_ << ounit(minMass_,GeV) 
     << BRminimum_ << intOrder_ << interpolators_;
}

void GenericWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> particle_ >> iunit(mass_,GeV) >> prefactor_ >> MEtype_ >> MEcode_ 
     >> iunit(MEmass1_,GeV) >> iunit(MEmass2_,GeV) >> MEcoupling_ >>modeOn_
     >> iunit(interMasses_,GeV) >> iunit(interWidths_,GeV)
     >> noOfEntries_ >> initialize_ >> output_ >> BRnorm_ >> twoBodyOnly_
     >> npoints_ >> decayModes_ >> decayTags_ >> iunit(minMass_,GeV)
     >> BRminimum_ >> intOrder_ >> interpolators_;
}

void GenericWidthGenerator::setParticle(string p) {
  if ( (particle_ = Repository::GetPtr<tPDPtr>(p)) ) return;
  particle_ = Repository::findParticle(StringUtils::basename(p));
  if ( ! particle_ ) 
    Throw<InterfaceException>() 
      << "Could not set Particle interface "
      << "for the object \"" << name()
      << "\". Particle \"" << StringUtils::basename(p) << "\" not found."
      << Exception::runerror;
}

string GenericWidthGenerator::getParticle() const {
  return particle_ ? particle_->fullName() : "";
}

void GenericWidthGenerator::Init() {

  static ClassDocumentation<GenericWidthGenerator> documentation
    ("The GenericWidthGenerator class is the base class for running widths");

  static Parameter<GenericWidthGenerator,string> interfaceParticle
    ("Particle",
     "The particle for which this is the width generator",
     0, "", true, false,
     &GenericWidthGenerator::setParticle, 
     &GenericWidthGenerator::getParticle);

  static Switch<GenericWidthGenerator,bool> interfaceInitialize
    ("Initialize",
     "Initialize the width using the particle data object",
     &GenericWidthGenerator::initialize_, false, false, false);
  static SwitchOption interfaceInitializeInitialization
    (interfaceInitialize,
     "Yes",
     "Do the initialization",
     true);
  static SwitchOption interfaceInitializeNoInitialization
    (interfaceInitialize,
     "No",
     "Don't do the initalization",
     false);

  static Switch<GenericWidthGenerator,bool> interfaceOutput
    ("Output",
     "Output the setup",
     &GenericWidthGenerator::output_, false, false, false);
  static SwitchOption interfaceOutputYes
    (interfaceOutput,
     "Yes",
     "Output the data",
     true);
  static SwitchOption interfaceOutputNo
    (interfaceOutput,
     "No",
     "Don't output the data",
     false);

  static ParVector<GenericWidthGenerator,int> interfacemetype
    ("MEtype",
     "The type of matrix element either 2-body from this class or higher from"
     " class inheriting from this",
     &GenericWidthGenerator::MEtype_,
     0, 0, 0, 0, 3, false, false, true);

  static ParVector<GenericWidthGenerator,int> interfacemecode
    ("MEcode",
     "The code of matrix element either 2-body from this class or higher from"
     " class inheriting from this",
     &GenericWidthGenerator::MEcode_,
     0, 0, 0, -1, 200, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceMinimumMasses
    ("MinimumMasses",
     "The minimum mass of the decay products",
     &GenericWidthGenerator::minMass_,
     GeV, 0, ZERO, ZERO,  1.E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,double> interfaceMEcoupling
    ("MEcoupling",
     "The coupling for a given ME",
     &GenericWidthGenerator::MEcoupling_,
     0, 0, 0, 0, 1.E12, false, false, true);

  static ParVector<GenericWidthGenerator,bool> interfaceModeOn
    ("ModeOn",
     "Is this mode included in the total width calculation",
     &GenericWidthGenerator::modeOn_,
     0, 0, 0, 0, 1, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceMEmass1
    ("MEmass1",
     "The mass for first particle in a two body mode",
     &GenericWidthGenerator::MEmass1_,
     GeV, 0, ZERO, ZERO,  1.E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceMEmass2
    ("MEmass2",
     "The mass for second particle in a two body mode",
     &GenericWidthGenerator::MEmass2_,
     GeV, 0, ZERO, ZERO,  1.E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceInterpolationMasses
    ("InterpolationMasses",
     "The masses for interpolation table",
     &GenericWidthGenerator::interMasses_,
     GeV, 0, ZERO, ZERO,  1E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,Energy> interfaceInterpolationWidths
    ("InterpolationWidths",
     "The widths for interpolation table",
     &GenericWidthGenerator::interWidths_,
     GeV, 0, ZERO, ZERO,  1E12*GeV, false, false, true);

  static ParVector<GenericWidthGenerator,int> interfacenoofenteries
    ("NumberofEntries",
     "The number of entries in the table after this mode",
     &GenericWidthGenerator::noOfEntries_,
     0, 0, 0, 0, 100000000, false, false, true);

  static Switch<GenericWidthGenerator,bool> interfaceBRNormalize
    ("BRNormalize",
     "Normalize the partial widths so that they have the value BR*Total Width"
     " for an on-shell particle",
     &GenericWidthGenerator::BRnorm_, false, false, false);
  static SwitchOption interfaceBRNormalizeNormalize
    (interfaceBRNormalize,
     "Yes",
     "Perform the normalization",
     true);
  static SwitchOption interfaceBRNormalizeNoNormalisation
    (interfaceBRNormalize,
     "No",
     "Do not perform the normalization",
     false);

  static Parameter<GenericWidthGenerator,double> interfaceBRMinimum
    ("BRMinimum",
     "Minimum branching ratio for inclusion in the running width calculation.",
     &GenericWidthGenerator::BRminimum_, 0.01, 0.0, 1.0,
     false, false, true);

  static Parameter<GenericWidthGenerator,double> interfacePrefactor
    ("Prefactor",
     "The prefactor to get the correct on-shell width",
     &GenericWidthGenerator::prefactor_, 1.0, 0., 1000.,
     false, false, false);

  static Parameter<GenericWidthGenerator,int> interfacePoints
    ("Points",
     "Number of points to use for interpolation tables when needed",
     &GenericWidthGenerator::npoints_, 50, 5, 1000,
     false, false, true);

  static ParVector<GenericWidthGenerator,string> interfaceDecayModes
    ("DecayModes",
     "The tags for the decay modes used in the width generator",
     &GenericWidthGenerator::decayTags_, -1, "", "", "",
     false, false, Interface::nolimits);

  static Parameter<GenericWidthGenerator,unsigned int> interfaceInterpolationOrder
    ("InterpolationOrder",
     "The interpolation order for the tables",
     &GenericWidthGenerator::intOrder_, 1, 1, 5,
     false, false, Interface::limited);

  static Switch<GenericWidthGenerator,bool> interfaceTwoBodyOnly
    ("TwoBodyOnly",
     "Only Use two-body modes for the calculation of the running "
     "width, higher multiplicity modes fixed partial width",
     &GenericWidthGenerator::twoBodyOnly_, false, false, false);
  static SwitchOption interfaceTwoBodyOnlyYes
    (interfaceTwoBodyOnly,
     "Yes",
     "Only include two-body modes",
     true);
  static SwitchOption interfaceTwoBodyOnlyNo
    (interfaceTwoBodyOnly,
     "No",
     "Include all modes",
     false);

}

Energy GenericWidthGenerator::width(const ParticleData &, Energy m) const {
  Energy gamma= ZERO;
  for(unsigned int ix =0;ix<MEcoupling_.size();++ix) {
    if(modeOn_[ix]) gamma +=partialWidth(ix,m);
  }
  return gamma*prefactor_;
}

void GenericWidthGenerator::doinit() {
  WidthGenerator::doinit();
  if(particle()->widthGenerator()!=this) return;
  // make sure the particle data object was initialized
  particle_->init();
  tDecayIntegratorPtr decayer;
  // mass of the decaying particle
  mass_ = particle_->mass();
  if(initialize_) {
    // the initial prefactor
    prefactor_=1.;
    // resize all the storage vectors
    MEtype_.clear();
    MEcode_.clear();
    MEmass1_.clear();
    MEmass2_.clear();
    MEcoupling_.clear(); 
    modeOn_.clear();
    minMass_.clear();
    interMasses_.clear();
    interWidths_.clear();
    noOfEntries_.clear();
    decayTags_.clear();
    // integrators that we may need
    WidthCalculatorBasePtr widthptr;
    // get the list of decay modes as a decay selector
    DecayMap modes=particle_->decaySelector();
    if ( Debug::level > 0 )
      Repository::cout() << "Width generator for "
			 << particle_->PDGName() << endl;
    DecayMap::const_iterator start=modes.begin();
    DecayMap::const_iterator end=modes.end();
    tPDPtr part1,part2;
    tGenericMassGeneratorPtr massgen1,massgen2;
    // loop over the decay modes to get the partial widths
    for(;start!=end;++start) {
      // the decay mode
      tcDMPtr mode=(*start).second;
      clock_t time = std::clock();
      if ( Debug::level > 1 ) {
	Repository::cout() << "Partial width " 
			   << left
			   << std::setw(40)
			   << mode->tag() << flush;
      }
      decayModes_.push_back(const_ptr_cast<DMPtr>(mode));
      decayTags_.push_back(decayModes_.back()->tag());
      ParticleMSet::const_iterator pit(mode->products().begin());
      // minimum mass for the decaymode
      Energy minmass = ZERO;
      for(;pit!=mode->products().end();++pit) {
	(**pit).init();
	minmass+=(**pit).massMin();
      }
      minMass_.push_back(minmass);
      pit=mode->products().begin();
      // its decayer
      decayer=dynamic_ptr_cast<tDecayIntegratorPtr>(mode->decayer());
      if(decayer) decayer->init();
      // if there's no decayer then set the partial width to the br times the
      // on-shell value
      if(!decayer) {
	MEtype_.push_back(0);
	MEcode_.push_back(0);
	MEcoupling_.push_back(mode->brat());
	MEmass1_.push_back(ZERO);
	MEmass2_.push_back(ZERO);
	noOfEntries_.push_back(interMasses_.size());
	modeOn_.push_back(mode->brat()>BRminimum_);
	setupMode(mode,decayer,MEtype_.size()-1);
      }
      else if(mode->products().size()==2) {
	// the outgoing particles
	ParticleMSet::const_iterator pit = mode->products().begin();
	part1=*pit;++pit;
	part2=*pit;
	// mass generators
	if( part1->stable() || part1->massGenerator())
	  massgen1=dynamic_ptr_cast<tGenericMassGeneratorPtr>(part1->massGenerator());
	else
	  massgen1=tGenericMassGeneratorPtr();
	if(part2->stable() || part2->massGenerator())
	  massgen2=dynamic_ptr_cast<tGenericMassGeneratorPtr>(part2->massGenerator());
	else
	  massgen2=tGenericMassGeneratorPtr();
	if(massgen1) massgen1->init();
	if(massgen2) massgen2->init();
	double coupling(0.);
	int mecode(-1);
	bool order(decayer->twoBodyMEcode(*mode,mecode,coupling));
	MEcode_.push_back(mecode);
	MEcoupling_.push_back(coupling);
	modeOn_.push_back(mode->brat()>BRminimum_);
	if(order) {
	  MEmass1_.push_back(part1->mass());
	  MEmass2_.push_back(part2->mass());
	}
	else {
	  MEmass1_.push_back(part2->mass());
	  MEmass2_.push_back(part1->mass());
	}
	// perform setup in the inheriting class
	setupMode(mode,decayer,MEcode_.size()-1);
	// both particles on shell
	if(!massgen1&&!massgen2) {
	  MEtype_.push_back(1);
	  noOfEntries_.push_back(interMasses_.size());
	  if(BRnorm_) {
	    if(mass_>MEmass1_[MEtype_.size()-1]+MEmass2_[MEtype_.size()-1]) {
	      Energy gamma(partial2BodyWidth(MEtype_.size()-1,mass_));
	      if(gamma==ZERO) {
		cerr << "Partial width for " << mode->tag()
		     << " is zero in GenericWidthGenerator::doinit().\n"
		     << "If doing BSM physics this is probably a problem with your input "
		     << "parameters.\n"
		     << "Zeroing mode\n";
		MEcoupling_.back() = 0.;
	      }
	      else {
		double ratio(mode->brat()*mode->parent()->width()/gamma);
		ratio=sqrt(ratio);
		MEcoupling_.back() *=ratio;
	      }
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
	      (new_ptr(TwoBodyAllOnCalculator(ttthis,MEcode_.size()-1,
					      MEmass1_[MEcode_.size()-1],
					      MEmass2_[MEcode_.size()-1])));
	    int ioff = ((part1->massGenerator()&&!order)||
			(part2->massGenerator()&&order)) ? 2 : 1;
	    if(massgen1)
	      widthptr=new_ptr(OneOffShellCalculator(ioff,twobody,massgen1,ZERO));
	    else
	      widthptr=new_ptr(OneOffShellCalculator(ioff,twobody,massgen2,ZERO));
	  }
	  else {
	    int ioff   = order ? 1 : 2;
	    int iother = order ? 2 : 1;
	    // create the width calculator
	    tGenericWidthGeneratorPtr 
	      ttthis(const_ptr_cast<tGenericWidthGeneratorPtr>(this));
	    // this is the both on-shell case
	    WidthCalculatorBasePtr twobody
	      (new_ptr(TwoBodyAllOnCalculator(ttthis,MEcode_.size()-1,
					      MEmass1_[MEcode_.size()-1],
					      MEmass2_[MEcode_.size()-1])));
	    // this is the first off-shell
	    WidthCalculatorBasePtr widthptr2=
	      new_ptr(OneOffShellCalculator(ioff,twobody,massgen1,ZERO));
	    widthptr=new_ptr(TwoOffShellCalculator(iother,widthptr2,massgen2,
						   ZERO,massgen1->lowerLimit()));
	  }
	  // set up the interpolation table
	  Energy test(part1->massMin()+part2->massMin());
	  Energy min(max(particle_->massMin(),test)),upp(particle_->massMax());
	  Energy step((upp-min)/(npoints_-1));
	  Energy moff(min);
	  Energy2 moff2;
	  // additional points to improve the interpolation
	  if(min==test) {
	    interMasses_.push_back(moff-2.*step);interWidths_.push_back(ZERO);
	    interMasses_.push_back(moff-   step);interWidths_.push_back(ZERO);
	    interMasses_.push_back(moff        );interWidths_.push_back(ZERO);
	    double fact(exp(0.1*log(1.+step/moff)));
	    for(unsigned int ix=0;ix<10;++ix) {
	      moff*=fact;
	      moff2=sqr(moff);
	      interMasses_.push_back(moff);
	      interWidths_.push_back(widthptr->partialWidth(moff2));
	    }
	    moff+=step;
	  }
	  else if(test>min-2.*step) {
	    interMasses_.push_back(moff-2.*step);interWidths_.push_back(ZERO);
	    interMasses_.push_back(test        );interWidths_.push_back(ZERO);
	  }
	  else {
	    interMasses_.push_back(moff-2.*step);
	    interWidths_.push_back(widthptr->partialWidth((moff-2.*step)*
							  (moff-2.*step)));
	    interMasses_.push_back(moff-   step);
	    interWidths_.push_back(widthptr->partialWidth((moff-   step)*
							  (moff-   step)));
	  }
	  for(; moff<upp+2.5*step;moff+=step) {
	    moff2=moff*moff;
	    interMasses_.push_back(moff);
	    interWidths_.push_back(widthptr->partialWidth(moff2));
	  }
	  if(BRnorm_) {
	    double ratio(1.);
	    if((massgen1&&massgen2&&
		mass_>massgen1->lowerLimit()+massgen2->lowerLimit())||
	       (massgen1&&!massgen2&&
		mass_>massgen1->lowerLimit()+part2->mass())||
	       (massgen2&&!massgen1&&
		mass_>massgen2->lowerLimit()+part1->mass())||
	       (!massgen1&&!massgen2&&
		mass_>part1->mass()+part2->mass())) {
	      Energy gamma(widthptr->partialWidth(mass_*mass_));
	      if(gamma==ZERO) {
		cerr << "Partial width for " << mode->tag()
		     << " is zero in GenericWidthGenerator::doinit()"
		     << " if doing BSM physics this is probably a problem with your input "
		     << "parameters.\n"
		     << "Zeroing mode\n";
		ratio = 0.;
	      }
	      else {
		ratio=mode->brat()*mode->parent()->width()/gamma;
	      }
	    }
	    MEcoupling_.back()=ratio;
	  }
	  else MEcoupling_.back()=1.;
	  MEtype_.push_back(2);
	  MEcode_.back()=0;
	  noOfEntries_.push_back(interMasses_.size());
	  interpolators_.resize(MEtype_.size());
	  // get the vectors we will need
	  vector<Energy>::iterator istart= interMasses_.begin();
	  if(MEtype_.size()>1){istart+=noOfEntries_[MEtype_.size()-2];}
	  vector<Energy>::iterator iend=interMasses_.end();
	  vector<Energy> masses(istart,iend);

	  istart= interWidths_.begin();
	  if(MEtype_.size()>1){istart+=noOfEntries_[MEtype_.size()-2];}
	  iend=interWidths_.end();
	  vector<Energy> widths(istart,iend);
	  interpolators_.back() = make_InterpolatorPtr(widths,masses,intOrder_);
	}
      }
      // higher multiplicities
      else {
	setupMode(mode,decayer,MEcode_.size());
	widthptr = twoBodyOnly_ ? WidthCalculatorBasePtr() : decayer->threeBodyMEIntegrator(*mode);
	if(!widthptr) {
	  MEtype_.push_back(0);
	  MEcode_.push_back(0);
	  MEcoupling_.push_back(mode->brat());
	  MEmass1_.push_back(ZERO);
	  MEmass2_.push_back(ZERO);
	  noOfEntries_.push_back(interMasses_.size());
	  modeOn_.push_back(mode->brat()>BRminimum_);
	}
	else {
	  Energy step((particle_->widthUpCut()+particle_->widthLoCut())/
		      (npoints_-1));
	  Energy moff(particle_->massMin()),upp(particle_->massMax());
	  for( ; moff<upp+0.5*step;moff+=step) {
	    Energy2 moff2=sqr(moff);
	    Energy wtemp=widthptr->partialWidth(moff2);
	    interMasses_.push_back(moff);
	    interWidths_.push_back(wtemp);
	  }
	  double coupling(1.);
	  if(BRnorm_) {
	    Energy gamma = widthptr->partialWidth(mass_*mass_);
	    if(gamma==ZERO) {
	      cerr << "Partial width for " << mode->tag()
		   << " is zero in GenericWidthGenerator::doinit()"
		   << " if doing BSM physics this is probably a problem with your input "
		   << "parameters.\n"
		   << "Zeroing mode\n";
	      coupling = 0.;
	    }
	    else {
	      coupling = mode->brat()*mode->parent()->width()/gamma;
	    }
	  }
	  MEtype_.push_back(2);
	  MEcode_.push_back(0);
	  MEcoupling_.push_back(coupling);
	  MEmass1_.push_back(ZERO);
	  MEmass2_.push_back(ZERO);
	  modeOn_.push_back(mode->brat()>BRminimum_);
	  noOfEntries_.push_back(interMasses_.size());
	  interpolators_.resize(MEtype_.size());
	  // get the vectors we will need
	  vector<Energy>::iterator istart( interMasses_.begin()),
	    iend(interMasses_.end());
	  if(MEtype_.size()>1){istart+=noOfEntries_[MEtype_.size()-2];}
	  vector<Energy> masses(istart,iend);
	  
	  istart= interWidths_.begin();
	  if(MEtype_.size()>1){istart+=noOfEntries_[MEtype_.size()-2];}
	  iend=interWidths_.end();
	  vector<Energy> widths(istart,iend);
	  interpolators_.back() = make_InterpolatorPtr(widths,masses,intOrder_);
	}
      }
      if ( Debug::level > 1 ) {
	double diff = double(std::clock()-time)/CLOCKS_PER_SEC;
	if ( diff > 0.2 )
	  Repository::cout() << ' ' << diff << " s";
	Repository::cout() << endl;
      }
    }
    // now check the overall normalisation of the running width
    Energy gamma = width(*particle_,mass_);
    if(gamma>ZERO) prefactor_ = particle_->width()/gamma;
    // output the info so it can be read back in
  }
  else {
    // get the decay modes from the tags
    if(decayTags_.size()!=0) {
      decayModes_.clear();
      for(unsigned int ix=0;ix<decayTags_.size();++ix) {
	decayModes_.push_back(CurrentGenerator::current().findDecayMode(decayTags_[ix]));
	if(!decayModes_.back()) 
	  generator()->log() << "Error in GenericWidthGenerator::doinit(). "
			     << "Failed to find DecayMode  for tag" 
			     << decayTags_[ix] << "\n";
      }
    }
    // otherwise just use the modes from the selector
    else {
      DecaySet modes(particle_->decayModes());
      DecaySet::const_iterator start(modes.begin()),end(modes.end());
      tcDMPtr mode;
      for(;start!=end;++start) {   
	decayModes_.push_back(const_ptr_cast<DMPtr>(*start));
      }
    }
    // set up the interpolators
    interpolators_.resize(MEtype_.size());
    vector<Energy>::iterator estart(interMasses_.begin()),eend;
    vector<Energy>::iterator wstart(interWidths_.begin()),wend;
    vector<Energy> masses,widths;
    for(unsigned int ix=0;ix<MEtype_.size();++ix) {
      eend=interMasses_.begin()+noOfEntries_[ix];
      wend=interWidths_.begin()+noOfEntries_[ix];
      if(MEtype_[ix]==2) {
	masses.assign(estart,eend);
	widths.assign(wstart,wend);
	interpolators_[ix]= make_InterpolatorPtr(widths,masses,intOrder_);
      }
      estart=eend;
      wstart=wend;
    }
  }
  // setup the partial widths in the decayers for normalization
  tDecayIntegratorPtr temp;
  for(unsigned int ix=0;ix<decayModes_.size();++ix) {
    if(!decayModes_[ix]) continue;
    decayModes_[ix]->init();
    decayer=dynamic_ptr_cast<tDecayIntegratorPtr>(decayModes_[ix]->decayer());
    if(!decayer) continue;
    decayer->init();
    if(particle_->widthGenerator() && 
       particle_->widthGenerator()==this ) decayer->setPartialWidth(*decayModes_[ix],ix);
  }
  if ( Debug::level > 29 ) {
  //  code to output plots
    string fname = CurrentGenerator::current().filename() + 
      string("-") + name() + string(".top");
    ofstream output(fname.c_str());
    Energy step = (particle_->massMax()-particle_->massMin())/100.;
    output << "SET FONT DUPLEX\n";
    output << "TITLE TOP \"Width for " << particle_->name() << "\"\n";
    output << "TITLE BOTTOM \"m/GeV\"\n";
    output << "TITLE LEFT \"G/GeV\"\n";
    output << "CASE       \"F    \"\n";
    output << "SET LIMITS X " 
  	 << (particle_->massMin()-10.*step)/GeV << " " 
  	 << particle_->massMax()/GeV << "\n";
    Energy upper(ZERO);
    for(Energy etest=particle_->massMin();etest<particle_->massMax();etest+=step) {
      Energy gamma=width(*particle_,etest);
      upper = max(gamma,upper);
      output << etest/GeV << "\t" << gamma/GeV << "\n";
    }
    output << "SET LIMITS Y 0. " << upper/GeV << "\n";
    output << "JOIN\n";
    output << (particle_->massMin()-9.*step)/GeV << "\t" 
  	 <<  upper*(MEcode_.size()+1)/(MEcode_.size()+2)/GeV << "\n";
    output << (particle_->massMin()-7.*step)/GeV << "\t" 
  	 <<  upper*(MEcode_.size()+1)/(MEcode_.size()+2)/GeV << "\n";
    output << "JOIN\n";
    output << "TITLE DATA " 
  	 << (particle_->massMin()-6.*step)/GeV << "\t" 
  	 <<  upper*(MEcode_.size()+1)/(MEcode_.size()+2)/GeV 
  	 << " \"total\"\n";
    for(unsigned int ix=0;ix<MEcode_.size();++ix) {
      for(Energy etest=particle_->massMin();etest<particle_->massMax();etest+=step) {
        output << etest/GeV << "\t" << partialWidth(ix,etest)*prefactor_/GeV << "\n";
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
      output << (particle_->massMin()-9.*step)/GeV << "\t" 
  	   <<  upper*(MEcode_.size()-ix)/(MEcode_.size()+2)/GeV << "\n";
      output << (particle_->massMin()-7.*step)/GeV << "\t" 
  	   <<  upper*(MEcode_.size()-ix)/(MEcode_.size()+2)/GeV << "\n"; 
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
  	   << (particle_->massMin()-6.*step)/GeV << "\t" 
  	   <<  upper*(MEcode_.size()-ix)/(MEcode_.size()+2)/GeV 
  	   << " \"" << decayTags_[ix] << "\"\n";
    }
}
}

void GenericWidthGenerator::dataBaseOutput(ofstream & output, bool header) {
  if(header) output << "update Width_Generators set parameters=\"";
  // prefactor and general switiches
  output << "newdef " << name() << ":Prefactor "   << prefactor_ << "\n";
  output << "newdef " << name() << ":BRNormalize " << BRnorm_    << "\n";
  output << "newdef " << name() << ":BRMinimum "   << BRminimum_ << "\n";
  output << "newdef " << name() << ":Points "      << npoints_   << "\n";
  output << "newdef " << name() << ":InterpolationOrder " << intOrder_ << "\n";
  // the type of the matrix element
  for(unsigned int ix=0;ix<MEtype_.size();++ix) {
    output << "insert " << name() << ":MEtype " << ix << " " 
	   << MEtype_[ix] << "\n";
  }
  // the code for thew two body matrix elements
  for(unsigned int ix=0;ix<MEcode_.size();++ix) {
    output << "insert " << name() << ":MEcode " 
	   << ix << " " << MEcode_[ix] << "\n";
  }
  // the coupling for trhe two body matrix elements
  for(unsigned int ix=0;ix<MEcoupling_.size();++ix) {
    output << "insert " << name() << ":MEcoupling " 
	   << ix << " " << MEcoupling_[ix] << "\n";
  }
  // use this mode for the running width
  for(unsigned int ix=0;ix<modeOn_.size();++ix) {
    output << "insert " << name() << ":ModeOn " 
	   << ix << " " << modeOn_[ix] << "\n";
  }
  // first outgoing mass
  for(unsigned int ix=0;ix<minMass_.size();++ix) {
    output << "insert " << name() << ":MinimumMasses " 
	   << ix << " " << minMass_[ix]/GeV << "\n";
  }
  // first outgoing mass
  for(unsigned int ix=0;ix<MEmass1_.size();++ix) {
    output << "insert " << name() << ":MEmass1 " 
	   << ix << " " << MEmass1_[ix]/GeV << "\n";
  }
  // second outgoing mass
  for(unsigned int ix=0;ix<MEmass2_.size();++ix) {
    output << "insert " << name() << ":MEmass2 " 
	   << ix << " " << MEmass2_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<decayModes_.size();++ix) {
    output << "insert " << name() << ":DecayModes "
	   << ix << " " << decayTags_[ix] << " \n";
  }
  // data for the interpolation tables
  std::streamsize curpre=output.precision();
  output.precision(curpre+2);
  for(unsigned int ix=0;ix<interMasses_.size();++ix) {
    output << "insert " << name() 
	   << ":InterpolationMasses " 
	   << ix << " " << interMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<interWidths_.size();++ix) {
    output << "insert " << name() 
	   << ":InterpolationWidths " 
	   << ix << " " << interWidths_[ix]/GeV << "\n";
  }
  output.precision(curpre);
  for(unsigned int ix=0;ix<noOfEntries_.size();++ix) {
    output << "insert " << name() 
	   << ":NumberofEntries " 
	   << ix << " " << noOfEntries_[ix] << "\n";
  }  
  if(header) output << "\n\" where BINARY ThePEGName=\"" << name() 
		    << "\";" << endl;
}

DecayMap GenericWidthGenerator::rate(const Particle & p) {
  // return default if not using running widths
  if(!particle_->variableRatio()) return p.data().decaySelector();
  // use the running widths to generate the branching ratios
  Energy scale(p.mass());
  DecayMap dm;
  Energy width = particle_->width();
  for(unsigned int ix=0;ix<decayModes_.size();++ix) {
    dm.insert(partialWidth(ix,scale)/width,
	      p.id()==particle_->id() ? 
	      decayModes_[ix] : decayModes_[ix]->CC());
  }
  return dm;
}

void GenericWidthGenerator::setupMode(tcDMPtr, tDecayIntegratorPtr,
				      unsigned int)
{}

Energy GenericWidthGenerator::partialWidth(int imode,Energy q) const {
  if(q<minMass_[imode]) return ZERO;
  Energy gamma;
  if(MEtype_[imode]==0) {
    gamma=MEcoupling_[imode]*particle_->width();
  }
  else if(MEtype_[imode]==1) {
    gamma=partial2BodyWidth(imode,q);
  }
  else if(MEtype_[imode]==2) {
    gamma=MEcoupling_[imode]*(*interpolators_[imode])(q);
  }
  else {
    throw Exception() << "Unknown type of mode " << MEtype_[imode] 
		      << "in GenericWidthGenerator::partialWidth()"
		      << Exception::runerror;
  }
  return max(gamma,ZERO);
}

void GenericWidthGenerator::dofinish() {
  if(output_) {
    string fname = CurrentGenerator::current().filename() + 
      string("-") + name() + string(".output");
    ofstream output(fname.c_str());
    dataBaseOutput(output,true);
  }
  WidthGenerator::dofinish();
}

void GenericWidthGenerator::rebind(const TranslationMap & trans) {
  particle_ = trans.translate(particle_);
  WidthGenerator::rebind(trans);
}

IVector GenericWidthGenerator::getReferences() {
  IVector ret = WidthGenerator::getReferences();
  ret.push_back(particle_);
  return ret;
}

Length GenericWidthGenerator::lifeTime(const ParticleData &, Energy m, Energy w) const {
  if(m<particle_->massMin()) w = width(*particle_,particle_->massMin());
  else if(m>particle_->massMax()) w = width(*particle_,particle_->massMax());
  return UseRandom::rndExp(hbarc/w);
}

Energy GenericWidthGenerator::partial2BodyWidth(int imode, Energy q,Energy m1,
						Energy m2) const {
  using Constants::pi;
  if(q<m1+m2) return ZERO;
  // calcluate the decay momentum
  Energy2 q2(q*q),m02(mass_*mass_),m12(m1*m1),m22(m2*m2),
    pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
  if(MEcode_[imode]==-1) return q/mass_*particle_->width();
  Energy  pcm(sqrt(pcm2));
  double gam(0.);
  switch(MEcode_[imode]) {
    // V -> P P
  case  0: gam = pcm2/6./q2;
    break;
    // V -> P V
  case  1: gam = pcm2/12./m02;
    break;
    // V -> f fbar
  case  2: gam = 1./12.*(q2*(2.*q2-m12-m22+6.*m1*m2)
			 -(m12-m22)*(m12-m22))/q2/q2;
    break;
    // P -> VV
  case  3: gam = 0.25*pcm2/m02;
    break;
    // A -> VP 
  case  4: gam = (2.*pcm2+3.*m12)/24./m02;
    break;
    // V -> VV
  case  5: gam = pcm2/3./q2*(1.+m12/q2+m22/q2);
    break;
    // S -> SS
  case  6: gam = 0.125/q2*m02;
    break;
    // T -> PP
  case  7: gam = pcm2*pcm2/60./q2/m02;
    break;
    // T -> VP
  case  8: gam = pcm2*pcm2/40./m02/m02;
    break;
    // T -> VV
  case  9: gam = 1./30./q2/q2/m02*
      (3.*q2*(8.*pcm2*pcm2+5.*(m12*m22+pcm2*(m12+m22)))
       -5.*(m12-m22)*(m12-m22)*pcm2);
    break;
    // P -> PV
  case 10: gam = 0.5*pcm2/m22;
    break;
    // P -> PT
  case 11: gam = sqr(pcm2)/12.*q2/m12/m12/m02;
    break;
    // S -> VV
  case 12: gam = 0.125*(2.*pcm2+3.*m12*m22/q2)/m02;
    break;
    // unknown
  default:
    throw Exception() << "Unknown type of mode " << MEcode_[imode] 
		      << " in GenericWidthGenerator::partial2BodyWidth() " 
		      << Exception::abortnow;
  }
  return gam*pcm*sqr(MEcoupling_[imode])/pi;
}

pair<Energy,Energy> 
GenericWidthGenerator::width(Energy m, const ParticleData & ) const {
  pair<Energy,Energy> gamma(make_pair(ZERO,ZERO));
  for(unsigned int ix=0;ix<decayModes_.size();++ix) {
    if(modeOn_[ix]) {
      Energy partial = partialWidth(ix,m);
      gamma.second += partial;
      if(decayModes_[ix]->on())
	gamma.first += partial;
    }
  }
  gamma.first  *= prefactor_;
  gamma.second *= prefactor_;
  return gamma;
}
