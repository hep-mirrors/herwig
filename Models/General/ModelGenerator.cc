// -*- C++ -*-
//
// ModelGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ModelGenerator class.
//

#include "ModelGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "BSMWidthGenerator.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "Herwig/Decay/DecayIntegrator.h"
#include "ThePEG/Repository/BaseRepository.h"

using namespace Herwig;

IBPtr ModelGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr ModelGenerator::fullclone() const {
  return new_ptr(*this);
}

void ModelGenerator::persistentOutput(PersistentOStream & os) const {
  os << hardProcessConstructors_ << _theDecayConstructor << particles_ 
     << offshell_ << Offsel_ << BRnorm_ << twoBodyOnly_ << howOffShell_
     << Npoints_ << Iorder_ << BWshape_ << brMin_ << decayOutput_;
}

void ModelGenerator::persistentInput(PersistentIStream & is, int) {
  is >> hardProcessConstructors_ >> _theDecayConstructor >> particles_
     >> offshell_ >> Offsel_ >> BRnorm_ >> twoBodyOnly_ >> howOffShell_
     >> Npoints_ >> Iorder_ >> BWshape_ >> brMin_ >> decayOutput_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<ModelGenerator,Interfaced>
describeThePEGModelGenerator("Herwig::ModelGenerator", "Herwig.so");


void ModelGenerator::Init() {

  static ClassDocumentation<ModelGenerator> documentation
    ("This class controls the the use of BSM physics.",
     "BSM physics was produced using the algorithm of "
     "\\cite{Gigg:2007cr,Gigg:2008yc}",
     "\\bibitem{Gigg:2007cr} M.~Gigg and P.~Richardson, \n"
     "Eur.\\ Phys.\\ J.\\  C {\\bf 51} (2007) 989.\n"
     "%%CITATION = EPHJA,C51,989;%%\n"
     " %\\cite{Gigg:2008yc}\n"
     "\\bibitem{Gigg:2008yc}\n"
     "  M.~A.~Gigg and P.~Richardson,\n"
     "  %``Simulation of Finite Width Effects in Physics Beyond the Standard Model,''\n"
     "  arXiv:0805.3037 [hep-ph].\n"
     "  %%CITATION = ARXIV:0805.3037;%%\n"
     );
 
  static RefVector<ModelGenerator,HardProcessConstructor> 
    interfaceHardProcessConstructors
    ("HardProcessConstructors",
     "The objects to construct hard processes",
     &ModelGenerator::hardProcessConstructors_, -1, 
     false, false, true, false, false);

  static Reference<ModelGenerator,Herwig::DecayConstructor> 
     interfaceDecayConstructor
     ("DecayConstructor",
      "Pointer to DecayConstructor helper class",
      &ModelGenerator::_theDecayConstructor, false, false, true, false);
  
  static RefVector<ModelGenerator,ThePEG::ParticleData> interfaceModelParticles
    ("DecayParticles",
     "ParticleData pointers to the particles requiring spin correlation "
     "decayers. If decay modes do not exist they will also be created.",
     &ModelGenerator::particles_, -1, false, false, true, false);
    
  static RefVector<ModelGenerator,ParticleData> interfaceOffshell
    ("Offshell",
     "The particles to treat as off-shell",
     &ModelGenerator::offshell_, -1, false, false, true, false);

  static Switch<ModelGenerator,int> interfaceWhichOffshell
    ("WhichOffshell",
     "A switch to determine which particles to create mass and width "
     "generators for.",
     &ModelGenerator::Offsel_, 0, false, false);
  static SwitchOption interfaceWhichOffshellSelected
    (interfaceWhichOffshell,
     "Selected",
     "Only create mass and width generators for the particles specified",
     0);
  static SwitchOption interfaceWhichOffshellAll
    (interfaceWhichOffshell,
     "All",
     "Treat all particles specified in the DecayParticles "
     "list as off-shell",
     1);
  
  static Switch<ModelGenerator,bool> interfaceBRNormalize
    ("BRNormalize",
     "Whether to normalize the partial widths to BR*total width for an "
     "on-shell particle",
     &ModelGenerator::BRnorm_, true, false, false);
  static SwitchOption interfaceBRNormalizeNormalize
    (interfaceBRNormalize,
     "Yes",
     "Normalize the partial widths",
     true);
  static SwitchOption interfaceBRNormalizeNoNormalize
    (interfaceBRNormalize,
     "No",
     "Do not normalize the partial widths",
     false);

  static Parameter<ModelGenerator,int> interfacePoints
    ("InterpolationPoints",
     "Number of points to use for interpolation tables when needed",
     &ModelGenerator::Npoints_, 10, 5, 1000,
     false, false, true);
  
  static Parameter<ModelGenerator,unsigned int> 
    interfaceInterpolationOrder
    ("InterpolationOrder", "The interpolation order for the tables",
     &ModelGenerator::Iorder_, 1, 1, 5,
     false, false, Interface::limited);

  static Switch<ModelGenerator,int> interfaceBreitWignerShape
    ("BreitWignerShape",
     "Controls the shape of the mass distribution generated",
     &ModelGenerator::BWshape_, 0, false, false);
  static SwitchOption interfaceBreitWignerShapeDefault
    (interfaceBreitWignerShape,
     "Default",
     "Running width with q in numerator and denominator width factor",
     0);
  static SwitchOption interfaceBreitWignerShapeFixedWidth
    (interfaceBreitWignerShape,
     "FixedWidth",
     "Use a fixed width",
     1);
  static SwitchOption interfaceBreitWignerShapeNoq
    (interfaceBreitWignerShape,
     "Noq",
     "Use M rather than q in the numerator and denominator width factor",
     2);
  static SwitchOption interfaceBreitWignerShapeNoNumerator
    (interfaceBreitWignerShape,
     "NoNumerator",
     "Neglect the numerator factors",
     3);

  static Switch<ModelGenerator,bool> interfaceTwoBodyOnly
    ("TwoBodyOnly",
     "Whether to use only two-body or all modes in the running width calculation",
     &ModelGenerator::twoBodyOnly_, false, false, false);
  static SwitchOption interfaceTwoBodyOnlyYes
    (interfaceTwoBodyOnly,
     "Yes",
     "Only use two-body modes",
     true);
  static SwitchOption interfaceTwoBodyOnlyNo
    (interfaceTwoBodyOnly,
     "No",
     "Use all modes",
     false);
  
  static Parameter<ModelGenerator,double> interfaceMinimumBR
    ("MinimumBR",
     "The minimum branching fraction to include",
     &ModelGenerator::brMin_, 1e-6, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<ModelGenerator,unsigned int> interfaceDecayOutput
    ("DecayOutput",
     "Option to control the output of the decay mode information",
     &ModelGenerator::decayOutput_, 1, false, false);
  static SwitchOption interfaceDecayOutputNone
    (interfaceDecayOutput,
     "None",
     "No output",
     0);
  static SwitchOption interfaceDecayOutputPlain
    (interfaceDecayOutput,
     "Plain",
     "Default plain text output",
     1);
  static SwitchOption interfaceDecayOutputSLHA
    (interfaceDecayOutput,
     "SLHA",
     "Output in the Susy Les Houches Accord format",
     2);

  static Parameter<ModelGenerator,double> interfaceMinimumWidthFraction
    ("MinimumWidthFraction",
     "Minimum fraction of the particle's mass the width can be"
     " for the off-shell treatment.",
     &ModelGenerator::minWidth_, 1e-6, 1e-15, 1.,
     false, false, Interface::limited);

  static Parameter<ModelGenerator,double> interfaceHowMuchOffShell
    ("HowMuchOffShell",
     "The multiple of the particle's width by which it is allowed to be off-shell",
     &ModelGenerator::howOffShell_, 5., 0.0, 100.,
     false, false, Interface::limited);

}

namespace {
  /// Helper function for sorting by mass
  inline bool massIsLess(tcPDPtr a, tcPDPtr b) {
    return a->mass() < b->mass();
  }

  // Helper function to find minimum possible mass of a particle
  inline Energy minimumMass(tcPDPtr parent) {
    Energy output(Constants::MaxEnergy);
    for(set<tDMPtr>::const_iterator dit = parent->decayModes().begin();
	dit != parent->decayModes().end(); ++dit) {
      Energy outMass(ZERO);
      for(unsigned int ix=0;ix<(**dit).orderedProducts().size();++ix) {
	outMass += (**dit).orderedProducts()[ix]->massMin();
      }
      output = min(output,outMass);
    }
    return output;
  }
}

void ModelGenerator::doinit() {
  useMe();
  Interfaced::doinit();
  // make sure the model is initialized
  Ptr<Herwig::StandardModel>::pointer model 
    = dynamic_ptr_cast<Ptr<Herwig::StandardModel>::pointer>(generator()->standardModel());
  model->init();
  // and the vertices
  for(size_t iv = 0; iv < model->numberOfVertices(); ++iv)
    model->vertex(iv)->init();
  // uniq and sort DecayParticles list by mass
  set<PDPtr> tmp(particles_.begin(),particles_.end());
  particles_.assign(tmp.begin(),tmp.end());
  sort(particles_.begin(),particles_.end(),massIsLess);
  //create decayers and decaymodes (if necessary)
  if( _theDecayConstructor ) {
    _theDecayConstructor->init();
    _theDecayConstructor->createDecayers(particles_,brMin_);
  }

  // write out decays with spin correlations
  ostream & os = CurrentGenerator::current().misc();
  ofstream ofs;
  if ( decayOutput_ > 1 ) {
    string filename 
      = CurrentGenerator::current().filename() + "-BR.spc";
    ofs.open(filename.c_str());
  }


  if(decayOutput_!=0) {
    if(decayOutput_==1) {
      os << "# The decay modes listed below will have spin\n"
	  << "# correlations included when they are generated.\n#\n#";
    }
    else {
      ofs << "#  Herwig decay tables in SUSY Les Houches accord format\n";
      ofs << "Block DCINFO                           # Program information\n";
      ofs << "1   Herwig          # Decay Calculator\n";
      ofs << "2   " << generator()->strategy()->versionstring() 
	  << "     # Version number\n";
    }
  }
  //create mass and width generators for the requested particles
  set<PDPtr> offShell;
  if( Offsel_ == 0 ) offShell = set<PDPtr>(offshell_.begin() ,offshell_.end() );
  else               offShell = set<PDPtr>(particles_.begin(),particles_.end());
  
  for(PDVector::iterator pit = particles_.begin(); 
      pit != particles_.end(); ++pit) {
    tPDPtr parent = *pit;
    // Check decays for ones where quarks cannot be put on constituent
    // mass-shell
    checkDecays(parent);
    parent->reset();
    parent->update();
    if( parent->CC() ) parent->CC()->synchronize();
    if( parent->decaySelector().empty() ) {
      parent->stable(true);
      parent->width(ZERO);
      parent->widthCut(ZERO);
      parent->massGenerator(tGenericMassGeneratorPtr());
      parent->widthGenerator(tGenericWidthGeneratorPtr());
    }
    else {
      if(parent->mass()*minWidth_>parent->width()) {
	parent->massGenerator(tGenericMassGeneratorPtr());
	parent->widthGenerator(tGenericWidthGeneratorPtr());
      }
      else {
	if( offShell.find(*pit) != offShell.end() ) {
	  createWidthGenerator(*pit);
	}
	else {
	  parent->massGenerator(tGenericMassGeneratorPtr());
	  parent->widthGenerator(tGenericWidthGeneratorPtr());
	}
      }

    }

    if( parent->massGenerator() ) {
      Energy minMass = minimumMass(parent);
      Energy offShellNess = howOffShell_*parent->width();
      if(minMass>parent->mass()-offShellNess) {
	offShellNess = parent->mass()-minMass;
      }
      parent->widthCut(offShellNess);

      parent->massGenerator()->reset();
      if(decayOutput_==1)
	os << "# " <<parent->PDGName() << " will be considered off-shell.\n#\n";
    }
    if( parent->widthGenerator() ) parent->widthGenerator()->reset();
  }
  // loop again to initialise mass and width generators
  // switch off modes and write output
  for(PDVector::iterator pit = particles_.begin();
      pit != particles_.end(); ++pit) {
    tPDPtr parent = *pit;
    if(parent->widthGenerator())
      parent->widthGenerator()->init();
    if(parent->massGenerator())
      parent->massGenerator()->init();
    // Now switch off the modes if needed
    for(DecaySet::const_iterator it=parent->decayModes().begin();
	it!=parent->decayModes().end();++it) {
      if( _theDecayConstructor->disableDecayMode((**it).tag()) )
	generator()->preinitInterface(*it, "OnOff", "set", "Off");
    }
    // output the modes if needed
    if( !parent->decaySelector().empty() ) {
      if ( decayOutput_ == 2 )
	writeDecayModes(ofs, parent);
      else
	writeDecayModes(os, parent);
    }
  }

  //Now construct hard processes given that we know which
  //objects have running widths
  for(unsigned int ix=0;ix<hardProcessConstructors_.size();++ix) {
    hardProcessConstructors_[ix]->init();
    hardProcessConstructors_[ix]->constructDiagrams();
  }
}

void ModelGenerator::checkDecays(PDPtr parent) {
  if( parent->stable() ) {
    if(parent->coloured())
      cerr << "Warning: No decays for coloured particle " << parent->PDGName() << "\n\n" 
	   << "have been calcluated in BSM model.\n"
	   << "This may cause problems in the hadronization phase.\n"
	   << "You may have forgotten to switch on the decay mode calculation using\n"
	   << "  set TwoBodyDC:CreateDecayModes Yes\n"
	   << "  set ThreeBodyDC:CreateDecayModes Yes\n"
	   << "  set WeakDecayConstructor:CreateDecayModes Yes\n"
	   << "or the decays of this particle are missing from your\n"
	   << "input spectrum and decay file in the SLHA format.\n\n";
    return;
  }
  DecaySet::iterator dit = parent->decayModes().begin();
  DecaySet::iterator dend = parent->decayModes().end();
  Energy oldwidth(parent->width()), newwidth(ZERO);
  bool rescalebrat(false);
  double brsum(0.);
  for(; dit != dend; ++dit ) {
    if( !(**dit).on() ) continue;
    Energy release((**dit).parent()->mass());
    tPDVector::const_iterator pit = (**dit).orderedProducts().begin();
    tPDVector::const_iterator pend =(**dit).orderedProducts().end();
    for( ; pit != pend; ++pit ) {
      release -= (**pit).constituentMass();
    }
    if( (**dit).brat() < brMin_ || release < ZERO ) {
      if( release < ZERO )
	cerr << "Warning: The shower cannot be generated using this decay " 
	     << (**dit).tag() << " because it is too close to threshold.\nIt "
	     << "will be switched off and the branching fractions of the "
	     << "remaining modes rescaled.\n";
      rescalebrat = true;
      generator()->preinitInterface(*dit, "OnOff", "set", "Off");
      generator()->preinitInterface(*dit, "BranchingRatio", 
				    "set", "0.0");
      DecayIntegratorPtr decayer = dynamic_ptr_cast<DecayIntegratorPtr>((**dit).decayer());
      if(decayer) {
      	generator()->preinitInterface(decayer->fullName(), "Initialize", "set","0");
      }
    }
    else {
      brsum += (**dit).brat();
      newwidth += (**dit).brat()*oldwidth;
    }
  }
  // if no modes left set stable
  if(newwidth==ZERO) {
    parent->stable(true);
    parent->width(ZERO);
    parent->widthCut(ZERO);
    parent->massGenerator(tGenericMassGeneratorPtr());
    parent->widthGenerator(tGenericWidthGeneratorPtr());
  }
  // otherwise rescale if needed
  else if( ( rescalebrat || abs(brsum - 1.) > 1e-12 ) && !parent->decayModes().empty()) {
    dit = parent->decayModes().begin();
    dend = parent->decayModes().end();
    double factor = oldwidth/newwidth;
    brsum = 0.;
    for( ; dit != dend; ++dit ) {
      if( !(**dit).on() ) continue;
      double newbrat = ((**dit).brat())*factor;
      brsum += newbrat;
      ostringstream brf;
      brf << setprecision(13) << newbrat;
      generator()->preinitInterface(*dit, "BranchingRatio",
				    "set", brf.str());
    }
    parent->width(newwidth);
    if( newwidth > ZERO ) parent->cTau(hbarc/newwidth);
  }
}

namespace {
  struct DecayModeOrdering {
    bool operator()(tcDMPtr m1, tcDMPtr m2) {
      if(m1->brat()!=m2->brat()) {
	return m1->brat()>m2->brat();
      }
      else {
	if(m1->products().size()==m2->products().size()) {
	  ParticleMSet::const_iterator it1=m1->products().begin();
	  ParticleMSet::const_iterator it2=m2->products().begin();
	  do {
	    if((**it1).id()!=(**it2).id()) {
	      return (**it1).id()>(**it2).id();
	    }
	    ++it1;
	    ++it2;
	  }
	  while(it1!=m1->products().end()&&
		it2!=m2->products().end());
	  assert(false);
	}
	else
	  return m1->products().size()<m2->products().size();
      }
      return false;
    }
  };
}

void ModelGenerator::writeDecayModes(ostream & os, tcPDPtr parent) const {
  if(decayOutput_==0) return;
  set<tcDMPtr,DecayModeOrdering> modes(parent->decayModes().begin(),
				       parent->decayModes().end());
  if(decayOutput_==1) {
    os << " Parent: " << parent->PDGName() << "  Mass (GeV): " 
       << parent->mass()/GeV << "  Total Width (GeV): " 
       << parent->width()/GeV << endl;
    os << std::left << std::setw(40) << '#' 
       << std::left << std::setw(20) << "Partial Width/GeV"
       << std::left << std::setw(20) << "BR" << "On/Off\n";
    for(set<tcDMPtr,DecayModeOrdering>::iterator dit=modes.begin();
	dit!=modes.end();++dit)
      os << std::left << std::setw(40) << (**dit).tag() 
	 << std::left << std::setw(20) << (**dit).brat()*parent->width()/GeV 
	 << std::left << std::setw(20)  << (**dit).brat()
	 << ((**dit).on() ? "On" : "Off" ) << '\n';
    os << "#\n#";
  }
  else if(decayOutput_==2) {
    os << "#    \t PDG \t Width\n";
    os << "DECAY\t" << parent->id() << "\t" << parent->width()/GeV << "\t # " << parent->PDGName() << "\n";
    for(set<tcDMPtr,DecayModeOrdering>::iterator dit=modes.begin();
	dit!=modes.end();++dit) {
      os << "\t" << std::left << std::setw(10) 
	 << (**dit).brat() << "\t" << (**dit).orderedProducts().size() 
	 << "\t";
      for(unsigned int ix=0;ix<(**dit).orderedProducts().size();++ix)
	os << std::right << std::setw(10)
	   << (**dit).orderedProducts()[ix]->id() ;
      for(unsigned int ix=(**dit).orderedProducts().size();ix<4;++ix)
	os << "\t";
      os << "# " << (**dit).tag() << "\n";
    }
  }
}

void ModelGenerator::createWidthGenerator(tPDPtr p) {
  string wn = p->fullName() + string("-WGen");
  string mn = p->fullName() + string("-MGen");
  GenericMassGeneratorPtr mgen = dynamic_ptr_cast<GenericMassGeneratorPtr>
    (generator()->preinitCreate("Herwig::GenericMassGenerator", mn));
  BSMWidthGeneratorPtr wgen = dynamic_ptr_cast<BSMWidthGeneratorPtr>
    (generator()->preinitCreate("Herwig::BSMWidthGenerator", wn));

  //set the particle interface
  mgen->particle(p);
  wgen->particle(p);

  //set the generator interfaces in the ParticleData object
  generator()->preinitInterface(p, "Mass_generator","set", mn);
  generator()->preinitInterface(p, "Width_generator","set", wn);
  //allow the branching fraction of this particle type to vary
  p->variableRatio(true);
  if( p->CC() ) p->CC()->variableRatio(true);
  
  //initialize the generators
  generator()->preinitInterface(mgen, "Initialize", "set", "Yes");
  generator()->preinitInterface(wgen, "Initialize", "set", "Yes");

  string norm = BRnorm_ ? "Yes" : "No";
  generator()->preinitInterface(wgen, "BRNormalize", "set", norm);
  string twob = twoBodyOnly_ ? "Yes" : "No";
  generator()->preinitInterface(wgen, "TwoBodyOnly", "set", twob);
  ostringstream os;
  os << Npoints_;
  generator()->preinitInterface(wgen, "Points", "set", os.str());
  os.str("");
  os << Iorder_;
  generator()->preinitInterface(wgen, "InterpolationOrder", "set",
				  os.str());
  os.str("");
  os << BWshape_;
  generator()->preinitInterface(mgen, "BreitWignerShape", "set", 
				  os.str());
}
