// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BSMModel class.
//

#include "BSMModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/MassGenerator.h"
#include "ThePEG/PDT/WidthGenerator.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;

namespace {
struct ParticleOrdering {
  bool operator()(tcPDPtr p1, tcPDPtr p2) {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};
}

BSMModel::BSMModel() : decayFile_(), readDecays_(true),
		       topModesFromFile_(false),
		       tolerance_(1e-6)
{}

void BSMModel::persistentOutput(PersistentOStream & os) const {
  os << decayFile_ << topModesFromFile_ << tolerance_;
}

void BSMModel::persistentInput(PersistentIStream & is, int) {
  is >> decayFile_ >> topModesFromFile_ >> tolerance_;
}

DescribeAbstractClass<BSMModel,Herwig::StandardModel>
describeHerwigBSMModel("Herwig::BSMModel", "Herwig.so");

void BSMModel::Init() {

  static ClassDocumentation<BSMModel> documentation
    ("The BSMModel class provides a base class for BSM models including the"
     " features to read decays in the SLHA format");

  static Parameter<BSMModel,string> interfaceDecayFileName
    ("DecayFileName",
     "Name of the file from which to read decays in the SLHA format",
     &BSMModel::decayFile_, "",
     false, false);

  static Switch<BSMModel,bool> interfaceTopModes
    ("TopModes",
     "Whether ro use the Herwig SM top decays or those from the SLHA file",
     &BSMModel::topModesFromFile_, false, false, false);
  static SwitchOption interfaceTopModesFile
    (interfaceTopModes,
     "File",
     "Take the modes from the files",
     true);
  static SwitchOption interfaceTopModesHerwig
    (interfaceTopModes,
     "Herwig",
     "Use the SM ones", false);

  static Parameter<BSMModel,double> interfaceBRTolerance
    ("BRTolerance",
     "Tolerance for the sum of branching ratios to be difference from one.",
     &BSMModel::tolerance_, 1e-6, 1e-8, 0.01,
     false, false, Interface::limited);

}

void BSMModel::doinit() {
  StandardModel::doinit();
  // check if need to read decays
  if(decayFile()==""||!readDecays_) return;
  decayRead();
}

void BSMModel::decayRead() {
  // read decays
  CFileLineReader cfile;
  cfile.open(decayFile_);
  if( !cfile ) throw SetupException() 
		 << "BSMModel::doinit - An error occurred in opening the "
		 << "decay file \"" << decayFile_ << "\"."
		 << Exception::runerror;
  //Before reading the spectrum/decay files the SM higgs 
  //decay modes, mass and width generators need to be turned off.
  PDPtr h0 = getParticleData(ParticleID::h0);
  h0->widthGenerator(WidthGeneratorPtr());
  h0->massGenerator(MassGenPtr());
  h0->width(ZERO);
  h0->stable(true);
  DecaySet::const_iterator dit = h0->decayModes().begin();
  DecaySet::const_iterator dend = h0->decayModes().end();
  for( ; dit != dend; ++dit ) {
    generator()->preinitInterface(*dit, "BranchingRatio", "set", "0.");
    generator()->preinitInterface(*dit, "OnOff", "set", "Off");
  }
  // if taking the top modes from the file
  // delete the SM stuff
  if(topModesFromFile_) {
    PDPtr top = getParticleData(ParticleID::t);
    top->widthGenerator(WidthGeneratorPtr());
    top->massGenerator(MassGenPtr());
    DecaySet::const_iterator dit = top->decayModes().begin();
    DecaySet::const_iterator dend = top->decayModes().end();
    for( ; dit != dend; ++dit ) {
      generator()->preinitInterface(*dit, "BranchingRatio", "set", "0.");
      generator()->preinitInterface(*dit, "OnOff", "set", "Off");
    }
  }
  // read first line and check if this is a Les Houches event file
  cfile.readline();
  bool lesHouches = cfile.find("<LesHouchesEvents");
  bool reading = !lesHouches;
  if(lesHouches) cfile.readline();
  // function pointer for putting all characters to lower case.
  int (*pf)(int) = tolower;
  while (true) {
    string line = cfile.getline();
    // check for start of slha block in SLHA files
    if(lesHouches && !reading) {
      if(line.find("<slha")==0) reading = true;
      if(!cfile.readline()) break;
      continue;
    }
    // ignore comment lines
    if(line[0] == '#') {
      if(!cfile.readline()) break;
      continue;
    }
    // make everything lower case
    transform(line.begin(), line.end(), line.begin(), pf);
    // start of a block
    if(line.find("decay") == 0) {
      readDecay(cfile, line);
      if(!cfile) break;
      continue;
    }
    else if( lesHouches && line.find("</slha") == 0 ) {
      break;
    }
    if(!cfile.readline()) break;
  }
}

void BSMModel::readDecay(CFileLineReader & cfile, 
			 string decay) const{
  // extract parent PDG code and width
  long parent(0);
  Energy width(ZERO);
  istringstream iss(decay);
  string dummy;
  iss >> dummy >> parent >> iunit(width, GeV);
  PDPtr inpart = getBSMParticleData(parent);
  // check this ain't a SM particle
  if(abs(parent)<=5||abs(parent)==23||abs(parent)==24||
     (abs(parent)>=11&&abs(parent)<=16))
    cerr << "BSMModel::readDecay() Resetting width of " 
	   << inpart->PDGName() << " using SLHA "
	   << "file,\nthis can affect parts of the Standard Model simulation and"
	   << " is strongly discouraged.\n";
  if(!topModesFromFile_&&abs(parent)==ParticleID::t) {
    cfile.readline();
    return;
  }
  if(!inpart)  throw SetupException() 
  		 << "BSMModel::readDecay() - "
  		 << "A ParticleData object with the PDG code "
  		 << parent << " does not exist. " 
  		 << Exception::runerror;
  inpart->width(width);
  if( width > ZERO ) inpart->cTau(hbarc/width);
  inpart->widthCut(5.*width);
  Energy inMass = inpart->mass();
  string prefix(inpart->name() + "->");
  double brsum(0.);
  unsigned int nmode = 0;
  while(cfile.readline()) {
    string line = cfile.getline();
    line = StringUtils::stripws(line);
    // skip comments
    if(line[0] == '#') continue;
    // reached the end
    if( line[0] == 'B' || line[0] == 'b' ||
	line[0] == 'D' || line[0] == 'd' ||
	line[0] == '<' ) {
      cfile.resetline();
      break;
    }
    // read the mode
    // get the branching ratio and no of decay products
    istringstream is(line);
    double brat(0.);
    unsigned int nda(0),npr(0);
    is >> brat >> nda;
    vector<tcPDPtr> products,bosons;
    Energy mout(ZERO),moutnoWZ(ZERO);
    string tag = prefix;
    multiset<tcPDPtr,ParticleOrdering> outgoing;
    int charge = -inpart->iCharge();
    while( true ) {
      long t;
      is >> t;
      if( is.fail() ) break; 
      if( t == abs(parent) ) {
  	throw SetupException() 
  	  << "An error occurred while read a decay of the " 
  	  << inpart->PDGName() << ". One of its products has the same PDG code "
  	  << "as the parent particle. Please check the SLHA file.\n"
  	  << Exception::runerror;
      }
      tcPDPtr p = getBSMParticleData(t);
      if( !p ) {
  	throw SetupException()
  	  << "BSMModel::readDecay() - An unknown PDG code has been encounterd "
  	  << "while reading a decay mode. ID: " << t
  	  << Exception::runerror;
      }
      charge += p->iCharge();
      ++npr;
      outgoing.insert(p);
      Energy mass =  p->mass();
      mout += mass;
      if(abs(p->id())==ParticleID::Wplus||p->id()==ParticleID::Z0) {
  	bosons.push_back(p);
      }
      else {
  	products.push_back(p);
  	moutnoWZ += mass;
      }
    }
    if( npr != nda ) {
      throw SetupException()
 	<< "BSMModel::readDecay - While reading a decay of the " 
 	<< inpart->PDGName() << " from an SLHA file, an inconsistency "
 	<< "between the number of decay products and the value in "
 	<< "the 'NDA' column was found. Please check if the spectrum "
 	<< "file is correct.\n"
 	<< Exception::warning;
    }

    // must be at least two decay products
    if(npr<=1) continue;

    // create the tag
    for(multiset<tcPDPtr,ParticleOrdering>::iterator it=outgoing.begin();
	it!=outgoing.end();++it)
      tag += (**it).name() + ",";
    tag.replace(tag.size() - 1, 1, ";");
    if(charge!=0) {
      cerr << "BSMModel::readDecay() "
	   << "Decay mode " << tag << " read from SLHA file does not conserve charge,"
	   << "\nare you really sure you want to do this?\n";
    }
    ++nmode;
    if(nmode==1) {
      if(abs(parent)<=5||abs(parent)==23||abs(parent)==24||
	 (abs(parent)>=11&&abs(parent)<=16)) {
	cerr << "BSMModel::readDecay() Resetting the decays of " 
	     << inpart->PDGName() << " using SLHA "
	     << "file,\nthis can affect parts of the Standard Model simulation,"
	     << " give unexpected results and"
	     << " is strongly discouraged.\n";
	cerr << "Switching off all the internal modes so only those from the SLHA file "
	     << "are used, this may have unintended consequences\n";
	for(DecaySet::iterator it=inpart->decayModes().begin();
	    it!=inpart->decayModes().end();++it) {
	  generator()->preinitInterface(*it, "OnOff", "set", "Off");
	}
	if(inpart->CC()) {
	  for(DecaySet::iterator it=inpart->CC()->decayModes().begin();
	      it!=inpart->CC()->decayModes().end();++it) {
	    generator()->preinitInterface(*it, "OnOff", "set", "Off");
	  }
	}
      }
    }
    // normal option
    if(mout<=inMass) {
      inpart->stable(false);
      brsum += brat;
      createDecayMode(tag, brat);
    }
    // no possible off-shell gauge bosons throw it away
    else if(bosons.empty() || bosons.size()>2 ||
	    moutnoWZ>=inMass) {
      cerr << "BSMModel::readDecay() "
	   << "The decay " << tag << " cannot proceed for on-shell "
	   << "particles, skipping it.\n";
    }
    else {
      Energy maxMass = inMass - moutnoWZ;
      string newTag = prefix;
      for(unsigned int ix=0;ix<products.size();++ix)
	newTag += products[ix]->name() + ",";
      if(bosons.size()==1) {
	cerr << "BSMModel::readDecay() "
	     << "The decay " << tag << " cannot proceed for on-shell\n"
	     << "particles, replacing gauge boson with its decay products\n";
	vector<pair<double,string> > modes = 
	  createWZDecayModes(newTag,brat,bosons[0],maxMass);
	for(unsigned int ix=0;ix<modes.size();++ix) {
	  modes[ix].second.replace(modes[ix].second.size() - 1, 1, ";");
	  createDecayMode(modes[ix].second,modes[ix].first);
	  brsum += modes[ix].first;
	}
      }
      else if(bosons.size()==2) {
	bool identical = bosons[0]->id()==bosons[1]->id();
	if(maxMass>bosons[0]->mass()&&maxMass>bosons[1]->mass()) {
	  cerr << "BSMModel::readDecay() "
	       << "The decay " << tag << " cannot proceed for on-shell\n"
	       << "particles, replacing one of the gauge bosons"
	       << " with its decay products\n";
	  unsigned int imax = identical ? 1 : 2;
	  if(imax==2) brat *= 0.5;
	  for(unsigned int ix=0;ix<imax;++ix) {
	    string newTag2 = newTag+bosons[ix]->name()+',';
	    unsigned int iother = ix==0 ? 1 : 0;
	    vector<pair<double,string> > modes = 
	      createWZDecayModes(newTag2,brat,bosons[iother],maxMass);
	    for(unsigned int ix=0;ix<modes.size();++ix) {
	      modes[ix].second.replace(modes[ix].second.size() - 1, 1, ";");
	      createDecayMode(modes[ix].second,modes[ix].first);
	      brsum += modes[ix].first;
	    }
	  }
	}
	else {
	  cerr << "BSMModel::readDecay() "
	       << "The decay " << tag << " cannot proceed for on-shell\n"
	       << "particles, and has too many off-shell gauge bosons,"
	       << " skipping it.\n";
	}
      }
      else {
	cerr << "BSMModel::readDecay() "
	     << "The decay " << tag << " cannot proceed for on-shell\n"
	     << "particles, and has too many outgoing gauge bosons skipping it.\n";
      }
    }
  }
  if( abs(brsum - 1.) > tolerance_ && nmode!=0 ) {
    cerr << "Warning: The total branching ratio for " << inpart->PDGName()
  	 << " from the spectrum file does not sum to 1.\nThe branching fractions"
  	 << " will be rescaled. Difference from 1 is " 
	 << abs(brsum - 1.) << "\n";
  }
  if(nmode>0) {
    inpart->update();
    inpart->reset();
    if(inpart->CC()) {
      inpart->CC()->update();
      inpart->CC()->reset();
    }
    if(inpart->massGenerator())
      inpart->massGenerator()->reset();
    if(inpart->widthGenerator())
      inpart->widthGenerator()->reset();
  }
}

void BSMModel::createDecayMode(string tag, double brat) const {
  tDMPtr dm = generator()->findDecayMode(tag);
  if(!dm) dm = generator()->preinitCreateDecayMode(tag);
  generator()->preinitInterface(dm, "OnOff", "set", "On");
  generator()->preinitInterface(dm, "Decayer", "set","/Herwig/Decays/Mambo");
  ostringstream brf;
  brf << setprecision(13)<< brat;
  generator()->preinitInterface(dm, "BranchingRatio","set", brf.str());
  if(dm->CC()) {
    generator()->preinitInterface(dm->CC(), "OnOff", "set", "On");
    generator()->preinitInterface(dm->CC(), "Decayer", "set","/Herwig/Decays/Mambo");
    generator()->preinitInterface(dm->CC(), "BranchingRatio","set", brf.str());
  }
}


vector<pair<double,string> > 
BSMModel::createWZDecayModes(string tag, double brat,
			     tcPDPtr boson, Energy maxMass) const {
  vector<pair<double,string> > modes;
  double sum(0.);
  for(DecaySet::const_iterator dit=boson->decayModes().begin();
      dit!=boson->decayModes().end();++dit) {
    tcDMPtr mode = *dit;
    if(!mode->on()) continue;
    string extra;
    Energy outMass(ZERO);
    for(ParticleMSet::const_iterator pit=mode->products().begin();
	pit!=mode->products().end();++pit) {
      extra += (**pit).name() + ",";
      outMass += (**pit).mass();
    }
    if(outMass<maxMass) {
      sum += mode->brat();
      modes.push_back(make_pair(mode->brat(),tag+extra));
    }
  }
  for(unsigned int ix=0;ix<modes.size();++ix) 
    modes[ix].first *= brat/sum;
  return modes;
}
