// -*- C++ -*-
//
// FxFxFileReader.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FxFxFileReader class.
//

#include "FxFxFileReader.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <sstream>
#include <iostream>

using namespace ThePEG;

FxFxFileReader::
FxFxFileReader(const FxFxFileReader & x)
  : FxFxReader(x), neve(x.neve), ieve(0),
    LHFVersion(x.LHFVersion), outsideBlock(x.outsideBlock),
    headerBlock(x.headerBlock), initComments(x.initComments),
    initAttributes(x.initAttributes), eventComments(x.eventComments),
    eventAttributes(x.eventAttributes),
    theFileName(x.theFileName), theQNumbers(x.theQNumbers),
    theIncludeFxFxTags(x.theIncludeFxFxTags),
    theIncludeCentral(x.theIncludeCentral),
    theDecayer(x.theDecayer) {}

FxFxFileReader::~FxFxFileReader() {}

IBPtr FxFxFileReader::clone() const {
  return new_ptr(*this);
}

IBPtr FxFxFileReader::fullclone() const {
  return new_ptr(*this);
}

bool FxFxFileReader::preInitialize() const {
  return true;
}

void FxFxFileReader::doinit() {
  FxFxReader::doinit();
  // are we using QNUMBERS
  if(!theQNumbers) return;
  // parse the header block and create 
  // any new particles needed in QNUMBERS blocks
  string block = headerBlock;
  string line  = "";
  bool readingSLHA = false;
  int (*pf)(int) = tolower;
  unsigned int newNumber(0);
  do {
    line  = StringUtils::car(block,"\r\n");
    block = StringUtils::cdr(block,"\r\n");
    if(line[0]=='#') continue;
    // are we reading the SLHA block
    if(readingSLHA) {
      // reached the end of slha block ?
      if(line.find("</slha") != string::npos) {
	readingSLHA = false;
	break;
      }
      // remove trailing comment from line
      vector<string> split = StringUtils::split(line,"#");
      // check for a qnumbers block
      transform(split[0].begin(), split[0].end(), split[0].begin(), pf);
      // if not contine
      if(split[0].find("block qnumbers")==string::npos)
	continue;
      // get name from comment
      string name;
      if(split.size()>=2) {
	name = StringUtils::stripws(split[1]);
      }
      else {
	++newNumber;
	ostringstream tname;
	tname << "NP" << newNumber;
	name = tname.str();
      }
      // extract the PDG code
      split = StringUtils::split(split[0]," ");
      istringstream is(split[2]);
      long PDGCode(0);
      is >> PDGCode;
      // get the charge, spin, colour and whether an antiparticle
      int charge(0),spin(0),colour(0),anti(0);
      for(unsigned int ix=0;ix<4;++ix) {
	line = StringUtils::car(block,"\r\n");
	block = StringUtils::cdr(block,"\r\n");
	int dummy[2];
	istringstream is(line);
	is >> dummy[0] >> dummy[1];
	switch (dummy[0]) {
	case 1:
	  charge = dummy[1];
	  break;
	case 2:
	  spin   = dummy[1];
	  break;
	case 3:
	  colour = dummy[1];
	  break;
	case 4:
	  anti   = dummy[1];
	  break;
	default:
	  assert(false);
	}
      }
      // check if particles already exist
      PDPair newParticle;
      newParticle.first = getParticleData(PDGCode);
      if(newParticle.first) Throw<SetupException>()
			      << "Particle with PDG code " << PDGCode 
			      << " whose creation was requested in a QNUMBERS Block"
			      << " already exists. Retaining the original particle" 
			      << Exception::warning;
      if(anti) {
	newParticle.second = getParticleData(-PDGCode);
	if(newParticle.second) Throw<SetupException>()
				 << "Anti-particle with PDG code " << -PDGCode 
				 << " whose creation was requested in a QNUMBERS Block"
				 << " already exists. Retaining the original particle" 
				 << Exception::warning;
	if(( newParticle.first  && !newParticle.second ) ||
	   ( newParticle.second && !newParticle.first  ) )
	  Throw<SetupException>()
	    << "Either particle or anti-particle with PDG code " << PDGCode 
	    << " whose creation was requested in a QNUMBERS Block"
	    << " already exists, but not both the particle and antiparticle. "
	    << " Something dodgy here stopping" 
	    << Exception::runerror;
      }
      // already exists continue
      if(newParticle.first) continue;
      // create the particles
      // particle with no anti particle
      if( anti == 0 ) {
	// construct the name
	if(name=="") {
	  ostringstream temp;
	  temp << PDGCode;
	  name = temp.str(); 
	}
	// create the ParticleData object
	newParticle.first = ParticleData::Create(PDGCode,name);
      }
      // particle anti-particle pair
      else {
	// construct the names
	string nameAnti;
	if(name=="") {
	  ostringstream temp;
	  temp << PDGCode;
	  name = temp.str(); 
	  ostringstream temp2;
	  temp << -PDGCode;
	  nameAnti = temp2.str(); 
	}
	else {
	  nameAnti=name;
	  for(string::iterator it=nameAnti.begin();it!=nameAnti.end();++it) {
	    if(*it=='+')      nameAnti.replace(it,it+1,"-");
	    else if(*it=='-') nameAnti.replace(it,it+1,"+");
	  }
	  if(nameAnti==name) nameAnti += "bar";
	}
	// create the ParticleData objects
	newParticle = ParticleData::Create(PDGCode,name,nameAnti);
      }
      // set the particle properties
      if(colour==1) colour = 0;
      newParticle.first->iColour(PDT::Colour(colour));
      newParticle.first->iSpin  (PDT::Spin  (spin  ));
      newParticle.first->iCharge(PDT::Charge(charge));
      // register it
      generator()->preinitRegister(newParticle.first,
				   "/Herwig/Particles/"+newParticle.first->PDGName());
      // set the antiparticle properties
      if(newParticle.second) {
	if(colour==3||colour==6) colour *= -1;
	charge = -charge;
	newParticle.second->iColour(PDT::Colour(colour));
	newParticle.second->iSpin  (PDT::Spin  (spin  ));
	newParticle.second->iCharge(PDT::Charge(charge));
	// register it
	generator()->preinitRegister(newParticle.second,
				     "/Herwig/Particles/"+newParticle.second->PDGName());
      }
    }
    // start of SLHA block ?
    else if(line.find("<slha") != string::npos) {
      readingSLHA = true;
    }
  }
  while(line!="");
  // now set any masses/decay modes
  block = headerBlock;
  line="";
  readingSLHA=false;
  bool ok=true;
  do {
    line = StringUtils::car(block,"\r\n");
    block = StringUtils::cdr(block,"\r\n");
    // are we reading the SLHA block
    if(readingSLHA) {
      // reached the end?
      if(line.find("</slha") == 0 ) {
	readingSLHA = false;
	break;
      }
      // make lower case
      transform(line.begin(),line.end(),line.begin(), pf);
      // found the mass block ?
      if(line.find("block mass")!=string::npos) {
	// read it
	line = StringUtils::car(block,"\r\n");
	// check not at end
	while(line[0] != 'D' && line[0] != 'B' &&
	      line[0] != 'd' && line[0] != 'b' &&
	      line    != "") {
	  // skip comment lines
	  if(line[0] == '#') {
	    block = StringUtils::cdr(block,"\r\n");
	    line = StringUtils::car(block,"\r\n");
	    continue;
	  }
	  // get the mass and PGD code
	  istringstream temp(line);
	  long id;
	  double mass;
	  temp >> id >> mass;
	  // skip resetting masses on SM particles
	  // as it can cause problems later on in event generation
	  if(abs(id)<=6 || (abs(id)>=11 && abs(id)<=16) ||
	     abs(id)==23 || abs(id)==24) {
	    // Throw<SetupException>() << "Standard model mass for PID " 
	    // 			    << id
	    // 			    << " will not be changed."
	    // 			    << Exception::warning;
	    block = StringUtils::cdr(block,"\r\n");
	    line = StringUtils::car(block,"\r\n");
	    continue;
	  }
	  // magnitude of mass for susy models
	  mass = abs(mass);
	  // set the mass
	  tPDPtr particle = getParticleData(id);
	  if(!particle) throw SetupException() 
			  << "FxFxFileReader::doinit() - Particle with PDG code not"
			  << id << " not found." << Exception::runerror;
	  const InterfaceBase * ifb = BaseRepository::FindInterface(particle,
								    "NominalMass");
	  ostringstream os;
	  os << mass;
	  ifb->exec(*particle, "set", os.str());
	  // read the next line
	  block = StringUtils::cdr(block,"\r\n");
	  line = StringUtils::car(block,"\r\n");
	}; 
      }
      // found a decay block
      else if(line.find("decay") == 0) {
	// get PGD code and width
	istringstream iss(line);
	string dummy;
	long parent(0);
	Energy width(ZERO);
	iss >> dummy >> parent >> iunit(width, GeV);
	// get the ParticleData object
	PDPtr inpart = getParticleData(parent);
	if(!inpart)  {
	  throw SetupException() 
	    << "FxFxFileReader::doinit() - A ParticleData object with the PDG code "
	    << parent << " does not exist. " << Exception::runerror;
	  return;
	}
	if ( abs(inpart->id()) == 6 || 
	     abs(inpart->id()) == 15 || 
	     abs(inpart->id()) == 23 || 
	     abs(inpart->id()) == 24 || 
	     abs(inpart->id()) == 25 ) {
	  Throw<SetupException>() << "\n"
	    "************************************************************************\n"
	    "* Your LHE file changes the width of " << inpart->PDGName() << ".\n"
	    "* This can cause serious problems in the event generation!\n"
	    "************************************************************************\n"
	    "\n" << Exception::warning;
	}
	else if (inpart->width() > ZERO && width <= ZERO) {
	  Throw<SetupException>() << "\n"
	    "************************************************************************\n"
	    "* Your LHE file zeroes the non-zero width of " << inpart->PDGName() << ".\n"
	    "* If " << inpart->PDGName() << " is a decaying SM particle,\n"
	    "*     this can cause serious problems in the event generation!\n"
	    "************************************************************************\n"
				  "\n" << Exception::warning;
	}
	// set the width
	inpart->width(width);
	if( width > ZERO ) {
	  inpart->cTau(hbarc/width);
	  inpart->widthCut(5.*width);
	  inpart->stable(false);
	}
	// construct prefix for DecayModes
	string prefix(inpart->name() + "->"), tag(prefix),line("");
	unsigned int nmode(0);
	// read any decay modes
	line = StringUtils::car(block,"\r\n");
	while(line[0] != 'D' && line[0] != 'B' &&
	      line[0] != 'd' && line[0] != 'b' && 
	      line[0] != '<' && line    != "") {
	  // skip comments
	  if(line[0] == '#') {
	    block = StringUtils::cdr(block,"\r\n");
	    line = StringUtils::car(block,"\r\n");
	    continue;
	  }
	  // read decay mode and construct the tag
	  istringstream is(line);
	  double brat(0.);
	  unsigned int nda(0),npr(0);
	  is >> brat >> nda;
	  while( true ) {
	    long t;
	    is >> t;
	    if( is.fail() ) break; 
	    if( t == abs(parent) )
	      throw SetupException() 
		<< "An error occurred while read a decay of the " 
		<< inpart->PDGName() << ". One of its products has the same PDG code "
		<< "as the parent particle in FxFxFileReader::doinit()."
		<< " Please check the Les Houches file.\n"
		<< Exception::runerror;
	    tcPDPtr p = getParticleData(t);
	    if( !p )
	      throw SetupException()
		<< "FxFxFileReader::doinit() -"
		<< " An unknown PDG code has been encounterd "
		<< "while reading a decay mode. ID: " << t
		<< Exception::runerror;
	    ++npr;
	    tag += p->name() + ",";
	  }
	  if( npr != nda )
	    throw SetupException()
	      << "FxFxFileReader::doinit() - While reading a decay of the " 
	      << inpart->PDGName() << " from an SLHA file, an inconsistency "
	      << "between the number of decay products and the value in "
	      << "the 'NDA' column was found. Please check if the spectrum "
	      << "file is correct.\n"
	      << Exception::warning;
	  // create the DecayMode
	  if( npr > 1 ) {
	    if( nmode==0 ) {
	      generator()->preinitInterface(inpart, "VariableRatio" , "set","false");
	      if(inpart->massGenerator()) {
		ok = false;
		Throw<SetupException>()
		  << inpart->PDGName() << " already has a WidthGenerator set"
		  << " this is incompatible with using QNUMBERS "
		  << "Use\n"
		  << "set " << inpart->fullName() << ":Width_generator NULL\n"
		  << "to fix this." << Exception::warning;
	      }
	      unsigned int ntemp=0;
	      for(DecaySet::const_iterator dit = inpart->decayModes().begin();
		  dit != inpart->decayModes().end(); ++dit ) {
		if((**dit).on()) ++ntemp;
	      }
	      if(ntemp!=0) {
		ok = false;
		Throw<SetupException>()
		  << inpart->PDGName() << " already has DecayModes"
		  << " this is incompatible with using QNUMBERS "
		  << "Use\n"
		  << "do " << inpart->fullName() << ":SelectDecayModes none\n" 
		  << " to fix this." << Exception::warning;
	      }
	    }
	    inpart->stable(false);
	    tag.replace(tag.size() - 1, 1, ";");
	    DMPtr dm = generator()->findDecayMode(tag);
	    if(!theDecayer) Throw<SetupException>()
			      << "FxFxFileReader::doinit() Decayer must be set using the "
			      << "FxFxFileReader:Decayer"
			      << " must be set to allow the creation of new"
			      << " decay modes."
			      << Exception::runerror;
	    if(!dm) {
	      dm = generator()->preinitCreateDecayMode(tag);
	      if(!dm)
		Throw<SetupException>()  
		  << "FxFxFileReader::doinit() - Needed to create "
		  << "new decaymode but one could not be created for the tag " 
		  << tag << Exception::warning;
	    }
	    generator()->preinitInterface(dm, "Decayer", "set",
					  theDecayer->fullName());
	    ostringstream br;
	    br << setprecision(13) << brat;
	    generator()->preinitInterface(dm, "BranchingRatio", "set", br.str());
	    generator()->preinitInterface(dm, "Active", "set", "Yes");
	    if(dm->CC()) {
	      generator()->preinitInterface(dm->CC(), "BranchingRatio", "set", br.str());
	      generator()->preinitInterface(dm->CC(), "Active", "set", "Yes");
	    }
	    ++nmode;
	  }
	  tag=prefix;
	  // read the next line
	  block = StringUtils::cdr(block,"\r\n");
	  line = StringUtils::car(block,"\r\n");
	};
	if(nmode>0) {
	  inpart->update();
	  if(inpart->CC())
	    inpart->CC()->update();
	}
      }
    }
    // start of SLHA block ?
    else if(line.find("<slha") != string::npos) {
      readingSLHA = true;
    }
  }
  while(line!="");
  if(!ok)
    throw SetupException() << "Problem reading QNUMBERS blocks in FxFxFileReader::doinit()"
			   << Exception::runerror;
}

void FxFxFileReader::initialize(FxFxEventHandler & eh) {
  FxFxReader::initialize(eh);
  if ( LHFVersion.empty() )
    Throw<FxFxFileError>()
      << "The file associated with '" << name() << "' does not contain a "
      << "proper formatted Les Houches event file. The events may not be "
      << "properly sampled." << Exception::warning;
}

//vector<string> FxFxFileReader::optWeightNamesFunc() { return optionalWeightsNames; }
vector<string> FxFxFileReader::optWeightsNamesFunc() { return optionalWeightsNames; }

void FxFxFileReader::open() {
  if ( filename().empty() )
    throw FxFxFileError()
      << "No Les Houches file name. "
      << "Use 'set " << name() << ":FileName'."
      << Exception::runerror;
  cfile.open(filename());
  if ( !cfile )
    throw FxFxFileError()
      << "The FxFxFileReader '" << name() << "' could not open the "
      << "event file called '" << theFileName << "'."
      << Exception::runerror;

  cfile.readline();
  if ( !cfile.find("<LesHouchesEvents") ) return;
  map<string,string> attributes =
    StringUtils::xmlAttributes("LesHouchesEvents", cfile.getline());
  LHFVersion = attributes["version"];
  //cout << LHFVersion << endl;
  if ( LHFVersion.empty() ) return;

  bool readingHeader = false;
  bool readingInit = false;
  headerBlock = "";
  
//  char (cwgtinfo_weights_info[250][15]);
  string hs;
//  int cwgtinfo_nn(0);

  // Loop over all lines until we hit the </init> tag.
  bool readingInitWeights = false, readingInitWeights_sc = false;
  string weightinfo;
  while ( cfile.readline() ) {
    
    // found the init block for multiple weights
    if(cfile.find("<initrwgt")) { /*cout << "reading weights" << endl;*/ readingInitWeights = true; }
    
     // found end of init block for multiple weights: end the while loop
    if(cfile.find("</initrwgt")) { readingInitWeights = false; readingInitWeights_sc = false; continue;}

    // found the end of init block
    if(cfile.find("</init")) { readingInit = false; break; } 

    /* read the weight information block 
     * optionalWeightNames will contain information about the weights
     * this will be appended to the weight information
     */ 
    if(readingInitWeights) {
      string scalename = "";
      if(cfile.find("<weightgroup")) {
	/* we found a weight group
	 * start reading the information 
	 * within it
	 */
	readingInitWeights_sc = true;
	weightinfo = cfile.getline();
	/* to make it shorter, erase some stuff
	 */
	string str_weightgroup = "<weightgroup";
	string str_arrow = ">";
	string str_newline = "\n";
	erase_substr(weightinfo, str_weightgroup);
	erase_substr(weightinfo, str_arrow);
	erase_substr(weightinfo, str_newline);
      }
      /* if we are reading a new weightgroup, go on 
       * until we find the end of it
       */
      if(readingInitWeights_sc && !cfile.find("</weightgroup") && !cfile.find("<weightgroup")) {
	hs = cfile.getline();
	//cout << "hs=" << hs << endl;
	//cout << "weightinfo= " << weightinfo << endl;
	//fix for potential new lines:
	if(!cfile.find("<weight") and !cfile.find("</weightgroup")) {
	  weightinfo = weightinfo + hs;
	  //cout << "weightinfo fixed= " << weightinfo << endl;
	  continue;
	}
	istringstream isc(hs);
	int ws = 0;
	/* get the name that will be used to identify the scale 
	 */
	do {
	  string sub; isc >> sub;
	  if(ws==1) { string str_arrow =  ">"; erase_substr(sub, str_arrow); scalename = sub; }
	  ++ws;
	} while (isc);
	/* now get the relevant information
	 * e.g. scales or PDF sets used
	 */
	string startDEL = ">"; //starting delimiter
	string stopDEL = "</weight>"; //end delimiter
	unsigned firstLim = hs.find(startDEL); //find start of delimiter
//	unsigned lastLim = hs.find(stopDEL); //find end of delimitr
	string scinfo = hs.substr(firstLim); //define the information for the scale
	erase_substr(scinfo,stopDEL);
	erase_substr(scinfo,startDEL);
        scinfo = StringUtils::stripws(scinfo);
	//cout << "scinfo = " << scinfo << endl;
	/* fill in the map 
	 * indicating the information to be appended to each scale
	 * i.e. scinfo for each scalname
	 */
	scalemap[scalename] = scinfo.c_str();
	string str_id = "id=";
	string str_prime = "'";
	erase_substr(scalename, str_id);
	erase_substr(scalename, str_prime);
	optionalWeightsNames.push_back(scalename);
      }
    }
   
    if ( cfile.find("<header") ) {
      // following lines to headerBlock until we hit the end of it.
      readingHeader = true;
      headerBlock = cfile.getline() + "\n";
    }
    if ( (cfile.find("<init ") && !cfile.find("<initrwgt")) || cfile.find("<init>") ) {
      //cout << "found init block" << endl;
      // We have hit the init block, so we should expect to find the
      // standard information in the following. But first check for
      // attributes.
      initAttributes = StringUtils::xmlAttributes("init", cfile.getline());
      readingInit = true;
      cfile.readline();
      if ( !( cfile >> heprup.IDBMUP.first >> heprup.IDBMUP.second
		    >> heprup.EBMUP.first >> heprup.EBMUP.second
	            >> heprup.PDFGUP.first >> heprup.PDFGUP.second
	            >> heprup.PDFSUP.first >> heprup.PDFSUP.second
		    >> heprup.IDWTUP >> heprup.NPRUP ) ) {
	heprup.NPRUP = -42;
	LHFVersion = "";
	return;
      }
      heprup.resize();
      for ( int i = 0; i < heprup.NPRUP; ++i ) {
	cfile.readline();
	if ( !( cfile >> heprup.XSECUP[i] >> heprup.XERRUP[i]
	              >> heprup.XMAXUP[i] >> heprup.LPRUP[i] ) ) {
	  heprup.NPRUP = -42;
	  LHFVersion = "";
	  return;
	}
      }
    }
    if ( cfile.find("</header") ) {
      readingHeader = false;
      headerBlock += cfile.getline() + "\n";
    }
    if ( readingHeader ) {
      /* We are in the process of reading the header block. Dump the
	 line to headerBlock.*/
      headerBlock += cfile.getline() + "\n";
    }
    if ( readingInit ) {
      // Here we found a comment line. Dump it to initComments.
      initComments += cfile.getline() + "\n";
    }
  }
  string central = "central";
  if (theIncludeCentral) optionalWeightsNames.push_back(central);

  //  cout << "reading init finished" << endl;
  if ( !cfile ) {
    heprup.NPRUP = -42;
    LHFVersion = "";
    return;
  }

}

bool FxFxFileReader::doReadEvent() {
  if ( !cfile ) return false;
  if ( LHFVersion.empty() ) return false;
  if ( heprup.NPRUP < 0 ) return false;
  eventComments = "";
  outsideBlock = "";
  hepeup.NUP = 0;
  hepeup.XPDWUP.first = hepeup.XPDWUP.second = 0.0;
  optionalWeights.clear();
  optionalWeightsTemp.clear();

  // Keep reading lines until we hit the next event or the end of
  // the event block. Save any inbetween lines. Exit if we didn't
  // find an event.
  while ( cfile.readline() && !cfile.find("<event") )
    outsideBlock += cfile.getline() + "\n";

  // We found an event. First scan for attributes.
  eventAttributes = StringUtils::xmlAttributes("event", cfile.getline());

  /* information necessary for FxFx merging:
   * the npLO and npNLO tags
   */
  istringstream ievat(cfile.getline());
  int we(0), npLO(-10), npNLO(-10);
  do {
    string sub; ievat >> sub;
    if(we==2) { npLO = atoi(sub.c_str()); }
    if(we==5) { npNLO = atoi(sub.c_str()); }
    ++we;
  } while (ievat);
  optionalnpLO = npLO;
  optionalnpNLO = npNLO;
  std::stringstream npstringstream;
  npstringstream << "np " << npLO << " " << npNLO;
  std::string npstrings = npstringstream.str();
  /* the FxFx merging information 
   * becomes part of the optionalWeights, labelled -999 
   * for future reference
   */

  if(theIncludeFxFxTags) optionalWeights[npstrings.c_str()] = -999;

  string ecomstring = "ecom";
  optionalWeights[ecomstring.c_str()] = heprup.EBMUP.first+heprup.EBMUP.second;
  
  if ( !cfile.readline()  ) return false;

  // The first line determines how many subsequent particle lines we
  // have.
  if ( !( cfile >> hepeup.NUP >> hepeup.IDPRUP >> hepeup.XWGTUP
	        >> hepeup.SCALUP >> hepeup.AQEDUP >> hepeup.AQCDUP ) )
    return false;
  hepeup.resize();
  // Read all particle lines.
  for ( int i = 0; i < hepeup.NUP; ++i ) {
    if ( !cfile.readline() ) return false;
    if ( !( cfile >> hepeup.IDUP[i] >> hepeup.ISTUP[i]
	          >> hepeup.MOTHUP[i].first >> hepeup.MOTHUP[i].second
         	  >> hepeup.ICOLUP[i].first >> hepeup.ICOLUP[i].second
	          >> hepeup.PUP[i][0] >> hepeup.PUP[i][1] >> hepeup.PUP[i][2]
	          >> hepeup.PUP[i][3] >> hepeup.PUP[i][4]
        	  >> hepeup.VTIMUP[i] >> hepeup.SPINUP[i] ) )
      return false;
    if(std::isnan(hepeup.PUP[i][0])||std::isnan(hepeup.PUP[i][1])||
       std::isnan(hepeup.PUP[i][2])||std::isnan(hepeup.PUP[i][3])||
       std::isnan(hepeup.PUP[i][4])) 
      throw Exception() 
	<< "nan's as momenta in Les Houches file "
	<< Exception::eventerror;
    if(hepeup.MOTHUP[i].first -1==i || 
       hepeup.MOTHUP[i].second-1==i) {
      throw Exception()
	<< "Particle has itself as a mother in Les Houches "
	<< "file, this is not allowed\n"
	<< Exception::eventerror;
    } 
  }

 // Now read any additional comments and named weights.
  // read until the end of rwgt is found
  bool readingWeights = false, readingaMCFast = false, readingMG5ClusInfo = false;
  while ( cfile.readline() && !cfile.find("</event>")) {
    
    if(cfile.find("</applgrid")) { readingaMCFast=false; } //aMCFast weights end
    if(cfile.find("</rwgt")) { readingWeights=false; } //optional weights end
    if(cfile.find("</clustering")) { readingMG5ClusInfo=false; } // mg5 mclustering scale info end
    
    /* reading of optional weights
     */
    if(readingWeights) { 
      if(!cfile.find("<wgt")) { continue; }
      istringstream iss(cfile.getline());
      int wi = 0;
      double weightValue(0);
      string weightName = "";
      // we need to put the actual weight value into a double
      do {
	string sub; iss >> sub;
	if(wi==1) { string str_arrow = ">" ; erase_substr(sub, str_arrow); weightName = sub; }
	if(wi==2) weightValue = atof(sub.c_str());
	++wi;
      } while (iss);
      // store the optional weights found in the temporary map
      //cout << "weightName, weightValue= " << weightName << ", " << weightValue << endl;
      optionalWeightsTemp[weightName] = weightValue; 
    }
    
    /* reading of aMCFast weights
     */
    if(readingaMCFast) {
      std::stringstream amcfstringstream;
      amcfstringstream << "aMCFast " << cfile.getline();
      std::string amcfstrings = amcfstringstream.str();
      string str_newline = "\n";
      erase_substr(amcfstrings,str_newline);
      optionalWeights[amcfstrings.c_str()] = -111; //for the aMCFast we give them a weight -111 for future identification
    }

    /* read additional MG5 Clustering information 
     * used in LO merging
     */
    if(readingMG5ClusInfo) {
      string hs = cfile.getline();
      string startDEL = "<clus scale="; //starting delimiter
      string stopDEL = "</clus>"; //end delimiter
      unsigned firstLim = hs.find(startDEL); //find start of delimiter
      //   unsigned lastLim = hs.find(stopDEL); //find end of delimitr
      string mg5clusinfo = hs.substr(firstLim); //define the information for the scale
      erase_substr(mg5clusinfo,stopDEL);
      erase_substr(mg5clusinfo,startDEL);
      string str_arrow = ">";
      erase_substr(mg5clusinfo,str_arrow);
      string str_quotation = "\"";
      erase_substr(mg5clusinfo,str_quotation);
      string str_newline= "\n";
      erase_substr(mg5clusinfo,str_newline);
      optionalWeights[mg5clusinfo.c_str()] = -222; //for the mg5 scale info weights we give them a weight -222 for future identification
    }
    
    //store MG5 clustering information 
    if(cfile.find("<scales")) {
      string hs = cfile.getline();
      string startDEL = "<scales"; //starting delimiter
      string stopDEL = "</scales>"; //end delimiter
      unsigned firstLim = hs.find(startDEL); //find start of delimiter
      //    unsigned lastLim = hs.find(stopDEL); //find end of delimitr
      string mg5scinfo = hs.substr(firstLim); //define the information for the scale
      erase_substr(mg5scinfo,stopDEL);
      erase_substr(mg5scinfo,startDEL);
      string str_arrow = ">";
      erase_substr(mg5scinfo,str_arrow);
      string str_newline= "\n";
      erase_substr(mg5scinfo,str_newline);
      optionalWeights[mg5scinfo.c_str()] = -333; //for the mg5 scale info weights we give them a weight -333 for future identification
    }

    //determine start of aMCFast weights
    if(cfile.find("<applgrid")) { readingaMCFast = true;}
    //determine start of optional weights
    if(cfile.find("<rwgt")) { readingWeights = true; }
    //determine start of MG5 clustering scale information
    if(cfile.find("<clustering")) { readingMG5ClusInfo = true;}
  }
  // loop over the optional weights and add the extra information as found in the init
  for (map<string,double>::const_iterator it=optionalWeightsTemp.begin(); it!=optionalWeightsTemp.end(); ++it){
    string itfirst = it->first;
    erase_substr(itfirst, "'");
    erase_substr(itfirst, "\"");
    for (map<string,string>::const_iterator it2=scalemap.begin(); it2!=scalemap.end(); ++it2){
      //find the scale id in the scale information and add this information
      string it2first = it2->first;
      erase_substr(it2first, "'");
      erase_substr(it2first, "\"");
      //cout << "itfirst, it2first = " << itfirst << "\t" << it2first << endl;
      if(itfirst==it2first) { 
        string info = it2->second;
	string str_newline = "\n";
	erase_substr(info, str_newline);
	//set the optional weights
	optionalWeights[info] = it->second;
      }
    }
  }
  /* additionally, we set the "central" scale
   * this is actually the default event weight 
   */
  string central = "central";
  if (theIncludeCentral) optionalWeights[central] = hepeup.XWGTUP;

  if ( !cfile ) return false;
  return true;

}

void FxFxFileReader::close() {
  cfile.close();
}

void FxFxFileReader::persistentOutput(PersistentOStream & os) const {
  os << neve << LHFVersion << outsideBlock << headerBlock << initComments
     << initAttributes << eventComments << eventAttributes << theFileName
     << theQNumbers << theIncludeFxFxTags << theIncludeCentral << theDecayer;
}

void FxFxFileReader::persistentInput(PersistentIStream & is, int) {
  is >> neve >> LHFVersion >> outsideBlock >> headerBlock >> initComments
     >> initAttributes >> eventComments >> eventAttributes >> theFileName
     >> theQNumbers >> theIncludeFxFxTags >> theIncludeCentral >> theDecayer;
  ieve = 0;
}

ClassDescription<FxFxFileReader>
FxFxFileReader::initFxFxFileReader;
// Definition of the static class description member.

void FxFxFileReader::Init() {

  static ClassDocumentation<FxFxFileReader> documentation
    ("ThePEG::FxFxFileReader is an base class to be used for objects "
     "which reads event files from matrix element generators. This class is "
     "able to read plain event files conforming to the Les Houches Event File "
     "accord.");

  static Parameter<FxFxFileReader,string> interfaceFileName
    ("FileName",
     "The name of a file containing events conforming to the Les Houches "
     "protocol to be read into ThePEG. A file name ending in "
     "<code>.gz</code> will be read from a pipe which uses "
     "<code>zcat</code>. If a file name ends in <code>|</code> the "
     "preceeding string is interpreted as a command, the output of which "
     "will be read through a pipe.",
     &FxFxFileReader::theFileName, "", false, false);

  interfaceFileName.fileType();
  interfaceFileName.rank(11);

  static Switch<FxFxFileReader,bool> interfaceQNumbers
    ("QNumbers",
     "Whether or not to read search for and read a QNUMBERS"
     " block in the header of the file.",
     &FxFxFileReader::theQNumbers, false, false, false);
  static SwitchOption interfaceQNumbersYes
    (interfaceQNumbers,
     "Yes",
     "Use QNUMBERS",
     true);
  static SwitchOption interfaceQNumbersNo
    (interfaceQNumbers,
     "No",
     "Don't use QNUMBERS",
     false);


  static Switch<FxFxFileReader,bool> interfaceIncludeFxFxTags
    ("IncludeFxFxTags",
     "Include FxFx tags",
     &FxFxFileReader::theIncludeFxFxTags, true, true, false);
  static SwitchOption interfaceIncludeFxFxTagsYes
    (interfaceIncludeFxFxTags,
     "Yes",
     "Use the FxFx tags",
     true);
  static SwitchOption interfaceIncludeFxFxTagsNo
    (interfaceIncludeFxFxTags,
     "No",
     "Don't use the FxFx tags",
     false);

   static Switch<FxFxFileReader,bool> interfaceIncludeCentral
    ("IncludeCentral",
     "Include definition of central weight",
     &FxFxFileReader::theIncludeCentral, false, true, false);
  static SwitchOption interfaceIncludeCentralYes
    (interfaceIncludeCentral,
     "Yes",
     "include definition of central weight",
     true);
  static SwitchOption interfaceIncludeCentralNo
    (interfaceIncludeCentral,
     "No",
     "Don't include definition of central weight",
     false);


  static Reference<FxFxFileReader,Decayer> interfaceDecayer
    ("Decayer",
     "Decayer to use for any decays read from the QNUMBERS Blocks",
     &FxFxFileReader::theDecayer, false, false, true, true, false);

}


void FxFxFileReader::erase_substr(std::string& subject, const std::string& search) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
      subject.erase( pos, search.length() );
    }
}
