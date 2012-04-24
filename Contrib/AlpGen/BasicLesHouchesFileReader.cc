// -*- C++ -*-
//
// BasicLesHouchesFileReader.cc is a part of Herwig++ - A multi-purpose
// Monte Carlo event generator.
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BasicLesHouchesFileReader class.
//
#include "BasicLesHouchesFileReader.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/NoPDF.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/EventRecord/TmpTransform.h"
#include "ThePEG/Utilities/UtilityBase.h"

using namespace Herwig;

BasicLesHouchesFileReader::
BasicLesHouchesFileReader(const BasicLesHouchesFileReader & x)
  : LesHouchesReader(x), neve(x.neve), ieve(0),
    LHFVersion(x.LHFVersion), outsideBlock(x.outsideBlock),
    headerBlock(x.headerBlock), initComments(x.initComments),
    initAttributes(x.initAttributes), eventComments(x.eventComments),
    eventAttributes(x.eventAttributes),
    theFileName(x.theFileName),overSampling_(x.overSampling_) {}

BasicLesHouchesFileReader::~BasicLesHouchesFileReader() {}

IBPtr BasicLesHouchesFileReader::clone() const {
  return new_ptr(*this);
}

IBPtr BasicLesHouchesFileReader::fullclone() const {
  return new_ptr(*this);
}

bool BasicLesHouchesFileReader::preInitialize() const {
  return true;
}

void BasicLesHouchesFileReader::doinit() {
  LesHouchesReader::doinit();
}

void BasicLesHouchesFileReader::initialize(LesHouchesEventHandler & eh) {
  LesHouchesReader::initialize(eh);
  if ( LHFVersion.empty() )
    Throw<LesHouchesFileError>()
      << "The file associated with '" << name() << "' does not contain a "
      << "proper formatted Les Houches event file. The events may not be "
      << "properly sampled." << Exception::warning;
}


long BasicLesHouchesFileReader::scan() {
  
  open();

  // Shall we write the events to a cache file for fast reading? If so
  // we write to a temporary file if the caches events should be
  // randomized.
  if ( cacheFileName().length() ) openWriteCacheFile();

  // Keep track of the number of events scanned.
  long neve = 0;
  long cuteve = 0;
  bool negw = false;

  // If the open() has not already gotten information about subprocesses
  // and cross sections we have to scan through the events.
  if ( !heprup.NPRUP || cacheFile() || abs(heprup.IDWTUP) != 1 ) { // why scan if IDWTUP != 1?

    HoldFlag<> isScanning(scanning);

    double oldsum = 0.0;
    vector<int> lprup;
    vector<double> newmax;
    vector<long> oldeve;
    vector<long> neweve;
    for ( int i = 0; ( maxScan() < 0 || i < maxScan() ) && readEvent(); ++i ) {
      if ( !checkPartonBin() ) Throw<LesHouchesInitError>()
	<< "Found event in LesHouchesReader '" << name()
	<< "' which cannot be handeled by the assigned PartonExtractor '"
	<< partonExtractor()->name() << "'." << Exception::runerror;
      vector<int>::iterator idit =
	find(lprup.begin(), lprup.end(), hepeup.IDPRUP);
      int id = lprup.size();
      if ( idit == lprup.end() ) {
	lprup.push_back(hepeup.IDPRUP);
	newmax.push_back(0.0);
	neweve.push_back(0);
	oldeve.push_back(0);
      } else {
	id = idit - lprup.begin();
      }
      ++neve;
      ++oldeve[id];
      oldsum += hepeup.XWGTUP;
      if ( cacheFile() ) {
 	if ( eventWeight() == 0.0 ) {
 	  ++cuteve;
	  continue;
 	}
	cacheEvent();
      }
      ++neweve[id];
      newmax[id] = max(newmax[id], abs(eventWeight()));
      if ( eventWeight() < 0.0 ) negw = true;
    }

    xSecWeights.resize(oldeve.size(), 1.0);
    for ( int i = 0, N = oldeve.size(); i < N; ++i )
      if ( oldeve[i] ) xSecWeights[i] = double(neweve[i])/double(oldeve[i]);

    if ( maxScan() < 0 || neve > NEvents() ) NEvents(neve - cuteve);

    if ( lprup.size() == heprup.LPRUP.size() ) {
      for ( int id = 0, N = lprup.size(); id < N; ++id ) {
	vector<int>::iterator idit =
	  find(heprup.LPRUP.begin(), heprup.LPRUP.end(), hepeup.IDPRUP);
	if ( idit == heprup.LPRUP.end() ) {
	  Throw<LesHouchesInitError>()
	    << "When scanning events, the LesHouschesReader '" << name()
	    << "' found undeclared processes."  << Exception::warning;
	  heprup.NPRUP = 0;
	  break;
	}
	int idh = idit - heprup.LPRUP.begin();
	heprup.XMAXUP[idh] = newmax[id];
      }	
    }
    if ( heprup.NPRUP == 0 ) {
      // No heprup block was supplied or something went wrong.
      heprup.NPRUP = lprup.size();
      heprup.LPRUP.resize(lprup.size());
      heprup.XMAXUP.resize(lprup.size());
      for ( int id = 0, N = lprup.size(); id < N; ++id ) {
	heprup.LPRUP[id] = lprup[id];
	heprup.XMAXUP[id] = newmax[id];
      }
    } else if ( abs(heprup.IDWTUP) != 1 ) {
      // Try to fix things if abs(heprup.IDWTUP) != 1.
      double sumxsec = 0.0;
      for ( int id = 0; id < heprup.NPRUP; ++id ) sumxsec += heprup.XSECUP[id];
      weightScale = picobarn*neve*sumxsec/oldsum;
    }
  }

  if ( cacheFile() ) closeCacheFile();

  if ( negw ) heprup.IDWTUP = min(-abs(heprup.IDWTUP), -1);

  return neve;

}

void BasicLesHouchesFileReader::open() {
  if ( filename().empty() )
    throw LesHouchesFileError()
      << "No Les Houches file name. "
      << "Use 'set " << name() << ":FileName'."
      << Exception::runerror;
  cfile.open(filename());
  if ( !cfile )
    throw LesHouchesFileError()
      << "The BasicLesHouchesFileReader '" << name() << "' could not open the "
      << "event file called '" << theFileName << "'."
      << Exception::runerror;

  cfile.readline();
  if ( !cfile.find("<LesHouchesEvents") ) return;
  map<string,string> attributes =
    StringUtils::xmlAttributes("LesHouchesEvents", cfile.getline());
  LHFVersion = attributes["version"];
  if ( LHFVersion.empty() ) return;

  bool readingHeader = false;
  bool readingInit = false;
  headerBlock = "";

  // Loop over all lines until we hit the </init> tag.
  while ( cfile.readline() && !cfile.find("</init>") ) {
    if ( cfile.find("<header") ) {
      // We have hit the header block, so we should dump this and all
      // following lines to headerBlock until we hit the end of it.
      readingHeader = true;
      headerBlock = cfile.getline() + "\n";
    }
    else if ( cfile.find("<init") ) {
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
    else if ( cfile.find("</header") ) {
      readingHeader = false;
      headerBlock += cfile.getline() + "\n";
    }
    else if ( readingHeader ) {
      // We are in the process of reading the header block. Dump the
	// line to headerBlock.
      headerBlock += cfile.getline() + "\n";
    }
    else if ( readingInit ) {
      // Here we found a comment line. Dump it to initComments.
      initComments += cfile.getline() + "\n";
    }
    else {
      // We found some other stuff outside the standard tags.
      outsideBlock += cfile.getline() + "\n";
    }
  }
  if ( !cfile ) {
    heprup.NPRUP = -42;
    LHFVersion = "";
    return;
  }

}
bool BasicLesHouchesFileReader::readEvent() {

  reset();

  if ( !doReadEvent() ) return false;

  // If we are just skipping event we do not need to reweight or do
  // anything fancy.
  if ( skipping ) return true;

  if ( cacheFile() && !scanning ) return true;

  // Reweight according to the re- and pre-weights objects in the
  // LesHouchesReader base class.
  lastweight = reweight();

  if ( !reweightPDF && !cutEarly() ) return true;
  // We should try to reweight the PDFs or make early cuts here.

  fillEvent();

  double x1 = incoming().first->momentum().plus()/
    beams().first->momentum().plus();

  if ( reweightPDF &&
       inPDF.first && outPDF.first && inPDF.first != outPDF.first ) {
    if ( hepeup.XPDWUP.first <= 0.0 )
      hepeup.XPDWUP.first =
	inPDF.first->xfx(inData.first, incoming().first->dataPtr(),
			 sqr(hepeup.SCALUP*GeV), x1);
    double xf = outPDF.first->xfx(inData.first, incoming().first->dataPtr(),
				  sqr(hepeup.SCALUP*GeV), x1);
    lastweight *= xf/hepeup.XPDWUP.first;
    hepeup.XPDWUP.first = xf;
  }

  double x2 = incoming().second->momentum().minus()/
    beams().second->momentum().minus();

  if ( reweightPDF &&
       inPDF.second && outPDF.second && inPDF.second != outPDF.second ) {
    if ( hepeup.XPDWUP.second <= 0.0 )
      hepeup.XPDWUP.second =
	inPDF.second->xfx(inData.second, incoming().second->dataPtr(),
			 sqr(hepeup.SCALUP*GeV), x2);
    double xf =
      outPDF.second->xfx(inData.second, incoming().second->dataPtr(),
			 sqr(hepeup.SCALUP*GeV), x2);
    lastweight *= xf/hepeup.XPDWUP.second;
    hepeup.XPDWUP.second = xf;
  }

  if ( cutEarly() ) {
    if ( !cuts().initSubProcess((incoming().first->momentum() +
				 incoming().second->momentum()).m2(),
				0.5*log(x1/x2)) ) lastweight = 0.0;
    tSubProPtr sub = getSubProcess();
    TmpTransform<tSubProPtr> tmp(sub, Utilities::getBoostToCM(sub->incoming()));
    if ( !cuts().passCuts(*sub) ) lastweight = 0.0;
  }

  return true;
}

double BasicLesHouchesFileReader::getEvent() {
  if ( cacheFile() ) {
    if (overSampling_) {
      if ( !uncacheEvent() ) reopen();
    } else {
      if ( !uncacheEvent() || stats.attempts()==NEvents() )
      throw LesHouchesReopenWarning()
	<< "More events requested than available in LesHouchesReader "
	<< name() << Exception::runerror;
    }
  } else {
    if (overSampling_) {
      if ( !readEvent() ) reopen();
    } else {
      if ( !readEvent() || stats.attempts()==NEvents() )
      throw LesHouchesReopenWarning()
	<< "More events requested than available in LesHouchesReader "
	<< name() << Exception::runerror;
    }
  }
  ++position;

  double max = maxWeights[hepeup.IDPRUP]*maxFactor;
  return max != 0.0? eventWeight()/max: 0.0;

}

void BasicLesHouchesFileReader::skip(long n) {
  HoldFlag<> skipflag(skipping);
  if(overSampling_) while ( n-- ) getEvent();
}

bool BasicLesHouchesFileReader::doReadEvent() {
  if ( !cfile ) return false;
  if ( LHFVersion.empty() ) return false;
  if ( heprup.NPRUP < 0 ) return false;
  eventComments = "";
  outsideBlock = "";
  hepeup.NUP = 0;
  hepeup.XPDWUP.first = hepeup.XPDWUP.second = 0.0;

  // Keep reading lines until we hit the next event or the end of
  // the event block. Save any inbetween lines. Exit if we didn't
  // find an event.
  while ( cfile.readline() && !cfile.find("<event") )
    outsideBlock += cfile.getline() + "\n";

  // We found an event. First scan for attributes.
  eventAttributes = StringUtils::xmlAttributes("event", cfile.getline());
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
    if(isnan(hepeup.PUP[i][0])||isnan(hepeup.PUP[i][1])||
       isnan(hepeup.PUP[i][2])||isnan(hepeup.PUP[i][3])||
       isnan(hepeup.PUP[i][4])) 
      throw Exception() 
	<< "nan's as momenta in Les Houches file "
	<< Exception::eventerror;
  }

  // Now read any additional comments.
  while ( cfile.readline() && !cfile.find("</event>") )
    eventComments += cfile.getline() + "\n";

  if ( !cfile ) return false;
  return true;

}

void BasicLesHouchesFileReader::close() {
  cfile.close();
}

void BasicLesHouchesFileReader::persistentOutput(PersistentOStream & os) const {
  os << neve << LHFVersion << outsideBlock << headerBlock << initComments
     << initAttributes << eventComments << eventAttributes << theFileName
     << overSampling_;
}

void BasicLesHouchesFileReader::persistentInput(PersistentIStream & is, int) {
  is >> neve >> LHFVersion >> outsideBlock >> headerBlock >> initComments
     >> initAttributes >> eventComments >> eventAttributes >> theFileName
     >> overSampling_;
  ieve = 0;
}

ClassDescription<BasicLesHouchesFileReader>
BasicLesHouchesFileReader::initBasicLesHouchesFileReader;
// Definition of the static class description member.

void BasicLesHouchesFileReader::Init() {

  static ClassDocumentation<BasicLesHouchesFileReader> documentation
    ("Herwig::BasicLesHouchesFileReader is an base class to be used for objects "
     "which reads event files from matrix element generators. This class is "
     "able to read plain event files conforming to the Les Houches Event File "
     "accord.");

  static Parameter<BasicLesHouchesFileReader,string> interfaceFileName
    ("FileName",
     "The name of a file containing events conforming to the Les Houches "
     "protocol to be read into ThePEG. A file name ending in "
     "<code>.gz</code> will be read from a pipe which uses "
     "<code>zcat</code>. If a file name ends in <code>|</code> the "
     "preceeding string is interpreted as a command, the output of which "
     "will be read through a pipe.",
     &BasicLesHouchesFileReader::theFileName, "", false, false);

  static Switch<BasicLesHouchesFileReader,bool> interfaceOverSampling
    ("OverSampling",
     "Allow / Forbid reading of LH events more than once by the "
     "LH reader, allowing / protecting against statistical problems.",
     &BasicLesHouchesFileReader::overSampling_, true, false, false);
  static SwitchOption AllowOverSampling
    (interfaceOverSampling,
     "AllowOverSampling",
     "The reader will read events in the file more than once if more "
     "events are needed to generate the requested number than that in "
     "the LH file.",
     true);
  static SwitchOption ForbidOverSampling
    (interfaceOverSampling,
     "ForbidOverSampling",
     "The reader will NOT read events in the file more than once if more "
     "events are needed to generate the requested number than that in "
     "the LH file - instead it will stop when all have been read.",
     false);

  interfaceFileName.fileType();
  interfaceFileName.rank(11);

}

