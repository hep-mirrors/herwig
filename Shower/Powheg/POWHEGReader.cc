// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the POWHEGReader class.
//

#include "POWHEGReader.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void POWHEGReader::persistentOutput(PersistentOStream & os) const {
}

void POWHEGReader::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<POWHEGReader> POWHEGReader::initPOWHEGReader;
// Definition of the static class description member.

void POWHEGReader::Init() {

  static ClassDocumentation<POWHEGReader> documentation
    ("There is no documentation for the POWHEGReader class");

}


void POWHEGReader::open() {
  LesHouchesFileReader::open();
}

bool POWHEGReader::doReadEvent() {
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
  }

  // Now read any additional comments.
  while ( cfile.readline() && !cfile.find("</event>") )
    eventComments += cfile.getline() + "\n";

  if ( !cfile ) return false;
  return true;

}

void POWHEGReader::close() {
  LesHouchesFileReader::close();
}
