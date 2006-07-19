// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GraphvizPlot class.
//

#include "GraphvizPlot.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Event.h"
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GraphvizPlot.tcc"
#endif
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/CLHEPWrap/GenEventConverter.h"

using namespace Herwig;
using namespace ThePEG;

namespace {
  const string header = "digraph test {\nrankdir=LR;\nranksep=1.5;\n";
}

void GraphvizPlot::analyze(tEventPtr event, long ieve, int loop, int state) {
  if (event->number() != 1) return;

  ostringstream filename;
  filename << _fileBaseName << '-' << event->number() << ".dot";
  ofstream hepmcdotfile(filename.str().c_str());
  
  hepmcdotfile << header 
	       << "node [width=0.1,height=0.1,shape=point,label=\"\"];\n";

  CLHEPMC::GenEvent * hepmc = 
    GenEventConverter::convert(*event);
  
  //  hepmc->print(hepmcfile);

  // loop over all vertices
  for (CLHEPMC::GenEvent::vertex_const_iterator 
	 it = hepmc->vertices_begin();
       it !=  hepmc->vertices_end();
       ++it) {

    // loop over incoming lines
    for (CLHEPMC::GenVertex::particles_in_const_iterator 
	   jt = (*it)->particles_in_const_begin() ;
	 jt != (*it)->particles_in_const_end() ;
	 ++jt) {

      if ((*jt)->production_vertex()) continue;
      
      hepmcdotfile << (*jt)->barcode() << " -> "
		   << (*it)->barcode() << " [label=\""
		   << generator()->
	getParticleData((*jt)->pdg_id())->PDGName()
		   << "\"]\n";
      
    }
    
    // loop over outgoing lines
    for (CLHEPMC::GenVertex::particles_out_const_iterator 
	   jt = (*it)->particles_out_const_begin() ;
	 jt != (*it)->particles_out_const_end() ;
	 ++jt) {
      hepmcdotfile << (*it)->barcode() << " -> ";

      if ((*jt)->end_vertex())
	hepmcdotfile << (*jt)->end_vertex()->barcode();
      else
	hepmcdotfile << (*jt)->barcode();
      
      hepmcdotfile << " [label=\"" 
		   << generator()->
	getParticleData((*jt)->pdg_id())->PDGName()
		   << "\"]\n";
    }
    
    if ((*it)->check_momentum_conservation() > 1.0 * MeV)
      hepmcdotfile << (*it)->barcode() 
		   << " [color=red,width=0.2,height=0.2]\n";
    
    hepmcdotfile << '\n';
  }
  hepmcdotfile << '}' << endl;
  delete hepmc;
  hepmcdotfile.close();
}

LorentzRotation GraphvizPlot::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}


void GraphvizPlot::analyze(tPPtr) {}

void GraphvizPlot::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void GraphvizPlot::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<GraphvizPlot> GraphvizPlot::initGraphvizPlot;
// Definition of the static class description member.

void GraphvizPlot::Init() {

  static ClassDocumentation<GraphvizPlot> documentation
    ("There is no documentation for the GraphvizPlot class");

  static Parameter<GraphvizPlot,string> interfaceFileName
    ("BaseName",
     "The base name of the output file. The event number will be added.",
     &GraphvizPlot::_fileBaseName, "graphviz",
     true, false);

}

