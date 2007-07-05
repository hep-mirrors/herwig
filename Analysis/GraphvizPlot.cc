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

#include "ThePEG/Vectors/HepMCConverter.h"
#include "HepMC/GenEvent.h"


using namespace Herwig;
using namespace ThePEG;

namespace {
  const string header = "digraph test {\nrankdir=LR;\nranksep=1.5;\n";
}

namespace ThePEG {
template<> 
struct HepMCTraits<HepMC::GenEvent> 
  : public HepMCTraitsBase<HepMC::GenEvent,
			   HepMC::GenParticle,
			   HepMC::GenVertex,
			   HepMC::Polarization> 
{};
}

void GraphvizPlot::analyze(tEventPtr event, long, int, int) {
  if (event->number() != _eventNumber) return;

  ostringstream filename;
  filename << _fileBaseName << '-' << event->number() << ".dot";
  ofstream hepmcdotfile(filename.str().c_str());
  
  hepmcdotfile << header 
	       << "node [width=0.1,height=0.1,shape=point,label=\"\"];\n";

  HepMC::GenEvent * hepmc = 
    HepMCConverter<HepMC::GenEvent>::convert(*event);
  
  //  hepmc->print(hepmcfile);

  // loop over all vertices
  for (HepMC::GenEvent::vertex_const_iterator 
	 it = hepmc->vertices_begin();
       it !=  hepmc->vertices_end();
       ++it) {

    // loop over incoming lines
    for (HepMC::GenVertex::particles_in_const_iterator 
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
    for (HepMC::GenVertex::particles_out_const_iterator 
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
    
    if ((*it)->check_momentum_conservation() > 1.0)
      hepmcdotfile << (*it)->barcode() 
		   << " [color=red,width=0.2,height=0.2]\n";
    
    hepmcdotfile << '\n';
  }
  hepmcdotfile << '}' << endl;
  delete hepmc;
  hepmcdotfile.close();
}

LorentzRotation GraphvizPlot::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}


void GraphvizPlot::analyze(tPPtr) {}

void GraphvizPlot::persistentOutput(PersistentOStream & os) const {
  os << _fileBaseName << _eventNumber;
}

void GraphvizPlot::persistentInput(PersistentIStream & is, int) {
  is >> _fileBaseName >> _eventNumber;
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

  static Parameter<GraphvizPlot,long> interfaceEventNumber
    ("EventNumber",
     "The number of the event that should be drawn.",
     &GraphvizPlot::_eventNumber, 1, 1, 1,
     false, false, Interface::lowerlim);
}

