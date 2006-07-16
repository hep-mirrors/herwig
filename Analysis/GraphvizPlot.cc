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

// namespace {
//   void dotLabel(std::ostream & os, tcPPtr p)
//   {
//     switch (abs(p->id())) {
//     case ParticleID::g: 
//       os << "[shape=point,width=0.25,height=0.25,color=orange,label=\"\"]"; 
//       break;
//     case ParticleID::gamma:
//       os << "[shape=point,width=0.15,height=0.15,color=yellow,label=\"\"]"; 
//       break;
//     case ParticleID::pi0:
//       os << "[shape=point,width=0.2,height=0.2,color=lightblue,label=\"\"]"; 
//       break; 
//     case ParticleID::piplus:
//       os << "[shape=point,width=0.2,height=0.2,color=blue,label=\"\"]"; 
//       break;
//     case ParticleID::K0:
//     case ParticleID::K_L0:
//     case ParticleID::K_S0:
//       os << "[shape=point,width=0.25,height=0.25,color=palegreen,label=\"\"]"; 
//       break;   
//     case ParticleID::Kplus:
//       os << "[shape=point,width=0.25,height=0.25,color=green,label=\"\"]"; 
//       break;
//     case ParticleID::pplus:
//       os << "[shape=point,width=0.45,height=0.45,color=red,label=\"\"]"; 
//       break;
//     case ParticleID::n0:
//       os << "[shape=point,width=0.45,height=0.45,color=brown,label=\"\"]"; 
//       break;
//     case 81:
//       os << "[shape=box,width=0.25,height=0.25,fillcolor=black,style=filled,label=\"\"]"; 
//       break;
//     case 82:
//       os << "[shape=box,width=0.35,height=0.35,fillcolor=purple,style=filled,label=\"\"]"; 
//       break;
//     default:
//       os << "[label=\"" << p->PDGName() << "\"]";
//     }
//   }
// }

namespace {
  const string header = "digraph test {\nrankdir=LR;\nranksep=0.1;\n";
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
  
//         set<tcPPtr> particles;
//         hw.eventGenerator()->currentEvent()->select(inserter(particles), ThePEG::AllSelector());

//         for(set<tcPPtr>::const_iterator it = particles.begin(); 
//             it != particles.end(); ++it) {
//           //if ((*it)->data().cTau() != 0.0)
//           //   cerr << (*it)->number() << ' ' << (*it)->PDGName() << ' ' 
// //               << (((*it)->birthStep() && (*it)->birthStep()->collision()) ? 'y' : 'n') 
// //               << "\n=============\n"
// //               << (*it)->vertex()/m << '\n' << (*it)->decayVertex()/m << "\n-\n" 
// //               << (*it)->labVertex()/m << '\n' << (*it)->labDecayVertex()/m << "\n-\n" 
// //               << (*it)->lifeLength()/m << ' ' << (*it)->lifeTime()/m << "\n\n";

//             dotfile << "\"" << (*it)->number() << "\" ";
//             dotLabel(dotfile,(*it));
//             //              << "\" [label=\"" << (*it)->number() <<" "
//             //    <<(*it)->PDGName() 
//               //<< "\\n"
//               //            <<(*it)->momentum()/GeV 
//             //    << "\""
//             dotfile << "\n";
//             if ((*it)->next()) {
//               dotfile << "\"" << (*it)->number() << "\" -> \""
//                       << (*it)->next()->number() << "\" ["
//                       << "style=bold, weight=32, dir=none]\n";

//             }




//             if ((*it)->children().empty() && !(*it)->next()) {
//               dotfile << "{rank=sink {\"" << (*it)->number() << "\"} }\n";
//             } else {
//               dotfile << "\"" << (*it)->number() << "\" -> { ";
//               for(ParticleVector::const_iterator jt=(*it)->children().begin();
//                   jt != (*it)->children().end(); ++jt) {
//                 dotfile << "\"" << (*jt)->number() << "\" ";
//               }
//               dotfile << " }";

//            unsigned int num = (*it)->children().size();
//            Lorentz5Momentum p = -(*it)->momentum();
//            for (unsigned int i=0; i < num; ++i) {
//              p += (*it)->children()[i]->momentum();
//              cerr << p << '\n';
// //   unsigned int num2 = (*it)->children()[i]->parents().size();
// //           for (unsigned int j=0; j < num2; ++j) {
// //             p -= (*it)->children()[i]->parents()[j]->momentum();
// //             cerr << p << '\n';
// //           }

//            }
              
//            if (p.mag() > 1000. || abs(p.e()) > 1000.0)
//              dotfile << " [color=red, style=bold]"; 

//              dotfile << "\n";
//            }
//        }
//      }


//      dotfile << '}' << endl;
//      dotfile.close();



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

