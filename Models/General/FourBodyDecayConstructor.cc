// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourBodyDecayConstructor class.
//

#include "FourBodyDecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Decay/General/GeneralFourBodyDecayer.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "DecayConstructor.h"
#include <queue>

using namespace Herwig;

FourBodyDecayConstructor::~FourBodyDecayConstructor() {}

IBPtr FourBodyDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr FourBodyDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void FourBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << removeOnShell_ << interopt_ << widthopt_ << minReleaseFraction_ 
     << maxBoson_ << maxList_ << excludedVector_ << excludedSet_;
}

void FourBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> removeOnShell_ >> interopt_ >> widthopt_ >> minReleaseFraction_ 
     >> maxBoson_ >> maxList_ >> excludedVector_ >> excludedSet_;
}

DescribeClass<FourBodyDecayConstructor,NBodyDecayConstructorBase>
describeFourBodyDecayConstructor("Herwig::FourBodyDecayConstructor","Herwig.so");

void FourBodyDecayConstructor::Init() {

  static ClassDocumentation<FourBodyDecayConstructor> documentation
    ("The FourBodyDecayConstructor class implements a small number"
     " of 4-body decays in general models");

  static Switch<FourBodyDecayConstructor,unsigned int> interfaceRemoveOnShell
    ("RemoveOnShell",
     "Remove on-shell diagrams as should be treated as a sequence of 1->2 decays",
     &FourBodyDecayConstructor::removeOnShell_, 1, false, false);
  static SwitchOption interfaceRemoveOnShellYes
    (interfaceRemoveOnShell,
     "Yes",
     "Remove the diagrams if neither the production of decay or the intermediate"
     " can happen",
     1);
  static SwitchOption interfaceRemoveOnShellNo
    (interfaceRemoveOnShell,
     "No",
     "Never remove the intermediate",
     0);
  static SwitchOption interfaceRemoveOnShellProduction
    (interfaceRemoveOnShell,
     "Production",
     "Remove the diagram if the on-shell production of the intermediate is allowed",
     2);

  static Switch<FourBodyDecayConstructor,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option for the treatment of the widths of the intermediates",
     &FourBodyDecayConstructor::widthopt_, 1, false, false);
  static SwitchOption interfaceWidthOptionFixed
    (interfaceWidthOption,
     "Fixed",
     "Use fixed widths",
     1);
  static SwitchOption interfaceWidthOptionRunning
    (interfaceWidthOption,
     "Running",
     "Use running widths",
     2);
  static SwitchOption interfaceWidthOptionZero
    (interfaceWidthOption,
     "Zero",
     "Set the widths to zero",
     3);

  static Switch<FourBodyDecayConstructor,unsigned int> interfaceIntermediateOption
    ("IntermediateOption",
     "Option for the inclusion of intermediates in the event",
     &FourBodyDecayConstructor::interopt_, 0, false, false);
  static SwitchOption interfaceIntermediateOptionAlways
    (interfaceIntermediateOption,
     "Always",
     "Always include the intermediates",
     1);
  static SwitchOption interfaceIntermediateOptionNever
    (interfaceIntermediateOption,
     "Never",
     "Never include the intermediates",
     2);
  static SwitchOption interfaceIntermediateOptionOnlyIfOnShell
    (interfaceIntermediateOption,
     "OnlyIfOnShell",
     "Only if there are on-shell diagrams",
     0);
  
  static Parameter<FourBodyDecayConstructor,double> interfaceMinReleaseFraction
    ("MinReleaseFraction",
     "The minimum energy release for a three-body decay, as a "
     "fraction of the parent mass.",
     &FourBodyDecayConstructor::minReleaseFraction_, 1e-3, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<FourBodyDecayConstructor,unsigned int> interfaceMaximumGaugeBosons
    ("MaximumGaugeBosons",
     "Maximum number of electroweak gauge bosons"
     " to be produced as decay products",
     &FourBodyDecayConstructor::maxBoson_, 1, false, false);
  static SwitchOption interfaceMaximumGaugeBosonsNone
    (interfaceMaximumGaugeBosons,
     "None",
     "Produce no W/Zs",
     0);
  static SwitchOption interfaceMaximumGaugeBosonsSingle
    (interfaceMaximumGaugeBosons,
     "Single",
     "Produce at most one W/Zs",
     1);
  static SwitchOption interfaceMaximumGaugeBosonsDouble
    (interfaceMaximumGaugeBosons,
     "Double",
     "Produce at most two W/Zs",
     2);
  static SwitchOption interfaceMaximumGaugeBosonsTriple
    (interfaceMaximumGaugeBosons,
     "Triple",
     "Produce at most three W/Zs",
     3);

  static Switch<FourBodyDecayConstructor,unsigned int> interfaceMaximumNewParticles
    ("MaximumNewParticles",
     "Maximum number of particles from the list of "
     "decaying particles to be allowed as decay products",
     &FourBodyDecayConstructor::maxList_, 0, false, false);
  static SwitchOption interfaceMaximumNewParticlesNone
    (interfaceMaximumNewParticles,
     "None",
     "No particles from the list",
     0);
  static SwitchOption interfaceMaximumNewParticlesSingle
    (interfaceMaximumNewParticles,
     "Single",
     "A single particle from the list",
     1);
  static SwitchOption interfaceMaximumNewParticlesDouble
    (interfaceMaximumNewParticles,
     "Double",
     "Two particles from the list",
     2);
  static SwitchOption interfaceMaximumNewParticlesTriple
    (interfaceMaximumNewParticles,
     "Triple",
     "Four particles from the list",
     3);

  static RefVector<FourBodyDecayConstructor,VertexBase> interfaceExcludedVertices
    ("ExcludedVertices",
     "Vertices which are not included in the three-body decayers",
     &FourBodyDecayConstructor::excludedVector_, -1, false, false, true, true, false);
}

void FourBodyDecayConstructor::doinit() {
  NBodyDecayConstructorBase::doinit();
  excludedSet_ = set<VertexBasePtr>(excludedVector_.begin(),
				    excludedVector_.end());
  if(removeOnShell_==0) 
    generator()->log() << "Warning: Including diagrams with on-shell "
		       << "intermediates in four-body BSM decays, this"
		       << " can lead to double counting and is not"
		       << " recommended unless you really know what you are doing\n"
		       << "This can be switched off using\n set "
		       << fullName() << ":RemoveOnShell Yes\n"; 
}

void FourBodyDecayConstructor::DecayList(const set<PDPtr> & particles) {
  if( particles.empty() ) return;
  // cast the StandardModel to the Hw++ one to get the vertices
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  model->init();
  unsigned int nv(model->numberOfVertices());
  // make sure vertices are initialized
  for(unsigned int i = 0; i < nv; ++i) model->vertex(i)->init();
  // loop over the particles and create the decayers
  for(set<PDPtr>::const_iterator ip=particles.begin();
      ip!=particles.end();++ip) {
    if((**ip).id()!=ParticleID::SUSY_tau_1minus) continue;
    // first create prototype 1->2 decays
    std::queue<PrototypeVertexPtr> prototypes;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      VertexBasePtr vertex = model->vertex(iv);
      tPDPtr parent = *ip;
      if(excludedSet_.find(vertex)!=excludedSet_.end()) continue;
      PrototypeVertex::createPrototypes(parent, vertex, prototypes,ZERO);
    }
    // then expand
    static const int nbody = 4;
    vector<vector<PrototypeVertexPtr> > modes;
    while(!prototypes.empty()) {
      PrototypeVertexPtr proto = prototypes.front();
      prototypes.pop();
      if(proto->npart==nbody && proto->incomingMass() > proto->outgoingMass()) {
 	if(!proto->checkExternal()) continue;
 	// check if first piece on-shell
	if(removeOnShell_|=0) {
	  Energy mass(ZERO);
	  for(OrderedVertices::const_iterator it = proto->outgoing.begin();
	      it!=proto->outgoing.end();++it) {
	    mass += it->first->mass();
	  }
	  if(mass<proto->incomingMass()) {
	    if(removeOnShell_==2)
	      continue;
	    else {
	      bool onShell=true;
	      for(OrderedVertices::const_iterator it = proto->outgoing.begin();
		  it!=proto->outgoing.end();++it) {
		if(it->second && it->second->incomingMass()<it->second->outgoingMass())
		  onShell=false;
	      }
	      if(onShell) continue;
	    }
	  }
	}
	// check if should be added to an existing decaymode
	bool added = false;
	for(unsigned int iy = 0; iy < modes.size(); ++iy) {
	  if(modes[iy][0]->sameDecay(*proto)) {
	    added = true;
	    bool already = false;
	    for(unsigned int iz = 0; iz < modes[iy].size(); ++iz) {
	      if( *modes[iy][iz] == *proto) {
		already = true;
		break;
	      }
	    }
	    if(!already) modes[iy].push_back(proto);
	    break;
	  }
	}
	if(!added) modes.push_back(vector<PrototypeVertexPtr>(1,proto));
      }
      if(proto->npart>=nbody) continue;
      // loop over all vertices
      for(unsigned int iv = 0; iv < nv; ++iv) {
	VertexBasePtr vertex = model->vertex(iv);
	if(excludedSet_.find(vertex)!=excludedSet_.end()) continue;
	PrototypeVertex::expandPrototypes(proto,vertex,prototypes);
      }
    }
    // now look at the decay modes
    for(vector<vector<PrototypeVertexPtr> >::iterator mit = modes.begin();
	mit!=modes.end();++mit) {
      // count the number of gauge bosons and particles from the list
      unsigned int nlist(0),nbos(0);
      for(OrderedParticles::const_iterator it=(*mit)[0]->outPart.begin();
	  it!=(*mit)[0]->outPart.end();++it) {
	if(abs((**it).id()) == ParticleID::Wplus ||
	   abs((**it).id()) == ParticleID::Z0) ++nbos;
	if(particles.find(*it)!=particles.end()) ++nlist;
      }
      // if do many ignore the mode
      if(nbos > maxBoson_ || nlist > maxList_) continue;
      // now create the decay mode
      bool inter=false;
      createDecayMode(*mit,inter);
    }

//   // loop over the outgoing particles
//   for(unsigned int ix=0;ix<2;++ix) {
//     tPDPtr dec   = proto.outgoing.first ;
//     tPDPtr other = proto.outgoing.second;
//     if(ix==1) swap(dec,other);
//     int id = dec->id();
//     if( !vertex->isIncoming(dec) ) continue;
//     tPDVector decaylist = vertex->search(list, dec);
//     tPDVector::size_type nd = decaylist.size();
//     for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
//       tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
//       if( pb->id() == id ) swap(pa, pb);
//       if( pc->id() == id ) swap(pa, pc);
//       //vertices are defined with all particles incoming
//       if( pb->CC() ) pb = pb->CC();
//       if( pc->CC() ) pc = pc->CC();
//       // create the three body diagram
//       TBDiagram diag(proto.incoming->id(), other->id(), 
// 		     make_pair(pb->id(),pc->id()));
//       diag.intermediate = pa;
//       diag.vertices   = make_pair(proto.vertex,vertex);
//       diag.colourFlow = vector<CFPair>(1,make_pair(1,1.));
//       diag.largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
//       decays.push_back(diag);
//     }
//   }
//   return decays; 



//   static Switch<FourBodyDecayConstructor,unsigned int> interfaceRemoveOnShell
//     ("RemoveOnShell",
//      "Remove on-shell diagrams as should be treated as a sequence of 1->2 decays",
//      &FourBodyDecayConstructor::removeOnShell_, 1, false, false);
//   static SwitchOption interfaceRemoveOnShellYes
//     (interfaceRemoveOnShell,
//      "Yes",
//      "Remove the diagrams if neither the production of decay or the intermediate"
//      " can happen",
//      1);
//   static SwitchOption interfaceRemoveOnShellNo
//     (interfaceRemoveOnShell,
//      "No",
//      "Never remove the intermediate",
//      0);
//   static SwitchOption interfaceRemoveOnShellProduction
//     (interfaceRemoveOnShell,
//      "Production",
//      "Remove the diagram if the on-shell production of the intermediate is allowed",
//      2);

//   static Switch<FourBodyDecayConstructor,unsigned int> interfaceWidthOption
//     ("WidthOption",
//      "Option for the treatment of the widths of the intermediates",
//      &FourBodyDecayConstructor::widthopt_, 1, false, false);
//   static SwitchOption interfaceWidthOptionFixed
//     (interfaceWidthOption,
//      "Fixed",
//      "Use fixed widths",
//      1);
//   static SwitchOption interfaceWidthOptionRunning
//     (interfaceWidthOption,
//      "Running",
//      "Use running widths",
//      2);
//   static SwitchOption interfaceWidthOptionZero
//     (interfaceWidthOption,
//      "Zero",
//      "Set the widths to zero",
//      3);

//   static Switch<FourBodyDecayConstructor,unsigned int> interfaceIntermediateOption
//     ("IntermediateOption",
//      "Option for the inclusion of intermediates in the event",
//      &FourBodyDecayConstructor::interopt_, 0, false, false);
//   static SwitchOption interfaceIntermediateOptionAlways
//     (interfaceIntermediateOption,
//      "Always",
//      "Always include the intermediates",
//      1);
//   static SwitchOption interfaceIntermediateOptionNever
//     (interfaceIntermediateOption,
//      "Never",
//      "Never include the intermediates",
//      2);
//   static SwitchOption interfaceIntermediateOptionOnlyIfOnShell
//     (interfaceIntermediateOption,
//      "OnlyIfOnShell",
//      "Only if there are on-shell diagrams",
//      0);
  
//   static Parameter<FourBodyDecayConstructor,double> interfaceMinReleaseFraction
//     ("MinReleaseFraction",
//      "The minimum energy release for a three-body decay, as a "
//      "fraction of the parent mass.",
//      &FourBodyDecayConstructor::minReleaseFraction_, 1e-3, 0.0, 1.0,
//      false, false, Interface::limited);

//   static RefVector<FourBodyDecayConstructor,VertexBase> interfaceExcludedVertices
//     ("ExcludedVertices",
//      "Vertices which are not included in the three-body decayers",
//      &FourBodyDecayConstructor::excludedVector_, -1, false, false, true, true, false);




    // finally make the four-body decay
  }
}

void FourBodyDecayConstructor::
createDecayMode(vector<PrototypeVertexPtr> & diagrams, bool inter) {
  cerr << "!!!!!!!!!!!!!!!!!! MODE !!!!!!!!!!!\n";
  cerr << "Number of diagrams " << diagrams.size() << "\n";
  cerr << diagrams[0]->incoming->PDGName() << " -> ";
  for(OrderedParticles::const_iterator it=diagrams[0]->outPart.begin();
      it!=diagrams[0]->outPart.end();++it)
    cerr << (**it).PDGName() << " ";
  cerr << "\n";
  for(unsigned int iy=0;iy<diagrams.size();++iy)
    cerr << "DIAGRAM " << iy << "\n" << *diagrams[iy] << "\n";
  // incoming particle  
  tPDPtr inpart = diagrams[0]->incoming;
  // outgoing particles
  OrderedParticles outgoing=diagrams[0]->outPart;
  // incoming particle is now unstable
  inpart->stable(false);
  // construct the tag for the decay mode
  string tag = inpart->name() + "->";
  for(OrderedParticles::const_iterator it = outgoing.begin();
      it != outgoing.end(); ++it) {
    if(it!=outgoing.begin()) tag += ",";
    tag += (**it).name();
  }
  tag += ";";
  tDMPtr dm = generator()->findDecayMode(tag);
  if( decayConstructor()->disableDecayMode(tag) ) {
    // If mode alread exists, ie has been read from file, 
    // disable it
    if( dm ) {
      generator()->preinitInterface(dm, "BranchingRatio", "set", "0.0");
      generator()->preinitInterface(dm, "OnOff", "set", "Off");
    }
    return;
  }
  if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
    cerr << "testing calling create " << tag << "\n";
    // create the decayer
    GeneralFourBodyDecayerPtr decayer = createDecayer(diagrams,inter);
    if(!decayer) return;
  }
  else if (dm) {
    assert(false);
  }
  //update CC mode if it exists
  if( inpart->CC() )
    inpart->CC()->synchronize();
}

GeneralFourBodyDecayerPtr 
FourBodyDecayConstructor::createDecayer(vector<PrototypeVertexPtr> & diagrams, 
					bool inter) const {
  if(diagrams.empty()) return GeneralFourBodyDecayerPtr();
  // extract the external particles for the process
  PDPtr incoming = diagrams[0]->incoming;
  // outgoing particles
  vector<PDPtr> outgoing(diagrams[0]->outPart.begin(),
			 diagrams[0]->outPart.end());
  // get the name for the object
  string objectname ("/Herwig/Decays/");
  string classname = DecayerClassName(incoming, diagrams[0]->outPart, objectname);
  if(classname=="") return GeneralFourBodyDecayerPtr();
  // create the object
  GeneralFourBodyDecayerPtr decayer = 
    dynamic_ptr_cast<GeneralFourBodyDecayerPtr>
    (generator()->preinitCreate(classname, objectname));
  // set up the decayer 
  cerr << "testing made the decayer\n";
  decayer->setDecayInfo(incoming,outgoing,diagrams);


//   // get the colour flows
//   unsigned int ncf(0);
//   pair<vector<DVector>, vector<DVector> > cfactors;
//   try {
//     cfactors = getColourFactors(incoming,outgoing,diagrams,ncf);
//   }
//   catch ( Veto ) { return GeneralThreeBodyDecayerPtr(); }
//   // set decayer options from base class
//   setDecayerInterfaces(objectname);
//   // set the width option
//   ostringstream value;
//   value << _widthopt;
//   generator()->preinitInterface(objectname, "WidthOption", "set", value.str());
//   // set the intermediates option
//   ostringstream value2;
//   value2 << inter;
//   generator()->preinitInterface(objectname, "GenerateIntermediates", "set", 
// 				value2.str());
//   // initialize the decayer
//   decayer->init();
//   // return the decayer
//   return decayer;


  assert(false);
}

string  FourBodyDecayConstructor::DecayerClassName(tcPDPtr incoming,
						   const OrderedParticles & outgoing, 
						   string & objname) const {
  string classname("Herwig::");
  // spins of the outgoing particles
  unsigned int ns(0),nf(0),nv(0);
  objname += incoming->PDGName() + "2";
  for(OrderedParticles::const_iterator it=outgoing.begin();
      it!=outgoing.end();++it) {
    if     ((**it).iSpin()==PDT::Spin0    ) ++ns;
    else if((**it).iSpin()==PDT::Spin1Half) ++nf;
    else if((**it).iSpin()==PDT::Spin1    ) ++nv;
    objname += (**it).PDGName();
  }
  objname   += "Decayer";
  if(incoming->iSpin()==PDT::Spin0) {
    if(nf==4) classname += "StoFFFFDecayer";
    else      classname  = "";
  }
  else {
    classname="";
  }
  return classname;
}
