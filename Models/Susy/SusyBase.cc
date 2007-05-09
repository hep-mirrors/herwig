// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SusyBase class.
//

#include "SusyBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/StringUtils.h"

using namespace Herwig;

void SusyBase::doinit() throw(InitException) {
  addVertex(theZSFSFVertex);
   addVertex(thePSFSFVertex);
   addVertex(theWSFSFVertex);
   addVertex(theNFSFVertex);
   addVertex(theGFSFVertex);
   addVertex(theH1SFSFVertex);
   addVertex(theH2SFSFVertex);
   addVertex(theH3SFSFVertex);
   addVertex(theCFSFVertex);
   addVertex(theGSFSFVertex);
   addVertex(theGGSQSQVertex);
   addVertex(theGSGSGVertex);
   addVertex(theNNZVertex);
   addVertex(theCCPVertex);
   addVertex(theCCZVertex);
   addVertex(theCNWVertex);
  StandardModel::doinit();
}

void SusyBase::persistentOutput(PersistentOStream & os) const {
  os << theMinPar << theHMix << theNMix << theUMix << theVMix << theStopMix 
     << theSbotMix << theStauMix << theAlpha << theAu << theAd << theAe
     << theZSFSFVertex << thePSFSFVertex << theWSFSFVertex 
     << theNFSFVertex << theGFSFVertex << theH1SFSFVertex 
     << theH2SFSFVertex << theH3SFSFVertex << theCFSFVertex 
     << theGSFSFVertex << theGGSQSQVertex << theGSGSGVertex 
     << theNNZVertex << theCCPVertex << theCCZVertex << theCNWVertex;
}

void SusyBase::persistentInput(PersistentIStream & is, int) {
  is >> theMinPar >> theHMix >> theNMix >> theUMix >> theVMix >> theStopMix 
     >> theSbotMix >> theStauMix >> theAlpha >> theAu >> theAd >> theAe
     >> theZSFSFVertex >> thePSFSFVertex >> theWSFSFVertex 
     >> theNFSFVertex >> theGFSFVertex >> theH1SFSFVertex 
     >> theH2SFSFVertex >> theH3SFSFVertex >> theCFSFVertex 
     >> theGSFSFVertex >> theGGSQSQVertex >> theGSGSGVertex 
     >> theNNZVertex >> theCCPVertex >> theCCZVertex >> theCNWVertex;
}

ClassDescription<SusyBase> SusyBase::initSusyBase;
// Definition of the static class description member.

void SusyBase::Init() {

  static ClassDocumentation<SusyBase> documentation
    ("This is the base class for any SUSY model.");

  static Reference<SusyBase,Helicity::VSSVertex> interfaceVertexZSS
    ("Vertex/ZSFSF",
     "Reference to Susy ZSSVertex",
     &SusyBase::theZSFSFVertex, false, false, true, false);

  static Reference<SusyBase,Helicity::VSSVertex> interfaceVertexPSS
    ("Vertex/PSFSF",
     "Reference to Susy P SF SF vertex",
     &SusyBase::thePSFSFVertex, false, false, true, false);
  
  static Reference<SusyBase,Helicity::VSSVertex> interfaceVertexWSS
    ("Vertex/WSFSF",
     "Reference to Susy W SF SF vertex",
     &SusyBase::theWSFSFVertex, false, false, true, false);

  
  static Reference<SusyBase,Helicity::FFSVertex> interfaceVertexNFSF
    ("Vertex/NFSF",
     "Reference to the neutralino-fermion-sfermion vertex",
     &SusyBase::theNFSFVertex, false, false, true, false);

  static Reference<SusyBase,Helicity::FFSVertex> interfaceVertexGFSF
    ("Vertex/GFSF",
     "Reference to the gluino-fermion-sfermion vertex",
     &SusyBase::theGFSFVertex, false, false, true, false);
  
  static Reference<SusyBase,Helicity::SSSVertex> interfaceVertexH1SFSF
    ("Vertex/H1SFSF",
     "Reference to the first higgs fermion-sfermion vertex",
     &SusyBase::theH1SFSFVertex, false, false, true, false);

  static Reference<SusyBase,Helicity::SSSVertex> interfaceVertexH2SFSF
  ("Vertex/H2SFSF",
     "Reference to the 2nd higgs fermion-sfermion vertex",
     &SusyBase::theH2SFSFVertex, false, false, true, false);

  static Reference<SusyBase,Helicity::SSSVertex> interfaceVertexH3SFSF
  ("Vertex/H3SFSF",
     "Reference to the 3nd higgs fermion-sfermion vertex",
     &SusyBase::theH3SFSFVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFSVertex> interfaceVertexCFSF
   ("Vertex/CFSF",
      "Reference to the chargino-fermion-sfermion vertex",
      &SusyBase::theCFSFVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::VSSVertex> interfaceVertexGSFSF
   ("Vertex/GSFSF",
      "Reference to the gluon-sfermion-sfermion vertex",
      &SusyBase::theGSFSFVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::VVSSVertex> interfaceVertexGGSS
     ("Vertex/GGSQSQ",
      "Reference to the gluon-gluon-squark-squark vertex",
      &SusyBase::theGGSQSQVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexGSGSG
     ("Vertex/GSGSG",
      "Reference to the gluon-gluino-gluino vertex",
      &SusyBase::theGSGSGVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexNNZ
    ("Vertex/NNZ",
     "Reference to Z-~chi_i0-~chi_i0 vertex",
     &SusyBase::theNNZVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexCCP
    ("Vertex/CCP",
     "Reference to ~chi_i+-~chi_i-photon vertex",
     &SusyBase::theCCPVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexCCZ
    ("Vertex/CCZ",
     "Reference to ~chi_i+-~chi_i-Z vertex",
     &SusyBase::theCCZVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexCNW
    ("Vertex/CNW",
     "Reference to ~chi_i+-chi_i0-W vertex",
     &SusyBase::theCNWVertex, false, false, true, false);

}

void SusyBase::readSetup(istream &is) throw(SetupException) {
  string name = dynamic_ptr_cast<istringstream*>(&is)->str();
  ifstream file(name.c_str());
  if(!file)
     throw SetupException() << "A problem occurred in opening the file: "
			    << name;
  string line;
  //function pointer for putting all characters to lower case.
  int (*pf)(int) = tolower;
  //Create a directory for the mixing matrices
  BaseRepository::CreateDirectory("/Defaults/MixingMatrices/");
  //Start to read in the file
  while(getline(file, line)) {
    if(line[0] == '#') continue;
    transform(line.begin(), line.end(), line.begin(), pf);
    if(line.find("alpha") != string::npos) {
      getline(file, line);
      istringstream iss(line);
      iss >> theAlpha;
      continue;      
    }
    if(line.find("block minpar") != string::npos)
      theMinPar = readBlock(file);
    if(line.find("block  hmix") != string::npos)
      theHMix = readBlock(file);
    if(line.find("block mass") != string::npos)
      theMasses = readBlock(file);
    if(line.find("block decay") != string::npos)
      readDecay(file, line);
    if(line.find("block nmix") != string::npos || 
       line.find("block umix") != string::npos ||
       line.find("block vmix") != string::npos || 
       line.find("block stopmix") != string::npos ||
       line.find("block sbotmix") != string::npos || 
       line.find("block staumix") != string::npos) {
      unsigned int size(0);
      string name = StringUtils::car(StringUtils::cdr(line), " #");
      vector<Complex> vals = readMatrix(file, size);
      createMixingMatrix(name, vals, size);
    }
    if(line.find("block au") != string::npos) {
      unsigned int dummy;
      theAu = readMatrix(file, dummy);
    }
    if(line.find("block ad") != string::npos) {
      unsigned int dummy;
      theAd = readMatrix(file, dummy);
    }
    if(line.find("block ae") != string::npos) {
      unsigned int dummy;
      theAe = readMatrix(file, dummy);
    }
  }
  //This has to be done afterwards because the mixing matrices
  //must have been read in first
  resetRepositoryMasses();
}

SusyBase::ParamMap
SusyBase::readBlock(ifstream & ifs) const throw(SetupException) {
  if(!ifs)
    throw SetupException() << "The input stream is in a bad state"
			   << Exception::runerror;
  string line;
  ParamMap store;
  while(getline(ifs, line)) {
    if(line[0] == '#') continue;
    istringstream is(line);
    long index;
    double value;   
    is >> index >> value;
    store.insert(make_pair(index, value));
    if(ifs.peek() == 'B' || ifs.peek() == 'D' ||
       ifs.peek() == 'b' || ifs.peek() == 'd' ||
       ifs.peek() == '#') break;
  }
  return store;
}

void SusyBase::readDecay(ifstream & ifs, 
			    string decay) const throw(SetupException){
  if(!ifs)
    throw SetupException() <<"The input stream is in a bad state";
  istringstream iss(decay);
  string dummy;
  long parent;
  double width;
  iss >> dummy >> parent >> width;
  string line;
  while(getline(ifs, line)) {
    if(line[0] == '#') continue;
    istringstream is(line);
    double brat(0.);
    unsigned int nda(0);
    is >> brat >> nda;
    vector<long> products(nda);
    vector<long>::size_type i(0);
    while(is && ++i < nda + 1)
      is >> products[i - 1];
    createDecayMode(parent, products, brat);
    if(ifs.peek() == 'B' || ifs.peek() == 'D' ||
       ifs.peek() == 'b' || ifs.peek() == 'd' || 
       ifs.peek() == '#') break;
  }
}

const vector<Complex>
SusyBase::readMatrix(ifstream & ifs, unsigned int & size) throw(SetupException) {
  if(!ifs)
    throw SetupException() << "The input stream is in a bad state";
  string line;
  unsigned int rowmax(0), colmax(0);
  vector<Complex> values;
  while(getline(ifs, line)) {
    if(line[0] == '#') continue;
    istringstream is(line);
    unsigned int index1, index2;
    double real(0.), imag(0.);   
    is >> index1 >> index2 >> real >> imag;
    values.push_back(Complex(real, imag));
    if(index1 > rowmax) rowmax = index1; 
    if(index2 > colmax) colmax = index2; 
    if(ifs.peek() == 'B' || ifs.peek() == 'D' ||
       ifs.peek() == 'b' || ifs.peek() == 'd' ||
       ifs.peek() == '#') {
      if(rowmax != colmax)
	throw SetupException() << "A square matrix has unequal row/col size!";
      size = rowmax;
      break;
    }
  }
  return values;
}

void SusyBase::createDecayMode(long parent, vector<long> products,
				  double brat) const {
  ostringstream cmd;
  cmd << "decaymode " + getParticleData(parent)->PDGName() + "->";
  vector<long>::size_type nda = products.size();
  for(vector<long>::size_type i = 0; i < nda; ++i) {
    cmd << getParticleData(products[i])->PDGName();
    if(i != nda - 1) cmd << ",";
    else cmd << "; ";
  }
  cmd << brat << " 1 /Defaults/Decays/Mambo";
  Repository::exec(cmd.str(), cerr);
}

void SusyBase::createMixingMatrix(string name, vector<Complex> & values,
				     unsigned int size) {
  string dir = "/Defaults/MixingMatrices/";
  MixingMatrixPtr mix = new_ptr(MixingMatrix(size));
  for(unsigned int r = 0; r < size; ++r)
    for(unsigned int c = 0; c < size; ++c)
      (*mix)(r, c) = values[size*r + c];
  if(name == "nmix") {
    theNMix = mix;
    vector<long> ids(4);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035; 
    (*mix).setIds(ids);
  }
  else if(name == "umix") {
    theUMix = mix;
    vector<long> ids(2);
    ids[0] = 1000024; ids[1] = 1000037; 
    (*mix).setIds(ids);
  }
  else if(name == "vmix") {
    theVMix = mix;
    vector<long> ids(2);
    ids[0] = 1000024; ids[1] = 1000037; 
    (*mix).setIds(ids);
  }
  else if(name == "stopmix") {
    theStopMix = mix;
    vector<long> ids(2);
    ids[0] = 1000006; ids[1] = 2000006; 
    (*mix).setIds(ids);
  }
  else if(name == "sbotmix") {
    theSbotMix = mix;
    vector<long> ids(2);
    ids[0] = 1000005; ids[1] = 2000005; 
    (*mix).setIds(ids);
  }
  else if(name == "staumix") {
    theStauMix = mix;
    vector<long> ids(2);
    ids[0] = 1000015; ids[1] = 2000015; 
    (*mix).setIds(ids);
  }
  else
    throw SetupException() << "Cannot find correct title for mixing matrix "
			   << name << Exception::runerror;
  BaseRepository::Register(mix, dir + name);
}

void SusyBase::resetRepositoryMasses() {
  for(ParamMap::iterator it = theMasses.begin(); it != theMasses.end();){
    long id = it->first;
    double mass = it->second;
    //a negative mass requires an adjustment to the associated mixing matrix 
    if(mass < 0.0) adjustMixingMatrix(id);
    ostringstream cmd;
    cmd << "set /Defaults/Particles/" 
	<< getParticleData(id)->PDGName() 
	<< ":NominalMass " << abs(it->second);
    if(!cmd) 
      throw SetupException() << "SusyBase::resetRepositoryMasses - "
			     << "The stream constructed to reset particle "
			     << "masses is not usuable. " 
			     << Exception::runerror;
    Repository::exec(cmd.str(), cerr);
    //no need to store all masses after this stage
    theMasses.erase(it++);
  }
}

void SusyBase::adjustMixingMatrix(long id) {
  //get correct mixing matrix
  switch(id) {
  case 1000022 :
  case 1000023 :
  case 1000025 :
  case 1000035 : 
    if(theNMix)
      theNMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The neutralino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    break;
  case 1000024 :
  case 1000037 : 
    if(theUMix)
      theUMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The U-Type chargino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    if(theVMix)
      theVMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The V-Type chargino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    break;
  case 1000006 :
  case 2000006 :
    if(theStopMix)
      theStopMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The stop mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  case 1000005 :
  case 2000005 :
    theSbotMix->adjustPhase(id);
if(theStopMix)
      theStopMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The stop mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  case 1000015 :
  case 2000015 :
    if(theStopMix)
      theStopMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The stop mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  default : 
    throw SetupException() << "SusyBase::adjustMixingMatrix - "
			   << "PDG code does not correspond to anything that "
			   << "has a mixing matrix associated with it"
			   << Exception::runerror;
  }
}
