// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MultiPartonExtractor class.
//

#include "MultiPartonExtractor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDF/RemnantHandler.h"
#include "ThePEG/PDF/NoPDF.h"

using namespace Herwig;

IBPtr MultiPartonExtractor::clone() const {
  return new_ptr(*this);
}

IBPtr MultiPartonExtractor::fullclone() const {
  return new_ptr(*this);
}

void MultiPartonExtractor::persistentOutput(PersistentOStream & os) const {
  os << firstPDF_ << secondPDF_;
}

void MultiPartonExtractor::persistentInput(PersistentIStream & is, int) {
  is >> firstPDF_ >> secondPDF_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MultiPartonExtractor,PartonExtractor>
describeHerwigMultiPartonExtractor("Herwig::MultiPartonExtractor", "HwPartonExtractor.so");

void MultiPartonExtractor::Init() {

  static ClassDocumentation<MultiPartonExtractor> documentation
    ("The MultiPartonExtractor class allows more control over the PDFs used"
     " than the default PartonExtractor of ThePEG");
  
  static RefVector<MultiPartonExtractor,PDFBase> interfaceFirstPDFs
    ("FirstPDFs",
     "The PDFs for the first beam",
     &MultiPartonExtractor::firstPDF_, -1, false, false, true, false, false);
  
  static RefVector<MultiPartonExtractor,PDFBase> interfaceSecondPDFs
    ("SecondPDFs",
     "The PDFs for the second beam",
     &MultiPartonExtractor::secondPDF_, -1, false, false, true, false, false);

}

PartonPairVec MultiPartonExtractor::
getPartons(Energy maxEnergy, const cPDPair & incoming,
	   const Cuts & kc) const {
  PartonPairVec result;
  PartonVector first;
  PDFCuts cuts1(kc, true, maxEnergy);
  PBPtr p1 =
    new_ptr(PartonBin(PDPtr(), PBPtr(), incoming.first, PDFPtr(), cuts1));
  std::deque<tcPDFPtr> pdfs(firstPDF_.begin(),firstPDF_.end());
  addPartons(p1, cuts1,  pdfs, first);
  PartonVector second;
  PDFCuts cuts2(kc, false, maxEnergy);
  pdfs = std::deque<tcPDFPtr>(secondPDF_.begin(),secondPDF_.end());
  PBPtr p2 =
    new_ptr(PartonBin(PDPtr(), PBPtr(), incoming.second, PDFPtr(), cuts2));
  addPartons(p2, cuts2, pdfs, second);
  for ( PartonVector::iterator it1 = first.begin();
  	it1 != first.end(); ++it1 )
    for ( PartonVector::iterator it2 = second.begin();
  	it2 != second.end(); ++it2 )
      result.push_back(PBPair(*it1, *it2));

  // We add the original parton bins as well to avoid them being
  // deleted.
  result.push_back(PBPair(p1, p2));
  return result;
}

void MultiPartonExtractor::
addPartons(tPBPtr incoming, const PDFCuts & cuts, std::deque<tcPDFPtr> pdfs,
	   PartonVector & pbins) const {
  tcPDFPtr pdf;
  if(!pdfs.empty()) {
    pdf = *pdfs.begin();
    pdfs.pop_front();
  }
  if(!pdf) pdf = getPDF(incoming->parton());
  if ( dynamic_ptr_cast<Ptr<NoPDF>::tcp>(pdf) ||
       incoming->parton() == incoming->particle() ) {
    pbins.push_back(incoming);
    return;
  }
  cPDVector partons = pdf->partons(incoming->parton());
  for ( int i = 0, N = partons.size(); i < N; ++i ) {
    PBPtr pb =
      new_ptr(PartonBin(incoming->parton(), incoming, partons[i], pdf, cuts));
    incoming->addOutgoing(pb);
    addPartons(pb, cuts, pdfs, pbins);
  }
}
