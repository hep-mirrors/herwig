#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"
#include "DC_Hjjj_CC_qBq.h"

using namespace Herwig;

IBPtr DC_Hjjj_CC_qBq::clone() const {
  return new_ptr(*this);
}
IBPtr DC_Hjjj_CC_qBq::fullclone() const {
  return new_ptr(*this);
}

vector<DiagPtr> DC_Hjjj_CC_qBq::getDiagrams() const {
  // get the particle data objects

  vector<DiagPtr> allPossibleDiagrams;
  for (int genU = 0; genU < 2; genU++) {
    for (int genL = 0; genL < 2; genL++) {
      PDPtr qu, qd, ku, kd;
      if (genU == 0){
        qu = theME->getParticleData(ParticleID::u);
        qd = theME->getParticleData(ParticleID::d);
      }
      else if (genU == 1){
        qu = theME->getParticleData(ParticleID::c);
        qd = theME->getParticleData(ParticleID::s);
      }
      else if (genU == 2){
        qu = theME->getParticleData(ParticleID::t);
        qd = theME->getParticleData(ParticleID::b);
      }
      if (genL == 0){
        ku = theME->getParticleData(ParticleID::u);
        kd = theME->getParticleData(ParticleID::d);
      }
      else if (genL == 1){
        ku = theME->getParticleData(ParticleID::c);
        kd = theME->getParticleData(ParticleID::s);
      }
      else if (genL == 2){
        ku = theME->getParticleData(ParticleID::t);
        kd = theME->getParticleData(ParticleID::b);
      }

      tcPDPtr quB = (*qu).CC();
      tcPDPtr qdB = (*qd).CC();

      tcPDPtr quC = SU2Helper::SU2CC(qu);
      tcPDPtr qdC = SU2Helper::SU2CC(qd);
      tcPDPtr quBC = SU2Helper::SU2CC(quB);
      tcPDPtr qdBC = SU2Helper::SU2CC(qdB);

      tcPDPtr kuB = (*ku).CC();
      tcPDPtr kdB = (*kd).CC();

      tcPDPtr kuC = SU2Helper::SU2CC(ku);
      tcPDPtr kdC = SU2Helper::SU2CC(kd);
      tcPDPtr kuBC = SU2Helper::SU2CC(kuB);
      tcPDPtr kdBC = SU2Helper::SU2CC(kdB);
      PDPtr Wplus = theME->getParticleData(ParticleID::Wplus); 
      PDPtr Wminus = theME->getParticleData(ParticleID::Wminus); 
      PDPtr g = theME->getParticleData(ParticleID::g); 
      PDPtr h0 = theME->getParticleData(ParticleID::h0); 
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Wminus, Wminus, kdB, kd, 1, qdBC, 3, kdC, 4, g, 2, h0, -1)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Wminus, Wminus, kd, 3, kdC, 1, qdBC, 5, kdC, 5, g, 2, h0, -2)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, qdB, Wminus, Wminus, kd, 2, qdBC, 4, kdC, 1, g, 3, h0, -3)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Wminus, Wminus, kd, 1, qdBC, 5, qdBC, 3, kdC, 5, g, 2, h0, -4)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Wplus, Wplus, kuB, ku, 1, quBC, 3, kuC, 4, g, 2, h0, -5)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Wplus, Wplus, ku, 3, kuC, 1, quBC, 5, kuC, 5, g, 2, h0, -6)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, quB, Wplus, Wplus, ku, 2, quBC, 4, kuC, 1, g, 3, h0, -7)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Wplus, Wplus, ku, 1, quBC, 5, quBC, 3, kuC, 5, g, 2, h0, -8)));
    }
  }
  return allPossibleDiagrams;
}



Selector<const ColourLines *> DC_Hjjj_CC_qBq::colourGeometries(tcDiagPtr diag) const {
  // colour lines for DC_Hjjj_CC_qBq

  static const ColourLines diag1[1] = { 
    ColourLines("-1 -6, 5 8, 7 -4 -8")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("-1 -6, 4 5 8, 7 -8")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("-1 -8, 5 7, -6 -2 8")
  }; 
  static const ColourLines diag4[1] = { 
    ColourLines("-1 -5 -8, 4 7, -6 8")
  }; 
  static const ColourLines diag5[1] = { 
    ColourLines("-1 -6, 5 8, 7 -4 -8")
  }; 
  static const ColourLines diag6[1] = { 
    ColourLines("-1 -6, 4 5 8, 7 -8")
  }; 
  static const ColourLines diag7[1] = { 
    ColourLines("-1 -8, 5 7, -6 -2 8")
  }; 
  static const ColourLines diag8[1] = { 
    ColourLines("-1 -5 -8, 4 7, -6 8")
  }; 

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( 1.0,  &(diag1[0]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( 1.0,  &(diag2[0]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( 1.0,  &(diag3[0]) );
  } 
  else if( diag->id() == -4 )  {
   sel.insert( 1.0,  &(diag4[0]) );
  } 
  else if( diag->id() == -5 )  {
   sel.insert( 1.0,  &(diag5[0]) );
  } 
  else if( diag->id() == -6 )  {
   sel.insert( 1.0,  &(diag6[0]) );
  } 
  else if( diag->id() == -7 )  {
   sel.insert( 1.0,  &(diag7[0]) );
  } 
  else if( diag->id() == -8 )  {
   sel.insert( 1.0,  &(diag8[0]) );
  } 
  return sel;
}
