SOURCES := Herwig++.cc

OBJECTS := $(SOURCES:.cc=.o)

PROGS   := $(SOURCES:.cc=.exe)

HWLIBS = -lHwHadronization -lHw -lHwDecay -lHwShower -lHwME

KTJETPATH = -I/usera/stephens/projects/ktjet

KTJETLIBS = -lHwKtJet -L$(KTJETPATH)/lib -lKtEvent

.PHONY: Lib

all: Herwig++ testKtJet

.o: .cc
	$(CXX) -c $(CXXFLAGS) $(COMPILER_FLAG) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(DEFS) $(INCLUDES) $(CXXEXT) $<

include ../Config/Makefile.common

CXXDEBUGFLAG = -g

SimpleLEP:
	setupHerwig.exe -init; setupHerwig.exe SimpleLEP.in

Lib:
	cd .. ; $(MAKE) -$(MAKEFLAGS) lib

current: $(PROGS) Herwig++

$(PROGS): ../lib/libHwHadronization.so ../lib/libHw.so

%.exe: %.o 
	export LD_RUN_PATH=$(RPATH); $(CXX) $(CXXEXT) $(LINKFLAGS) $(CXXDEBUGFLAG) -o $@ $< -L$(SRCDIR)/lib $(HWLIBS) $(LIBS)

Herwig++.o: Herwig++.cc
	$(CXX) -c $(CXXFLAGS) $(COMPILER_FLAG) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(DEFS) $(INCLUDES) $(CXXEXT) $<

Herwig++: Herwig++.o
	export LD_RUN_PATH=$(RPATH); $(CXX) $(CXXEXT) $(LINKFLAGS) $(CXXDEBUGFLAG) -o $@ $< -L$(SRCDIR)/lib $(HWLIBS) $(LIBS)

testKtJet: testKtJet.o
	export LD_RUN_PATH=$(RPATH); $(CXX) $(CXXEXT) $(LINKFLAGS) $(CXXDEBUGFLAG) -o $@ $< -L$(SRCDIR)/lib $(HWLIBS) $(LIBS) $(KTJETLIBS)

testKtJet.o: testKtJet.cc
	$(CXX) -c $(CXXFLAGS) $(COMPILER_FLAG) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(DEFS) $(INCLUDES) $(CXXEXT) $(KTJETPATH) $<  
# DO NOT DELETE THIS LINE -- g++dep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.


Herwig++.o: Herwig++.cc ../../Herwig++/Utilities/HerwigRun.h\
	 ../../Pythia7/Repository/EventGenerator.h \
	 ../../Pythia7/Config/Pythia7.h \
	 ../../Pythia7/Config/config.h \
	 ../../Pythia7/Pointer/Ptr.h \
	 ../../Pythia7/Pointer/RCPtr.fh \
	 ../../Pythia7/Pointer/PtrTraits.h \
	 ../../Pythia7/Pointer/RCPtr.h \
	 ../../Pythia7/Pointer/ReferenceCounted.h \
	 ../../Pythia7/Pointer/ReferenceCounted.icc \
	 ../../Pythia7/Pointer/RCPtr.icc \
	 ../../Pythia7/Utilities/Rebinder.fh \
	 ../../Pythia7/Utilities/Interval.fh \
	 ../../Pythia7/Utilities/ClassDescription.fh \
	 ../../Pythia7/Interface/InterfaceBase.fh \
	 ../../Pythia7/Persistency/PersistentOStream.fh \
	 ../../Pythia7/Persistency/PersistentIStream.fh \
	 ../../Pythia7/Config/Units.h \
	 ../../Pythia7/CLHEPWrap/Lorentz5Vector.fh \
	 ../../Pythia7/CLHEPWrap/LorentzVector.fh \
	 ../../Pythia7/CLHEPWrap/ThreeVector.fh \
	 ../../Pythia7/CLHEPWrap/Transverse.fh \
	 ../../Pythia7/CLHEPWrap/SystemOfUnits.h \
	 ../../Pythia7/CLHEPWrap/CLHEP.h \
	 ../../Pythia7/Utilities/UnitIO.h \
	 ../../Pythia7/Config/Constants.h \
	 ../../Pythia7/Config/std.h \
	 ../../Pythia7/Config/Pointers.h \
	 ../../Pythia7/Config/Containers.h \
	 ../../Pythia7/Utilities/FixedSizeAllocator.h \
	 ../../Pythia7/Utilities/Debug.h \
	 ../../Pythia7/Utilities/Debug.icc \
	 ../../Pythia7/Utilities/FixedSizeAllocator.icc \
	 ../../Pythia7/Misc/Named.h \
	 ../../Pythia7/Misc/Named.icc \
	 ../../Pythia7/Repository/EventGenerator.fh \
	 ../../Pythia7/Repository/RandomGenerator.h \
	 ../../Pythia7/Interface/Interfaced.h \
	 ../../Pythia7/Interface/InterfacedBase.h \
	 ../../Pythia7/Utilities/ClassDescription.h \
	 ../../Pythia7/Utilities/ClassTraits.h \
	 ../../Pythia7/Utilities/DescriptionList.h \
	 ../../Pythia7/Utilities/DescriptionList.icc \
	 ../../Pythia7/Utilities/ClassDescription.icc \
	 ../../Pythia7/Utilities/ClassDescription.tcc \
	 ../../Pythia7/Utilities/HoldFlag.h \
	 ../../Pythia7/Utilities/HoldFlag.icc \
	 ../../Pythia7/Interface/InterfacedBase.xh \
	 ../../Pythia7/Interface/InterfaceBase.h \
	 ../../Pythia7/Interface/InterfaceBase.xh \
	 ../../Pythia7/Misc/Exception.h \
	 ../../Pythia7/Misc/Exception.fh \
	 ../../Pythia7/Utilities/sstream.h \
	 ../../Pythia7/Misc/Exception.icc \
	 ../../Pythia7/Interface/InterfaceBase.icc \
	 ../../Pythia7/Interface/InterfacedBase.icc \
	 ../../Pythia7/Interface/Interfaced.icc \
	 ../../Pythia7/CLHEPWrap/RandomEngine.fh \
	 ../../Pythia7/Repository/RandomGenerator.icc \
	 ../../Pythia7/Repository/UseRandom.h \
	 ../../Pythia7/Repository/UseRandom.icc \
	 ../../Pythia7/Repository/CurrentGenerator.fh \
	 ../../Pythia7/Repository/EventGenerator.xh \
	 ../../Pythia7/Repository/EventGenerator.icc \
	 ../../Pythia7/Repository/EventGenerator.tcc \
	 ../../Pythia7/Persistency/PersistentIStream.h \
	 ../../Pythia7/Persistency/InputDescription.h \
	 ../../Pythia7/Persistency/InputDescription.icc \
	 ../../Pythia7/Persistency/PersistentIStream.icc \
	 ../../Pythia7/Persistency/PersistentIStream.tcc \
	 ../../Pythia7/PDT/StandardMatchers.h \
	 ../../Pythia7/PDT/Matcher.h \
	 ../../Pythia7/PDT/MatcherBase.h \
	 ../../Pythia7/PDT/ParticleData.h \
	 ../../Pythia7/PDT/PDT.h \
	 ../../Pythia7/PDT/PDT.icc \
	 ../../Pythia7/CLHEPWrap/PhysicalConstants.h \
	 ../../Pythia7/CLHEPWrap/LorentzVector.h \
	 ../../Pythia7/CLHEPWrap/Transverse.h \
	 ../../Pythia7/CLHEPWrap/Lorentz5Vector.h \
	 ../../Pythia7/Utilities/Math.h \
	 ../../Pythia7/Utilities/Math.icc \
	 ../../Pythia7/Utilities/Direction.h \
	 ../../Pythia7/Utilities/Direction.xh \
	 ../../Pythia7/Utilities/Direction.icc \
	 ../../Pythia7/CLHEPWrap/Lorentz5Vector.icc \
	 ../../Pythia7/CLHEPWrap/Transverse.icc \
	 ../../Pythia7/Utilities/Selector.h \
	 ../../Pythia7/Utilities/Selector.icc \
	 ../../Pythia7/Utilities/Selector.tcc \
	 ../../Pythia7/PDT/WidthGenerator.fh \
	 ../../Pythia7/PDT/MassGenerator.fh \
	 ../../Pythia7/PDT/DecayMode.fh \
	 ../../Pythia7/PDT/ParticleData.icc \
	 ../../Pythia7/EventRecord/Particle.h \
	 ../../Pythia7/EventRecord/EventConfig.h \
	 ../../Pythia7/CLHEPWrap/LorentzRotation.h \
	 ../../Pythia7/CLHEPWrap/LorentzRotation.fh \
	 ../../Pythia7/EventRecord/ColourBase.h \
	 ../../Pythia7/EventRecord/EventInfoBase.h \
	 ../../Pythia7/EventRecord/EventInfoBase.icc \
	 ../../Pythia7/EventRecord/ColourLine.h \
	 ../../Pythia7/EventRecord/ColourSinglet.h \
	 ../../Pythia7/EventRecord/ColourSinglet.icc \
	 ../../Pythia7/EventRecord/ColourLine.icc \
	 ../../Pythia7/EventRecord/ColourBase.icc \
	 ../../Pythia7/EventRecord/SpinBase.h \
	 ../../Pythia7/EventRecord/SpinBase.icc \
	 ../../Pythia7/EventRecord/Particle.icc \
	 ../../Pythia7/EventRecord/Particle.tcc \
	 ../../Pythia7/PDT/MatcherBase.icc \
	 ../../Pythia7/PDT/Matcher.icc \
	 ../../Pythia7/PDT/Matcher.tcc \
	 ../../Pythia7/PDT/EnumParticles.h \
	 ../../Pythia7/PDT/PYDECYDummy.h \
	 ../../Pythia7/PDT/Decayer.h \
	 ../../Pythia7/PDT/Decayer.fh \
	 ../../Pythia7/Handlers/HandlerBase.h \
	 ../../Pythia7/Handlers/HandlerBase.icc \
	 ../../Pythia7/MatrixElement/Amplitude.h \
	 ../../Pythia7/Config/Complex.h \
	 ../../Pythia7/MatrixElement/Amplitude.icc \
	 ../../Pythia7/PDT/Decayer.icc \
	 ../../Pythia7/PDT/PYDECYDummy.icc \
	 ../../Herwig++/Utilities/HwDebug.h \
	 ../../Herwig++/Utilities/SmplHist.h \
	 ../../Herwig++/Utilities/SmplStat.h \
	 ../../Pythia7/Utilities/Timer.h \
	 ../../Pythia7/Utilities/Timer.icc \
	 ../../Pythia7/Utilities/DynamicLoader.h \
	 ../../Pythia7/Utilities/DynamicLoader.icc \
	 ../../Pythia7/EventRecord/Event.h \
	 ../../Pythia7/EventRecord/StandardSelectors.h \
	 ../../Pythia7/EventRecord/SelectorBase.h \
	 ../../Pythia7/EventRecord/SelectorBase.icc \
	 ../../Pythia7/EventRecord/ParticleTraits.h \
	 ../../Pythia7/EventRecord/Event.icc \
	 ../../Pythia7/EventRecord/Collision.h \
	 ../../Pythia7/EventRecord/Collision.icc \
	 ../../Pythia7/EventRecord/Step.h \
	 ../../Pythia7/EventRecord/Step.icc \
	 ../../Pythia7/EventRecord/SubProcess.h \
	 ../../Pythia7/EventRecord/SubProcess.icc \
	 ../../Pythia7/EventRecord/SubProcess.tcc \
	 ../../Pythia7/EventRecord/Step.tcc \
	 ../../Pythia7/EventRecord/Collision.tcc \
	 ../../Pythia7/EventRecord/Event.tcc \
	 ../../Pythia7/Repository/Repository.h \
	 ../../Pythia7/Repository/BaseRepository.h \
	 ../../Pythia7/Repository/BaseRepository.xh \
	 ../../Pythia7/Interface/ClassDocumentation.fh \
	 ../../Pythia7/Repository/BaseRepository.icc \
	 ../../Pythia7/Repository/BaseRepository.tcc \
	 ../../Pythia7/PDT/DecayMode.h \
	 ../../Pythia7/PDT/DecayMode.icc \
	 ../../Pythia7/Repository/Repository.icc \
	 ../../Pythia7/Repository/Repository.tcc

# IF YOU PUT ANYTHING HERE IT WILL GO AWAY
