build:
	ln -s ../../.. Herwig++
	g++ -c CutsInterface.cc
	g++ -c Resolution.cc
	rm -f Herwig++

install:
	if [ ! -d "$(HERWIGPATH)/include/Herwig++/Shower/CKKW/GenerationCuts" ]; then mkdir $(HERWIGPATH)/include/Herwig++/Shower/CKKW/GenerationCuts ; fi
	cp *.h $(HERWIGPATH)/include/Herwig++/Shower/CKKW/GenerationCuts
	cp *.o $(HERWIGPATH)/lib/Herwig++

