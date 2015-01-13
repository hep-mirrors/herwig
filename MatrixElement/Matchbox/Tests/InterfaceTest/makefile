HEPMC=/opt/HepMC/2.06.08
HERWIG=/...../bin/Herwig++
BUILDIN = Buildin
INTERFACES = MadGraph  OpenLoops Buildin GoSam


define dep1
echo '$(var)';\
if [ -f '$(var).hepmc' ];then  rm -f $(var).hepmc;fi ;\
cd $(var);\
if [ -f '$(var).run' ]; \
then $(HERWIG) run $(var).run -N100 &> $(var).run.log ;\
fi ;\
cd ..;
endef

define dep2
echo '$1 trying to produce events.';\
if [ ! -d '$1' ]; then mkdir $1;fi;\
cd $1;\
cp ../global.in $1.in;\
sed -i 's/INTERFACE/$1/g' $1.in;\
if [ '$1' = '$(BUILDIN)' ];then sed -i 's/BUILDIN/#/g' $1.in; else sed -i 's/BUILDIN//g' $1.in; fi;\
if [ -f '$1.run' ]; then rm  $1.run;fi;\
$(HERWIG) read $1.in &> $1.read.log ;\
$(HERWIG) run $1.run -N10 &> $1.run.log;\
cd ..;\
if [ -f '$1.hepmc' ]; then echo '$1 produced events';fi;
endef


define dep3
rm -rf $1
endef


all:
	g++ compare-interfaces.cc\
	 -I $(HEPMC)/include/\
	 -L $(HEPMC)/lib64/ -lHepMC\
	 -o compare-interfaces

check:
	@$(foreach var,$(INTERFACES), $(call dep2,$(var)))
	@echo ''
	@echo ''
	@echo 'run make compare to compare the interfaces per PS-point.';

compare:
	@$(foreach var,$(INTERFACES),$(call dep1,$(var)))
	@sleep 1
	./compare-interfaces; 
	

clean_hepmc:
	@rm -f *.hepmc;

clean:
	@if [ -f 'compare-interfaces' ]; then rm  compare-interfaces;fi;
	@rm -f *.hepmc;
	@$(foreach var,$(INTERFACES), $(call dep3,$(var)))
	
settruth:
	cp ${TRUTH} Truth.hepmc
