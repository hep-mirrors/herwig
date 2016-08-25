#########################################################
#							#
#  usage: 					      	#
#							#
#    make NAME=LEP JOBS=4 LO=3 NLO=0 POINTS=10000  	#
#							#
#    possible NAME: LEP, LHC-Z, LHC-W			#
#							#
#########################################################

NAME=LEP
JOBS=4
NLO=0
LO=2
POINTS=10000


all: $(NAME)

build$(NAME):
	../../Herwig build Merging/$(LO)LO-$(NLO)NLO-$(NAME).in -y$(JOBS)
	ls
integrate$(NAME): build$(NAME)
	echo "$(shell  for i in $$(seq 0 $$(( $(JOBS)-1 )) );do ../../Herwig integrate --jobid=$$i $(LO)LO-$(NLO)NLO-$(NAME).run & done)"
	ls
run$(NAME): integrate$(NAME)
	../../Herwig run $(LO)LO-$(NLO)NLO-$(NAME).run -N$(POINTS) -j$(JOBS)
	yodamerge $(LO)LO-$(NLO)NLO-$(NAME)-*yoda > $(LO)LO-$(NLO)NLO-$(NAME).yoda

mergeyodas$(NAME): run$(NAME)
	yodamerge $(LO)LO-$(NLO)NLO-$(NAME)-*yoda > $(LO)LO-$(NLO)NLO-$(NAME).yoda
	mkdir -p $(LO)LO-$(NLO)NLO-$(NAME)
	mkdir -p tars

tar$(NAME): mergeyodas$(NAME)
	echo "$(shell mv $$(ls|grep "$(LO)LO-$(NLO)NLO-$(NAME)"| grep "\." |grep -v $(LO)LO-$(NLO)NLO-$(NAME).yoda)  $(LO)LO-$(NLO)NLO-$(NAME))"
	tar -czvf $(LO)LO-$(NLO)NLO-$(NAME).tar.gz $(LO)LO-$(NLO)NLO-$(NAME)
	mv $(LO)LO-$(NLO)NLO-$(NAME).tar.gz tars

$(NAME): build$(NAME) integrate$(NAME) run$(NAME) mergeyodas$(NAME) tar$(NAME)
	rm -fr $(LO)LO-$(NLO)NLO-$(NAME)


plots:
	rivet-mkhtml *LEP.yoda  -o LEP
	rivet-mkhtml *LHC-Z.yoda -o LHC-Z
	rivet-mkhtml *LHC-W.yoda -o LHC-W

