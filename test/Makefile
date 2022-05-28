RM=rm -f
MPI=mpirun --machinefile ~/.mpi-openmpi-cat  --np 8
COMMON_OPT=--DE-dither=0.2 --DE-crossover-prob=0.9 -R 1

DATA=tl.data folded-multi.data folded3-multi.data folded-refl.data \
    folded-refl-multi.data folded-bc.data hb9cv.data
VRFY=$(DATA:.data=.vrfy)
NEC=$(DATA:.data=.nec)
NECVRFY=$(DATA:.data=.necvrfy)

# Commands and arguments for easier rules
CMD_tl=transmission_line
ARG_tl=-c belden_8295 -z 50-500j -f 28.85 -v optimize
CMD_folded-multi=folded_antenna
ARG_folded-multi=--force-forward --multi optimize
CMD_folded3-multi=folded_3ele_antenna
ARG_folded3-multi=--multi optimize
CMD_folded-refl=folded_bigrefl_antenna
ARG_folded-refl=--relax-swr --force-horizontal --popsize 30 \
    --large-refldist optimize
CMD_folded-refl-multi=folded_bigrefl_antenna
ARG_folded-refl-multi=--relax-swr --force-horizontal --large-refldist \
    --multi --popsize=30 optimize
CMD_folded-bc=folded_bc_antenna
ARG_folded-bc=optimize
CMD_hb9cv=hb9cv_antenna
ARG_hb9cv=--popsize=30 --multi optimize

all: $(NECVRFY) $(VRFY)

%.vrfy: %.data
	$(MPI) $(CMD_$*) $(COMMON_OPT) $(ARG_$*) > $@
	diff $< $@

%.necvrfy: %.data %.nec
	$(CMD_$*) $$(tac $< | sed -n '/Title/{g;p;q};h') necout \
	| diff $(word 2,$^) - > $@

clean:
	$(RM) *.vrfy
