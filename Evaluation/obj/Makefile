.SUFFIXES : .o .f90

VPATH = $(SRCGPF):

# Emulator file
seqmpi = me_emulator_$(LES_ARC).f90

SRCSF90 = me_util.f90 me_gp.f90 me_optim.f90 me_training.f90 $(seqmpi) me_ifs_sample.f90 me_ifs_day.f90 me_main.f90

OBJS =	$(SRCSF90:.f90=.o)

all: $(OUT)

$(OUT): $(OBJS)
	$(F90) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

%.o: %.f90 
	$(F90) -c $(FFLAGS) $<

#
# Dependencies
# --------------------
me_util.o:
me_gp.o: me_util.o
me_optim.o: me_gp.o
me_training.o: me_gp.o me_optim.o
$(seqmpi:.f90=.o): me_gp.o
me_ifs_sample.o: $(seqmpi:.f90=.o)
me_ifs_day.o: $(seqmpi:.f90=.o)
me_main.o: me_training.o me_ifs_sample.o me_ifs_day.o
