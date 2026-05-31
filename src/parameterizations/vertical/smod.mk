FC = mpifort
FCFLAGS = -r8
FCFLAGS += -I../../framework
FCFLAGS += -I../../../config_src/memory/dynamic_symmetric
FCFLAGS += -I../../../../../ocean_only/build_arm
FCFLAGS += -I../../../../../shared/fms/build_arm
FCFLAGS += -I/scratch4/GFDL/gfdloceans/Marshall.Ward/aarch64/include

all: MOM_energetic_PBL.o MOM_energetic_PBL_smod.o

MOM_energetic_PBL.o: MOM_energetic_PBL.F90
	$(FC) $(FCFLAGS) -c $<

MOM_energetic_PBL_smod.o: MOM_energetic_PBL_smod.F90 mom_energetic_pbl.mod
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f MOM_energetic_PBL.o MOM_energetic_PBL_smod.o \
	  mom_energetic_pbl.mod mom_energetic_pbl-mom_energetic_pbl_impl.mod
