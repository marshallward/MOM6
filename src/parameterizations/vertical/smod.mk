FC = mpifort
FCFLAGS = -r8
FCFLAGS += -I../../framework
FCFLAGS += -I../../../config_src/memory/dynamic_symmetric
FCFLAGS += -I../../../../../ocean_only/build_arm
FCFLAGS += -I../../../../../shared/fms/build_arm
FCFLAGS += -I/scratch4/GFDL/gfdloceans/Marshall.Ward/aarch64/include

all: MOM_energetic_PBL_int.o

MOM_energetic_PBL_int.o: MOM_energetic_PBL_int.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f MOM_energetic_PBL_int.o mom_energetic_pbl.mod
