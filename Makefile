#IN1=M_DEF_$(MACH)
IN1=M_DEF_Mac2

include $(IN1) 

POLYLIBo = -L$(LIB_O) -lpolyo
POLYLIBg = -L$(LIB_D) -lpolyg

#CHEMIELIB = -L$(LIB_O) -lchemie
CHEMIELIB = -L$(LIB_O) 

OBJo =  $(POLYLIBo) $(CHEMIELIB)
OBJg =  $(POLYLIBg) $(CHEMIELIB)

all: CHEMIE CHEMIE_g
       
CHEMIE:  PolyDo  
	$(LINK) $(OPTL) $(KINCL_O) -o CHEMIE \
        MAIN_CHEMKIN.f90 $(OP) $(OBJo) $(METIS) $(MUMPS) $(UMF) $(LAPACK) $(SCOTCH) $(SCALAP) $(PARMETIS) $(NETCDF) $(LMPI) $(CL) -lblas -llapack  ;

CHEMIE_g:  PolyDg  
	$(LINK) $(OPTL) $(KINCL_D) -o CHEMIE_g \
        MAIN_CHEMKIN.f90 $(OP) $(OBJg) $(METIS) $(MUMPS) $(UMF) $(LAPACK) $(SCOTCH)  $(SCALAP) $(PARMETIS) $(NETCDF) $(LMPI) $(CL) -lblas -llapack  ;

PolyDo: 
	@make -f Make_src "IN2=$(IN1)" "LIB2=$(LIB_O)" "OPT2=$(OPT_O)" "POLY=polyo" "KINCL=$(KINCL_O)"

PolyDg: 
	@make -f Make_src "IN2=$(IN1)" "LIB2=$(LIB_D)" "OPT2=$(OPT_D)" "POLY=polyg" "KINCL=$(KINCL_D)"

clean:
	rm -f *.o *.mod LIB*/*

RTOL = 1.0D-5
ATOL = 1.0D-7

RowM1 = TSRosW2P.fort
RowM2 = Ros34PW3.fort
RowM3 = TSRosWRodas3.fort

Par_INORG:
	mpirun -np 1 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np1
	mpirun -np 2 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np2
	mpirun -np 3 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np3
	mpirun -np 4 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np4
	mpirun -np 5 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np5
	mpirun -np 6 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np6
	mpirun -np 7 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np7
	mpirun -np 8 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np8
	mpirun -np 9 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np9
	mpirun -np 10 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np10
	mpirun -np 11 ./CHEMIE RUN/INORGcl.run | tee Parallel_Test1/Parallel_INORG_np11

Par_MCMC:
	mpirun -np 1 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np1
	mpirun -np 2 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np2
	mpirun -np 3 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np3
	mpirun -np 4 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np4
	mpirun -np 5 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np5
	mpirun -np 6 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np6
	mpirun -np 7 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np7
	mpirun -np 8 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np8
	mpirun -np 9 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np9
	mpirun -np 10 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np10
	mpirun -np 11 ./CHEMIE RUN/MCMcl.run | tee Parallel_Test1/Parallel_MCM+CAPRAM_np11
	

MCMcl_ROW:
	./CHEMIE RUN/MCMcl.run  1.0d-2 1.0d-3 $(RowM1) | tee OUTm/ROWcl_2Norm_M1/MCMex_2Norm_M1_01
	./CHEMIE RUN/MCMcl.run  1.0d-3 1.0d-5 $(RowM1) | tee OUTm/ROWcl_2Norm_M1/MCMex_2Norm_M1_02
	./CHEMIE RUN/MCMcl.run  1.0d-4 1.0d-5 $(RowM1) | tee OUTm/ROWcl_2Norm_M1/MCMex_2Norm_M1_03
	./CHEMIE RUN/MCMcl.run  1.0d-5 1.0d-7 $(RowM1) | tee OUTm/ROWcl_2Norm_M1/MCMex_2Norm_M1_04
	./CHEMIE RUN/MCMcl.run  1.0d-6 1.0d-7 $(RowM1) | tee OUTm/ROWcl_2Norm_M1/MCMex_2Norm_M1_05


MCMcl_ROW1:
	./CHEMIE RUN/MCMcl_Ros3Pw_03.run | tee MCMcl_Ro3Pw_03
	./CHEMIE RUN/MCMcl_Ros3Pw_04.run | tee MCMcl_Ro3Pw_04

INORGcl:  
	./CHEMIE RUN/INORGcl1.run
	./CHEMIE RUN/INORGcl2.run
	./CHEMIE RUN/INORGcl3.run
	./CHEMIE RUN/INORGcl4.run
	./CHEMIE RUN/INORGcl5.run
	./CHEMIE RUN/INORGcl6.run
	./CHEMIE RUN/INORGcl7.run
	./CHEMIE RUN/INORGcl8.run
	./CHEMIE RUN/INORGcl9.run
	./CHEMIE RUN/INORGcl10.run
	./CHEMIE RUN/INORGcl11.run

INORGex:  
	./CHEMIE RUN/INORGex6.run
	./CHEMIE RUN/INORGex7.run
	./CHEMIE RUN/INORGex8.run
	./CHEMIE RUN/INORGex9.run
	./CHEMIE RUN/INORGex10.run
	./CHEMIE RUN/INORGex11.run

MCMCAPcl:
	./CHEMIE RUN/MCM+CAPRAMcl4.run
	./CHEMIE RUN/MCM+CAPRAMcl5.run
	./CHEMIE RUN/MCM+CAPRAMcl6.run
	./CHEMIE RUN/MCM+CAPRAMcl7.run
	./CHEMIE RUN/MCM+CAPRAMcl8.run
	./CHEMIE RUN/MCM+CAPRAMcl9.run
	./CHEMIE RUN/MCM+CAPRAMcl10.run
	./CHEMIE RUN/MCM+CAPRAMcl11.run

MCMCAPex:
	./CHEMIE RUN/MCM+CAPRAMex1.run
	./CHEMIE RUN/MCM+CAPRAMex2.run
	./CHEMIE RUN/MCM+CAPRAMex3.run
	./CHEMIE RUN/MCM+CAPRAMex4.run
	./CHEMIE RUN/MCM+CAPRAMex5.run
	./CHEMIE RUN/MCM+CAPRAMex6.run
	./CHEMIE RUN/MCM+CAPRAMex7.run
	./CHEMIE RUN/MCM+CAPRAMex8.run
	./CHEMIE RUN/MCM+CAPRAMex9.run
	./CHEMIE RUN/MCM+CAPRAMex10.run
	./CHEMIE RUN/MCM+CAPRAMex11.run

Stratoex:
	./CHEMIE RUN/SmallStratoKPPex1.run
	./CHEMIE RUN/SmallStratoKPPex2.run
	./CHEMIE RUN/SmallStratoKPPex3.run
	./CHEMIE RUN/SmallStratoKPPex4.run
	./CHEMIE RUN/SmallStratoKPPex5.run
	./CHEMIE RUN/SmallStratoKPPex6.run
	./CHEMIE RUN/SmallStratoKPPex7.run
	./CHEMIE RUN/SmallStratoKPPex8.run
	./CHEMIE RUN/SmallStratoKPPex9.run
	./CHEMIE RUN/SmallStratoKPPex10.run
	./CHEMIE RUN/SmallStratoKPPex11.run

Stratocl:
	./CHEMIE RUN/SmallStratoKPPcl1.run
	./CHEMIE RUN/SmallStratoKPPcl2.run
	./CHEMIE RUN/SmallStratoKPPcl3.run
	./CHEMIE RUN/SmallStratoKPPcl4.run
	./CHEMIE RUN/SmallStratoKPPcl5.run
	./CHEMIE RUN/SmallStratoKPPcl6.run
	./CHEMIE RUN/SmallStratoKPPcl7.run
	./CHEMIE RUN/SmallStratoKPPcl8.run
	./CHEMIE RUN/SmallStratoKPPcl9.run
	./CHEMIE RUN/SmallStratoKPPcl10.run
	./CHEMIE RUN/SmallStratoKPPcl11.run
