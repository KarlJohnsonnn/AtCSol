IN1=M_DEF_Mac

include $(IN1) 

OBJo = -L$(LIB_O) -lchemieo
OBJg = -L$(LIB_D) -lchemieg

all: AtCSol AtCSol_dbg
       
AtCSol:  optimized  
	$(LINK) $(OPT_O) $(KINCL_O) -o AtCSol.exe \
        AtCSol.f90 SOLVER/opkdmain.o SOLVER/opkda1.o SOLVER/opkda2.o  \
				$(OP) $(OBJo) $(LAPACK) $(NETCDF) $(LMPI) $(CL) $(COARRAY)  ;

AtCSol_dbg:  debug
	$(LINK) $(OPT_D) $(KINCL_D) -o AtCSol_dbg.exe \
        AtCSol.f90 SOLVER/opkdmain.o SOLVER/opkda1.o SOLVER/opkda2.o  \
				$(OP) $(OBJg) $(LAPACK) $(NETCDF) $(LMPI) $(CL) $(COARRAY) ;

optimized: 
	cd SOLVER/ ; \
	gfortran -c $(OPT_O) opkda1.f opkda2.f opkdmain.f ; \
	cd .. ; \
	make -f Make_src "IN2=$(IN1)" "LIB2=$(LIB_O)" "OPT2=$(OPT_O)" "CHEM=chemieo" "KINCL=$(KINCL_O)"

debug: 
	cd SOLVER/ ; \
	gfortran -c $(OPT_D) opkda1.f opkda2.f opkdmain.f ; \
	cd .. ; \
	make -f Make_src "IN2=$(IN1)" "LIB2=$(LIB_D)" "OPT2=$(OPT_D)" "CHEM=chemieg" "KINCL=$(KINCL_D)"

clean:
	rm -f *.o *.mod LIB*/* *.exe

distclean:
	rm -f *.exe

test:
	./AtCSol.exe RUN/SmallStratoKPP.run
	./AtCSol.exe RUN/RACM+C24.run
	./AtCSol.exe RUN/MCM+CAPRAM.run
	./AtCSol.exe RUN/ERC_nheptane.run
	./AtCSol.exe RUN/LLNL_nHeptane.run
	./AtCSol.exe RUN/LLNL_MD.run
