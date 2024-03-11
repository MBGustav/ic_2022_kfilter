# Declaring Compilers
CXX=gcc #icx-cc
CFLAGS=-lm -D_PROFILING
ICC = icpx
ICC_FLAGS = -fsycl -xHost -qopenmp -Wno-writable-strings


# Declaring Folders
REPORT_DIRECTORY= report
SRCDIR= include
INCDIR= include


# Project name and other files
EXEC = main
BENCHMARK= bench


# Declaring paths for LIBS
ONEAPI_ROOT?=opt/intel/oneapi

APM=$(ONEAPI_ROOT)/advisor/latest/pythonapi/advisor/
MKL_COPTS = -DMKL_ILP64  -I"${MKLROOT}/include"
MKL_LIBS = -L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lsycl -lOpenCL -lpthread -lm -ldl
DPCPP_OPTS = $(MKL_COPTS) -fsycl-device-code-split=per_kernel $(MKL_LIBS) -fopenmp-targets=spir64

# Here we use Ocl/sycl with device declaration
MKL_PAR = $(MKL_COPTS) -fsycl-device-code-split=per_kernel $(MKL_LIBS)
# And here we use gemm cblas declaration
MKL_SEQ = $(MKL_COPTS) -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl 

C_FILES = $(wildcard ./$(SRCDIR)/*.c)
H_FILES = $(wildcard ./$(INCDIR)/*.h)


.PHONY:all
all: $(C_FILES) $(H_FILES) $(EXEC)
 $(EXEC): $(C_FILES) $(H_FILES)
	$(ICC) $(ICC_FLAGS) $(CFLAGS)  -I ./$(INCDIR) $(C_FILES) -lm $(EXEC).cpp -o $(EXEC) $(DPCPP_OPTS)


.PHONY: UnitTest
UnitTest: $(C_FILES) $(H_FILES)
	$(ICC) UnitTest.cpp -fsycl -I ./include ./include/*.c* -o UnitTest $(MKL_PAR)

.PHONY: CheckResult
CheckResult: $(C_FILES) $(H_FILES)	
	$(CXX) 	$(CFLAGS) -D_DEBUG -pg -I  ./$(INCDIR) $(C_FILES) -lm CheckResult.c -o CheckResult

.PHONY: OutputCompare
OutputCompare : plot_output.py
	python3 plot_output.py


# OMP, Serial and MKL Implementation Test using different parameters
.PHONY: testbench
testbench:
	$(ICC) $(ICC_FLAGS) testbench-omp.cpp $(DPCPP_OPTS) -I ./$(INCDIR) $(C_FILES) -o  testbench-omp_nopad 
	$(ICC) $(ICC_FLAGS) testbench-omp.cpp $(DPCPP_OPTS) -I ./$(INCDIR) $(C_FILES) -o  testbench-omp_pad -DPAD_LD
	$(ICC) $(ICC_FLAGS) testbench-oneMKL.cpp -I ./$(INCDIR) $(C_FILES) $(DPCPP_OPTS) -o testbench-oneMKL_nopad  -DMKL_ONEAPI    $(MKL_LIBS)
	$(ICC) $(ICC_FLAGS) testbench-oneMKL.cpp -I ./$(INCDIR) $(C_FILES) $(DPCPP_OPTS) -o testbench-oneMKL -DMKL_ONEAPI  -DPAD_LD $(MKL_LIBS)
	$(ICC) $(ICC_FLAGS) testbench-cblas.cpp $(DPCPP_OPTS) -I ./$(INCDIR) $(C_FILES) -o  testbench-omp_cblas -DPAD_LD
	$(ICC) $(ICC_FLAGS) testbench-cblas.cpp $(DPCPP_OPTS) -I ./$(INCDIR) $(C_FILES) -o  testbench-omp_cblas
	$(ICC) $(ICC_FLAGS) testbench-oneAPI.cpp -I ./$(INCDIR) $(C_FILES) $(DPCPP_OPTS) -o testbench-oneAPI_nopad  -DMKL_ONEAPI    $(MKL_LIBS)
	$(ICC) $(ICC_FLAGS) testbench-oneAPI.cpp -I ./$(INCDIR) $(C_FILES) $(DPCPP_OPTS) -o testbench-oneAPI -DMKL_ONEAPI  -DPAD_LD $(MKL_LIBS)

	
	$(ICC) testbench-serial.cpp -I ./$(INCDIR) $(C_FILES) -o  testbench-serial 
	


.PHONY: TimeExecution_gnu #GNU Profiler (only works with gcc)
TimeExecution_gnu: $(EXEC)
	mkdir -p $(REPORT_DIRECTORY)
	gcc $(CFLAGS) -D_PROFILING -pg -I  ./$(INCDIR) $(C_FILES) -lm $(EXEC).c -o $(EXEC) 1> /dev/null
	./$(EXEC)
	gprof ./$(EXEC) gmon.out -b > $(REPORT_DIRECTORY)/report.txt
	gprof2dot $(REPORT_DIRECTORY)/report.txt > $(REPORT_DIRECTORY)/report.dot
	dot -Tpng -o $(REPORT_DIRECTORY)/report.png $(REPORT_DIRECTORY)/report.dot
	python3 graphic.py fun_time.out
	rm gmon.out


.PHONY: clean
clean:
	rm -rf $(EXEC) test *.out UnitTest report/* $(BENCHMARK) $(BENCHMARK)_par $(BENCHMARK)_seq 
