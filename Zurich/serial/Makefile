CXX=gcc #icx-cc
CFLAGS=  -lm #-std=c11 -DNDEBUG

NSYS=nsys profile
NSYSFLAGS=--stats=true --force-overwrite=true

ONEAPI_ROOT=opt/intel/oneapi
ADV=advisor -collect=survey
APM=$(ONEAPI_ROOT)/advisor/latest/pythonapi/advisor/

PROFILE_DIR:= report

SRCDIR := include
INCDIR := include

C_FILES = $(wildcard ./$(SRCDIR)/*.c)
H_FILES = $(wildcard ./$(INCDIR)/*.h)

EXEC = main

.PHONY:all
all: main
 $(EXEC): $(C_FILES) $(H_FILES)
	$(CXX) 	-pg $(CFLAGS) -I ./$(INCDIR) $(C_FILES) -lm $(EXEC).c -o $(EXEC)


.PHONY: UnitTest
UnitTest: $(C_FILES) $(H_FILES)
	$(CXX) 	$(CFLAGS) -D_DEBUG -pg -I  ./$(INCDIR) $(C_FILES) -lm UnitTest.c -o UnitTest

.PHONY: CheckResult
CheckResult: $(C_FILES) $(H_FILES)	
	$(CXX) 	$(CFLAGS) -D_DEBUG -pg -I  ./$(INCDIR) $(C_FILES) -lm CheckResult.c -o CheckResult

.PHONY: OutputCompare
OutputCompare : plot_output.py
	python3 plot_output.py


.PHONY: TimeExecution #GNU Profiler
TimeExecution: $(EXEC)
	mkdir -p $(PROFILE_DIR)
	gcc $(CFLAGS) -D_PROFILING -pg -I  ./$(INCDIR) $(C_FILES) -lm $(EXEC).c -o $(EXEC) 1> /dev/null
	./$(EXEC)
	gprof ./$(EXEC) gmon.out -b > $(PROFILE_DIR)/report.txt
	gprof2dot $(PROFILE_DIR)/report.txt > $(PROFILE_DIR)/report.dot
	dot -Tpng -o $(PROFILE_DIR)/report.png $(PROFILE_DIR)/report.dot
	python3 graphic.py fun_time.out
	rm gmon.out

# .PHONY: adv_profile
# adv_profile: $(EXEC)
# 	mkdir -p $(PROFILE_DIR)
# 	$(ADV) $(APM)/collect.py advisor_project --config gen9 -- ./$(EXEC)
# 	$(ADV) $(APM)/analyze.py advisor_project --config gen9 --out-dir ./analyze



.PHONY: clean
clean:
	rm -rf $(EXEC) test *.out UnitTest report/*
 
