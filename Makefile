
-include Makefile.conf

AMREX_TARGET ?= 
CC ?= mpicxx -cxx=g++
MPI_LIB ?= -lgfortran -lmpich

RESET              = \033[0m
B_ON               = \033[1m
FG_RED             = \033[31m
FG_DIM             = \033[2m
FG_LIGHTRED        = \033[91m
FG_LIGHTGRAY       = \033[90m
FG_GREEN           = \033[32m
FG_LIGHTGREEN      = \033[92m
FG_YELLOW          = \033[33m
FG_LIGHTYELLOW     = \033[93m
FG_BLUE            = \033[34m
FG_LIGHTBLUE       = \033[94m
FG_CYAN            = \033[36m
FG_MAGENTA         = \033[35m




METADATA_GITHASH  = $(shell git log --pretty=format:'%H' -n 1)
METADATA_USER     = $(shell whoami)
METADATA_PLATFORM = $(shell hostname)
METADATA_COMPILER = $(COMP)
METADATA_DATE     = $(shell date +%x)
METADATA_TIME     = $(shell date +%H:%M:%S)
BUILD_DIR         = ${shell pwd}

METADATA_FLAGS = -DMETADATA_GITHASH=\"$(METADATA_GITHASH)\" -DMETADATA_USER=\"$(METADATA_USER)\" -DMETADATA_PLATFORM=\"$(METADATA_PLATFORM)\" -DMETADATA_COMPILER=\"$(METADATA_COMPILER)\" -DMETADATA_DATE=\"$(METADATA_DATE)\" -DMETADATA_TIME=\"$(METADATA_TIME)\" -DBUILD_DIR=\"${BUILD_DIR}\" $(if ${MEME}, -DMEME)


CXX_COMPILE_FLAGS += -Winline -Wpedantic -Wextra -Wall  -std=c++11 $(METADATA_FLAGS)
ifeq ($(DEBUG),TRUE)
 CXX_COMPILE_FLAGS += -ggdb -g3
else 
 CXX_COMPILE_FLAGS += -O3
endif

LINKER_FLAGS += -Bsymbolic-functions

INCLUDE = $(if ${EIGEN}, -isystem ${EIGEN})  $(if ${AMREX}, -isystem ${AMREX}/include/) -I./src/ $(for pth in ${CPLUS_INCLUDE_PATH}; do echo -I"$pth"; done)
LIB     = -L${AMREX}/lib/ -lamrex 

HDR_ALL = $(shell find src/ -name *.H)
HDR_TEST = $(shell find src/ -name *Test.H)
HDR = $(filter-out $(HDR_TEST),$(HDR_ALL))
SRC = $(shell find src/ -mindepth 2  -name "*.cpp" )
SRC_F = $(shell find src/ -mindepth 2  -name "*.F90" )
SRC_MAIN = $(shell find src/ -maxdepth 1  -name "*.cc" )
EXE = $(subst src/,bin/, $(SRC_MAIN:.cc=-$(POSTFIX))) 
OBJ = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC:.cpp=.cpp.o)) 
DEP = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC:.cpp=.cpp.d)) $(subst src/,obj/obj-$(POSTFIX)/, $(SRC_MAIN:.cc=.cc.d))
OBJ_MAIN = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC_MAIN:.cpp=.cc.o))
OBJ_F = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC_F:.F90=.F90.o))



.SECONDARY: 

default: $(DEP) $(EXE)
	@printf "$(B_ON)$(FG_GREEN)DONE $(RESET)\n" 


python: $(OBJ)
	@printf "$(B_ON)$(FG_MAGENTA)PYTHON  $(RESET)    Compiling library\n" 
	@$(CC) -x c++ -c py/alamo.cpy -fPIC -o py/alamo.cpy.o ${INCLUDE} ${PYTHON_INCLUDE} ${CXX_COMPILE_FLAGS} 
	@$(CC) -shared -Wl,-soname,alamo.so -o alamo.so py/alamo.cpy.o ${OBJ} ${LIB} ${MPI_LIB} $(PYTHON_LIB) 

tidy:
	@printf "$(B_ON)$(FG_RED)TIDYING  $(RESET)\n" 
	rm -f Backtrace*
	rm -f amrex.build.log
	
clean: tidy
	@printf "$(B_ON)$(FG_RED)CLEANING  $(RESET)\n" 
	find src/ -name "*.o" -exec rm {} \;
	rm -f bin/*
	rm -rf obj
	rm -f Backtrace*
	rm -rf docs/build docs/doxygen docs/html docs/latex
	rm -f amrex.build.log

realclean: clean
	@printf "$(B_ON)$(FG_RED)CLEANING AMREX $(RESET)\n" 
	-make -C amrex realclean
	@printf "$(B_ON)$(FG_RED)CLEANING OLD CONFIGURATIONS $(RESET)\n" 
	rm -f Makefile.conf Makefile.amrex.conf


info:
	@printf "$(B_ON)$(FG_BLUE)Compiler version information$(RESET)\n"
	@$(CC) --version

-include Makefile.amrex.conf

bin/%: bin/%-$(POSTFIX) ;

bin/%-$(POSTFIX): ${OBJ_F} ${OBJ} obj/obj-$(POSTFIX)/%.cc.o 
	@printf "$(B_ON)$(FG_BLUE)LINKING$(RESET)     $@ \n" 
	@mkdir -p bin/
	@$(CC) -o $@ $^ ${LIB}  ${MPI_LIB}  ${LINKER_FLAGS}

obj/obj-$(POSTFIX)/test.cc.o: src/test.cc ${AMREX_TARGET}
	@printf "$(B_ON)$(FG_YELLOW)COMPILING$(RESET)   $< \n" 
	@mkdir -p $(dir $@)
	@$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} 

obj/obj-$(POSTFIX)/%.cc.o: src/%.cc ${AMREX_TARGET}
	@printf "$(B_ON)$(FG_YELLOW)COMPILING$(RESET)   $< \n" 
	@mkdir -p $(dir $@)
	@$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} 

obj/obj-$(POSTFIX)/%.cpp.o: 
	@printf "$(B_ON)$(FG_YELLOW)COMPILING$(RESET)   $< \n" 
	@mkdir -p $(dir $@)
	@$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} 

obj/obj-$(POSTFIX)/%.cpp.d: src/%.cpp  ${AMREX_TARGET}
	@printf "$(B_ON)$(FG_LIGHTGRAY)DEPENDENCY$(RESET)  $< \n" 
	@mkdir -p $(dir $@)
	@g++ -I./src/ $< ${INCLUDE} ${CXX_COMPILE_FLAGS} -MM -MT $(@:.cpp.d=.cpp.o) -MF $@

obj/obj-$(POSTFIX)/%.cc.d: src/%.cc ${AMREX_TARGET}
	@printf "$(B_ON)$(FG_LIGHTGRAY)DEPENDENCY$(RESET)  $< \n" 
	@mkdir -p $(dir $@)
	@g++ -I./src/ $< ${INCLUDE} ${CXX_COMPILE_FLAGS} -MM -MT $(@:.cc.d=.cc.o) -MF $@

obj/obj-$(POSTFIX)/IO/WriteMetaData.cpp.o: .FORCE ${AMREX_TARGET}
	@printf "$(B_ON)$(FG_LIGHTYELLOW)$(FG_DIM)COMPILING$(RESET)   ${subst obj/obj-$(POSTFIX)/,src/,${@:.cpp.o=.cpp}} \n" 
	@mkdir -p $(dir $@)
	@$(CC) -c ${subst obj/obj-$(POSTFIX)/,src/,${@:.cpp.o=.cpp}} -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} 

.PHONY: .FORCE

FORT_INCL = $(shell for i in ${CPLUS_INCLUDE_PATH//:/ }; do echo -I'$i'; done)

obj/obj-$(POSTFIX)/%.F90.o: src/%.F90 
	@printf "$(B_ON)$(FG_YELLOW)COMPILING  $(RESET)$<\n" 
	@mkdir -p $(dir $@)
	mpif90 -c $< -o $@  -I${subst :, -I,$(CPLUS_INCLUDE_PATH)}
	rm *.mod -rf

docs: docs/doxygen/index.html docs/build/html/index.html .FORCE 
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Done\n" 

docs/doxygen/index.html: $(SRC) $(SRC_F) $(SRC_MAIN) $(HDR_ALL)
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Generating doxygen files\n" 	
	@cd docs && doxygen 
docs/build/html/index.html: $(shell find docs/source/ -type f) Readme.md
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Generating sphinx\n" 	
	@make -C docs html > /dev/null


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
ifneq ($(MAKECMDGOALS),info)
ifneq ($(MAKECMDGOALS),help)
ifneq ($(MAKECMDGOALS),docs)
-include $(DEP)
endif
endif
endif
endif
endif

