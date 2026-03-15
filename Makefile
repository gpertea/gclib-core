SEARCHDIRS := -I.

CXX   := $(if $(CXX),$(CXX),g++)
LINKER  := $(if $(LINKER),$(LINKER),g++)
LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)
LIBS := -lz

BASEFLAGS  := -Wall -Wextra -std=c++11 ${SEARCHDIRS} -D_FILE_OFFSET_BITS=64 \
 -D_LARGEFILE_SOURCE -D_REENTRANT -fno-strict-aliasing \
 -fno-exceptions -fno-rtti
STRICT_CHECK_FLAGS := -Wconversion -Wsign-conversion

ifdef STRICT_WARNINGS
BASEFLAGS += $(STRICT_CHECK_FLAGS)
endif

GCCV8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCV8)" "1"
 BASEFLAGS += -Wno-class-memaccess
endif

CXXFLAGS := $(if $(CXXFLAGS),$(BASEFLAGS) $(CXXFLAGS),$(BASEFLAGS))

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  ifneq (,$(findstring static,$(MAKECMDGOALS)))
    LDFLAGS += -static-libstdc++ -static-libgcc
  endif
  CXXFLAGS := -O3 -DNDEBUG $(CXXFLAGS)
else
  ifneq (,$(filter %memcheck %memdebug, $(MAKECMDGOALS)))
     GCCVER49 := $(shell expr `${CXX} -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address
     GCCVER5 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 5)
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
     endif
     CXXFLAGS += -fno-common -fstack-protector
     LIBS := -lasan -lubsan -ldl $(LIBS)
  else
     CXXFLAGS += -g -O0 -DDEBUG -D_DEBUG -DGDEBUG
  endif
endif

RM = rm -f

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

OBJS := GBase.o GArgs.o GFaSeqGet.o gdna.o codons.o gff.o GStr.o GFastaIndex.o

.PHONY : all release debug memcheck memdebug test tests large-tests strict-coords clean

all release debug memcheck memdebug: gclib-test

gclib-test: $(OBJS) gclib-test.o
	${LINKER} ${LDFLAGS} -o $@ $(OBJS) gclib-test.o ${LIBS}

test tests: gclib-test
	@./run_tests.sh

large-tests: gclib-test
	@./run_large_tests.sh

strict-coords:
	@for src in GFaSeqGet.cpp GFastaIndex.cpp gff.cpp gclib-test.cpp; do \
		echo "Checking $$src with $(STRICT_CHECK_FLAGS)"; \
		${CXX} ${BASEFLAGS} $(STRICT_CHECK_FLAGS) -fsyntax-only $$src; \
	done

clean:
	@${RM} gclib-test gclib-test.o $(OBJS)
	@rm -rf tests/wrk
	@${RM} core.*
