CXX=g++

MARCH = native
CXXFLAGS =	-std=c++17  -Wall -DNDEBUG -O3 -g -ffast-math -funsafe-math-optimizations -fno-math-errno \
			-pipe -faligned-new  \
			-mavx2 -march=$(MARCH) -mfma -ftree-vectorize -msse -msse2 -msse3 \
			-Wno-unused-variable -Wno-strict-aliasing 

LDFLAGS =	-O3 -g -lm -lstdc++fs
LDFLAGS_MAIN = -lhts -lzstd
LDFLAGS_FF = -larb -lflint -lpthread -lgsl -lgslcblas

.SUFFIXES:.c .cpp .o
.PHONY:all clean

%.o: %.cpp	
	$(CXX) $(CXXFLAGS) -o $@ -c $<

all: bin/gamma_smc

bin/gamma_smc: src/gamma_smc.o
	mkdir -p bin
	$(CXX) -o $@ $< $(LDFLAGS) $(LDFLAGS_MAIN) 

bin/generate_canonical_flow_field: src/generate_canonical_flow_field.o
	mkdir -p bin
	$(CXX) -o $@ $< $(LDFLAGS) $(LDFLAGS_FF)

clean:
	rm -f src/*.o
	rm -f bin/*

src/generate_canonical_flow_field.o: 
src/gamma_smc.o: src/cxxopts.hpp src/io.h src/common.h src/flow_field.h src/gamma_smc.h src/data_processor.h src/sys.h src/screenoutput.h src/indicators.h



