CXX=g++

CXXFLAGS =	-std=c++14  -Wall -DNDEBUG -O3 -g -ffast-math -funsafe-math-optimizations -fno-math-errno \
			-mavx2 -mfma -march=native -ftree-vectorize -msse -msse2 -msse3 -pipe -faligned-new  \
			-Wno-unused-variable -Wno-strict-aliasing 

LDFLAGS =	-O3 -g -lm -lboost_iostreams -lboost_program_options -lboost_system -lboost_filesystem
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
src/gamma_smc.o: src/io.h src/common.h src/flow_field.h src/gamma_smc.h src/data_processor.h src/sys.h src/screenoutput.h src/indicators.h



