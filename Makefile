CXX=g++

CXXFLAGS=	-std=c++14  -Wall -DNDEBUG -O3 -g -ffast-math -funsafe-math-optimizations -fno-math-errno \
			-mavx2 -mfma -march=native -ftree-vectorize -msse -msse2 -msse3 -pipe -faligned-new  \
			-Wno-unused-variable -Wno-strict-aliasing 

LDFLAGS=	-lboost_program_options -O3 -g -larb -lflint -lpthread -lgsl -lgslcblas -lm \
			-lhts -lboost_iostreams

.SUFFIXES:.c .cpp .o
.PHONY:all clean

%.o: %.cpp	
	$(CXX) $(CXXFLAGS) -o $@ -c $<

all: bin/gamma_smc bin/generate_canonical_flow_field

bin/%: src/%.o
	mkdir -p bin
	$(CXX) -o $@ $(LDFLAGS) $<

clean:
	rm -f src/*.o
	rm -f bin/*

src/generate_canonical_flow_field.o: 
src/gamma_smc.o: src/io.h src/common.h src/flow_field.h src/gamma_smc.h src/data_processor.h src/sys.h



