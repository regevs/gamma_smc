CXX=g++

CXXFLAGS=-std=c++14  -Wall -DNDEBUG -O3  -g -ffast-math -funsafe-math-optimizations -fno-math-errno -mavx2 -mfma -march=native -ftree-vectorize -msse -msse2 -msse3 -pipe -faligned-new -I/home/rs2145/rs2145/software/arb/include -Wno-unused-variable -Wno-strict-aliasing -Wl,-rpath=/home/rs2145/.conda/envs/mamba/envs/snakemake/lib
LDFLAGS=-lboost_program_options -O3 -g -L/home/rs2145/rs2145/software/arb/lib/ -L/home/rs2145/rs2145/software/mpfr-4.1.0/lib  -L/home/rs2145/.conda/envs/mamba/envs/snakemake/lib/  -larb -lflint -lpthread -lgsl -lgslcblas -lm -lhts -lboost_iostreams

# bin/gamma_smc_tests
all: bin/gamma_smc bin/generate_canonical_flow_field

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

bin/%: src/%.o src/io.h src/common.h src/flow_field.h src/gamma_smc.h src/data_processor.h 
	mkdir -p bin
	$(CXX) -o $@ $(LDFLAGS) $<

clean:
	rm -f src/*.o
	rm -f bin/*



