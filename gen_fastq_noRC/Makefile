.PHONY: gen_fastq_noRC

all: gen_fastq_noRC

CXX = g++
CXX_FLAGS += -std=c++11 -O3 -DNDEBUG -flto -fwhole-program
CXX_FLAGS += -m64 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -static

gen_fastq_noRC:
	$(CXX) $(CXX_FLAGS) -o $@ gen_fastq_noRC.cpp
	strip $@

clean:
	-rm -f gen_fastq_noRC
