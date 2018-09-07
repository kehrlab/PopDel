# Use version 4.9 of g++

# Set dependencies, s.t. make detects changes in header files.
%.o: %.cc %.h
	g++ -c $^

CXX=g++
CC=$(CXX)

# Set this to include SeqAn libraries, either system wide
# or download into current folder and set to .
SEQAN_LIB=.

CXXFLAGS+=-I$(SEQAN_LIB) -DSEQAN_HAS_ZLIB=1 -std=c++14 -DSEQAN_DISABLE_VERSION_CHECK
LDLIBS=-lz -lpthread

DATE=on $(shell git log --pretty=format:"%cd" --date=iso | cut -f 1,2 -d " " | head -n 1)
VERSION=1.0-$(shell git log --pretty=format:"%h" --date=iso | head -n 1)
CXXFLAGS+=-DDATE=\""$(DATE)"\" -DVERSION=\""$(VERSION)"\"

# Enable warnings
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result


HEADERS=parse_popdel.h insert_histogram_popdel.h workflow_popdel.h utils_popdel.h popdel_view_parameter_parsing.h
HEADERS+=popdel_profile/*.h
HEADERS+=popdel_call/*.h

all: CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
all: popdel

profiling: CXXFLAGS+=-g -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1
profiling: popdel

debug: CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1
debug: popdel

popdel: popdel.o

popdel.o: popdel.cpp $(HEADERS)

clean:
	rm -f *.o popdel
