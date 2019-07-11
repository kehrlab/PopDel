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

DATE=on 2019-07-11
VERSION=1.1.1
CXXFLAGS+=-DDATE=\""$(DATE)"\" -DVERSION=\""$(VERSION)"\"

# Enable warnings
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result


HEADERS=parse_popdel.h insert_histogram_popdel.h workflow_popdel.h utils_popdel.h popdel_view_parameter_parsing.h
HEADERS+=popdel_profile/*.h
HEADERS+=popdel_call/*.h

.PHONY: all
all: CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
all: popdel

.PHONY: profiling
profiling: CXXFLAGS+=-g -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1
profiling: popdel

.PHONY: debug
debug: CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1
debug: popdel

PREFIX = /usr/local
.PHONY: install
install: CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
install: popdel
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp $< $(DESTDIR)$(PREFIX)/bin/popdel

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/popdel

popdel: popdel.o

popdel.o: popdel.cpp $(HEADERS)

.PHONY: clean
clean:
	rm -f *.o popdel
