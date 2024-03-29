# sources contains a list of the source files (i.e., the .cpp files)
sources := $(patsubst %.cpp,%.cpp,$(wildcard *.cpp))
# there is one .o file for each .cpp file
objects := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
# this is used to automatically generate the dependencies
dependencies := $(patsubst %.cpp,%.d,$(wildcard *.cpp))

globals = globals.h

pflags = $(includedirs) -pg -O3
oflags = $(includedirs) -O3 -Wall
gflags = $(includedirs) -Wall -ggdb

cflags = $(gflags) # the default setting

ifeq ($(optimize),true)
	NDEBUG = true
	cflags = $(oflags)
endif

ifdef profile
	cflags = $(pflags)
endif

CXX = g++ -std=c++11

all: coxeter #clean

coxeter: $(objects)
	$(CXX) $(cflags) -o coxeter $(objects)

clean:
	rm -f $(objects)

%.o:%.cpp
	$(CXX) -c $(cflags) $*.cpp

# dependencies --- these were generated automatically by make depend on my
# system; they are explicitly copied for portability. Only local dependencies
# are considered. If you add new #include directives you should add the
# corresponding dependencies here; the best way if your compiler supports the
# -MM option is probably to simply say make depend > tmp, then paste in the
# contents of tmp in lieu of the dependencies listed here.

%.d:%.cpp
	@$(CXX) -MM $*.cpp
depend: $(dependencies)

affine.o: affine.cpp affine.h globals.h coxgroup.h coxtypes.h constants.h \
 io.h list.h memory.h list.hpp error.h files.h hecke.h interface.h \
 automata.h bits.h containers.h bitmap.h minroots.h graph.h type.h \
 dotval.h transducer.h schubert.h sl_list.h hecke.hpp polynomials.h \
 polynomials.hpp invkl.h klsupport.h search.h search.hpp kl.h uneqkl.h \
 wgraph.h files.hpp cells.h
automata.o: automata.cpp automata.h globals.h bits.h containers.h \
 memory.h constants.h bitmap.h io.h list.h list.hpp error.h
bitmap.o: bitmap.cpp bitmap.h constants.h globals.h containers.h memory.h \
 bits.h io.h list.h list.hpp error.h
bits.o: bits.cpp bits.h globals.h containers.h memory.h constants.h \
 bitmap.h io.h list.h list.hpp error.h
cells.o: cells.cpp cells.h globals.h bits.h containers.h memory.h \
 constants.h bitmap.h io.h list.h list.hpp error.h kl.h coxtypes.h \
 klsupport.h polynomials.h polynomials.hpp schubert.h graph.h type.h \
 sl_list.h interface.h automata.h minroots.h dotval.h transducer.h \
 hecke.h hecke.hpp uneqkl.h wgraph.h stack.h stack.hpp
commands.o: commands.cpp commands.h globals.h dictionary.h sl_list.h \
 containers.h memory.h constants.h io.h list.h list.hpp error.h \
 dictionary.hpp coxgroup.h coxtypes.h files.h hecke.h interface.h \
 automata.h bits.h bitmap.h minroots.h graph.h type.h dotval.h \
 transducer.h schubert.h hecke.hpp polynomials.h polynomials.hpp invkl.h \
 klsupport.h search.h search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h \
 directories.h fcoxgroup.h help.h interactive.h special.h typeA.h
constants.o: constants.cpp constants.h globals.h
coxgroup.o: coxgroup.cpp coxgroup.h globals.h coxtypes.h constants.h io.h \
 list.h memory.h list.hpp error.h files.h hecke.h interface.h automata.h \
 bits.h containers.h bitmap.h minroots.h graph.h type.h dotval.h \
 transducer.h schubert.h sl_list.h hecke.hpp polynomials.h \
 polynomials.hpp invkl.h klsupport.h search.h search.hpp kl.h uneqkl.h \
 wgraph.h files.hpp cells.h
coxtypes.o: coxtypes.cpp coxtypes.h globals.h constants.h io.h list.h \
 memory.h list.hpp error.h
error.o: error.cpp error.h globals.h coxtypes.h constants.h io.h list.h \
 memory.h list.hpp directories.h graph.h bits.h containers.h bitmap.h \
 type.h interactive.h interface.h automata.h minroots.h dotval.h \
 transducer.h kl.h klsupport.h polynomials.h polynomials.hpp schubert.h \
 sl_list.h hecke.h hecke.hpp version.h
fcoxgroup.o: fcoxgroup.cpp fcoxgroup.h globals.h coxgroup.h coxtypes.h \
 constants.h io.h list.h memory.h list.hpp error.h files.h hecke.h \
 interface.h automata.h bits.h containers.h bitmap.h minroots.h graph.h \
 type.h dotval.h transducer.h schubert.h sl_list.h hecke.hpp \
 polynomials.h polynomials.hpp invkl.h klsupport.h search.h search.hpp \
 kl.h uneqkl.h wgraph.h files.hpp cells.h
files.o: files.cpp files.h globals.h hecke.h list.h memory.h constants.h \
 list.hpp error.h interface.h automata.h bits.h containers.h bitmap.h \
 io.h coxtypes.h minroots.h graph.h type.h dotval.h transducer.h \
 schubert.h sl_list.h hecke.hpp polynomials.h polynomials.hpp invkl.h \
 klsupport.h search.h search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h \
 directories.h posets.h version.h
general.o: general.cpp general.h globals.h coxgroup.h coxtypes.h \
 constants.h io.h list.h memory.h list.hpp error.h files.h hecke.h \
 interface.h automata.h bits.h containers.h bitmap.h minroots.h graph.h \
 type.h dotval.h transducer.h schubert.h sl_list.h hecke.hpp \
 polynomials.h polynomials.hpp invkl.h klsupport.h search.h search.hpp \
 kl.h uneqkl.h wgraph.h files.hpp cells.h
graph.o: graph.cpp graph.h globals.h list.h memory.h constants.h list.hpp \
 error.h bits.h containers.h bitmap.h io.h coxtypes.h type.h sl_list.h \
 directories.h interactive.h interface.h automata.h minroots.h dotval.h \
 transducer.h
help.o: help.cpp help.h globals.h commands.h dictionary.h sl_list.h \
 containers.h memory.h constants.h io.h list.h list.hpp error.h \
 dictionary.hpp coxgroup.h coxtypes.h files.h hecke.h interface.h \
 automata.h bits.h bitmap.h minroots.h graph.h type.h dotval.h \
 transducer.h schubert.h hecke.hpp polynomials.h polynomials.hpp invkl.h \
 klsupport.h search.h search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h \
 directories.h
interactive.o: interactive.cpp interactive.h globals.h bits.h \
 containers.h memory.h constants.h bitmap.h io.h list.h list.hpp error.h \
 coxtypes.h graph.h type.h interface.h automata.h minroots.h dotval.h \
 transducer.h affine.h coxgroup.h files.h hecke.h schubert.h sl_list.h \
 hecke.hpp polynomials.h polynomials.hpp invkl.h klsupport.h search.h \
 search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h directories.h \
 fcoxgroup.h general.h typeA.h
interface.o: interface.cpp interface.h globals.h automata.h bits.h \
 containers.h memory.h constants.h bitmap.h io.h list.h list.hpp error.h \
 coxtypes.h minroots.h graph.h type.h dotval.h transducer.h
invkl.o: invkl.cpp invkl.h globals.h containers.h memory.h constants.h \
 coxtypes.h io.h list.h list.hpp error.h klsupport.h bitmap.h \
 polynomials.h polynomials.hpp schubert.h graph.h bits.h type.h sl_list.h \
 interface.h automata.h minroots.h dotval.h transducer.h hecke.h \
 hecke.hpp search.h search.hpp
io.o: io.cpp io.h globals.h list.h memory.h constants.h list.hpp error.h
kl.o: kl.cpp kl.h globals.h coxtypes.h constants.h io.h list.h memory.h \
 list.hpp error.h containers.h klsupport.h bitmap.h polynomials.h \
 polynomials.hpp schubert.h graph.h bits.h type.h sl_list.h interface.h \
 automata.h minroots.h dotval.h transducer.h hecke.h hecke.hpp iterator.h
klsupport.o: klsupport.cpp klsupport.h globals.h coxtypes.h constants.h \
 io.h list.h memory.h list.hpp error.h bitmap.h containers.h \
 polynomials.h polynomials.hpp schubert.h graph.h bits.h type.h sl_list.h \
 interface.h automata.h minroots.h dotval.h transducer.h
main.o: main.cpp constants.h globals.h commands.h dictionary.h sl_list.h \
 containers.h memory.h io.h list.h list.hpp error.h dictionary.hpp \
 coxgroup.h coxtypes.h files.h hecke.h interface.h automata.h bits.h \
 bitmap.h minroots.h graph.h type.h dotval.h transducer.h schubert.h \
 hecke.hpp polynomials.h polynomials.hpp invkl.h klsupport.h search.h \
 search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h version.h
memory.o: memory.cpp memory.h globals.h constants.h error.h
minroots.o: minroots.cpp minroots.h globals.h graph.h list.h memory.h \
 constants.h list.hpp error.h bits.h containers.h bitmap.h io.h \
 coxtypes.h type.h dotval.h
posets.o: posets.cpp posets.h globals.h bits.h containers.h memory.h \
 constants.h bitmap.h io.h list.h list.hpp error.h wgraph.h interface.h \
 automata.h coxtypes.h minroots.h graph.h type.h dotval.h transducer.h \
 sl_list.h
schubert.o: schubert.cpp schubert.h globals.h coxtypes.h constants.h io.h \
 list.h memory.h list.hpp error.h graph.h bits.h containers.h bitmap.h \
 type.h sl_list.h interface.h automata.h minroots.h dotval.h transducer.h
special.o: special.cpp special.h globals.h commands.h dictionary.h \
 sl_list.h containers.h memory.h constants.h io.h list.h list.hpp error.h \
 dictionary.hpp coxgroup.h coxtypes.h files.h hecke.h interface.h \
 automata.h bits.h bitmap.h minroots.h graph.h type.h dotval.h \
 transducer.h schubert.h hecke.hpp polynomials.h polynomials.hpp invkl.h \
 klsupport.h search.h search.hpp kl.h uneqkl.h wgraph.h files.hpp cells.h \
 directories.h interactive.h
transducer.o: transducer.cpp transducer.h globals.h coxtypes.h \
 constants.h io.h list.h memory.h list.hpp error.h graph.h bits.h \
 containers.h bitmap.h type.h
typeA.o: typeA.cpp typeA.h globals.h fcoxgroup.h coxgroup.h coxtypes.h \
 constants.h io.h list.h memory.h list.hpp error.h files.h hecke.h \
 interface.h automata.h bits.h containers.h bitmap.h minroots.h graph.h \
 type.h dotval.h transducer.h schubert.h sl_list.h hecke.hpp \
 polynomials.h polynomials.hpp invkl.h klsupport.h search.h search.hpp \
 kl.h uneqkl.h wgraph.h files.hpp cells.h
type.o: type.cpp type.h globals.h io.h list.h memory.h constants.h \
 list.hpp error.h
uneqkl.o: uneqkl.cpp uneqkl.h globals.h containers.h memory.h constants.h \
 coxtypes.h io.h list.h list.hpp error.h hecke.h interface.h automata.h \
 bits.h bitmap.h minroots.h graph.h type.h dotval.h transducer.h \
 schubert.h sl_list.h hecke.hpp polynomials.h polynomials.hpp klsupport.h \
 interactive.h
wgraph.o: wgraph.cpp wgraph.h globals.h containers.h memory.h constants.h \
 list.h list.hpp error.h bits.h bitmap.h io.h interface.h automata.h \
 coxtypes.h minroots.h graph.h type.h dotval.h transducer.h sl_list.h
