/*
  This is interactive.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef INTERACTIVE_H  /* guard against multiple inclusions */
#define INTERACTIVE_H

#include "globals.h"

/******** type declarations *************************************************/

namespace interactive {
  class OutputFile;
};

#include "bits.h"
#include "coxtypes.h"
#include "graph.h"
#include "interface.h"
#include "transducer.h"
#include "type.h"

/******** function declarations **********************************************/

namespace interactive {
  coxgroup::CoxGroup* allocCoxGroup();
  coxgroup::CoxGroup* allocCoxGroup(const type::Type& x);
  void changeOrdering(coxgroup::CoxGroup *W, bits::Permutation& order);
  coxgroup::CoxGroup* coxeterGroup(const type::Type& x, const coxtypes::Rank& l);
  int endOfLine(FILE *f);
  const type::Type& getType();
  graph::CoxEntry getCoxEntry(const coxtypes::Rank& i, const coxtypes::Rank& j);
  coxtypes::CoxArr& getCoxArr(transducer::Transducer& T) /* not implemented */;
  coxtypes::CoxNbr& getCoxNbr(transducer::Transducer& T) /* not implemented */;
  const coxtypes::CoxWord& getCoxWord(coxgroup::CoxGroup *W);
  coxtypes::Generator getGenerator(coxgroup::CoxGroup *W);
  coxtypes::Generator getGenerator(coxgroup::CoxGroup *W, const bits::Lflags& f);
  void getLength
    (list::List<coxtypes::Length>& L, const graph::CoxGraph& G,
     const interface::Interface& I);
  coxtypes::Rank getRank(const type::Type& type);
  void printInterface(FILE* file, const interface::GroupEltInterface& GI,
			const bits::Permutation& a);
  void printInterface(FILE* file, const interface::GroupEltInterface& GI,
		      const interface::GroupEltInterface& WI,
		      const bits::Permutation& a);
  void printMatrix(FILE *file, const coxgroup::CoxGroup* W);
  void printOrdering(FILE* file, const coxgroup::CoxGroup* W);
  void printRepresentation(FILE *file, const coxgroup::CoxGroup* W);
  graph::CoxEntry readCoxEntry
    (const coxtypes::Rank& i, const coxtypes::Rank& j, FILE *inputfile);
  bool yesNo();
};

/* type definitions */

namespace interactive {

class OutputFile {
 private:
  FILE* d_file;
 public:
  OutputFile();
  ~OutputFile();
  FILE* f() {return d_file;}
};

};

#endif
