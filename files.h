/*
  This is files.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef FILES_H  /* guard against multiple inclusions */
#define FILES_H

#include "globals.h"

/******** type declarations *************************************************/

namespace files {

  enum Header { bettiH, basisH, closureH, dufloH, extremalsH, ihBettiH,
		lCOrderH, lCellsH, lCellWGraphsH, lWGraphH, lrCOrderH,
		lrCellsH, lrCellWGraphsH, lrWGraphH, rCOrderH, rCellsH,
		rCellWGraphsH, rWGraphH, slocusH, sstratificationH,
		numHeaders};

  struct AddHeckeTraits;
  struct HeckeTraits;
  struct OutputTraits;
  struct PolynomialTraits;
  struct PosetTraits;
  struct PartitionTraits;
  struct WgraphTraits;
};

/******** function definitions **********************************************/

#include "hecke.h"
#include "invkl.h"
#include "kl.h"
#include "uneqkl.h"
#include "wgraph.h"

namespace files {
  // do _not_ use namespace kl! creates conflicts in coxgroup
};

namespace files {
template <class C>
  void appendCoefficient(std::string& str, const C& c,
			PolynomialTraits& traits);
template <class E>
  void appendExponent(std::string& str, const E& e, PolynomialTraits& traits);
template <class M>
  void appendHeckeMonomial(std::string& str, const M& m, const schubert::SchubertContext& p,
			   const interface::Interface& I, HeckeTraits& hTraits,
			   PolynomialTraits& pTraits, const coxtypes::Length& l);
void appendHomology(std::string& str, const schubert::Homology& h, OutputTraits& traits);
template <class C>
  void appendMonomial(std::string& str, const C& c, const Ulong& e,
		      PolynomialTraits& traits,
		      const Ulong& d = 1, const long& m = 0);
void appendModifier(std::string& str, const Ulong& d, const long& m,
		    PolynomialTraits& traits);
template <class M>
  void appendMuMark(std::string& str, const M& m, const schubert::SchubertContext& p,
		    const coxtypes::Length& l, HeckeTraits& traits);
template <class P>
  void appendPolynomial(std::string& str, const P& p,
			PolynomialTraits& traits,
			const Ulong& d = 1, const long& m = 0);
void appendSeparator(std::string& str, const Ulong& n, HeckeTraits& traits);
template <class KL>
  void makeWGraph
    (wgraph::WGraph& X, const list::List<coxtypes::CoxNbr>& c,
     const Lflags& f, KL& kl);
void minReps
  (list::List<coxtypes::CoxNbr>& min, const bits::Partition& pi,
   schubert::NFCompare& c);
void pad(std::string& str, const Ulong& n, HeckeTraits& traits);
template<class H>
  void printAsBasisElt(FILE* file, const H& h, const schubert::SchubertContext& p,
		       interface::Interface& I, OutputTraits& traits);
void printBetti(FILE* file, const coxtypes::CoxNbr& y, const schubert::SchubertContext& p,
		OutputTraits& traits);
void printCellOrder(FILE* file, const wgraph::OrientedGraph& X,
		    const schubert::SchubertContext& p, const interface::Interface& I,
		    PosetTraits& traits);
void printCoatoms
  (FILE* file, const coxtypes::CoxNbr& y, const schubert::SchubertContext& p,
		  const interface::Interface& I, OutputTraits& traits);
template <class KL>
  void printClosure(FILE* file, const coxtypes::CoxNbr& y, KL& kl, const interface::Interface& I,
		    OutputTraits& traits);
template <class C>
  void printCoefficient(FILE* file, const C& c,
			PolynomialTraits& traits);
void printDescents(FILE* file, const Lflags& df, const Lflags& f,
		   const interface::Interface& I, WgraphTraits& traits);
template <class KL>
  void printDuflo(FILE* file, const list::List<coxtypes::CoxNbr>& d, const bits::Partition& pi,
		  KL& kl, const interface::Interface& I, OutputTraits& traits);
void printEltData(FILE* file, const coxtypes::CoxNbr& y, const schubert::SchubertContext& p,
		  const interface::Interface& I, OutputTraits& traits);
template <class E>
  void printExponent(FILE* file, const E& e, PolynomialTraits& traits);
template <class KL>
  void printExtremals(FILE* file, const coxtypes::CoxNbr& y, const KL& kl,
		      const interface::Interface& I, OutputTraits& traits);
void printHeader(FILE* file, const Header& header, OutputTraits& traits);
template <class H>
  void printHeckeElt(FILE* file, const H& h, const schubert::SchubertContext& p,
		     const interface::Interface& I, OutputTraits& traits,
		     const coxtypes::Length& l = coxtypes::undef_length);
template <class H>
  void printHeckeElt(FILE* file, const H& h, const bits::Permutation& a,
		     const schubert::SchubertContext& p, const interface::Interface& I,
		     HeckeTraits& hTraits,
		     PolynomialTraits& pTraits,
		     const coxtypes::Length& l = coxtypes::undef_length);
void printHomology(FILE* file, const schubert::Homology& h, OutputTraits& traits);
template <class KL>
  void printIHBetti(FILE* file, const coxtypes::CoxNbr& y, KL& kl, OutputTraits& traits);
template <class KL>
  void printLCOrder(FILE* file, KL& kl, const interface::Interface& I,
		    OutputTraits& traits);
template <class KL>
  void printLCells(FILE* file, const bits::Partition& lp, KL& kl, const interface::Interface& I,
		   OutputTraits& traits);
template <class KL>
  void printLCellWGraphs(FILE* file, const bits::Partition& lp, KL& kl,
			 const interface::Interface& I, OutputTraits& traits);
template <class KL>
  void printLRCOrder(FILE* file, KL& kl, const interface::Interface& I,
		     OutputTraits& traits);
template <class KL>
  void printLRCells(FILE* file, const bits::Partition& lp, KL& kl,
		    const interface::Interface& I, OutputTraits& traits);
template <class KL>
  void printLRCellWGraphs(FILE* file, const bits::Partition& lp, KL& kl,
			  const interface::Interface& I, OutputTraits& traits);
template <class KL>
  void printLRWGraph(FILE* file, KL& kl, const interface::Interface& I,
		     OutputTraits& traits);
template <class KL>
  void printLWGraph(FILE* file, KL& kl, const interface::Interface& I,
		    OutputTraits& traits);
template <class C>
  void printMonomial(FILE* file, const C& c, const Ulong& e,
		     PolynomialTraits& traits,
		     const Ulong& d = 1, const long& m = 0);
void printModifier(FILE* file, const Ulong& d, const long& m,
		   PolynomialTraits& traits);
template <class M>
  void printMuMark(FILE* file, const M& m, const schubert::SchubertContext& p,
		   const coxtypes::Length& l, HeckeTraits& traits);
void printPartition(FILE* file, const bits::Partition& pi, const schubert::SchubertContext& p,
		    const interface::Interface& I, PartitionTraits& traits);
template <class P>
  void printPolynomial(FILE* file, const P& p, PolynomialTraits& traits,
		       const Ulong& d = 1, const long& m = 0);
template <class KL>
  void printRCOrder(FILE* file, KL& kl, const interface::Interface& I,
		    OutputTraits& traits);
template <class KL>
  void printRCells(FILE* file, const bits::Partition& lp, KL& kl, const interface::Interface& I,
		   OutputTraits& traits);
template <class KL>
  void printRCellWGraphs(FILE* file, const bits::Partition& lp, KL& kl,
			 const interface::Interface& I, OutputTraits& traits);
template <class KL>
  void printRWGraph(FILE* file, KL& kl, const interface::Interface& I,
		    OutputTraits& traits);
void printSeparator(FILE* file, const Ulong& n, HeckeTraits& traits);
template <class KL>
  void printSingularLocus(FILE* file, const coxtypes::CoxNbr& y, KL& kl,
			  const interface::Interface& I, OutputTraits& traits);
template <class KL>
  void printSingularStratification(FILE* file, const coxtypes::CoxNbr& y, KL& kl,
				   const interface::Interface& I, OutputTraits& traits);
void printWGraph(FILE* file, const wgraph::WGraph& X, const Lflags& f,
		 const interface::Interface& I, WgraphTraits& traits);
template <class KL>
  void printWGraphList(FILE* file, const bits::Partition& pi, const Lflags& f,
		       const interface::Interface& I, KL& kl, OutputTraits& traits);
template <class H>
  bool setTwoSided(const H& h, const bits::Permutation& a,
		   const schubert::SchubertContext& p,
		   const interface::Interface& I, HeckeTraits& hTraits,
		   PolynomialTraits& pTraits,
		   const coxtypes::Length& l = coxtypes::undef_length);
void sortLists(list::List<list::List<coxtypes::CoxNbr> >& lc, schubert::NFCompare& nfc,
	       bits::Permutation& a);
void writeClasses
  (list::List<list::List<coxtypes::CoxNbr> >& lc, const bits::Partition& pi);
};

/******** type definitions **************************************************/

namespace files {

struct PolynomialTraits {
  std::string prefix;
  std::string postfix;
  std::string indeterminate;
  std::string sqrtIndeterminate;
  std::string posSeparator;
  std::string negSeparator;
  std::string product;
  std::string exponent;
  std::string expPrefix;
  std::string expPostfix;
  std::string zeroPol;
  std::string one;
  std::string negOne;
  std::string modifierPrefix;
  std::string modifierPostfix;
  std::string modifierSeparator;
  bool printExponent;
  bool printModifier;
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(PolynomialTraits));}
  PolynomialTraits(io::Pretty);
  PolynomialTraits(io::Terse);
  PolynomialTraits(io::GAP);
  ~PolynomialTraits();
};

struct HeckeTraits {
  std::string prefix;
  std::string postfix;
  std::string evenSeparator;
  std::string oddSeparator;
  std::string monomialPrefix;
  std::string monomialPostfix;
  std::string monomialSeparator;
  std::string muMark;
  std::string hyphens;
  Ulong lineSize;
  Ulong indent;
  Ulong evenWidth;
  Ulong oddWidth;
  char padChar;
  bool doShift;
  bool reversePrint;
  bool twoSided;
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(HeckeTraits));}
  HeckeTraits(const interface::Interface& I, io::Pretty);
  HeckeTraits(const interface::Interface& I, io::Terse);
  HeckeTraits(const interface::Interface& I, io::GAP);
  virtual ~HeckeTraits();
};

struct AddHeckeTraits:public HeckeTraits { // Hecke traits for additive output
  interface::GroupEltInterface* eltTraits;
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(AddHeckeTraits));}
  AddHeckeTraits(const interface::Interface& I, io::Pretty);
  AddHeckeTraits(const interface::Interface& I, io::Terse);
  AddHeckeTraits(const interface::Interface& I, io::GAP);
  ~AddHeckeTraits();
};

struct PartitionTraits {
  std::string prefix;
  std::string postfix;
  std::string separator;
  std::string classPrefix;
  std::string classPostfix;
  std::string classSeparator;
  std::string classNumberPrefix;
  std::string classNumberPostfix;
  bool printClassNumber;
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(PartitionTraits));}
  PartitionTraits(io::Pretty);
  PartitionTraits(io::Terse);
  PartitionTraits(io::GAP);
  ~PartitionTraits();
};

struct PosetTraits {
  std::string prefix;
  std::string postfix;
  std::string separator;
  std::string edgePrefix;
  std::string edgePostfix;
  std::string edgeSeparator;
  std::string nodePrefix;
  std::string nodePostfix;
  Ulong nodeShift;
  bool printNode;
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(PosetTraits));}
  PosetTraits(io::Pretty);
  PosetTraits(io::Terse);
  PosetTraits(io::GAP);
  ~PosetTraits();
};

struct WgraphTraits {
  std::string prefix;
  std::string postfix;
  std::string separator;
  std::string edgeListPrefix;
  std::string edgeListPostfix;
  std::string edgeListSeparator;
  std::string edgePrefix;
  std::string edgePostfix;
  std::string edgeSeparator;
  std::string nodePrefix;
  std::string nodePostfix;
  std::string nodeSeparator;
  std::string nodeNumberPrefix;
  std::string nodeNumberPostfix;
  Ulong nodeShift;
  int padSize;
  bool hasPadding;
  bool printNodeNumber;
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(WgraphTraits));}
  WgraphTraits(io::Pretty);
  WgraphTraits(io::Terse);
  WgraphTraits(io::GAP);
  ~WgraphTraits();
};

struct OutputTraits {
// strings
  std::string version_string;
  std::string type_string;
  // header file names
  std::string header[numHeaders];
  std::string prefix[numHeaders];
  std::string postfix[numHeaders];
  bool hasHeader[numHeaders];
  // prettyfying strings for printouts
  std::string closureSeparator1;
  std::string closureSeparator2;
  std::string closureSeparator3;
  std::string closureSeparator4;
  std::string closureSeparator5;
  std::string closureSeparator6;
  std::string eltList;
  std::string singularLocus;
  std::string singularStratification;
  std::string emptySingularLocus;
  std::string emptySingularStratification;
  // list formatting
  std::string bettiPrefix;
  std::string bettiPostfix;
  std::string bettiSeparator;
  std::string bettiRankPrefix;
  std::string bettiRankPostfix;
  std::string cellNumberPrefix;
  std::string cellNumberPostfix;
  std::string closureSizePrefix;
  std::string closureSizePostfix;
  std::string coatomPrefix;
  std::string coatomPostfix;
  std::string coatomSeparator;
  std::string compCountPrefix;
  std::string compCountPostfix;
  std::string dufloPrefix;
  std::string dufloPostfix;
  std::string dufloSeparator;
  std::string dufloListPrefix;
  std::string dufloListPostfix;
  std::string dufloListSeparator;
  std::string dufloNumberPrefix;
  std::string dufloNumberPostfix;
  std::string eltNumberPrefix;
  std::string eltNumberPostfix;
  std::string eltListPrefix;
  std::string eltListPostfix;
  std::string eltListSeparator;
  std::string eltPrefix;
  std::string eltPostfix;
  std::string eltDataPrefix;
  std::string eltDataPostfix;
  std::string graphListPrefix;
  std::string graphListPostfix;
  std::string graphListSeparator;
  std::string lDescentPrefix;
  std::string lDescentPostfix;
  std::string rDescentPrefix;
  std::string rDescentPostfix;
  std::string lengthPrefix;
  std::string lengthPostfix;
  std::string close_string;
  std::string bettiHyphens;
  Ulong lineSize;
// traits for the output of a polynomial
  PolynomialTraits polTraits;
// traits for the output of a Hecke element
  HeckeTraits heckeTraits;
  AddHeckeTraits addHeckeTraits;
// traits for the output of a partition
  PartitionTraits partitionTraits;
// traits for the output of a W-graph
  WgraphTraits wgraphTraits;
// traits for the output of a poset
  PosetTraits posetTraits;
// flags
  bool printBettiRank;
  bool printCellNumber;
  bool printClosureSize;
  bool printCoatoms;
  bool printCompCount;
  bool printDufloNumber;
  bool printEltDescents;
  bool printElt;
  bool printEltData;
  bool printEltNumber;
  bool printLength;
  bool printType;
  bool printVersion;
  bool hasBettiPadding;
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void* operator new(size_t size, void* ptr) {return ptr;}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(OutputTraits));}
  void operator delete(void* p1, void* p2) {};
  OutputTraits(const graph::CoxGraph& G, const interface::Interface& I, io::Pretty);
  OutputTraits(const graph::CoxGraph& G, const interface::Interface& I, io::Terse);
  OutputTraits(const graph::CoxGraph& G, const interface::Interface& I, io::GAP);
  ~OutputTraits();
// manipulators
  void setBasisTraits(HeckeTraits& hTraits);
  void setDefaultTraits(HeckeTraits& hTraits);
};

};

/******** inline definitions *************************************************/

#include "files.hpp"

#endif
