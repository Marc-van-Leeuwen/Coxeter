/*
  This is files.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "files.h"

#include <algorithm>

#include "cells.h"
#include "directories.h"
#include "posets.h"
#include "version.h"

namespace files {
  using namespace directories;
  using namespace wgraph;
}

/****************************************************************************

  This module contains code for the handling of output to files. It is our
  strong belief that the usefulness of a program depends greatly on the
  quality and flexibility of its output --- unfortunately there is no uniform
  measure of quality : a human reader will want a compact and well-formatted
  file, whereas a computer will want a file that is easy to parse, and a
  computer algebra system will want a file that is in its own format.

  At each point in time, the group has a component d_fileoutput, which is
  set to the current "style" of output (typically we envision pretty, terse,
  gap, tex and latex as desirable styles). This component contains the
  virtual function table corresponding to the desired output.

  The functions outputting to files are the following :

    - printClosure(file,h,p,I) : prints information about one schubert
      variety; the list of extremal pairs is contained in h;

*****************************************************************************/

namespace {
  using namespace files;
  Ulong maxLength(const schubert::Homology& h);
};

namespace {
  using namespace files;

  void makeVersionString(std::string& vstr, const std::string& str);
  void makeTypeString
    (std::string& tstr, const std::string& str, const graph::CoxGraph& G);
};

/****************************************************************************

        Chapter I -- The OutputTraits class.

  This class holds all the little details which govern output of data to
  files. Although this may seem a lot of detail, it seems that not much less
  can be done is one is to automatically produce output suitable for direct
  input into programs such as TeX or GAP.

  Mostly, constructors are provided for the various output tags :

    - OutputTraits(Pretty);
    - OutputTraits(Terse);
    - OutputTraits(GAP);

*****************************************************************************/

namespace files {

OutputTraits::OutputTraits(const graph::CoxGraph& G, const interface::Interface& I, io::Pretty)
  :version_string(""),
   type_string(""),
   closureSeparator1("P_{x,y} for x extremal w.r.t. y:\n\n"),
   closureSeparator2(""),
   closureSeparator3("rational singular locus:\n\n"),
   closureSeparator4("rational singular stratification:\n\n"),
   closureSeparator5("betti numbers:\n\n"),
   closureSeparator6("IH betti numbers:\n\n"),
   eltList("context :\n\n"),
   singularLocus("singular locus :\n\n"),
   singularStratification("rational singular stratification :\n\n"),
   emptySingularLocus("rational singular locus is empty\n"),
   emptySingularStratification("rational singular stratification is empty\n"),
   bettiPrefix(""),
   bettiPostfix(""),
   bettiSeparator(""),
   bettiRankPrefix("h["),
   bettiRankPostfix("] = "),
   cellNumberPrefix(""),
   cellNumberPostfix(" : "),
   closureSizePrefix("size : "),
   closureSizePostfix(""),
   coatomPrefix("coatoms :\n\n"),
   coatomPostfix("\n"),
   coatomSeparator("\n"),
   compCountPrefix("components : "),
   compCountPostfix(""),
   dufloPrefix(""),
   dufloPostfix(""),
   dufloSeparator(":"),
   dufloListPrefix(""),
   dufloListPostfix(""),
   dufloListSeparator("\n"),
   dufloNumberPrefix(""),
   dufloNumberPostfix(" : "),
   eltNumberPrefix(""),
   eltNumberPostfix(":"),
   eltListPrefix(""),
   eltListPostfix(""),
   eltListSeparator("\n"),
   eltPrefix("y = "),
   eltPostfix(""),
   eltDataPrefix(""),
   eltDataPostfix("\n"),
   graphListPrefix(""),
   graphListPostfix(""),
   graphListSeparator("\n\n"),
   lDescentPrefix(" L:"),
   lDescentPostfix(""),
   rDescentPrefix(" R:"),
   rDescentPostfix(""),
   lengthPrefix(" length "),
   lengthPostfix(""),
   close_string("\n"),
   bettiHyphens("h"),
   lineSize(io::LINESIZE),
   polTraits(io::Pretty()),
   heckeTraits(I,io::Pretty()),
   addHeckeTraits(I,io::Pretty()),
   partitionTraits(io::Pretty()),
   wgraphTraits(io::Pretty()),
   posetTraits(io::Pretty()),
   printBettiRank(true),
   printCellNumber(true),
   printClosureSize(true),
   printCoatoms(true),
   printCompCount(true),
   printDufloNumber(true),
   printEltDescents(true),
   printElt(true),
   printEltData(true),
   printEltNumber(true),
   printLength(true),
   printType(false),
   printVersion(false),
   hasBettiPadding(true)

{
  prefix[basisH] = "";
  prefix[bettiH] = "";
  prefix[closureH] = "";
  prefix[dufloH] = "";
  prefix[extremalsH] = "";
  prefix[ihBettiH] = "";
  prefix[lCOrderH] = "";
  prefix[lCellsH] = "";
  prefix[lCellWGraphsH] = "";
  prefix[lWGraphH] = "graph :\n\n";
  prefix[lrCOrderH] = "";
  prefix[lrCellsH] = "";
  prefix[lrCellWGraphsH] = "";
  prefix[lrWGraphH] = "graph :\n\n";
  prefix[rCOrderH] = "";
  prefix[rCellsH] = "";
  prefix[rCellWGraphsH] = "";
  prefix[rWGraphH] = "graph :\n\n";
  prefix[slocusH] = "";
  prefix[sstratificationH] = "";

  postfix[basisH] = "\n";
  postfix[bettiH] = "\n";
  postfix[closureH] = "\n";
  postfix[dufloH] = "\n";
  postfix[extremalsH] = "\n";
  postfix[ihBettiH] = "\n";
  postfix[lCOrderH] = "\n";
  postfix[lCellsH] = "\n";
  postfix[lCellWGraphsH] = "\n";
  postfix[lWGraphH] = "\n";
  postfix[lrCOrderH] = "\n";
  postfix[lrCellsH] = "\n";
  postfix[lrCellWGraphsH] = "\n";
  postfix[lrWGraphH] = "\n";
  postfix[rCOrderH] = "\n";
  postfix[rCellsH] = "\n";
  postfix[rCellWGraphsH] = "\n";
  postfix[rWGraphH] = "\n";
  postfix[slocusH] = "\n";
  postfix[sstratificationH] = "\n";

  hasHeader[bettiH] = false;
  hasHeader[basisH] = false;
  hasHeader[closureH] = false;
  hasHeader[dufloH] = false;
  hasHeader[extremalsH] = false;
  hasHeader[ihBettiH] = false;
  hasHeader[lCOrderH] = false;
  hasHeader[lCellsH] = false;
  hasHeader[lCellWGraphsH] = false;
  hasHeader[lWGraphH] = false;
  hasHeader[lrCOrderH] = false;
  hasHeader[lrCellsH] = false;
  hasHeader[lrCellWGraphsH] = false;
  hasHeader[lrWGraphH] = false;
  hasHeader[rCOrderH] = false;
  hasHeader[rCellsH] = false;
  hasHeader[rCellWGraphsH] = false;
  hasHeader[rWGraphH] = false;
  hasHeader[slocusH] = false;
  hasHeader[sstratificationH] = false;
}

OutputTraits::OutputTraits(const graph::CoxGraph& G, const interface::Interface& I, io::Terse)
  :version_string(""),
   type_string(""),
   closureSeparator1("# extremal pairs\n"),
   closureSeparator2(""),
   closureSeparator3("# rational singular locus\n"),
   closureSeparator4("# rational singular stratification\n"),
   closureSeparator5("# betti numbers\n"),
   closureSeparator6("# IH betti numbers\n"),
   eltList("# context enumeration\n"),
   singularLocus("# rational singular locus\n"),
   singularStratification("# rational singular stratification\n"),
   emptySingularLocus("# rational singular locus is empty"),
   emptySingularStratification("# rational singular stratification is empty"),
   bettiPrefix(""),
   bettiPostfix(""),
   bettiSeparator(","),
   dufloPrefix(""),
   dufloPostfix(""),
   dufloSeparator(":"),
   dufloListPrefix(""),
   dufloListPostfix(""),
   dufloListSeparator("\n"),
   eltListPrefix(""),
   eltListPostfix(""),
   eltListSeparator("\n"),
   eltPrefix(""),
   eltPostfix(""),
   eltDataPrefix("# the element y\n"),
   eltDataPostfix(""),
   graphListPrefix(""),
   graphListPostfix(""),
   graphListSeparator("\n#\n"),
   close_string(""),
   polTraits(io::Terse()),
   heckeTraits(I,io::Terse()),
   addHeckeTraits(I,io::Terse()),
   partitionTraits(io::Terse()),
   wgraphTraits(io::Terse()),
   posetTraits(io::Terse()),
   printBettiRank(false),
   printCellNumber(false),
   printClosureSize(false),
   printCoatoms(false),
   printCompCount(false),
   printDufloNumber(false),
   printEltDescents(false),
   printElt(true),
   printEltData(true),
   printEltNumber(false),
   printLength(false),
   printType(true),
   printVersion(true),
   hasBettiPadding(false)

{
  prefix[basisH] = "";
  prefix[bettiH] = "";
  prefix[closureH] = "";
  prefix[dufloH] = "";
  prefix[extremalsH] = "";
  prefix[ihBettiH] = "";
  prefix[lCOrderH] = "";
  prefix[lCellsH] = "";
  prefix[lCellWGraphsH] = "";
  prefix[lWGraphH] = "# graph\n";
  prefix[lrCOrderH] = "";
  prefix[lrCellsH] = "";
  prefix[lrCellWGraphsH] = "";
  prefix[lrWGraphH] = "# graph\n";
  prefix[rCOrderH] = "";
  prefix[rCellsH] = "";
  prefix[rCellWGraphsH] = "";
  prefix[rWGraphH] = "# graph\n";
  prefix[slocusH] = "";
  prefix[sstratificationH] = "";

  postfix[basisH] = "";
  postfix[bettiH] = "";
  postfix[closureH] = "";
  postfix[dufloH] = "";
  postfix[extremalsH] = "";
  postfix[ihBettiH] = "";
  postfix[lCOrderH] = "";
  postfix[lCellsH] = "";
  postfix[lCellWGraphsH] = "";
  postfix[lWGraphH] = "";
  postfix[lrCOrderH] = "";
  postfix[lrCellsH] = "";
  postfix[lrCellWGraphsH] = "";
  postfix[lrWGraphH] = "";
  postfix[rCOrderH] = "";
  postfix[rCellsH] = "";
  postfix[rCellWGraphsH] = "";
  postfix[rWGraphH] = "";
  postfix[slocusH] = "";
  postfix[sstratificationH] = "";

  header[basisH] = "terse_basis";
  header[closureH] = "terse_closure";
  header[dufloH] = "terse_duflo";
  header[extremalsH] = "terse_extremals";
  header[lCOrderH] = "terse_lcorder";
  header[lCellsH] = "terse_lcells";
  header[lCellWGraphsH] = "terse_lcellwgraphs";
  header[lWGraphH] = "terse_lwgraph";
  header[lrCOrderH] = "terse_lrcorder";
  header[lrCellsH] = "terse_lrcells";
  header[lrCellWGraphsH] = "terse_lrcellwgraphs";
  header[lrWGraphH] = "terse_lrwgraph";
  header[rCOrderH] = "terse_rcorder";
  header[rCellsH] = "terse_rcells";
  header[rCellWGraphsH] = "terse_rcellwgraphs";
  header[rWGraphH] = "terse_rwgraph";
  header[slocusH] = "terse_slocus";
  header[sstratificationH] = "terse_sstratification";

  hasHeader[bettiH] = false;
  hasHeader[basisH] = true;
  hasHeader[closureH] = true;
  hasHeader[dufloH] = true;
  hasHeader[extremalsH] = true;
  hasHeader[ihBettiH] = false;
  hasHeader[lCOrderH] = true;
  hasHeader[lCellsH] = true;
  hasHeader[lCellWGraphsH] = true;
  hasHeader[lWGraphH] = true;
  hasHeader[lrCOrderH] = true;
  hasHeader[lrCellsH] = true;
  hasHeader[lrCellWGraphsH] = true;
  hasHeader[lrWGraphH] = true;
  hasHeader[rCOrderH] = true;
  hasHeader[rCellsH] = true;
  hasHeader[rCellWGraphsH] = true;
  hasHeader[rWGraphH] = true;
  hasHeader[slocusH] = true;
  hasHeader[sstratificationH] = true;

  makeVersionString(version_string,"#");
  makeTypeString(type_string,"#",G);
}

OutputTraits::OutputTraits(const graph::CoxGraph& G, const interface::Interface& I, io::GAP)
  :version_string(""),
   type_string(""),
   closureSeparator1(""),
   closureSeparator2(""),
   closureSeparator3(""),
   closureSeparator4(""),
   closureSeparator5(""),
   closureSeparator6(""),
   eltList("coxeter_contextEnumeration:="),
   singularLocus("coxeter_slocus:="),
   singularStratification("coxeter_sstratification:="),
   emptySingularLocus("coxeter_slocus:=[];"),
   emptySingularStratification("coxeter_sstratification:=[];"),
   bettiPrefix("["),
   bettiPostfix("]"),
   bettiSeparator(","),
   dufloPrefix("["),
   dufloPostfix("]"),
   dufloSeparator(","),
   dufloListPrefix("[\n"),
   dufloListPostfix("]"),
   dufloListSeparator(",\n"),
   eltListPrefix("[\n"),
   eltListPostfix("]"),
   eltListSeparator(",\n"),
   eltPrefix("coxeter_currentElement:="),
   eltPostfix(";"),
   eltDataPrefix(""),
   eltDataPostfix(""),
   graphListPrefix("[\n"),
   graphListPostfix("]"),
   graphListSeparator(",\n"),
   close_string(";"),
   polTraits(io::GAP()),
   heckeTraits(I,io::GAP()),
   addHeckeTraits(I,io::GAP()),
   partitionTraits(io::GAP()),
   wgraphTraits(io::GAP()),
   posetTraits(io::GAP()),
   printBettiRank(false),
   printCellNumber(false),
   printClosureSize(false),
   printCoatoms(false),
   printCompCount(false),
   printDufloNumber(false),
   printEltDescents(false),
   printElt(true),
   printEltData(true),
   printEltNumber(false),
   printLength(false),
   printType(true),
   printVersion(true),
   hasBettiPadding(false)

{
  prefix[basisH] = "coxeter_cbasis:=";
  prefix[bettiH] = "coxeter_betti:=";
  prefix[closureH] = "";
  prefix[dufloH] = "coxeter_duflo:=";
  prefix[extremalsH] = "coxeter_criticalPairs:=";
  prefix[ihBettiH] = "coxeter_ihbetti:=";
  prefix[lCOrderH] = "coxeter_lcorder:=";
  prefix[lCellsH] = "coxeter_lcells:=";
  prefix[lCellWGraphsH] = "coxeter_lcwgraphs:=";
  prefix[lWGraphH] = "coxeter_lwgraph:=";
  prefix[lrCOrderH] = "coxeter_lrcorder:=";
  prefix[lrCellsH] = "coxeter_lrcells:=";
  prefix[lrCellWGraphsH] = "coxeter_lrcwgraphs:=";
  prefix[lrWGraphH] = "coxeter_lrwgraph:=";
  prefix[rCOrderH] = "coxeter_rcorder:=";
  prefix[rCellsH] = "coxeter_rcells:=";
  prefix[rCellWGraphsH] = "coxeter_rcwgraphs:=";
  prefix[rWGraphH] = "coxeter_rwgraph:=";
  prefix[slocusH] = "coxeter_slocus:=";
  prefix[sstratificationH] = "coxeter_sstratification:=";

  postfix[basisH] = ";";
  postfix[bettiH] = ";";
  postfix[closureH] = "";
  postfix[dufloH] = ";";
  postfix[extremalsH] = ";";
  postfix[ihBettiH] = ";";
  postfix[lCOrderH] = ";";
  postfix[lCellsH] = ";";
  postfix[lCellWGraphsH] = ";";
  postfix[lWGraphH] = ";";
  postfix[lrCOrderH] = ";";
  postfix[lrCellsH] = ";";
  postfix[lrCellWGraphsH] = ";";
  postfix[lrWGraphH] = ";";
  postfix[rCOrderH] = ";";
  postfix[rCellsH] = ";";
  postfix[rCellWGraphsH] = ";";
  postfix[rWGraphH] = ";";
  postfix[slocusH] = ";";
  postfix[sstratificationH] = ";";

  header[basisH] = "GAPbasis";
  header[closureH] = "GAPclosure";
  header[dufloH] = "GAPduflo";
  header[extremalsH] = "GAPextremals";
  header[lCOrderH] = "GAPlcorder";
  header[lCellsH] = "GAPlcells";
  header[lCellWGraphsH] = "GAPlcellwgraphs";
  header[lWGraphH] = "GAPlwgraph";
  header[lrCOrderH] = "GAPlrcorder";
  header[lrCellsH] = "GAPlrcells";
  header[lrCellWGraphsH] = "GAPlrcellwgraphs";
  header[lrWGraphH] = "GAPlrwgraph";
  header[rCOrderH] = "GAPrcorder";
  header[rCellsH] = "GAPrcells";
  header[rCellWGraphsH] = "GAPrcellwgraphs";
  header[rWGraphH] = "GAPrwgraph";
  header[slocusH] = "GAPslocus";
  header[sstratificationH] = "GAPsstratification";

  hasHeader[bettiH] = false;
  hasHeader[basisH] = true;
  hasHeader[closureH] = true;
  hasHeader[dufloH] = true;
  hasHeader[extremalsH] = true;
  hasHeader[ihBettiH] = false;
  hasHeader[lCOrderH] = true;
  hasHeader[lCellsH] = true;
  hasHeader[lCellWGraphsH] = true;
  hasHeader[lWGraphH] = true;
  hasHeader[lrCOrderH] = true;
  hasHeader[lrCellsH] = true;
  hasHeader[lrCellWGraphsH] = true;
  hasHeader[lrWGraphH] = true;
  hasHeader[rCOrderH] = true;
  hasHeader[rCellsH] = true;
  hasHeader[rCellWGraphsH] = true;
  hasHeader[rWGraphH] = true;
  hasHeader[slocusH] = true;
  hasHeader[sstratificationH] = true;

  makeVersionString(version_string,"##");
  makeTypeString(type_string,"##",G);
}

OutputTraits::~OutputTraits()

{}

}

/****************************************************************************

        Chapter II -- The HeckeTraits class.

  This class holds the part of the OutputTraits which is relevant for
  output of Hecke algebra elements.

  Mostly, constructors are provided for the various output tags :

    - HeckeTraits(Pretty);
    - HeckeTraits(Terse);
    - HeckeTraits(GAP);

  Also for the derived class AddHeckeTraits, which is used for printing out
  Hecke algebra elements in vector form.

    - AddHeckeTraits(Pretty);
    - AddHeckeTraits(Terse);
    - AddHeckeTraits(GAP);

*****************************************************************************/

namespace files {

HeckeTraits::HeckeTraits(const interface::Interface& I, io::Pretty)
  :prefix(""),
   postfix(""),
   evenSeparator(""),
   oddSeparator("\n"),
   monomialPrefix(""),
   monomialPostfix(""),
   monomialSeparator(" : "),
   muMark(" *"),
   hyphens("+"),
   lineSize(io::LINESIZE),
   indent(4),
   evenWidth(io::HALFLINESIZE),
   oddWidth(0),
   padChar(' '),
   doShift(false),
   reversePrint(false),
   twoSided(true)

{}

HeckeTraits::HeckeTraits(const interface::Interface& I, io::Terse)
  :prefix(""),
   postfix(""),
   evenSeparator(""),
   oddSeparator("\n"),
   monomialPrefix(""),
   monomialPostfix(""),
   monomialSeparator(":"),
   muMark(""),
   lineSize(0),
   evenWidth(0),
   oddWidth(0),
   padChar(' '),
   doShift(false),
   reversePrint(false),
   twoSided(false)

{}

HeckeTraits::HeckeTraits(const interface::Interface& I, io::GAP)
  :prefix("[\n"),
   postfix("]"),
   evenSeparator(""),
   oddSeparator(",\n"),
   monomialPrefix("["),
   monomialPostfix("]"),
   monomialSeparator(","),
   muMark(""),
   lineSize(0),
   evenWidth(0),
   oddWidth(0),
   doShift(false),
   reversePrint(false),
   twoSided(false)

{}

HeckeTraits::~HeckeTraits()

{}

AddHeckeTraits::AddHeckeTraits(const interface::Interface& I, io::Pretty)
  :HeckeTraits(I,io::Pretty())

{
  eltTraits = new interface::GroupEltInterface(I.outInterface());
}

AddHeckeTraits::AddHeckeTraits(const interface::Interface& I, io::Terse)
  :HeckeTraits(I,io::Terse())

{
  eltTraits = new interface::GroupEltInterface(I.outInterface());
  doShift = true;
}

AddHeckeTraits::AddHeckeTraits(const interface::Interface& I, io::GAP)
  :HeckeTraits(I,io::GAP())

{
  eltTraits = new interface::GroupEltInterface(I.outInterface());

  prefix = "";
  postfix = "";
  oddSeparator = "+";
  monomialPrefix = "(";
  monomialPostfix = ")";
  monomialSeparator = ")*t(";
  reversePrint = true;
  doShift = true;
  eltTraits->prefix = "";
  eltTraits->postfix = "";
}

AddHeckeTraits::~AddHeckeTraits()

{
  delete eltTraits;
}

};

/****************************************************************************

        Chapter III -- The PartitionTraits class.

  This class holds the part of the OutputTraits which is relevant for
  output of partitions.

  Mostly, constructors are provided for the various output tags :

    - PartitionTraits(Pretty);
    - PartitionTraits(Terse);
    - PartitionTraits(GAP);

*****************************************************************************/

namespace files {

PartitionTraits::PartitionTraits(io::Pretty)
  :prefix(""),
   postfix(""),
   separator("\n"),
   classPrefix("{"),
   classPostfix("}"),
   classSeparator(","),
   classNumberPrefix(""),
   classNumberPostfix(" : "),
   printClassNumber(true)

{}

PartitionTraits::PartitionTraits(io::Terse)
  :prefix(""),
   postfix(""),
   separator("\n"),
   classPrefix(""),
   classPostfix(""),
   classSeparator(","),
   classNumberPrefix(""),
   classNumberPostfix(""),
   printClassNumber(false)

{}

PartitionTraits::PartitionTraits(io::GAP)
  :prefix("[\n"),
   postfix("]"),
   separator(",\n"),
   classPrefix("["),
   classPostfix("]"),
   classSeparator(","),
   classNumberPrefix(""),
   classNumberPostfix(""),
   printClassNumber(false)

{}

PartitionTraits::~PartitionTraits()

{}

};

/****************************************************************************

        Chapter IV -- The PolynomialTraits class.

  This class holds the part of the OutputTraits which is relevant for
  output of polynomials
  Mostly, constructors are provided for the various output tags :

    - PolynomialTraits(Pretty);
    - PolynomialTraits(Terse);
    - PolynomialTraits(GAP);

*****************************************************************************/

namespace files {

PolynomialTraits::PolynomialTraits(io::Pretty)
  :prefix(""),
   postfix(""),
   indeterminate("q"),
   sqrtIndeterminate("u"),
   posSeparator("+"),
   negSeparator(""),
   product(""),
   exponent("^"),
   expPrefix(""),
   expPostfix(""),
   zeroPol("0"),
   one(""),
   negOne("-"),
   modifierPrefix(""),
   modifierPostfix(""),
   modifierSeparator(""),
   printExponent(true),
   printModifier(false)

{}

PolynomialTraits::PolynomialTraits(io::Terse)
  :prefix("["),
   postfix("]"),
   indeterminate(""),
   sqrtIndeterminate(""),
   posSeparator(","),
   negSeparator(","),
   product(""),
   exponent(""),
   expPrefix(""),
   expPostfix(""),
   zeroPol("[]"),
   one("1"),
   negOne("-1"),
   modifierPrefix("("),
   modifierPostfix(")"),
   modifierSeparator(","),
   printExponent(false),
   printModifier(true)

{}

PolynomialTraits::PolynomialTraits(io::GAP)
  :prefix(""),
   postfix(""),
   indeterminate("q"),
   sqrtIndeterminate("u"),
   posSeparator("+"),
   negSeparator(""),
   product(""),
   exponent("^"),
   expPrefix(""),
   expPostfix(""),
   zeroPol("0"),
   one(""),
   negOne("-"),
   modifierPrefix(""),
   modifierPostfix(""),
   modifierSeparator(""),
   printExponent(true),
   printModifier(false)

{}

PolynomialTraits::~PolynomialTraits()

{}

};

/****************************************************************************

        Chapter V --- The PosetTraits class.

  This class holds the part of the OutputTraits which is relevant for the
  printout of abstract posets.

  Mostly, constructors are provided for the various tags :

    - PosetTraits(Pretty);
    - PosetTraits(Terse);
    - PosetTraits(GAP);

*****************************************************************************/

namespace files {

PosetTraits::PosetTraits(io::Pretty)
  :prefix(""),
   postfix(""),
   separator("\n"),
   edgePrefix(""),
   edgePostfix(""),
   edgeSeparator(","),
   nodePrefix(""),
   nodePostfix(" : "),
   nodeShift(0),
   printNode(true)

{}

PosetTraits::PosetTraits(io::Terse)
  :prefix(""),
   postfix(""),
   separator("\n"),
   edgePrefix(""),
   edgePostfix(""),
   edgeSeparator(","),
   nodePrefix(""),
   nodePostfix(""),
   nodeShift(0),
   printNode(false)

{}

PosetTraits::PosetTraits(io::GAP)
  :prefix("[\n"),
   postfix("]"),
   separator(",\n"),
   edgePrefix("["),
   edgePostfix("]"),
   edgeSeparator(","),
   nodePrefix(""),
   nodePostfix(""),
   nodeShift(1),
   printNode(false)

{}

PosetTraits::~PosetTraits()

{}

};

/****************************************************************************

        Chapter VI -- The WgraphTraits class.

  This class holds the part of the OutputTraits which is relevant for
  output of W-graphs.

  Mostly, constructors are provided for the various output tags :

    - WgraphTraits(Pretty);
    - WgraphTraits(Terse);
    - WgraphTraits(GAP);

*****************************************************************************/

namespace files {

WgraphTraits::WgraphTraits(io::Pretty)
  :prefix(""),
   postfix(""),
   separator("\n"),
   edgeListPrefix("{"),
   edgeListPostfix("}"),
   edgeListSeparator(","),
   edgePrefix("("),
   edgePostfix(")"),
   edgeSeparator(","),
   nodePrefix(""),
   nodePostfix(""),
   nodeSeparator(":"),
   nodeNumberPrefix(""),
   nodeNumberPostfix(":"),
   nodeShift(0),
   hasPadding(true),
   printNodeNumber(true)

{}

WgraphTraits::WgraphTraits(io::Terse)
  :prefix(""),
   postfix(""),
   separator("\n"),
   edgeListPrefix("{"),
   edgeListPostfix("}"),
   edgeListSeparator(","),
   edgePrefix("("),
   edgePostfix(")"),
   edgeSeparator(","),
   nodePrefix(""),
   nodePostfix(""),
   nodeSeparator(":"),
   nodeShift(0),
   padSize(0),
   hasPadding(false),
   printNodeNumber(false)

{}

WgraphTraits::WgraphTraits(io::GAP)
  :prefix("[\n"),
   postfix("]"),
   separator(",\n"),
   edgeListPrefix("["),
   edgeListPostfix("]"),
   edgeListSeparator(","),
   edgePrefix("["),
   edgePostfix("]"),
   edgeSeparator(","),
   nodePrefix("["),
   nodePostfix("]"),
   nodeSeparator(","),
   nodeShift(0),
   hasPadding(false),
   printNodeNumber(false)

{}

WgraphTraits::~WgraphTraits()

{}

}

/****************************************************************************

        Chapter VII --- Output functions.

  This section contains the definition of the output functions declared
  in files.h, and which are not templates :

    - appendHomology(str,h,traits) : appends a homology vector;
    - printBetti(file,y,p,traits) : prints the betti numbers of [e,y];
    - printCellOrder(file,X,p,order,traits) : prints the set of cells in X as
      an abstract poset;
    - printCoatoms(file,y,p,traits) : prints data about the element y;
    - printDescents(file,df,f,traits) : prints df as a descent set;
    - printEltData(file,y,p,traits) : prints data about the element y;
    - printHeader(file,header,traits) : prints the header;
    - printHomology(file,h,traits) : prints a homology vector;
    - printPartition(file,pi,p,order,traits);
    - printWGraph(file,X,f,traits) : prints a W-graph;

*****************************************************************************/

namespace files {


// Append the homology vector to str.
void appendHomology
  (std::string& str, const schubert::Homology& h, OutputTraits& traits)
{
  Ulong initLength = str.length();
  Ulong maxWidth = maxLength(h); // maximum width of output

  str.append(traits.bettiPrefix);

  for (Ulong j = 0; j < h.size(); ++j) {
    if (traits.printBettiRank) {
      str.append(traits.bettiRankPrefix);
      io::append(str,j);
      str.append(traits.bettiRankPostfix);
    }
    io::append(str,static_cast<Ulong>(h[j]));
    if (traits.hasBettiPadding)
      io::pad(str,(j+1)*(maxWidth+1)+initLength);
    if (j+1 < h.size()) // there is more to come
      str.append(traits.bettiSeparator);
  }

  str.append(traits.bettiPostfix);

  return;
}

void printBetti(FILE* file, const coxtypes::CoxNbr& y,
		const schubert::SchubertContext& p,
		OutputTraits& traits)

{
  schubert::Homology h = betti(y,p);

  io::print(file,traits.prefix[bettiH]);
  printHomology(file,h,traits);
  io::print(file,traits.postfix[bettiH]);
  fprintf(file,"\n");

  return;
}


/*
  Print to |file| the poset of cells defined by the graph |X|. The interface
  argument is needed to determine the shortlex ordering of the elements.
*/
void printCellOrder(FILE* file, const OrientedGraph& X,
		    const schubert::SchubertContext& p,
		    const interface::Interface& I,
		    PosetTraits& traits)
{
  OrientedGraph P(0);
  bits::Partition pi = X.cells(&P);
  posets::Poset Q(P);
  OrientedGraph H(0);
  Q.hasseDiagram(H);

  list::List<list::List<coxtypes::CoxNbr> > lc(0);
  writeClasses(lc,pi);

  schubert::NFCompare nfc(p,I.order());
  bits::Permutation a(0);
  sortLists(lc,nfc,a);
  H.permute(a.inverse());

  io::print(file,traits.prefix);

  for (Ulong j = 0; j < pi.class_count(); ++j) {
    if (traits.printNode) {
      io::print(file,traits.nodePrefix);
      fprintf(file,"%lu",j+traits.nodeShift);
      io::print(file,traits.nodePostfix);
    }
    const EdgeList& e = H.edge(j);
    io::print(file,traits.edgePrefix);
    for (Ulong i = 0; i < e.size(); ++i) {
      fprintf(file,"%lu",e[i]+traits.nodeShift);
      if (i+1 < e.size()) /* there is more to come */
	io::print(file,traits.edgeSeparator);
    }
    io::print(file,traits.edgePostfix);

    if (j+1 < pi.class_count()) // there is more to come
      io::print(file,traits.separator);
  }

  io::print(file,traits.postfix);
}

void printCoatoms(FILE* file, const coxtypes::CoxNbr& y, const schubert::SchubertContext& p,
		  const interface::Interface& I, OutputTraits& traits)

/*
  Prints out the coatoms of y.
*/

{
  const schubert::CoxNbrList& c = p.hasse(y);

  io::print(file,traits.coatomPrefix);

  for (Ulong j = 0; j < c.size(); ++j) {
    print(p,file,c[j],I);
    if (j+1 < c.size()) // there is more to come
      io::print(file,traits.coatomSeparator);
  }

  io::print(file,traits.coatomPostfix);

  return;
}

void printDescents(FILE* file, const Lflags& df, const Lflags& f,
		   const interface::Interface& I, WgraphTraits& traits)

{
  if (!(f&1)) { // left descents
    print(file,df,I);
  }
  else if ((f >> I.rank()) == 0) { // right descents
    print(file,df,I);
  }
  else { // two-sided descents
    printTwosided(file,df,I);
  }

  return;
}

void printEltData(FILE* file, const coxtypes::CoxNbr& y, const schubert::SchubertContext& p,
		  const interface::Interface& I, OutputTraits& traits)

/*
  This function prints out some data about the element y :
    - normal form;
    - left and right descent sets;
    - length;
*/

{
  io::print(file,traits.eltDataPrefix);

  if (traits.printElt) {
    io::print(file,traits.eltPrefix);
    print(p,file,y,I);
    io::print(file,traits.eltPostfix);
  }

  if (traits.printEltDescents) {
    io::print(file,traits.lDescentPrefix);
    print(file,p.ldescent(y),I);
    io::print(file,traits.lDescentPostfix);
    io::print(file,traits.rDescentPrefix);
    print(file,p.rdescent(y),I);
    io::print(file,traits.rDescentPostfix);
  }

  if (traits.printLength) {
    io::print(file,traits.lengthPrefix);
    fprintf(file,"%lu",static_cast<Ulong>(p.length(y)));
    io::print(file,traits.lengthPostfix);
  }

  io::print(file,traits.eltDataPostfix);

  return;
}

void printHeader(FILE* file, const Header& header, OutputTraits& traits)

/*
  Prints the header part of the file. It has already been checked if the
  header has to be printed.
*/

{
  using namespace directories;

  // preliminaries

  if (traits.printVersion)
    io::print(file,traits.version_string);
  if (traits.printType)
    io::print(file,traits.type_string);
  if (traits.hasHeader[header])
    io::printFile(file,traits.header[header].c_str(),HEADER_DIR);

  return;
}

void printHomology(FILE* file, const schubert::Homology& h, OutputTraits& traits)

{
  std::string buf;

  appendHomology(buf,h,traits);

  if (traits.lineSize)
    io::foldLine(file,buf,traits.lineSize,0,traits.bettiHyphens.c_str());
  else
    io::print(file,buf);

  if (traits.printClosureSize) {
    fprintf(file,"\n\n");
    Ulong size = 0;
    for (Ulong j = 0; j < h.size(); ++j)
      size += h[j];
    io::print(file,traits.closureSizePrefix);
    fprintf(file,"%lu",size);
    io::print(file,traits.closureSizePostfix);
  }

  return;
}


// Print the classes of the partition, according to the given traits.
void printPartition(FILE* file, const bits::Partition& pi,
		    const schubert::SchubertContext& p,
		    const interface::Interface& I, PartitionTraits& traits)
{
  list::List<list::List<coxtypes::CoxNbr> > lc(0);
  writeClasses(lc,pi);

  schubert::NFCompare nfc(p,I.order());
  bits::Permutation a(0);
  sortLists(lc,nfc,a);

  int d = io::digits(lc.size()-1,10);
  io::print(file,traits.prefix);

  for (Ulong j = 0; j < lc.size(); ++j) {
    list::List<coxtypes::CoxNbr> l = lc[a[j]];
    if (traits.printClassNumber) {
      io::print(file,traits.classNumberPrefix);
      fprintf(file,"%*lu",d,j);
      io::print(file,traits.classNumberPostfix);
    }
    io::print(file,traits.classPrefix);
    for (Ulong i = 0; i < l.size(); ++i) {
      print(p,file,l[i],I);
      if (i+1 < l.size()) /* there is more to come */
	io::print(file,traits.classSeparator);
    }
    io::print(file,traits.classPostfix);
    if (j+1 < lc.size()) // there is more to come
      io::print(file,traits.separator);
  }

  io::print(file,traits.postfix);

  return;
}


/*
  Print out a W-graph. We represent a W-graph as a list of pairs [descent
  set,edge-list]. The edge-list itself is a set of pairs [dest,mu], where dest
  is the end-point of the edge, and mu is the mu-coefficient.
*/

void printWGraph(FILE* file, const WGraph& X, const Lflags& f,
		 const interface::Interface& I, WgraphTraits& traits)
{
  int d = io::digits(X.size()-1,10);

  io::print(file,traits.prefix);

  for (Ulong j = 0; j < X.size(); ++j) {
    if (traits.printNodeNumber) {
      io::print(file,traits.nodeNumberPrefix);
      fprintf(file,"%*lu",d,j);
      io::print(file,traits.nodeNumberPostfix);
    }
    io::print(file,traits.nodePrefix);
    printDescents(file,X.descent(j),f,I,traits);
    io::print(file,traits.nodeSeparator);
    const CoeffList& mu = X.coeffList(j);
    const EdgeList& e = X.edge(j);
    io::print(file,traits.edgeListPrefix);
    for (Ulong i = 0; i < e.size(); ++i) {
      io::print(file,traits.edgePrefix);
      fprintf(file,"%lu",static_cast<Ulong>(e[i]));
      io::print(file,traits.edgeSeparator);
      fprintf(file,"%ld",static_cast<long>(mu[i]));
      io::print(file,traits.edgePostfix);
      if (i+1 < e.size()) // there is more to come
	io::print(file,traits.edgeListSeparator);
    }
    io::print(file,traits.edgeListPostfix);
    io::print(file,traits.nodePostfix);
     if (j+1 < X.size()) { // there is more to come
      io::print(file,traits.separator);
      if (traits.hasPadding)
	fprintf(file,"%*s",traits.padSize,"");
     }
  }

  io::print(file,traits.postfix);

  return;
}

};

/****************************************************************************

       Chapter VIII --- Helper functions

  This section defines some helper functions declared in files.h:

    - appendModifier(str,d,m,traits) : does the printing of the modifier;
    - appendSeparator(str,n,traits) : does the printing of the separator;
    - minimal_class_reps(pi,c) : minimal representatives of the classes in |pi|;
    - printModifier(str,d,m,traits) : does the printing of the modifier;
    - printSeparator(file,n,traits) : does the printing of the separator;
    - setTwoSided(traits,f) : sets *traits.twoSided to f;
    - sortLists(lc,nfc,a) : sorts the various lists in lc;
    - writeClasses(lc,pi) : writes out the classes of pi in lc;

  and some others private to the present module:

    - makeVersion_string(vstr,str) : makes the version string with str as
      prefix;
    - makeTypeString(tstr,str) : makes the type string with str as prefix;

*****************************************************************************/

namespace files {

void appendModifier(std::string& str, const Ulong& d, const long& m,
		    PolynomialTraits& traits)

/*
  Appends the polynomial modifier to the string.
*/

{
  str.append(traits.modifierPrefix);
  io::append(str,static_cast<Ulong>(d));
  str.append(traits.modifierSeparator);
  io::append(str,static_cast<long>(m));
  str.append(traits.modifierPostfix);

  return;
}

void appendSeparator(std::string& str, const Ulong& n, HeckeTraits& traits)

{
  if (traits.twoSided) { // need to look at parity of n
    if (n%2)
      str.append(traits.oddSeparator);
    else
      str.append(traits.evenSeparator);
  }
  else
    str.append(traits.oddSeparator);

  return;
}


/*
  Extract the minimal representative from each class in pi, for the comparison
  functor c, and writes it in min, in the order of the partition numbers.
*/
schubert::CoxNbrList minimal_class_reps
  (const bits::Partition& pi, schubert::NFCompare& c)
{
  schubert::CoxNbrList result;
  result.reserve(pi.class_count());
  for (bits::Partition::iterator pit = pi.begin(); pit; ++pit)
  {
    auto pi_class = *pit;
    coxtypes::CoxNbr x = *std::min_element(pi_class.begin(),pi_class.end(),c);
    result.push_back(x);
  }
  return result;
}

void pad(std::string& str, const Ulong& n, HeckeTraits& traits)

/*
  Adds padding to the string, to make its length at least equal to
  traits.even(odd)Width, according to the parity of n. If the
  width is zero, no padding is added.
*/

{
  if (!traits.twoSided)
    return;

  if (n%2 && traits.oddWidth) {
    Ulong m = str.length();
    if (m < traits.oddWidth) { // odd padding
      for (Ulong j = m; j < traits.oddWidth; ++j)
	str.push_back(traits.padChar);
    }
  }

  if (!(n%2) && traits.evenWidth) {
    Ulong m = str.length();
    if (m < traits.evenWidth) { // odd padding
      for (Ulong j = m; j < traits.evenWidth; ++j)
	str.push_back(traits.padChar);
    }
  }

  return;
}

void printModifier(FILE* file, const Ulong& d, const long& m,
		   PolynomialTraits& traits)

/*
  Prints the modifier parameters (d,m).
*/

{
  io::print(file,traits.modifierPrefix);
  fprintf(file,"%lu",d);
  io::print(file,traits.modifierSeparator);
  fprintf(file,"%ld",m);
  io::print(file,traits.modifierPostfix);

  return;
}

void printSeparator(FILE* file, const Ulong& n, HeckeTraits& traits)

{
  if (traits.twoSided) { // need to look at parity of n
    if (n%2)
      io::print(file,traits.oddSeparator);
    else
      io::print(file,traits.evenSeparator);
  }
  else
    io::print(file,traits.oddSeparator);

  return;
}

void sortLists(list::List<list::List<coxtypes::CoxNbr> >& lc, schubert::NFCompare& nfc,
	       bits::Permutation& a)

/*
  This function sorts the various rows in lc, and puts in a the permutation
  which sorts ld according to the first term in each row.
*/

{
  list::List<coxtypes::CoxNbr> lfirst(0);
  lfirst.setSize(lc.size());

  for (Ulong j = 0; j < lc.size(); ++j) {
    lc[j].sort(nfc);
    lfirst[j] = lc[j][0];
  }

  sortI(lfirst,nfc,a);

  return;
}


/*
  Write out in lc the various classes of the partition pi.
*/
void writeClasses
  (list::List<list::List<coxtypes::CoxNbr> >& lc, const bits::Partition& pi)
{
  lc.setSize(pi.class_count());
  Ulong j = 0;

  for (bits::Partition::iterator pit=pi.begin(); pit; ++pit) {
    new(&lc[j]) list::List<coxtypes::CoxNbr>((*pit).begin(),(*pit).end());
    ++j;
  }

  return;
}

};

namespace {

void makeVersionString(std::string& vstr, const std::string& str)

{
  using namespace version;

  vstr.append(str);
  vstr.push_back('\n');
  vstr.append(str);
  vstr.append(" This file has been created by ");
  vstr.append(NAME);
  vstr.append(" version ");
  vstr.append(VERSION);
  vstr.push_back('\n');
}

void makeTypeString
  (std::string& tstr, const std::string& str, const graph::CoxGraph& G)

{
  tstr.append(str);
  tstr.push_back('\n');
  tstr.append(str);
  tstr.append(" Group type is ");
  tstr.append(G.type().name());
  io::append(tstr,G.rank());
  tstr.push_back('\n');
}

};

namespace {


// Return the maximal length of an entry in h.
Ulong maxLength(const schubert::Homology& h)
{
  std::string buf;

  Ulong maxl = 0;

  for (Ulong j = 0; j < h.size(); ++j)
  {
    buf.clear();
    buf.append("h[");
    io::append(buf,j);
    buf.append("] = ");
    io::append(buf,h[j]);
    if (maxl < buf.size())
      maxl = buf.size();
  }

  return maxl;
}

};
