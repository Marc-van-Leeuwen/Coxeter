/*
  This is cells.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef CELLS_H  /* guard against multiple inclusions */
#define CELLS_H

#include "globals.h"

#include "bits.h"
#include "kl.h"
#include "uneqkl.h"
#include "wgraph.h"

/******** function declarations **********************************************/

namespace cells {
  coxtypes::CoxNbr checkClasses(const bits::Partition& pi,
				const schubert::SchubertContext& p);
  void lCells(bits::Partition& pi, kl::KLContext& kl);
  void rCells(bits::Partition& pi, kl::KLContext& kl);
  void lrCells(bits::Partition& pi, kl::KLContext& kl);
  void lDescentPartition(bits::Partition& pi,
			 const schubert::SchubertContext& p);
  void rDescentPartition(bits::Partition& pi,
			 const schubert::SchubertContext& p);
  void lStringEquiv(bits::Partition& pi, const schubert::SchubertContext& p);
  void lStringEquiv(bits::Partition& pi,
		    const bits::SubSet& q, const schubert::SchubertContext& p);
  void rStringEquiv(bits::Partition& pi, const schubert::SchubertContext& p);
  void rStringEquiv(bits::Partition& pi,
		    const bits::SubSet& q, const schubert::SchubertContext& p);
  void lrStringEquiv(bits::Partition& pi, const schubert::SchubertContext& p);
  void lrStringEquiv(bits::Partition& pi,
		     const bits::SubSet& q, const schubert::SchubertContext& p);
  void lGeneralizedTau(bits::Partition& pi, schubert::SchubertContext& p);
  void rGeneralizedTau(bits::Partition& pi, schubert::SchubertContext& p);

  template<char> // one of 'l', 'r', 'b'
    wgraph::OrientedGraph graph(kl::KLContext& kl);
  template<char> // one of 'l', 'r', 'b'
    wgraph::OrientedGraph graph(uneqkl::KLContext& kl);
  template<char> // one of 'l', 'r', 'b'
    wgraph::WGraph W_graph(kl::KLContext& kl);
  template<char> // one of 'l', 'r', 'b'
    wgraph::WGraph W_graph(const bits::SubSet& q, kl::KLContext& kl);

  void printCellPartition(FILE* file, const kl::KLContext& kl);

}; // |namespace cells|

#endif
