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
  template<char> // one of 'l', 'r'
    bits::Partition descent_partition(const schubert::SchubertContext& p);
  template<char> // one of 'l', 'r'
    bits::Partition generalized_tau(schubert::SchubertContext& p);

  template<char> // one of 'l', 'r', 'b'
    bits::Partition cells(kl::KLContext& kl);

  template<char> // one of 'l', 'r'
    bits::Partition string_equiv(const schubert::SchubertContext& p);
  template<char> // one of 'l', 'r'
    bits::Partition string_equiv
      (const bits::SubSet& q, const schubert::SchubertContext& p);

  template<char> // one of 'l', 'r', 'b'
    wgraph::OrientedGraph graph(kl::KLContext& kl);
  template<char> // one of 'l', 'r', 'b'
    wgraph::OrientedGraph graph(uneqkl::KLContext& kl);
  template<char> // one of 'l', 'r', 'b'
    wgraph::WGraph W_graph(kl::KLContext& kl);
  template<char> // one of 'l', 'r', 'b'
    wgraph::WGraph W_graph(const bits::SubSet& q, kl::KLContext& kl);

  coxtypes::CoxNbr checkClasses(const bits::Partition& pi,
				const schubert::SchubertContext& p);
  void printCellPartition(FILE* file, const kl::KLContext& kl);

}; // |namespace cells|

#endif
