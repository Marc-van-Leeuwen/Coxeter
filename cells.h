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
  coxtypes::CoxNbr checkClasses(const bits::Partition& pi, const schubert::SchubertContext& p);
  void lCells(bits::Partition& pi, kl::KLContext& kl);
  void rCells(bits::Partition& pi, kl::KLContext& kl);
  void lrCells(bits::Partition& pi, kl::KLContext& kl);
  void lDescentPartition(bits::Partition& pi, const schubert::SchubertContext& p);
  void rDescentPartition(bits::Partition& pi, const schubert::SchubertContext& p);
  void lStringEquiv(bits::Partition& pi, const schubert::SchubertContext& p);
  void lStringEquiv(bits::Partition& pi, const bits::SubSet& q, const schubert::SchubertContext& p);
  void rStringEquiv(bits::Partition& pi, const schubert::SchubertContext& p);
  void rStringEquiv(bits::Partition& pi, const bits::SubSet& q, const schubert::SchubertContext& p);
  void lrStringEquiv(bits::Partition& pi, const schubert::SchubertContext& p);
  void lrStringEquiv(bits::Partition& pi, const bits::SubSet& q, const schubert::SchubertContext& p);
  void lGeneralizedTau(bits::Partition& pi, const schubert::SchubertContext& p);
  void rGeneralizedTau(bits::Partition& pi, const schubert::SchubertContext& p);
  void lGraph(wgraph::OrientedGraph& X, kl::KLContext& kl);
  void lrGraph(wgraph::OrientedGraph& X, kl::KLContext& kl);
  void rGraph(wgraph::OrientedGraph& X, kl::KLContext& kl);
  void lGraph(wgraph::OrientedGraph& X, uneqkl::KLContext& kl);
  void lrGraph(wgraph::OrientedGraph& X, uneqkl::KLContext& kl);
  void rGraph(wgraph::OrientedGraph& X, uneqkl::KLContext& kl);
  void lWGraph(wgraph::WGraph& X, kl::KLContext& kl);
  void lWGraph(wgraph::WGraph& X, const bits::SubSet& q, kl::KLContext& kl);
  void rWGraph(wgraph::WGraph& X, kl::KLContext& kl);
  void rWGraph(wgraph::WGraph& X, const bits::SubSet& q, kl::KLContext& kl);
  void lrWGraph(wgraph::WGraph& X, kl::KLContext& kl);
  void lrWGraph(wgraph::WGraph& X, const bits::SubSet& q, kl::KLContext& kl);
  void printCellPartition(FILE* file, const kl::KLContext& kl);
};

#endif
