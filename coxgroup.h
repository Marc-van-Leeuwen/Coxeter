/*
  This is coxgroup.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

/**************************************************************************

 This file presents the main type in this program, namely the CoxGroup
 class. This is the base class for the hierarchy of Coxeter group
 classes, ranging from the most general ones considered in this program
 (rank <= 255, coefficients of Coxeter matrix <= COXENTRY_MAX) to
 the most special, the class SmallCoxGroup.

 We have adhered to the philosophy that all non-leaf classes in the
 hierarchy should be abstract, marked '(a)' below; the leaf classes are
 concrete, marked '(c)' below. So this file contains only the root
 of the hierarchy, as an abstract class.

 The layout of the Coxeter hierarchy as considered in this program
 is as follows :

   CoxGroup(a)
     FiniteCoxGroup(a)
       FiniteBigRankCoxGroup(a)
         GeneralFBRCoxGroup(c)
       FiniteMedRankCoxGroup(a)
         GeneralFMRCoxGroup(c)
         FiniteSmallRankCoxGroup(a)
	   GeneralFSRCoxGroup(c)
	   SmallCoxGroup(a)
	     GeneralSCoxGroup(c)
       TypeACoxGroup(a)
         TypeABigRankCoxGroup(c)
         TypeAMedRankCoxGroup(c)
           TypeASmallRankCoxGroup(c)
             TypeASmallCoxGroup(c)
     AffineCoxGroup(a)
       AffineBigRankCoxGroup(a)
         GeneralABRCoxGroup(c)
       AffineMedRankCoxGroup(a)
         GeneralAMRCoxGroup(c)
         AffineSmallRankCoxGroup(a)
           GeneralASRCoxGroup(c)
     GeneralCoxGroup(a)
       BigRankCoxGroup(a)
         GeneralBRCoxGroup(c)
       MedRankCoxGroup(a)
         GeneralMRCoxGroup(c)
         SmallRankCoxGroup(a)
           GeneralSRCoxGroup(c)

 **************************************************************************/

#ifndef COXGROUP_H  /* guarantee single inclusion */
#define COXGROUP_H

#include "globals.h"
#include "coxtypes.h"
#include "files.h"
#include "graph.h"
#include "hecke.h"
#include "interface.h"
#include "invkl.h"
#include "kl.h"
#include "klsupport.h"
#include "minroots.h"
#include "transducer.h"
#include "uneqkl.h"

/******** type definitions **************************************************/

class coxgroup::CoxGroup { // has been declared in coxtypes.h

 protected:

  graph::CoxGraph d_graph;
  minroots::MinTable d_mintable;
  klsupport::KLSupport d_klsupport;
  std::unique_ptr<kl::KLContext> d_kl;
  std::unique_ptr<invkl::KLContext> d_invkl;
  std::unique_ptr<uneqkl::KLContext> d_uneqkl;
  std::unique_ptr<interface::Interface> d_interface; // ptr to maybe derived
  files::OutputTraits d_outputTraits;

  struct CoxHelper; // predeclare, defined in implementation part
  friend CoxHelper;  /* provides helper functions */
  std::unique_ptr<CoxHelper> d_help; // pointer level to hide implementation

 public:

  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(CoxGroup));}

/******** Chapter 0 : CoxGroup objects *************************************/

  CoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~CoxGroup();

  minroots::MinTable& mintable() { return d_mintable; }
  klsupport::KLSupport& klsupport() { return d_klsupport; }
  kl::KLContext& kl() { activateKL(); return *d_kl; }
  invkl::KLContext& invkl() { activateIKL(); return *d_invkl; }
  uneqkl::KLContext& uneqkl() { activateUEKL(); return *d_uneqkl; }
  virtual interface::Interface& interface() { return *d_interface; }
  virtual files::OutputTraits& outputTraits() { return d_outputTraits; }

  const graph::CoxGraph& graph() const { return d_graph; }
  virtual const interface::Interface& interface() const { return *d_interface; }
  const minroots::MinTable& mintable() const { return d_mintable; }
  const klsupport::KLSupport& klsupport() const { return d_klsupport; }
  const kl::KLContext& kl() const { return *d_kl; } // assumes |d_kl!=nullptr|
  const uneqkl::KLContext& uneqkl() const { return *d_uneqkl; }
  const schubert::SchubertContext& schubert() const
    { return d_klsupport.schubert(); }

  virtual const files::OutputTraits& outputTraits() const
    { return d_outputTraits; }

  void activateKL();
  void activateIKL();
  void activateUEKL();

  /* graph data */

  graph::CoxEntry M(coxtypes::Generator s, coxtypes::Generator t) const;/* inlined */
  coxtypes::Rank rank() const;                                   /* inlined */
  const type::Type& type() const;                                /* inlined */
  virtual coxtypes::CoxSize order() const = 0;
  virtual bool isFullContext() const;                            /* inlined */

/******** Chapter I : Elementary operations ********************************/

  /* word operations */

  virtual int insert(coxtypes::CoxWord& g, const coxtypes::Generator& s) const;      /* inlined */
  virtual const coxtypes::CoxWord& inverse(coxtypes::CoxWord& g) const;              /* inlined */
  virtual const coxtypes::CoxWord& normalForm(coxtypes::CoxWord& g) const;           /* inlined */
  virtual const coxtypes::CoxWord& power(coxtypes::CoxWord& g, const Ulong& m) const;
                                                                 /* inlined */
  virtual int prod(coxtypes::CoxWord& g, const coxtypes::Generator& s) const;        /* inlined */
  virtual int prod(coxtypes::CoxWord& g, const coxtypes::CoxWord& h) const;          /* inlined */
  virtual const coxtypes::CoxWord& reduced(coxtypes::CoxWord& g, coxtypes::CoxWord& h) const;  /* inlined */

  /* descent sets */

  virtual Lflags descent(const coxtypes::CoxWord& g) const;                /* inlined */
  virtual GenSet ldescent(const coxtypes::CoxWord& g) const;               /* inlined */
  virtual GenSet rdescent(const coxtypes::CoxWord& g) const;               /* inlined */
  bool isDescent(const coxtypes::CoxWord& g, const coxtypes::Generator& s) const;

/******** Chapter II : Schubert context **************************************/

  virtual coxtypes::CoxNbr contextNumber(const coxtypes::CoxWord& g) const;           /* inlined */
  coxtypes::CoxNbr contextSize() const;                                     /* inlined */
  coxtypes::Length length(const coxtypes::CoxNbr& x) const;                           /* inlined */

  virtual coxtypes::CoxNbr extendContext(const coxtypes::CoxWord& g);
  virtual void permute(const bits::Permutation& a);

  virtual Lflags descent(const coxtypes::CoxNbr& x) const;                  /* inlined */
  virtual GenSet ldescent(const coxtypes::CoxNbr& x) const;                 /* inlined */
  virtual GenSet rdescent(const coxtypes::CoxNbr& x) const;                 /* inlined */

  virtual coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& x) const;                  /* inlined */
  virtual int prod(coxtypes::CoxNbr& x, const coxtypes::Generator& s) const;
  virtual int lprod(coxtypes::CoxNbr& x, const coxtypes::Generator& s) const;         /* inlined */
  virtual int prod(coxtypes::CoxWord& g, const coxtypes::CoxNbr& x) const;
  virtual int prod(coxtypes::CoxNbr& x, const coxtypes::CoxWord& g) const;

  virtual const containers::vector<coxtypes::CoxNbr>& extrList
    (const coxtypes::CoxNbr& x) const  { return d_klsupport.extrList(x); }

/******** Chapter III : Bruhat ordering **************************************/

  virtual void coatoms(list::List<coxtypes::CoxWord>& c, const coxtypes::CoxWord& g) const;
  virtual const schubert::CoxNbrList& coatoms(const coxtypes::CoxNbr& x) const
    { return schubert().hasse(x); }
  virtual void extractClosure(bits::BitMap& b, const coxtypes::CoxNbr& x) const;  /* inlined */
  virtual bool inOrder(const coxtypes::CoxWord& h, const coxtypes::CoxWord& g) const; /* inlined */
  virtual bool inOrder(list::List<coxtypes::Length>& a, const coxtypes::CoxWord& h, const coxtypes::CoxWord& g)
    const;                                                        /* inlined */
  virtual bool inOrder(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y) const;   /* inlined */
  virtual bool isDihedral(const coxtypes::CoxWord& g) const;

/******** Chapter IV : Kazhdan-Lusztig functions *****************************/

  virtual void cBasis(kl::HeckeElt& h, const coxtypes::CoxNbr& y);
  virtual void fillKL();
  virtual void fillMu();
  virtual const kl::KLPol& klPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
  virtual void klRow(kl::HeckeElt& h, const coxtypes::CoxNbr& y);
  virtual klsupport::KLCoeff mu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);

  virtual void fillIKL();
  virtual void fillIMu();
  virtual const invkl::KLPol& invklPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
  virtual void invklRow(invkl::HeckeElt& h, const coxtypes::CoxNbr& y);

  virtual void fillUEKL();
  virtual void fillUEMu();
  virtual const uneqkl::KLPol& uneqklPol
   (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
  const uneqkl::MuPol uneqmu
   (const coxtypes::Generator& s,
    const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
  virtual void uneqcBasis(uneqkl::HeckeElt& h, const coxtypes::CoxNbr& y);
  virtual void uneqklRow(uneqkl::HeckeElt& h, const coxtypes::CoxNbr& y);

/******** Chapter V : I/O ***************************************************/

  /* elementary i/o functions */

  const bits::Permutation& ordering() const;                           /* inlined */

  std::string& append(std::string& str, const coxtypes::Generator& s) const;         /* inlined */
  std::string& append(std::string& str, const coxtypes::CoxWord& g) const;           /* inlined */
  std::string& append(std::string& str, const GenSet& f) const;            /* inlined */

  void printSymbol(FILE* file, const coxtypes::Generator& s) const;        /* inlined */
  void print(FILE* file, const coxtypes::CoxWord& g) const;                /* inlined */
  void print(FILE* file, const coxtypes::CoxNbr& x) const;                 /* inlined */
  void printFlags(FILE* file, const GenSet& f) const;            /* inlined */

  void parse(interface::ParseInterface& P) const;
  virtual bool parseGroupElement(interface::ParseInterface& P) const;
  bool parseBeginGroup(interface::ParseInterface& P) const;
  bool parseContextNumber(interface::ParseInterface& P) const;
  bool parseEndGroup(interface::ParseInterface& P) const;
  virtual bool parseModifier(interface::ParseInterface& P) const;
  virtual void modify
    (interface::ParseInterface& P, const interface::Token& tok) const;

  /* modifying the interface */

  template<class C> void setOutputTraits(C);

  void setInPostfix(const std::string& a);                       /* inlined */
  void setInPrefix(const std::string& a);                        /* inlined */
  void setInSeparator(const std::string& a);                     /* inlined */
  void setInSymbol(const coxtypes::Generator& s, const std::string& a);         /* inlined */
  void setOrdering(const bits::Permutation& order);                    /* inlined */
  void setOutPostfix(const std::string& a);                      /* inlined */
  void setOutPrefix(const std::string& a);                       /* inlined */
  void setOutSeparator(const std::string& a);                    /* inlined */
  void setOutSymbol(const coxtypes::Generator& s, const std::string& a);        /* inlined */

  template <class H> void printHeckeElt(FILE* file, const H& h); /* inlined */
};

/******** Inline implementations ******************************************/

namespace coxgroup {

/* Chapter 0 */


inline graph::CoxEntry CoxGroup::M
  (coxtypes::Generator s, coxtypes::Generator t) const
 {return(graph().M(s,t));}
inline coxtypes::Rank CoxGroup::rank() const {return(graph().rank());}
inline const type::Type& CoxGroup::type() const {return graph().type();}

inline bool CoxGroup::isFullContext() const {return false;}

/* Chapter I */

inline int CoxGroup::insert(coxtypes::CoxWord& g, const coxtypes::Generator& s) const
 {return mintable().insert(g,s,ordering());}
inline const coxtypes::CoxWord& CoxGroup::inverse(coxtypes::CoxWord& g) const
 {return mintable().inverse(g);}
inline const coxtypes::CoxWord& CoxGroup::normalForm(coxtypes::CoxWord& g) const
 {return mintable().normalForm(g,interface().order());}

inline const coxtypes::CoxWord& CoxGroup::power(coxtypes::CoxWord& g, const Ulong& m) const
 {return mintable().power(g,m);}
inline int CoxGroup::prod(coxtypes::CoxWord& g, const coxtypes::Generator& s) const
 {return mintable().prod(g,s);}
inline int CoxGroup::prod(coxtypes::CoxWord& g, const coxtypes::CoxWord& h) const
 {return mintable().prod(g,h);}
inline const coxtypes::CoxWord& CoxGroup::reduced(coxtypes::CoxWord& g, coxtypes::CoxWord& h) const
 {return mintable().reduced(g,h);}

inline Lflags CoxGroup::descent(const coxtypes::CoxWord& g) const
 {return mintable().descent(g);}
inline GenSet CoxGroup::ldescent(const coxtypes::CoxWord& g) const
 {return mintable().ldescent(g);}
inline GenSet CoxGroup::rdescent(const coxtypes::CoxWord& g) const
 {return mintable().rdescent(g);}

/* Chapter II */

inline coxtypes::CoxNbr CoxGroup::contextNumber(const coxtypes::CoxWord& g) const
 {return schubert().contextNumber(g);}
inline coxtypes::CoxNbr CoxGroup::contextSize() const
 {return d_klsupport.size();}
inline coxtypes::Length CoxGroup::length(const coxtypes::CoxNbr& x) const
 {return d_klsupport.length(x);}

inline Lflags CoxGroup::descent(const coxtypes::CoxNbr& x) const
  {return schubert().descent(x);}
inline GenSet CoxGroup::ldescent(const coxtypes::CoxNbr& x) const
 {return schubert().ldescent(x);}
inline GenSet CoxGroup::rdescent(const coxtypes::CoxNbr& x) const
 {return schubert().rdescent(x);}

inline coxtypes::CoxNbr CoxGroup::inverse(const coxtypes::CoxNbr& x) const
 {return d_klsupport.inverse(x);}
inline int CoxGroup::lprod(coxtypes::CoxNbr& x, const coxtypes::Generator& s) const
 {return prod(x,s+rank());}


/* Chapter III */

inline void CoxGroup::extractClosure(bits::BitMap& b, const coxtypes::CoxNbr& x) const
 {return schubert().extractClosure(b,x);}
inline bool CoxGroup::inOrder(const coxtypes::CoxWord& g, const coxtypes::CoxWord& h) const
 {return mintable().inOrder(g,h);}
inline bool CoxGroup::inOrder(list::List<coxtypes::Length>& a,
			      const coxtypes::CoxWord& g,
			      const coxtypes::CoxWord& h) const
 {return mintable().inOrder(a,g,h);}
inline bool CoxGroup::inOrder(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y) const
 {return schubert().inOrder(x,y);}

/* Chapter V */

inline const bits::Permutation& CoxGroup::ordering() const
  {return interface().order();}

inline std::string& CoxGroup::append(std::string& str, const coxtypes::Generator& s)
  const {return appendSymbol(str,s,interface());}
inline std::string& CoxGroup::append(std::string& str, const coxtypes::CoxWord& g) const
 {return interface::append(str,g,interface());}
inline std::string& CoxGroup::append(std::string& str, const GenSet& f) const
 {return interface::append(str,f,interface());}

inline void CoxGroup::printSymbol(FILE* file, const coxtypes::Generator& s)
  const {return interface::printSymbol(file,s,interface());}
inline void CoxGroup::print(FILE* file, const coxtypes::CoxWord& g) const
 {return interface().print(file,g);}
inline void CoxGroup::print(FILE* file, const coxtypes::CoxNbr& x) const
 {return schubert().print(file,x,interface());}
inline void CoxGroup::printFlags(FILE* file, const GenSet& f) const
 {return interface::print(file,f,interface());}

inline void CoxGroup::setInPostfix(const std::string& a)
  {interface().setInPostfix(a);}
inline void CoxGroup::setInPrefix(const std::string& a)
  {interface().setInPrefix(a);}
inline void CoxGroup::setInSeparator(const std::string& a)
  {interface().setInSeparator(a);}
inline void CoxGroup::setInSymbol(const coxtypes::Generator& s, const std::string& a)
  {interface().setInSymbol(s,a);}
inline void CoxGroup::setOrdering(const bits::Permutation& order)
  {interface().setOrder(order);}
inline void CoxGroup::setOutPostfix(const std::string& a)
  {interface().setOutPostfix(a);}
inline void CoxGroup::setOutPrefix(const std::string& a)
  {interface().setOutPrefix(a);}
inline void CoxGroup::setOutSeparator(const std::string& a)
  {interface().setOutSeparator(a);}
inline void CoxGroup::setOutSymbol
  (const coxtypes::Generator& s, const std::string& a)
  {interface().setOutSymbol(s,a);}

template <class H>
inline void CoxGroup::printHeckeElt(FILE* file, const H& h)
  {files::printHeckeElt(file,h,schubert(),outputTraits());}



/******** template definitions ***********************************************/

template<class C> void CoxGroup::setOutputTraits(C)
 { d_outputTraits = files::OutputTraits(graph(),interface(),C()); }

 }; // |namespace coxgroup|

#endif
