/*
  This is uneqkl.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  Coxeter version 3.1 Copyright (C) 2024 Marc van Leeuwen
  See file LICENSE for full copyright notice
*/

#include "uneqkl.h"
#include "sl_list.h"
#include "interactive.h"

namespace uneqkl {

namespace {
  void muSubtraction(KLPol& p, const MuPol& mp, const KLPol& q,
		     const Ulong& d, const long& m);
  KLPol positivePart(const KLPol& q, const Ulong& d, const long& m);
}; // |namespace|

/*
			      The |KLHelper| class

  This is a subclass that holds a non-|const| reference to the main |KLContext|
  that created it and vice versa; its methods can therefore freely operate on
  the private members of |KLContext|. This relation is dangerous, as the |const|
  correctness is completely left to the discipline of the implementations below:
  a |const| method of |KLContext| can freely call a non-|const| method of the
  helper class (through is non-|const| pointer to it), and vice versa. This is
  avoided however, in large part by providing methods in the helper class that
  just call |KLContext| methods with the same name and |const|-ness.

 */
struct KLContext::KLHelper
{
// data
  klsupport::KLSupport& d_klsupport; // unowned, |CoxGroup| owns it

  containers::vector<std::unique_ptr<KLRow> > d_klList;
  containers::vector<MuTable> d_muTable; // indexed by |s|
  containers::vector<coxtypes::Length> d_L; // lengths of generators
  containers::vector<coxtypes::Length> d_length; // lengths of context elements
  containers::bag<KLPol> KL_pool; // this owns the |KLPol| values
  containers::bag<MuPol> mu_pool; // this owns the |MuPol| values
  KLStats d_stats;

// constructors and destructors
  void* operator new(size_t size) { return memory::arena().alloc(size); }
  void operator delete(void* ptr)
    { return memory::arena().free(ptr,sizeof(KLHelper)); }
  KLHelper
    (klsupport::KLSupport& kls,
     const graph::CoxGraph& G, const interface::Interface& I);

// methods
// relay methods
  coxtypes::Rank rank() { return d_klsupport.rank(); }
  const schubert::SchubertContext& schubert() const
    { return d_klsupport.schubert();}
  coxtypes::Generator last(coxtypes::CoxNbr x) const
    { return d_klsupport.last(x); }
  const klsupport::ExtrRow& extrList(coxtypes::CoxNbr y) const
    { return d_klsupport.extrList(y); }
  bool isExtrAllocated(coxtypes::CoxNbr y) const
    { return d_klsupport.isExtrAllocated(y);}
  void ensure_extr_row_exists(coxtypes::CoxNbr y)
    { d_klsupport.allocExtrRow(y); }
  coxtypes::CoxNbr inverse (coxtypes::CoxNbr y) const
    { return d_klsupport.inverse(y);}

  Ulong gen_length(const coxtypes::Generator& s) const { return d_L[s]; }
  Ulong length(const coxtypes::CoxNbr& x) const { return d_length[x]; }
  Ulong size() const { return d_klList.size(); }
  KLStats& stats() { return d_stats; }

  KLRow& klList(coxtypes::CoxNbr y) { return *d_klList[y]; }
  MuTable& muTable(coxtypes::Generator s) { return d_muTable[s]; }
  MuRow& muList(coxtypes::Generator s, coxtypes::CoxNbr y)
    { return *d_muTable[s][y]; }
  bool row_needs_creation(coxtypes::CoxNbr y) const
    { return d_klList[y]==nullptr; }

// modifiers
  const MuPol* to_MuPol(const KLPol& p); // transform |p| to |MuPol| and record
  MuRow mu_row_frame(coxtypes::Generator s, coxtypes::CoxNbr y);
  // ensure a row exists to store polynomials |mu(s,x,y)| in
  void create_mu_row(coxtypes::Generator s, coxtypes::CoxNbr y);
  // whether the row corresponding to y in muTable[s] has been completely filled.
  bool mu_is_complete(coxtypes::Generator s, coxtypes::CoxNbr y);
  // transfer a complete row of computed |MuPol| values for $(y,s)$ to mu-table
  void store_row(const MuRow& row, // already transformed and recorded |MuPol|s
		 coxtypes::Generator s, coxtypes::CoxNbr y);
  void ensure_KL_row(coxtypes::CoxNbr y); // prepare KL row for mu computation
  void compute_mu_row // compute and pass through |store_row| a row of |mu|s
    (coxtypes::Generator s, coxtypes::CoxNbr y);

  const MuPol* compute_mu // bare computation of $\mu(s,x,y)$ (via |to_MuPol|)
    (coxtypes::Generator s, coxtypes::CoxNbr x, coxtypes::CoxNbr y);


  // KL computations
  bool row_is_complete(coxtypes::CoxNbr y); // whether already done for |y|
  // ensure that |d_klsupport| contains all extremal lists that can arise
  void prepareRowComputation (coxtypes::CoxNbr y, coxtypes::Generator s);
  containers::vector<KLPol> initial_polys
   (coxtypes::CoxNbr y, coxtypes::Generator s); // row of values $P(xs,ys)$
  // ensure a row exists to store polynomials |P(x,y)| in
  void create_KL_row(coxtypes::CoxNbr y);
  void add_second_terms(containers::vector<KLPol>& pols,
		  coxtypes::CoxNbr y,
		  coxtypes::Generator s);
  void mu_correct_row(containers::vector<KLPol>& pols,
		      coxtypes::Generator s, coxtypes::CoxNbr y);
  void store_KL_row(containers::vector<KLPol>& pols, coxtypes::CoxNbr y);
  void compute_KL_row(coxtypes::CoxNbr y); // do all computation for a whole row

  const KLPol* compute_KL_pol // bare computation of $P(x,y)$, with storing
    (coxtypes::CoxNbr x, coxtypes::CoxNbr y,
     coxtypes::Generator s = coxtypes::undef_generator);
  void mu_correct_KL_pol
    (KLPol& pol, coxtypes::Generator s, coxtypes::CoxNbr x, coxtypes::CoxNbr y);

  // direct service to |KLContext| methods
  void permute(const bits::Permutation& a);
  void grow(Ulong prev, Ulong n);
  void shrink(const Ulong& n);

  void fill_KL_table ();
  void fill_mu_table ();
  const KLPol& klPol (coxtypes::CoxNbr x, coxtypes::CoxNbr y);
  const MuPol mu
    (coxtypes::Generator s, coxtypes::CoxNbr x, coxtypes::CoxNbr y);
  void move_KL_row_to_inverse(coxtypes::CoxNbr x);
  void apply_inverse_permutation
    (coxtypes::CoxNbr y, const bits::Permutation& a);
  HeckeElt KL_row_as_HeckeElt(coxtypes::CoxNbr y);

private:
  bool no_mu_yet (const coxtypes::Generator& s, const coxtypes::CoxNbr& y) const
    { return d_muTable[s][y] == nullptr; }

}; // |struct KLContext::KLHelper|

/*
  This module contains code for the computation of kl polynomials with
  unequal parameters.

  Even though the general setup is similar to that for the ordinary kl
  polynomials, there are several important technical differences :

    - the polynomials are now Laurent polynomials in the indeterminate q;
      morevoer they need not have positive coefficients.
    - the mu-coefficients are no longer integers, but Laurent polynomials
      themselves; there is one family of mu-coefficients for each conjugacy
      class of generators (equivalently, for each connected component of the
      Coxeter graph, after suppression of the links with an even or infinite
      label.)

  The reference for this construction are Lusztig's course notes.

  The basic datum is for each generator s an element v_s in Z[q^{1/2},q^{-1/2}]
  , a positive integral power of v = q^{1/2}, constant under conjugacy in W.
  Then the Hecke algebra is defined by generators t_s, s in S, and relations :
  (t_s-v_s)(t_s+v_s^{-1}) = 0, plus the braid relations (these t_s are related
  to the "ordinary" T_s by t_s = v_s^{-1}T_s). This implies that t_s.t_w is
  t_{sw} when sw > w, and t_{sw} - a_s.t_w otherwise, where a_s = v_s^{-1}-v_s.
  Also we have t_s^{-1} = t_s + a_s for each s in S.

  Then there is an A-antilinear ring involution on H defined as usual by
  t_x^{-1} = t_{x^{-1}}^{-1}. This defines a formal R-function on the group W,
  hence a Kazhdan-Lusztig basis c_x, x in W. As usual, c_x is the unique
  self-adjoint element in H such that c_x = t_x mod H_{<0}, where H_{<0} is
  the direct sum of the v^{-1}Z[v^{-1}]t_x. It turns out that c_y is of
  the form :

      c_y = t_y + sum_{x<y}p(x,y)t_x

  where p(x,y) is in v^{-1}Z[v^{-1}] for all x < y. It is easy to see that
  c_e = 1, c_s = t_s + v_s^{-1} for all s in S. From this, one can try to
  construct the c_y inductively as in the equal parameter case, as follows.

  Let y > e in W, and let s in S such that sy < y. Then the element c_sc_{sy}
  is known; it is equal to :

      (t_s + v_s^{-1})(t_{ys} + sum_{x<ys}p(x,ys)t_x)

  which is t_y + sum_{z<sy,sz<z}v^s.p(z,ys).t_z mod H_{<0}. The positive part
  of v_s.p(z,sy) is a polynomial in v of degree at most deg(v_s)-1 (whereas
  in the classical case, this degree is at most zero).

  We try to correct this by subtracting an expression of the form

      sum_{z<sy,sz<z}mu^s(z,sy)c_z

  This will be selfadjoint if the mu^s(z,sy) are. The problem is that there
  will be an interaction between the various mu-corrections. To correct just
  the z-coefficent, we could symmetrize the positive part of v_s.p(z,sy);
  but then multiplying c_z with that would create new problems at levels
  below z. So the best thing here is to write what we get at the level of
  the polynomials. For x < y we have :

  p(x,y) = p(sx,sy) if sx < x, x not comparable to sy
         = p(sx,sy)+v_s.p(x,sy)-sum_{x<=z<sy,sz<z}mu^s(z,sy)p(x,z) if sx<x<=sy
	 = v_s^{-1}.p(x,sy) if sx > x, sx not comparable to sy
	 = v_s^{-1}.p(x,sy)+p(sx,sy)-sum_{x<=z<sy,sz<z}mu^s(z,sy)p(x,z)
	     if x < sx <= sy

  so we need in any case the condition (from the second line) :

  v_s.p(x,sy) = sum_{x<=z<sy,sz<z}mu^s(z,sy)p(x,z) mod A_{<0} if sx < x <= sy.
  If we define the mu^s(x,y) by induction on the l(y)-l(x), we see that this
  translates into : mu^s(x,sy) known mod A_{<0}, and hence mu^s(x,sy) known.

  It turns out that these mu^s(x,ys) then also verify the last formula. As
  a corollary, it turns out that the computation of the K-L polynomials, even
  in this case, reduces to the situation of extremal pairs.

  However, the main difference is that the mu-coefficents bear no obvious
  relation to the K-L polynomials. They have to be computed through an
  independent (although related) recursion, and for them it is not so
  clear that there is much extremality reduction.

  We need mu^s(x,y) for x < y, sx < x, sy > y. It is not clear to me to what
  extent the mu^s depend only on the conjugacy class of s. For one thing,
  the domain is not the same; but on the intersection of domains ? Anyway,
  for now I'll make one table for each s.
*/

/*
  This part of the program is of a more exploratory nature than the rest.
*/

/*****************************************************************************

        Chapter I -- The KLContext class

  Just like the ordinary K-L polynomials, the KLContext class holds the
  polynomial tables and the mu-lists, and provides the functions to
  access and output them. Since we restrict ourselves to positive length
  functions, the polynomials can be reduced to extremal pairs; hence
  the polynomial tables are synchronized with the extrList tables in
  the klsupport part (which is shared among all K-L tables.)

  The following functions are provided :

   - constructors and destructors :

     - KLContext(klsupport::KLSupport& kls, CoxGraph& G, Interface& I);
     - ~KLContext();

   - accessors not already inlined :

   - manipulators :

     - applyInverse(const coxtypes::CoxNbr& x) : auxiliary to permute;
     - fillKL() : fills the full K-L table;
     - fillMu() : fills the full mu-tables;
     - klPol(CoxNbr x, CoxNbr y) : maybe computes and returns P_{x,y};
     - row(Permutation& a, CoxNbr y) :
         fills the row for y and returns an appropriate sort of it in |a|;
     - permute(const Permutation&) : applies a permutation to the context;
     - revertSize(const Ulong&) : reverts the size to a previous value;
     - setSize(const Ulong&) : sets the size to a larger value;

 *****************************************************************************/


KLContext::KLContext
  (klsupport::KLSupport& kls,
   const graph::CoxGraph& G, const interface::Interface& I)
  : d_help(nullptr)
{
  d_help = new KLHelper(kls,G,I);
}

KLContext::~KLContext()
{
}

/******** accessors **********************************************************/

const klsupport::KLSupport& KLContext::klsupport() const
  { return d_help->d_klsupport; }
Ulong KLContext::size() const { return d_help->size(); }
const KLRow& KLContext::klList(const coxtypes::CoxNbr& y) const
  { return *d_help->d_klList[y]; }
const MuRow& KLContext::muList
  (const coxtypes::Generator& s, const coxtypes::CoxNbr& y) const
  { return *d_help->d_muTable[s][y]; }


/******** manipulators *******************************************************/

/*
  This function adjusts the size of the context to a context of size n.
*/
void KLContext::setSize(const Ulong& n)
{
  d_help->grow(size(),n);
}

void KLContext::KLHelper::grow(Ulong prev, Ulong n)
{
  error::CATCH_MEMORY_OVERFLOW = true;

  try {
    d_klList.resize(n);

    for (coxtypes::Generator s = 0; s < d_muTable.size(); ++s)
      d_muTable[s].resize(n);
  }
  catch (...) {
    goto revert;
  }

  d_length.reserve(n);
  if (error::ERRNO)
    goto revert;

  error::CATCH_MEMORY_OVERFLOW = false;

  // fill in new Lengths
  for (coxtypes::CoxNbr x = prev; x < n; ++x) {
    coxtypes::Generator s = last(x);
    coxtypes::CoxNbr xs = schubert().shift(x,s);
    assert(xs<x);
    d_length.push_back(d_length[xs] + gen_length(s));
  }

  return;

 revert:
  error::CATCH_MEMORY_OVERFLOW = false;
  shrink(prev);
} // |KLContext::KLHelper::grow|

void KLContext::fillKL() { d_help->fill_KL_table(); }

void KLContext::KLHelper::fill_KL_table ()
{
  for (coxtypes::CoxNbr y = 0; y < size(); ++y)
    if (inverse(y) >= y and not row_is_complete(y))
      compute_KL_row(y);
}

void KLContext::fillMu() { d_help->fill_mu_table(); }

void KLContext::KLHelper::fill_mu_table ()
{
  for (coxtypes::Generator s = 0; s < rank(); ++s)
    for (coxtypes::CoxNbr y = 0; y < size(); ++y)
      if (not schubert().isDescent(y,s) and not mu_is_complete(s,y))
	compute_mu_row(s,y);
}


/*
  This function returns the Kazhdan-Lusztig polynomial P_{x,y}. It is
  assumed that the condition x <= y has already been checked, and that
  x and y are valid context numbers.
*/
const KLPol& KLContext::klPol (coxtypes::CoxNbr x, coxtypes::CoxNbr y)
{ return d_help->klPol(x,y); }

const KLPol& KLContext::KLHelper::klPol (coxtypes::CoxNbr x, coxtypes::CoxNbr y)
{ // we modify |x| and |y| as local variables, but $P(x,y)$ remains the same

  const schubert::SchubertContext& p = schubert();

  // go to inverses if this helps
  if (inverse(y) < y) {
    y = inverse(y);
    x = inverse(x);
  }

  // ensure that |extrList[y]| is allocated
  if (row_needs_creation(y)) {
    create_KL_row(y);
    if (error::ERRNO)
      return errorPol();
  }

  // move up |x| to make it extremal for |y|
  x = p.maximize(x,p.descent(y));

  // find |x| in extrList[y]
  const auto& eL = extrList(y);
  Ulong m = std::lower_bound(eL.begin(),eL.end(),x)-eL.begin();

  // look up to see if $P(x,y)$ is already stored
  const KLPol* pol = (*d_klList[y])[m];

  // if not, compute it explicitly
  if (pol == nullptr) { /* we have to compute the polynomial */
    pol = compute_KL_pol(x,y);
    if (error::ERRNO)
      return errorPol();
  }

  return *pol;
}


/*
  This function returns mu^s_{x,y}, filling it in if necessary. It is
  assumed that the conditions x < y, ys > y, xs < x have already been
  checked.
*/
const MuPol KLContext::mu(const coxtypes::Generator& s,
			  const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y)
{ return d_help->mu(s,x,y); }

const MuPol KLContext::KLHelper::mu
  (coxtypes::Generator s, coxtypes::CoxNbr x, coxtypes::CoxNbr y)
{
  // either use existing row or start one from scratch
  if (no_mu_yet(s,y))
    create_mu_row(s,y);
  const MuRow& mu_row = muList(s,y);

  // find |x| in that row
  const MuData mx(x,nullptr);
  auto it = std::lower_bound(mu_row.begin(),mu_row.end(),mx);

  // absent |x|s certainly have $\mu(s,x,y)=0$
  if (it == mu_row.end() or it->x!=x)
    return MuPol::zero();

  // get location for |mu(s,x,y)|
  const MuPol* mp = it->pol;

  // if that mu polynomial is not yet computed, do so now
  if (mp == nullptr) {
    mp = compute_mu(s,x,y);
    if (error::ERRNO)
      return errorMuPol();
  }

  return *mp; // return found or coputed mu polynomial
} // |KLContext::KLHelper::mu(s,x,y)|

/*
  This function makes sure that the row corresponding to y in the K-L table
  is entirely filled, and returns the corresponding data, sorted in
  the order of increasing context numbers.
*/
void KLContext::row(HeckeElt& h, const coxtypes::CoxNbr& y)
{ h = d_help->KL_row_as_HeckeElt(y); }

HeckeElt KLContext::KLHelper::KL_row_as_HeckeElt(coxtypes::CoxNbr y)
{
  HeckeElt h;
  if (not row_is_complete(y)) {
    d_klsupport.allocRowComputation(y);
    if (error::ERRNO)
      goto error_exit;
    compute_KL_row(y);
    if (error::ERRNO)
      goto error_exit;
  }

  { if (y <= inverse(y)) {
      const klsupport::ExtrRow& e = extrList(y);
      h.reserve(e.size());
      const KLRow& klr = klList(y);
      for (Ulong j = 0; j < e.size(); ++j) {
	h.emplace_back(e[j],klr[j]);
      }
    }
    else { /* go over to inverses */
      coxtypes::CoxNbr yi = inverse(y);
      const klsupport::ExtrRow& e = extrList(yi);
      h.reserve(e.size());
      const KLRow& klr = klList(yi);
      for (Ulong j = 0; j < e.size(); ++j) {
	h.emplace_back(inverse(e[j]),klr[j]);
      }
      std::sort(h.begin(),h.end()); // make sure list is ordered
    }
  }

  return h;

 error_exit:
  error::Error(error::ERRNO);
  error::ERRNO = error::ERROR_WARNING;
  return h;
} // |KLContext::row|



/*
  Reverts the sizes of the lists to size n. This is meant to be used
  only immediately after a failing context extension, to preserve the
  consistency of the various list sizes. In particular, it will fail
  miserably if a premutation has taken place in between.
*/
void KLContext::revertSize(const Ulong& n)
{
  d_help->shrink(n);
}

void KLContext::KLHelper::shrink(const Ulong& n)
{
  d_klList.resize(n);

  for (coxtypes::Generator s = 0; s < d_muTable.size(); ++s) {
    MuTable& t = d_muTable[s];
    t.resize(n); // drop any extension of |d_muTable[s]|
  }
  d_length.resize(n);
}



/*
  Exchange rows for |x| and |inverse(x)| in klList. It is assumed that the row
  for |inverse(x)| is allocated.
*/
void KLContext::applyInverse(const coxtypes::CoxNbr& x)
{ d_help->move_KL_row_to_inverse(x);
}

void KLContext::KLHelper::move_KL_row_to_inverse(coxtypes::CoxNbr x)
{
  coxtypes::CoxNbr xi = inverse(x);
  d_klList[x] = std::move(d_klList[xi]);
}

void KLContext::applyIPermutation
  (const coxtypes::CoxNbr& y, const bits::Permutation& a)
{ d_help->apply_inverse_permutation(y,a); }

void KLContext::KLHelper::apply_inverse_permutation
  (coxtypes::CoxNbr y, const bits::Permutation& a)
{ right_permute(*d_klList[y],a); }


/*
  Apply the permutation |a| to the context. See the |permute| function of
  |klsupport::KLSupport| for a detailed explanation.
*/
void KLContext::permute(const bits::Permutation& a)
{ d_help->permute(a); }

void KLContext::KLHelper::permute(const bits::Permutation& a)
{
  // permute values inside each row of |muTable|
  for (coxtypes::Generator s = 0; s < d_muTable.size(); ++s) {
    MuTable& t = d_muTable[s];
    for (coxtypes::CoxNbr y = 0; y < size(); ++y) {
      if (no_mu_yet(s,y))
	continue;
      MuRow& row = *t[y];
      for (Ulong j = 0; j < row.size(); ++j)
	row[j].x = a[row[j].x];
      std::sort(row.begin(),row.end());
    }
  }

  // permute ranges
  bits::BitMap seen(a.size());

  for (coxtypes::CoxNbr x = 0; x < size(); ++x)
  {
    if (seen.getBit(x))
      continue;
    // now |x| starts an as yet unseen nontrivial cycle in |a|
    containers::vector<MuRowPtr> mu_buf(d_muTable.size()); // Coxeter rank

    for (coxtypes::CoxNbr y = a[x]; y != x; y = a[y])
    {
      std::swap(d_length[y],d_length[x]);

      // back up values for y
      auto kl_buf = std::move(d_klList[y]);
      for (coxtypes::Generator s = 0; s < d_muTable.size(); ++s) {
	MuTable& t = d_muTable[s];
	mu_buf[s] = std::move(t[y]);
      }
      // move values from |x| in |y|
      d_klList[y] = std::move(d_klList[x]);
      for (coxtypes::Generator s = 0; s < d_muTable.size(); ++s) {
	MuTable& t = d_muTable[s];
	t[y] = std::move(t[x]);
      }

      /* store backed up values in x */
      d_klList[x] = std::move(kl_buf);
      for (coxtypes::Generator s = 0; s < d_muTable.size(); ++s) {
	MuTable& t = d_muTable[s];
	t[x] = std::move(mu_buf[s]);
      }
      seen.setBit(y);
    }  // |for(y)|

    seen.setBit(x);
  } // |for(x)|
} // |KLContext::KLHelper::permute|



/*****************************************************************************

        Chapter II -- The KLHelper class.

  The purpose of the KLHelper class is to hide from the public eye a number
  of helper functions, used in the construction and maintenance of the
  K-L context. This unclutters kl.h quite a bit.

  The following functions are defined :

    - create_KL_row(y) : allocates one row in the K-L table;
    - mu_row_frame(s,y) : a basic skeleton for a mu-row for y and s;
    - create_mu_row(s,y) : creates row |muTable(s)[y]| with x's but no mu-polys
    - row_is_complete(y) : whether row for y (or inverse(y)) is incomplete
      in the K-L table is incompletely filled;
    - compute_KL_row(y) : fill out the row for y or inverse(y);
    - compute_KL_pol(x,y) : compute P_{x,y} and store it into the table
    - pol = initial_polys(y,s) : auxiliary to |compute_KL_row|;
    - mu_correct_row(pols,s,y) : auxiliary to |compute_KL_row|;
    - prepareRowComputation(y,s) : auxiliary to |compute_KL_row|;
    - add_second_terms(pol,y,s) : ausiliary to |compute_KL_row|;
    - store_KL_row(pols,y) : auxiliary to |compute_KL_row|;

 *****************************************************************************/
/*
  This (unique) constructor gets the lengths interactively from the user. This
  makes it possible that an error is set during the construction. Fortunately we
  can check this before any memory is gotten from the heap, so that automatic
  destruction of the components on exit will be satisfactory in that case.
*/
KLContext::KLHelper::KLHelper
  (klsupport::KLSupport& kls,
   const graph::CoxGraph& G, const interface::Interface& I)
: d_klsupport(kls)
, d_klList(kls.size())
, d_muTable()
, d_L(2*kls.rank()) // the size expected by |interactive::getLength|
, d_length{0} // start with a single length entry 0
, KL_pool() // empty set of |KLPol|
, mu_pool() // empty set of |MuPol|
, d_stats()
{
  interactive::getLength(d_L,G,I);

  if (error::ERRNO) /* error code is ABORT */
    return;

  d_length.reserve(kls.size()); // |size()| also works, but is no different

  for (coxtypes::CoxNbr x = 1; x < size(); ++x) {
    coxtypes::Generator s = last(x);
    coxtypes::CoxNbr xs = schubert().shift(x,s);
    d_length.push_back(d_length[xs] + d_L[s]);
  }

  d_klList[0].reset(new KLRow(1,KL_pool.find(one())));

  d_muTable.reserve(rank());

  for (coxtypes::Generator s = 0; s < rank(); ++s) {
    d_muTable.emplace_back(kls.size()); // create a |MuTable| of |kls->size()|
    MuTable& t = d_muTable.back();
    t[0].reset(new MuRow); // create one entry (empty list) in initial slot
  }


  // moved from |KLContext| constructor
  d_stats.klrows++;
  d_stats.klnodes++;
  d_stats.klcomputed++;
} // |KLContext::KLHelper::KLHelper|


/*
  Allocate one previously unallocated row in the K-L table. It is assumed
  that y <= inverse(y). Allocate the corresponding extremal row if necessary.
*/
void KLContext::KLHelper::create_KL_row(coxtypes::CoxNbr y)
{
  ensure_extr_row_exists(y);

  Ulong n_extremals = extrList(y).size();

  d_klList[y].reset(new KLRow(n_extremals,nullptr));

  stats().klnodes += n_extremals;
  stats().klrows++;
}


/*
  Return an empty |MuRow| with sufficient |x| values for $(y,s)$, and no |MuPol|
  values yet. Not all |x| values will necessarily have a nonzero |mu|.

  For unequal parameters, we don't try to be particularly smart. The row
  for y contains one entry for each x < y s.t. xs < x (it is already
  implicit that ys > y, or the row would not even be allocated.)
*/
MuRow KLContext::KLHelper::mu_row_frame
  (coxtypes::Generator s, coxtypes::CoxNbr y)
{
  bits::BitMap b(0);
  schubert().extractClosure(b,y);
  b &= schubert().downset(s);

  MuRow row;

  for (coxtypes::CoxNbr x : b)
    row.emplace_back(x,nullptr);

  return row;
}

// Create an empty |MuRow| ready to store values |mu(x,y,s)| for various |x|.
void KLContext::KLHelper::create_mu_row
  (coxtypes::Generator s, coxtypes::CoxNbr y)
{
  muTable(s)[y].reset(new MuRow(mu_row_frame(s,y)));
  d_stats.munodes += muList(s,y).size();
  d_stats.murows++;
}



// Whether the row corresponding to |y| in the K-L table is completely filled.
// Actual the check is for |inverse(y)| if |inverse(y) < y|.
bool KLContext::KLHelper::row_is_complete(coxtypes::CoxNbr y)
{
  if (inverse(y) < y)
    y = inverse(y);

  if (row_needs_creation(y))
    return false;

  for (auto entry : klList(y))
    if (entry==nullptr)
      return false;

  return true;
}


// Whether the row corresponding to y in muTable[s] has been completely filled.
// Contrary to the KL row, no symmetry with respect to |inverse| can be used.
bool KLContext::KLHelper::mu_is_complete
  (coxtypes::Generator s, coxtypes::CoxNbr y)
{
  const MuTable& t = muTable(s);

  if (t[y] == nullptr)
    return false;

  const MuRow& mr = *t[y];

  for (Ulong j = 0; j < mr.size(); ++j) {
    if (mr[j].pol == nullptr)
      return false;
  }

  return true;
}


// Make sure that the K-L row for |y| is available.
void KLContext::KLHelper::ensure_KL_row(coxtypes::CoxNbr y)
{
  if (not row_is_complete(y)) {
    d_klsupport.allocRowComputation(y);
    if (error::ERRNO)
      goto abort;
    compute_KL_row(y);
    if (error::ERRNO)
      goto abort;
  }
  return;

 abort:
  error::Error(error::ERRNO);
  error::ERRNO = error::ERROR_WARNING;
}


/*
  This function computes a single polynomial in the K-L table (as opposed to
  |compute_KL_row|, which fills a whole row.) It isn't particularly optimized for
  speed; it is not a good idea to fill large parts of the table by repeated
  calls to compute_KL_pol.

  It is assumed that |inOrder(x,y)| (Bruhat order), that |y <= inverse(y)|, that
  |x| is extremal for |y|, and that the row for |y| in the K-L table is
  allocated.
*/
const KLPol* KLContext::KLHelper::compute_KL_pol
  (coxtypes::CoxNbr x, coxtypes::CoxNbr y, coxtypes::Generator d_s)
{
  const schubert::SchubertContext& p = schubert();

  coxtypes::Generator s = d_s;

  /* If d_s is undef_coxnbr, we compute the polynomial using descent by last
     term in normal form */

  if (s == coxtypes::undef_generator)
    s = last(y);

  coxtypes::CoxNbr ys = p.shift(y,s);
  coxtypes::CoxNbr xs = p.shift(x,s);

  /* check if x is comparable to ys */

  if (!p.inOrder(x,ys)) { /* return the answer immediately */
    stats().klcomputed++;
    const auto& eL = extrList(y);
    Ulong m = std::lower_bound(eL.begin(),eL.end(),x)-eL.begin();
    klList(y)[m] = &klPol(xs,ys);
    return klList(y)[m];
  }

  /* get workspace */

  error::CATCH_MEMORY_OVERFLOW = true;

  KLPol pol; // workspace

  /* initialize the workspace to P_{xs,ys} */

  const KLPol& p_xsys = klPol(xs,ys);
  if (error::ERRNO)
    goto abort;
  pol = p_xsys;

  /* add q.P_{x,ys} */

  {
    const KLPol& p_xys = klPol(x,ys);
    if (error::ERRNO)
      goto abort;
    pol.add(p_xys,gen_length(s));
    if (error::ERRNO)
      goto abort;
  }

  /* subtract correction terms */

  mu_correct_KL_pol(pol,s,x,y);
  if (error::ERRNO)
    goto abort;

  /* find address of polynomial */

  {
    const KLPol* p_xy = KL_pool.find(pol);
    if (error::ERRNO)
      goto abort;
    const auto& eL = extrList(y);
    Ulong m = std::lower_bound(eL.begin(),eL.end(),x)-eL.begin();
    klList(y)[m] = p_xy;

    /* return workspace and exit */

    error::CATCH_MEMORY_OVERFLOW = false;

    stats().klcomputed++;
    return p_xy;
  }

 abort: /* an error occurred */

  error::CATCH_MEMORY_OVERFLOW = false;
  if (error::ERRNO != error::MEMORY_WARNING)
    error::ERRNO = error::KL_FAIL;
  return nullptr;
}


/*
  This function fills one row in the K-L table entirely. This can be done
  rather more efficiently than computing each polynomial individually :
  in particular, most of the closure computations can be "factored" for the
  whole row at a time.

  It is assumed that checkKLRow(y) returns false. The row which is actually
  filled is the one for the smaller of (y,inverse(y)).
*/
void KLContext::KLHelper::compute_KL_row (coxtypes::CoxNbr y)
{
  if (inverse(y) < y) /* fill in the row for inverse(y) */
    y = inverse(y);

  { // group so that jumps to |abort| do not "cross" variable declarations
    if (row_needs_creation(y))
      create_KL_row(y);
    if (error::ERRNO)
      goto abort;

    const coxtypes::Generator s = last(y); // choose any proper |s| for |y|

    // make sure the necessary terms are available
    prepareRowComputation(y,s);
    if (error::ERRNO)
      goto abort;

    // initialize row of polynomials |pol| with $P_{xs,ys}$
    auto pols = initial_polys(y,s);
    if (error::ERRNO)
      goto abort;

    // add $q.P_{x,ys}$ to each $P(x,y)$ under construction with $x\leq ys$
    add_second_terms(pols,y,s);
    if (error::ERRNO)
      goto abort;

    // subtract correcting terms
    mu_correct_row(pols,s,y);
    if (error::ERRNO)
      goto abort;

    // copy results to row |y| of the table
    store_KL_row(pols,y);
    if (error::ERRNO)
      goto abort;

    return;
  }
 abort:
  error::Error(error::ERRNO);
  error::ERRNO = error::ERROR_WARNING;
  return;
}


/*
  Compute the mu-polynomial for $\mu(s,x,y)$, recording it in |kl.mu_pool|. It
  is assumed that $x < y$, $sy > y$, $sx < x$. This mu-polynomial is the unique
  symmetric Laurent polynomial whose positive part is the same as that of:

    q^{length(s)/2}P(x,y) - \sum_{x<z<y,zs<s}\mu(s,z,y)P(x,z)

  Hence it can be computed inductively if we assume that $P(x,y)$ is already
  known, and also the $\mu(s,z,y)$ for $z$ in $]x,y[$ (this is the case in the
  way our K-L computation is set up; when we need $\mu(s,x,y)$, the
  corresponding $P(x,y)$ is known.) Here the $P$'s are the actual Laurent
  polynomials, in $\Z[u]$ where $u^{-2}=q$; so some shifting is required when we
  read them from our |klPol|'s.

  Returns |nullptr| in case of error, otherwise a pointer into the
  |d_muTable| of our |KLContext|.
*/
const MuPol* KLContext::KLHelper::compute_mu
  (coxtypes::Generator s,
   coxtypes::CoxNbr x, coxtypes::CoxNbr y)
{
  KLPol pos_mu;
  MuRow& mu_row = muList(s,y);

  // initialize |pol| with the value $u^(L(x)-L(y)+L(s))P_{x,y}(u^2)$
  const KLPol& pol = klPol(x,y);
  if (error::ERRNO) {
    error::Error(error::MU_FAIL,x,y);
    error::ERRNO = error::ERROR_WARNING;
    return nullptr;
  }

  // extract the posisitve part, with respect ot degree $L(x)-L(y)+L(s)$
  pos_mu = positivePart(pol,2,length(x)-length(y)+gen_length(s));

  MuData mx(x,nullptr);
  auto loc_x =
    std::lower_bound(mu_row.begin(),mu_row.end(),mx); // search cannot fail

  // subtract correcting terms

  const schubert::SchubertContext& p = schubert();
  for (auto it = std::next(loc_x); it!=mu_row.end(); ++it)
  {
    coxtypes::CoxNbr z = it->x;
    if (not p.inOrder(x,z))
      continue;
    const KLPol& pol = klPol(x,z);
    if (error::ERRNO) {
      error::Error(error::MU_FAIL,x,y);
      error::ERRNO = error::ERROR_WARNING;
      return 0;
    }
    const MuPol mq = mu(s,z,y);
    if (not mq.isZero())
      muSubtraction(pos_mu,mq,pol,2,length(x)-length(z));
    if (error::ERRNO) {
      error::Error(error::MU_FAIL,x,y);
      error::ERRNO = error::ERROR_WARNING;
      return 0;
    }
  }

  // convert |pos_mu| to a mu-polynomial, then store and return it
  return loc_x->pol = to_MuPol(pos_mu);
} // |compute_mu|


/*
  This function fills one row in the mu-list for s. Recall that mu(s,x,y) is
  defined whenever x < y, sy > y, sx < x, and is a Laurent polynomial in
  $u=q^{1/2}$, which is symmetric w.r.t. the involution that sends $q->q^-1$.
  See |compute_mu| above for the formula defining |mu|.

  In order to avoid huge amounts of calls to the expensive method |inOrder|, we
  proceed as in |compute_KL_row|, computing the full row at a time (setting all
  their initial expressions first, then subtracting correction terms). This
  appears to be a bit more difficult than for the K-L-polynomials, because in
  the recursion we need mu's for the _same_ value of $y$, but as it turns out,
  when we need a $\mu(s,z,y)$, it is already fully computed, provided we proceed
  with the correcting terms in decreasing order.

  A number of K-L polynomials are needed in the process; it turns out that
  when one is needed, usually many will be for the same value of z; so
  again we fill the whole K-L row for z when a polynomial is needed. The
  problem with this is that it may trigger recursive calls to |compute_mu_row|
  [...]

  [It appears the implicit recursion happens at the call |ensure_KL_row(y)|.
  The problem that Fokko refers to can be solved simply by using non-|static|
  local variables |mu_row| end |pos_mu| that will automatically create a fresh
  instance in any recursive call. Fokko laboriously implemented stacks in
  |static| variables instead, which really seems a misguided solution. MvL]

*/
void KLContext::KLHelper::compute_mu_row
  (coxtypes::Generator s, coxtypes::CoxNbr y)
{
  ensure_KL_row(y); // ensure that calls |klPol(x,y)| will directly return a result
  MuRow mu_row = mu_row_frame(s,y); // construct a |MuRow| into |mu_row|
  containers::vector<KLPol> pos_mu(mu_row.size());

  coxtypes::CoxNbr x;

  // initialize |pos_mu((x))| with value from $u^(L(x)-L(y)+L(s))P_{x,y}[u^2]$

  for (Ulong j = 0; j < mu_row.size(); ++j)
  {
    x = mu_row[j].x;
    const KLPol& pol = klPol(x,y);
    if (error::ERRNO) // this should not happen in typical usage
      goto abort;
    pos_mu[j] = positivePart(pol,2,length(x)-length(y)+gen_length(s));
  }

  /* we run through mu_row in decreasing order; for each z in the list,
     we subtract from the correction term corresponding to z from mu(x,y)
     for each x < z.
  */

  for (Ulong j = mu_row.size(); j-->0;)
  { // this loop modifies |pos_mu[i]| for certain |i<j|
    mu_row[j].pol = to_MuPol(pos_mu[j]); // now consolidate |pos_mu[j]|
    d_stats.mucomputed++;
    if (mu_row[j].pol->isZero()) {
      d_stats.muzero++;
      continue;
    }

    /* subtract correcting terms */

    coxtypes::CoxNbr z = mu_row[j].x;
    ensure_KL_row(z); // ensure that |klPol(x,z)| will work
    if (error::ERRNO)
      goto abort;
    bits::BitMap b(0);
    schubert().extractClosure(b,z); // set |b| to interval $[e,z]$
    b &= schubert().downset(s);
    b.clearBit(z);
    bits::BitMap::Iterator b_end = b.end();

    Ulong i = 0;

    for (bits::BitMap::Iterator k = b.begin(); k != b_end; ++k) {
      x = *k;
      while (mu_row[i].x != x) // advance |i| to entry for |x|
	++i;
      const KLPol& pol = klPol(x,z);
      if (error::ERRNO)
	goto abort;
      muSubtraction(pos_mu[i],mu_row[j].pol[0],pol,2,
		    length(x)-length(z));
      if (error::ERRNO)
	goto abort;
      ++i;
    } // |for(k)| and |x|
  } // |for(j)| and |z|

  store_row(mu_row,s,y);

  return;

 abort:
  error::Error(error::MU_FAIL,x,y);
  error::ERRNO = error::ERROR_WARNING;

}


/*
  This function sets pol to a row of one polynomial for each x in klList(y),
  and initializes the corresponding pol[j] to klPol(xs,ys).

  It is assumed that |prepareRowComputation| has been called for y and s,
  so that the row for ys is available.
*/
containers::vector<KLPol>
  KLContext::KLHelper::initial_polys(coxtypes::CoxNbr y,
				     coxtypes::Generator s)
{
  const schubert::SchubertContext& p = schubert();
  const klsupport::ExtrRow& e = extrList(y);
  containers::vector<KLPol> pol;

  try {
    pol.reserve(e.size());
  }
  catch(...) {
    error::Error(error::MEMORY_WARNING);
    error::ERRNO = error::ERROR_WARNING;
    return pol;
  }

  /* initialize with values P_{xs,ys} */

  {
    coxtypes::CoxNbr ys = p.rshift(y,s);

    for (Ulong j = 0; j < e.size(); ++j) {
      coxtypes::CoxNbr xs = p.shift(e[j],s);
      pol.push_back(klPol(xs,ys)); // no error can occur here
    }
  }

  return pol;
}




/*
  This function carries out the "mu-correction" step in the computation of
  the row of K-L polynomials. It is assumed that pol holds one polynomial
  for each x <= y extremal w.r.t. y, that we are applying recursion w.r.t.
  s, and that the polynomials have been initialized to P_{xs,ys}
  +q^gen_length(s)P_{x,ys}.

  We have to subtract terms q^{(L(y)-L(z))/2}P_{x,z}mu(z,ys); from results
  of Lusztig it is known that these are actually polynomials in q as well.

  We minimize the number of calls to the Bruhat order functions by proceeding
  as follows : for each z in muList(s,ys), we extract the interval [e,z],
  extremalize w.r.t. y, and subtract the appropriate term from P_{x,y}.
*/
void KLContext::KLHelper::mu_correct_row
  (containers::vector<KLPol>& pols,
   coxtypes::Generator s, coxtypes::CoxNbr y)
{
  const schubert::SchubertContext& p = schubert();
  const klsupport::ExtrRow& e = extrList(y);
  coxtypes::CoxNbr ys = p.rshift(y,s);
  const MuRow& mu_row = muList(s,ys);

  for (Ulong j = 0; j < mu_row.size(); ++j)
  {
    const MuPol& mu_pol = *mu_row[j].pol;
    if (mu_pol.isZero())
      continue;

    coxtypes::CoxNbr z = mu_row[j].x;
    bits::BitMap b(size());
    p.extractClosure(b,z);
    schubert::select_maxima_for(p,b,p.descent(y));

    Ulong i = 0;
    bits::BitMap::Iterator b_end = b.end();

    for (bits::BitMap::Iterator k = b.begin(); k != b_end; ++k)
    {
      coxtypes::CoxNbr x = *k;
      while (e[i] < x)
	++i;
      Ulong h = length(y)-length(z);
      pols[i].subtract(klPol(x,z),mu_pol,h);
      if (error::ERRNO)
      {
	error::Error(error::ERRNO,this,x,y);
	error::ERRNO = error::ERROR_WARNING;
	return;
      }
    } // |for(k)|
  } // |for(j)|
} // |mu_correct_row|


/*
  This function carries out the "mu-correction" step in the computation of
  a single K-L polynomial |pol|. It is assumed that |pol| has been
  initialized to P_{xs,ys}+q^gen_length(s)P_{x,ys}.

  We have to subtract terms q^{(L(y)-L(z))/2}P_{x,z}mu(z,ys); from results
  of Lusztig it is known that these are actually polynomials in q as well.
*/
void KLContext::KLHelper::mu_correct_KL_pol
 (KLPol& pol, coxtypes::Generator s, coxtypes::CoxNbr x, coxtypes::CoxNbr y)
{
  const schubert::SchubertContext& p = schubert();
  coxtypes::CoxNbr ys = p.rshift(y,s);
  if (no_mu_yet(s,ys)) {
    create_mu_row(s,ys);
    if (error::ERRNO)
      goto abort;
  }

  {
    const MuRow& mu_row = muList(s,ys);

    for (Ulong j = 0; j < mu_row.size(); ++j) {

      coxtypes::CoxNbr z = mu_row[j].x;
      if (!p.inOrder(x,z))
	continue;

      const MuPol mp = mu(s,z,ys);
      if (mp.isZero())
	continue;

      Ulong h = length(y)-length(z);
      const KLPol& kl_pol = klPol(x,z);
      if (error::ERRNO)
	goto abort;
      pol.subtract(kl_pol,mp,h);
      if (error::ERRNO)
	goto abort;
    }
  }

  return;


 abort:
  error::Error(error::UEMU_FAIL,x,y);
  error::ERRNO = error::ERROR_WARNING;

} // |mu_correct_KL_pol|


/*
  This auxiliary function to |compute_KL_row| makes sure that the necessary terms
  for the filling of the row of |y| in the K-L table are available. This ensures
  that there will be no recursive calls to |compute_KL_row| once the actual
  computation starts.

  At this point it is assumed that y <= inverse(y).
*/

void KLContext::KLHelper::prepareRowComputation(coxtypes::CoxNbr y,
						coxtypes::Generator s)
{
  const schubert::SchubertContext& p = schubert();

  coxtypes::CoxNbr ys = p.rshift(y,s);

  // ensure all the KL polynomials $P(x,ys)$ are computed and stored
  if (not row_is_complete(ys)) {
    compute_KL_row(ys);
    if (error::ERRNO)
      goto abort;
  }

 { // ensure all the mu polynomials $\mu(s,x,ys)$ are computed and stored
    if (not mu_is_complete(s,ys)) {
      compute_mu_row(s,ys);
      if (error::ERRNO)
	goto abort;
    }

    const MuRow& mu_row = muList(s,ys);

    for (Ulong j = 0; j < mu_row.size(); ++j) {
      if (mu_row[j].pol->isZero())
	continue;
      coxtypes::CoxNbr z = mu_row[j].x;
      if (not row_is_complete(z)) { // ensure extremal lists in descent path |z|
	d_klsupport.allocRowComputation(z); // are all created in |klsupport()|
	if (error::ERRNO)
	  goto abort;
	compute_KL_row(z);
	if (error::ERRNO)
	  goto abort;
      }
    }
  }

  return;

 abort:
  error::Error(error::ERRNO);
  error::ERRNO = error::ERROR_WARNING;
  return;
}


/*
  This function adds the "second term", which is q^gen_length(s).P_{x,ys}, to the
  polynomials in |pol|. It is assumed that y <= inverse(y).

  In order to avoid calls to inOrder, we proceed as follows : we extract
  [e,ys], we extremalize it w.r.t. the descent set of y, and run through
  it to make the correction; this makes us run exactly through those
  x in the extremal list of y which are <= ys.
*/
void KLContext::KLHelper::add_second_terms
  (containers::vector<KLPol>& pol, coxtypes::CoxNbr y, coxtypes::Generator s)
{
  const schubert::SchubertContext& p = schubert();
  bits::BitMap b(size());
  coxtypes::CoxNbr ys = p.rshift(y,s);

  p.extractClosure(b,ys);
  schubert::select_maxima_for(p,b,p.descent(y));

  Ulong i = 0;
  const klsupport::ExtrRow& e = extrList(y);

  for (coxtypes::CoxNbr x : b)
  {
    while(e[i] < x) // advance in |extrList(y)| to synchronize
      ++i;
    assert(e[i]==x);

    pol[i].add(klPol(x,ys),gen_length(s));
    if (error::ERRNO) {
      error::Error(error::ERRNO,this,x,y);
      error::ERRNO = error::ERROR_WARNING;
      return;
    }
    ++i; // now we can advance one step unconditionally
  }
}


/*
  This function writes the polynomials from the list pol to klList(y);
  more precisely, it finds their adresses in |KL_pool|, and writes those
  to klList(y).

  It is assumed that y <= inverse(y).

  The only error that can occur here is memory overflow because of the
  allocation for new polynomials in |KL_pool|. In that case, the error
  is treated, and error::ERROR_WARNING is set.
*/
void KLContext::KLHelper::store_KL_row
  ( containers::vector<KLPol>& pols, coxtypes::CoxNbr y)
{
  KLRow& kl_row = klList(y);

  for (Ulong j = 0; j < kl_row.size(); ++j)
    if (kl_row[j] != nullptr)
      assert(pols[j] == *kl_row[j]); // check match with previously known value
    else
    {
      const KLPol* q = KL_pool.find(pols[j]);
      if (q == nullptr) { /* an error occurred */
	error::Error(error::ERRNO);
	error::ERRNO = error::ERROR_WARNING;
	return;
      }
      kl_row[j] = q; // store newly computed (and recorded) polynomial
      stats().klcomputed++;
    } // |else| and |for|
}

const MuPol* KLContext::KLHelper::to_MuPol(const KLPol& p)
{
  if (p.isZero())
    return mu_pool.find(MuPol::zero());

  // now we must symmetrize |p| before looking up in |t|
  MuPol mp(p.deg(),-p.deg());
  mp[0] = p[0];

  for (Ulong j = 1; j <= p.deg(); ++j)
  {
    mp[-static_cast<long>(j)] = p[j];
    mp[j] = p[j];
  }

  return mu_pool.find(mp);
} // |KLContext::KLHelper::to_MuPol|

/*
  This function copies |row| elements to the corresponding row in the mu-table
  for |s| at |y|, omitting the zero terms.
*/
void KLContext::KLHelper::store_row(const MuRow& row,
				    coxtypes::Generator s,
				    coxtypes::CoxNbr y)
{
  containers::sl_list<MuData> result;

  // collect nonzero terms
  for (Ulong j = 0; j < row.size(); ++j)
    if (not row[j].pol->isZero())
      result.push_back(row[j]);

  //stash away in |muTable|
  muTable(s)[y].reset(new MuRow(result.to_vector()));
} // |store_row|


/*****************************************************************************

        Chapter III -- The KLPol class.

  The KLPol class is derived form Polynomial<SKLcoeff>, because we
  want to re-define the arithmetic operations so that overflow is carefully
  checked. This makes them expensive, but arithmetic is only used when the
  polynomials are defined, and there we have to check anyway.

  The following functions are defined :

    - add(p,n) : adds p shifted by n to the current polynomial;
    - subtract(p,mp,n) : subtracts from the current polynomial the product of
      p and the MuPol mp, shifted by n;

 *****************************************************************************/

KLPol& KLPol::add(const KLPol& p, const long& n)

/*
  Increments the polynomial by p shifted by n, i.e. X^n.p, while checking
  that the coefficients of the result remain within bounds.

  NOTE : a correct implementation would check beforehand the size of the
  result, so as not to waste memory; we are content with setting the size
  to the correct value after the fact. This doesn't matter as this function
  will be used only on temporaries.
*/

{
  /* set degree and valuation of the result */

  if (deg() < p.deg()+n) {
    setDeg(p.deg()+n);
  }

  for (polynomials::Degree j = 0; j <= p.deg(); ++j) {
    klsupport::safeAdd((*this)[j+n],p[j]);
    if (error::ERRNO)
      return *this;
  }

  reduceDeg();

  return *this;
}

KLPol& KLPol::subtract(const KLPol& p, const MuPol& mp, const Ulong& n)

/*
  This function subtracts from the current polynomial the polynomial p*mu
  shifted by q^{n/2}. Here mp is a MuPol, i.e., a Laurent polynomial in
  q^{1/2}; it is assumed that n is such that mp*q^{n/2} is a polynomial in
  q.

  It is known that for unequal parameters, negative coefficients can occur
  in K-L polynomials. So we only check for overflow during the computation,
  and set the error KLCOEFF_OVERFLOW or KLCOEFF_UNDERFLOW accordingly.
*/

{
  KLPol q(0);
  q.setDeg((mp.deg()+n)/2);

  for (long j = mp.val(); j <= mp.deg(); ++j) {
    if (mp[j] == 0)
      continue;
    /* if we get here, n + j is even */
    q[(n+j)/2] = mp[j];
  }

  /* compute the product and check for overflow */

  for (Ulong i = 0; i <= q.deg(); ++i) {
    if (q[i] == 0)
      continue;
    for (Ulong j = 0; j <= p.deg(); ++j) {
      klsupport::SKLcoeff a = p[j];
      klsupport::safeMultiply(a,q[i]);
      if (error::ERRNO)
	return *this;
      if (isZero() || (i+j) > deg())
	setDeg(i+j);
      klsupport::safeAdd(v[i+j],-a);
      if (error::ERRNO)
	return *this;
    }
  }

 reduceDeg();
 return *this;
}

/****************************************************************************

        Chapter IV -- Kazhdan-Lustig bases.

  This section defines functions returning Kazhdan-Lusztig bases. Note that
  the coefficients of these bases should actually be Laurent polynomials;
  however we (perhaps mistakenly) leave it for now to the output functions
  to do the shifting that is required; this saves us from introducing
  a new type at this point.

  The following functions are defined :

    - cBasis(h,y,kl) : also called sometimes the C'-basis; in our opinion,
      the right choice of basis;

 ****************************************************************************/


void cBasis(HeckeElt& h, const coxtypes::CoxNbr& y, KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  bits::BitMap b(0);
  p.extractClosure(b,y);

  bits::BitMap::Iterator b_end = b.end();
  h.clear();

  for (bits::BitMap::Iterator x = b.begin(); x != b_end; ++x) {
    const KLPol& pol = kl.klPol(*x,y);
    h.emplace_back(*x,&pol); // add a |HeckeMonomial<KLPol>|
  }
}

/*****************************************************************************

        Chapter V -- Utilities.

  This section defines some utility functions for this module :

    - errorPol() : returns an error value;
    - one() : returns the K-L polynomial 1;
    - positivePart(q,d,m) : returns the positive part of q with u^d
      substituted and shifted by m;
    - zero() : returns the Laurent polynomial 0;

 *****************************************************************************/


const MuPol& errorMuPol()
{
  static MuPol p(klsupport::SKLCOEFF_MIN-1,MuPol::const_tag());
    /* cannot be a legal polynomial */
  return p;
}

const KLPol& errorPol()
{
  static KLPol p(klsupport::SKLCOEFF_MIN-1,KLPol::const_tag());
    /* cannot be a legal polynomial */
  return p;
}

const KLPol& one()
{
  static KLPol p(1,KLPol::const_tag());
  return p;
}


namespace {


/*
  This function is an auxiliary to |compute_mu|. It subtracts from |p|, which is
  destined to hold the positive degree part of a mu-polynomial, the positive
  degree part of $mp*r$, where $r=q[X:=u^d]*u^m$ (a Laurent polynomial in $u$)

  Forwards an error if there is overflow or underflow of the coefficients.

  NOTE : it is assumed that mp is non-zero!
*/
void muSubtraction(KLPol& p, const MuPol& mp, const KLPol& q,
		   const Ulong& d, const long& m)
{
  assert(not mp.isZero());
  MuPol r(d*q.deg()+m,m);

  for (Ulong j = 0; j <= q.deg(); ++j)
    r[static_cast<long>(d*j)+m] = q[j];

  for (long j = mp.val(); j <= mp.deg(); ++j)
    if (mp[j]!=0)
      for (long i = r.val(); i <= r.deg(); ++i)
	if (i+j >= 0) // ignore negative degree comtributions
        {
	  klsupport::SKLcoeff a = mp[j];
	  klsupport::safeMultiply(a,r[i]);
	  if (error::ERRNO)
	    return;
	  if (p.isZero() or i+j > static_cast<long>(p.deg())) // if required
	    p.setDeg(i+j); // zero-extend the coefficient vector
	  klsupport::safeAdd(p[i+j],-a);
	  if (error::ERRNO)
	    return;
	}

  p.reduceDeg();
}


/*
  Return the positive part (i.e. the part with positive degree) of the
  Laurent polynomial $q[X:=u^d]*u^m$ (in all calls |d==2|)
*/
KLPol positivePart(const KLPol& q, const Ulong& d, const long& m)
{
  KLPol p; // start out with zero

  /* compute degree of result */

  long h = q.deg()*d + m;

  if (h < 0)
    return p;

  p.setDeg(h); // resize
  p.setZero(h+1); // zero-fill

  for (polynomials::Degree j = q.deg()+1; h>=0 and j-->0; h-=d)
    // run index |h| through |p| decreasingly with steps |d|
    p[h] = q[j];
  return p;
}


// This function symmetrizes |p| and returns its address within the tree |t|


}; // |namespace|

}; // |namespace uneqkl|
