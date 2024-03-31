/*
  This is coxgroup.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "coxgroup.h"

#include "error.h"

/*************************************************************************

        Chapter I -- The CoxGroup class.

  a description of the class should go here ...

  This section defines the functions not already inlined or abstract :

    - CoxGroup(x,l) : constructs a CoxGroup of type x and rank l;
    - ~CoxGroup() : destructor;

  accessors :

    - coatoms(c,g) : puts the coatoms of g in c;
    - isDescent(g,s) : tells if s is a descent of g;
    - isDihedral(g) : tells if g is a dihedral element in the group;
    - modify(P,tok) : accessory to parse;
    - parse(P) : the parsing function (should perhaps become virtual);
    - parseBeginGroup(P) : parses the begin-group token in an expression;
    - parseEndGroup(P) : parses the end-group token in an expression;
    - parseGroupElement(P) : parses a group element;
    - parseModifier(P) : parses a modifier;
    - prod(x,s) : increments x by s;
    - prod(g,x) : increments g by x;
    - prod(x,g) : increments x by g;

  manipulators :

    - activateKL() : activates the K-L context;
    - activateIKL() : activates the inverse K-L context;
    - activateUEKL() : activates the unequal-parameter K-L context;
    - cBasis(h,y) : returns in h the data for the basis element c_y;
    - extendContext(g) : extends the active contexts to accomodate g;
    - fillKL() : fills the full K-L table;
    - fillMu() : fills the full mu-table;
    - klPol(x,y) : returns the ordinary K-L polynomial P_{x,y};
    - klRow(h,y) : returns the data for row y in kl-table;
    - mu(x,y) : returns the ordinary mu-coefficient mu(x,y);
    - permute(a) : permutes the context according to a;
    - setOutputStyle(C) : sets the output style;
    - sortContext() : sorts the context;
    - uneqcBasis(h,y) : returns in h the data for the unequal-parameter
      basis element c_y;

 *************************************************************************/

namespace coxgroup {

  class CoxGroup::CoxHelper {
  private:
    CoxGroup* d_W;
  public:
    void* operator new(size_t size) {return memory::arena().alloc(size);}
    void operator delete(void* ptr)
      {return memory::arena().free(ptr,sizeof(CoxHelper));}
    CoxHelper(CoxGroup* W);
    ~CoxHelper ();
    void sortContext();
    void checkInverses();
  };

};

namespace coxgroup {

/*
  Constructor for the abstract CoxGroup class. Does the basic initializations in
  the following order (the order is important said Fokko, but the code does not
  reflect what he wrote):

  - the Coxeter graph, which is where the Coxeter matrix gets constructed;
  - the interface;
  - the minroot table;
  - the Kazhdan-Lusztig context;
*/
CoxGroup::CoxGroup(const type::Type& x, const coxtypes::Rank& l)
  : d_graph(x,l)
  , d_mintable(d_graph)
  , d_klsupport(std::unique_ptr<schubert::SchubertContext>
		 (new schubert::SchubertContext(d_graph))
	       )
  , d_kl(nullptr)
  , d_invkl(nullptr)
  , d_uneqkl(nullptr)
  , d_interface(new interface::Interface(x,l))
  , d_outputTraits(d_graph,interface(),io::Pretty())
  , d_help(new CoxHelper(this))
{}


CoxGroup::~CoxGroup() {} // out-of-line, implicitly destructs |CoxHelper|

/******** accessors **********************************************************/

void CoxGroup::coatoms
  (list::List<coxtypes::CoxWord>& c, const coxtypes::CoxWord& g) const

/*
  This is a simple-minded "local" function that puts a list of reduced
  expressions for the coatoms of g in c.
*/

{
  c.setSize(0);

  for (Ulong j = 0; j < g.length(); ++j) {
    coxtypes::CoxWord h(0);
    Ulong i = 0;
    for (; i < j; ++i)
      h.append(g[i]);
    ++i;
    for (; i < g.length(); ++i) {
      coxtypes::Generator s = g[i]-1;
      int d = prod(h,s);
      if (d == -1)
	goto next;
    }
    // if we get here, h is a coatom
    c.append(h);
  next:
    continue;
  }

  return;
}


// Whether |s| is a descent of |g|.
bool CoxGroup::isDescent
  (const coxtypes::CoxWord& g, const coxtypes::Generator& s) const
{
  return (descent(g) & constants::eq_mask[s])!=0;
}

// Whether |g| is a dihedral element.
bool CoxGroup::isDihedral(const coxtypes::CoxWord& g) const
{
  if (g.length() < 3)
    return true;

  coxtypes::CoxLetter s = g[0];
  coxtypes::CoxLetter t = g[1];

  for (Ulong j = 2; j < g.length(); ++j) {
    if (j%2) { // g[j] should be t
      if (g[j] != t)
	return false;
    }
    else { // g[j] should be s
      if (g[j] != s)
	return false;
    }
  }

  return true;
}


/*
  Execute the modification indicated by |tok|, which is assumed to be of
  type-modifier_type. It is possible that further characters may have to
  be read from str.

  In the case of a general coxeter group, only two modifies are allowed :
  ! and ^
*/
void CoxGroup::modify
  (interface::ParseInterface& P, const interface::Token& tok) const
{
  if (interface::isInverse(tok)) {
    inverse(P.c);
  }

  if (interface::isPower(tok)) {
    Ulong m = readCoxNbr(P,ULONG_MAX);
    CoxGroup::power(P.c,m);
  }
}

void CoxGroup::parse(interface::ParseInterface& P) const

/*
  This function parses a group element from the line, starting
  at position r, and increments the coxtypes::CoxWord g with it. We have tried
  to include a number of convenient features, without overdoing it.

  The parser reads tokens from the string, skipping over any leading
  blank space at the beginning of each read (this means that tokens
  are not allowed to have leading blank spaces.) It always reads off
  the longest meaningful token. The following tokens are accepted :

    - a coxtypes::CoxWord prefix (postfix,separator);
    - a generator symbol;
    - a * (finite groups only) : represents the longest element;
    - a # : indicates that a number is to be read off;
    - a % : also indicates a number;
    - a ^ : indicates exponentiation;
    - a ! : indicates inversion;
    - a begin (end) group token;

  The symbols for prefix, postfix, separator, generators and grouping
  are provided by the interface, and already entered in the symbolTree.
  The grammar of expressions is as follows :

    expression -> expression elmt_expr
    elmt_expr -> (expression)[modifier] | group_elt[modifier]
    modifier -> empty | modifier ! | modifier * | modifier exp
    exp -> ^number
    group_elt -> coxword | #number | %number

  To parse an expression from a string, we try to parse a group element,
  followed by a modifier (= string of modifier symbols), and iterate
  this until we come up with an empty parse. Then we try to parse a
  begin-group token; if successful, we parse an expression, an end-group
  token, and a modifier.

  The parsing is managed by a stack of CoxWords; the topmost element serves
  to parse the current elementary expression, the next-to-topmost as an
  accumulator to hold the current string of elementary expressions.

  The return value is the number of characters read from the string.
*/

{
  for (;;) {
    if (parseGroupElement(P)) {
      if (error::ERRNO)
	return;
      continue;
    }
    if (parseBeginGroup(P)) { /* enter new nesting level */
      continue;
    }
    if (parseEndGroup(P)) { /* exit current nesting level */
      continue;
    }

    /* if we get to this point there is nothing more we can parse off */

    break;
  }

  if (P.nestlevel) { /* nesting error */
    error::ERRNO = error::PARSE_ERROR;
    return;
  }

  /* flush the current group element */

  prod(P.a[0],P.c);
  P.c.reset();

  return;

}

bool CoxGroup::parseBeginGroup(interface::ParseInterface& P) const

/*
  Tries to parse a begingroup token off P; in case of success, advances
  P.offset and changes the nesting level. This means simply that we
  increase the nestlevel, and reset it.
*/

{
  interface::Token tok = 0;
  const interface::Interface& I = interface();

  Ulong p = I.getToken(P,tok);

  if (p == 0)
    return false;

  if (not interface::isBeginGroup(tok))
    return false;

  P.nestlevel++;
  P.a.setSize(P.nestlevel+1);
  P.a[P.nestlevel].reset();
  P.offset += p;

  return true;
}

bool CoxGroup::parseContextNumber(interface::ParseInterface& P) const

/*
  Tries to parse a ContextNumber from P. This is a '%' character, followed
  by an integer which has to lie in the range [0,N[, where N is the current
  size of the enumerated part of the group.
*/

{
  const interface::Interface& I = interface();

  interface::Token tok = 0;
  Ulong p = I.getToken(P,tok);

  if (p == 0)
    return false;

  if (not interface::isContextNbr(tok))
    return false;

  // if we get to this point, we must read a valid integer

  P.offset += p;
  coxtypes::CoxNbr x = interface::readCoxNbr(P,contextSize());

  if (x == coxtypes::undef_coxnbr) { //error
    P.offset -= p;
    error::Error(error::CONTEXTNBR_OVERFLOW,contextSize());
    error::ERRNO = error::PARSE_ERROR;
  }
  else // x is valid
    prod(P.c,x);

  return true;
}

bool CoxGroup::parseEndGroup(interface::ParseInterface& P) const

/*
  Tries to parse an endgroup token; in case of success, reduces the nestlevel
  after doing the necessary bookkeeping.
*/

{
  interface::Token tok = 0;
  const interface::Interface& I = interface();

  Ulong p = I.getToken(P,tok);

  if (p == 0)
    return false;

  if (not interface::isEndGroup(tok))
    return false;

  if (P.nestlevel == 0) { /* error */
    error::ERRNO = error::PARSE_ERROR;
    return true;
  }

  // make the completed group the current group element

  P.c = P.a[P.nestlevel];
  P.nestlevel--;
  P.offset += p;

  // look for modifiers

  while (parseModifier(P)) {
    if (error::ERRNO)
      return true;
  }

  // flush modified group into accumulator

  prod(P.a[P.nestlevel],P.c);
  P.c.reset();

  return true;
}

bool CoxGroup::parseGroupElement(interface::ParseInterface& P) const

/*
  This function parses a group element from the string. A group element
  is one of (a) a coxword (b) a context number
  followed by a (possibly empty) string of modifiers. The modifiers
  are all treated as unary postfix operators (so that, for instance,
  g!^2* means (((g)!)^2)*).
*/

{
  Ulong r = P.offset;

  if (parseContextNumber(P)) { // the next token is a ContextNumber
    if (error::ERRNO) // parse error
      return true;
    else
      goto modify;
  }

  // if we get to this point, we have to read a coxtypes::CoxWord

  {
    interface().parseCoxWord(P,mintable());

    if (error::ERRNO) { // no coxtypes::CoxWord could be parsed
    if (P.offset == r) { // nothing was parsed
      error::ERRNO = 0;
      return false;
    }
    else // parse error
      return true;
    }
  }

 modify:

  // if we get to this point, a group element was successfully read

  while (parseModifier(P)) {
    if (error::ERRNO)
      return true;
  }

  // flush the current group element

  prod(P.a[P.nestlevel],P.c);
  P.c.reset();

  if (P.offset == r) // nothing was read; c is unchanged
    return false;
  else
    return true;
}

bool CoxGroup::parseModifier(interface::ParseInterface& P) const

/*
  This function parses a modifier from P.str at P.offset, and acts upon
  it accordingly : in case of success, it applies the modifier to P.c,
  and advances the offset.

  This is the default implementation, which doesn't allow the * modifier;
  this is accepted only for finite groups.
*/

{
  interface::Token tok = 0;
  const interface::Interface& I = interface();

  Ulong p = I.getToken(P,tok);

  if (p == 0)
    return false;

  if (not interface::isModifier(tok))
    return false;

  if (interface::isLongest(tok)) { /* error */
    error::ERRNO = error::PARSE_ERROR;
    return true;
  }

  P.offset += p;
  modify(P,tok);

  return true;
}

int CoxGroup::prod(coxtypes::CoxNbr& x, const coxtypes::Generator& s) const

/*
  This function increments x by right multiplication with s (i.e., it could
  have been written as x *= s). Returns +1 if the length goes up, -1 if the
  length goes down. Values of rank <= s < 2*rank correspond to right products.
*/

{
  coxtypes::CoxNbr x_old = x;
  x = schubert().shift(x,s);

  if (x > x_old)
    return 1;
  else
    return -1;
}

int CoxGroup::prod(coxtypes::CoxNbr& x, const coxtypes::CoxWord& g) const

/*
  Multiplies x consecutively by the terms in g. Stops at the first undefined
  operation. Returns the length increase.
*/

{
  int l = 0;

  for (Ulong j = 0; j < g.length(); ++j) {
    l += prod(x,g[j]-1);
    if (x == coxtypes::undef_coxnbr)
      break;
  }

  return l;
}

int CoxGroup::prod(coxtypes::CoxWord& g, const coxtypes::CoxNbr& d_x) const

/*
  Multiplies g by the terms in x. Returns the length increase.
*/

{
  int l = 0;
  coxtypes::CoxNbr x = d_x;

  while(x) {
    coxtypes::Generator s = constants::firstBit(ldescent(x));
    l += prod(g,s);
    prod(x,s+rank());
  }

  return l;
}

/******** manipulators ******************************************************/


/*
  This function activates the ordinary K-L context if it isn't already
  active.

  A memory error could happen in the process, but is not caught (i.e.
  CATCH_MEMORY_ERROR is not set); the program will simply exit printing
  the memory status.
*/
void CoxGroup::activateKL()
{
  if (d_kl == nullptr)
    d_kl.reset(new kl::KLContext(&d_klsupport));
}


/*
  This function activates the inverse K-L context if it isn't already
  active.
*/

void CoxGroup::activateIKL()
{
  if (d_invkl == nullptr)
    d_invkl.reset(new invkl::KLContext(&d_klsupport));
}


/*
  This function activates the unequal-parameter K-L context if it isn't already
  active.

  Forwards the error ABORT in case of failure (this means that there was a
  problem while getting the lengths from the user.)
*/
void CoxGroup::activateUEKL()
{
  if (d_uneqkl == nullptr) {
    d_uneqkl.reset(new uneqkl::KLContext(d_klsupport,graph(),interface()));
    if (error::ERRNO) {
      error::Error(error::ERRNO);
      d_uneqkl.reset();
    }
  }

  return;
}


/*
  Put in h the data of the full row for y in the K-L context corresponding
  to y, sorted in the order of the current normal forms. Activates the
  context if necessary.
*/
void CoxGroup::cBasis(kl::HeckeElt& h, const coxtypes::CoxNbr& y)
{
  activateKL();
  h = kl::cBasis(y,*d_kl);
  return;
}


/*
  This function extends the active contexts to acccomodate g. An active
  context is one that has a non-zero pointer.

  Currently there are three K-L contexts : ordinary, unequal parameter,
  and inverse. Parabolic stuff should be numbered independently.

  This is the place where we try to cope gently with a memory extension
  error; we wish to leave things as they were if the extension fails.

  Forwards the error ERROR_WARNING in case of failure.
*/
coxtypes::CoxNbr CoxGroup::extendContext(const coxtypes::CoxWord& g)
{
  Ulong prev_size = contextSize();
  coxtypes::CoxNbr x = d_klsupport.extendContext(g);

  if (error::ERRNO) {
    goto revert;
  }

  if (d_kl) {
    d_kl->setSize(contextSize());
  if (error::ERRNO)
    goto revert;
  }
  if (d_uneqkl) {
    d_uneqkl->setSize(contextSize());
  if (error::ERRNO)
    goto revert;
  }
  if (d_invkl) {
    d_invkl->setSize(contextSize());
  if (error::ERRNO)
    goto revert;
  }

  return x;

 revert:
  d_klsupport.revertSize(prev_size);
  if (d_kl)
    d_kl->revertSize(prev_size);
  if (d_uneqkl)
    d_uneqkl->revertSize(prev_size);
  if (d_invkl)
    d_invkl->revertSize(prev_size);
  error::ERRNO = error::ERROR_WARNING;
  return coxtypes::undef_coxnbr;
} // |extendContext|

coxtypes::CoxNbr CoxGroup::extend_context(const coxtypes::Cox_word& g)
{
  Ulong prev_size = contextSize();
  coxtypes::CoxNbr x = d_klsupport.extend_context(g);

  if (error::ERRNO) {
    goto revert;
  }

  if (d_kl) {
    d_kl->setSize(contextSize());
  if (error::ERRNO)
    goto revert;
  }
  if (d_uneqkl) {
    d_uneqkl->setSize(contextSize());
  if (error::ERRNO)
    goto revert;
  }
  if (d_invkl) {
    d_invkl->setSize(contextSize());
  if (error::ERRNO)
    goto revert;
  }

  return x;

 revert:
  d_klsupport.revertSize(prev_size);
  if (d_kl)
    d_kl->revertSize(prev_size);
  if (d_uneqkl)
    d_uneqkl->revertSize(prev_size);
  if (d_invkl)
    d_invkl->revertSize(prev_size);
  error::ERRNO = error::ERROR_WARNING;
  return coxtypes::undef_coxnbr;
}



/*
  This function fills the whole inverse K-L table up to the size of the current
  schubert context, for inverse kazhdan-lusztig polynomials, avtrer having
  activated the context if necessary.
*/
void CoxGroup::fillIKL()
{
  activateIKL();

  d_invkl->fillKL();
  return;
}


/*
  This function fills the whole inverse mu-table up to the size of the current
  schubert context, for the inverse kazhdan-lusztig polynomials, after having
  activated the context if necessary.
*/
void CoxGroup::fillIMu()
{
  activateIKL();

  d_invkl->fillMu();
  return;
}


/*
  This function fills the whole K-L table up to the size of the current
  schubert context, after having activated the context if necessary.
*/
void CoxGroup::fillKL()
{
  activateKL();

  d_kl->fillKL();
  return;
}

void CoxGroup::fillMu()

/*
  This function fills the whole mu-table up to the size of the current
  schubert context, after having activated the context if necessary.
*/

{
  activateKL();

  d_kl->fillMu();
  return;
}

void CoxGroup::fillUEKL()

/*
  This function fills the whole unequal-parameter K-L table up to the size
  of the current schubert context, after having activated the context if
  necessary.
*/

{
  activateUEKL();

  d_uneqkl->fillKL();
  return;
}

void CoxGroup::fillUEMu()

/*
  This function fills the whole unequal-parameter mu-tables up to the size of
  the current schubert context, after having activated the context if
  necessary.
*/

{
  activateUEKL();

  d_uneqkl->fillMu();
  return;
}

const invkl::KLPol& CoxGroup::invklPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y)

/*
  Returns the inverse K-L polynomial Q_{x,y}, after activating the context
  if necessary.
*/

{
  activateIKL();

  return d_invkl->klPol(x,y);
}

void CoxGroup::invklRow(invkl::HeckeElt& h, const coxtypes::CoxNbr& y)

/*
  Puts in h the data of the full row for y in the inverse K-L context
  corresponding to y, sorted in the order of the current normal forms.
  Activates the context if necessary.
*/

{
  activateIKL();
  d_invkl->row(h,y);
  return;
}

const kl::KLPol& CoxGroup::klPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y)

/*
  Returns the ordinary K-L polynomial P_{x,y}, after activating the context
  if necessary.
*/

{
  activateKL();

  return d_kl->klPol(x,y);
}

void CoxGroup::klRow(kl::HeckeElt& h, const coxtypes::CoxNbr& y)

/*
  Puts in h the data of the full row for y in the K-L context corresponding
  to y, sorted in the order of the current normal forms. Activates the
  context if necessary.
*/

{
  activateKL();
  d_kl->row(h,y);
  return;
}

klsupport::KLCoeff CoxGroup::mu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y)

/*
  Returns the ordinary mu-coefficent mu(x,y), after activating the context
  if necessary.
*/

{
  activateKL();

  return d_kl->mu(x,y);
}


/*
  This function permutes all the active contexts w.r.t. the permutation a.
  The idea is that we _renumber_ the context; a is the permutation of the
  range [0,size[ which gives for each number x the new number of the same
  group element. So in terms of group elements the correspondence is that
  new[a(x)] = old[x].

  Applying the permutation to group-valued functions is trivial : just permute
  the values using a. Applying it to functions with range in the group, we
  need to compose on the right with the inverse permutation : new_f(x) =
  old_f(a^{-1}(x)). For functions from group to group, we need to do both.

  What we have to do here is :

    - permute the schubert context itself;
    - permute inverse;
    - permute last;
    - permute involution;
    - permute the extrList;
    - permute the various kl-contexts;

  Note that extrList should be seen as a table of enumerated subsets of the
  group. The various kllists are tables of sequences of polynomials,
  enumerated in accordance with the extrList. The mulists are also tables
  of enumerated subsets of the group, together with additional data. A
  further requirement is that these enumerations be increasing; so we
  will furthermore have to sort each extrrow, and the rows that should
  be compatible with it, and each mu-row.

  Finally, we have imposed on ourselves the burden of writing only one
  of the pairs (y,y_inverse), viz. the one with the smaller index. So we
  need to maintain that requirement as well.
*/
void CoxGroup::permute(const bits::Permutation& a)
{
  d_klsupport.permute(a);

  if (d_kl)
    d_kl->permute(a);
  if (d_invkl)
    d_invkl->permute(a);
  if (d_uneqkl)
    d_uneqkl->permute(a);

  d_help->checkInverses();
  d_help->sortContext();

  return;
}


/*
  Put into h the data of the full row for y in the K-L context corresponding
  to y, sorted in the order of the current normal forms. Activates the
  context if necessary.
*/
void CoxGroup::uneqcBasis(uneqkl::HeckeElt& h, const coxtypes::CoxNbr& y)
{
  activateUEKL();
  h = uneqkl::cBasis(y,*d_uneqkl);
  return;
}

const uneqkl::KLPol& CoxGroup::uneqklPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y)

/*
  Returns the unequal-parameter K-L polynomial P_{x,y}, after activating the
  context if necessary.
*/

{
  activateUEKL();

  return d_uneqkl->klPol(x,y);
}


/*
  Returns the unequal-parameter mu-polynomial mu_{s,x,y}, after activating
  the context if necessary.
*/
const uneqkl::MuPol CoxGroup::uneqmu
  (const coxtypes::Generator& s,
   const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y)
{
  activateUEKL();

  return d_uneqkl->mu(s,x,y);
}

void CoxGroup::uneqklRow(uneqkl::HeckeElt& h, const coxtypes::CoxNbr& y)

/*
  Puts in e_row and kl_row the data of the full row of the unequal-parameter
  K-L context corresponding to y, sorted in the short-lex order of the
  current normal forms. Activates the context if necessary.
*/

{
  activateUEKL();

  d_uneqkl->row(h,y);
  return;
}

};

/*****************************************************************************

        Chapter II -- The CoxHelper class.

  This class provides some private helper functions for managing a CoxGroup
  structure, that we don't want in the coxgroup.h file. The following
  functions are provided :

   - CoxHelper(W);
   _ ~CoxHelper();

 *****************************************************************************/

namespace coxgroup {

CoxGroup::CoxHelper::CoxHelper(CoxGroup* W):d_W(W)

{}

CoxGroup::CoxHelper::~CoxHelper()

{}

void CoxGroup::CoxHelper::sortContext()

/*
  This function is an auxiliary to permute; it takes care of putting
  things in increasin context number order after the permutation.
*/

{
  klsupport::KLSupport* kls = &d_W->d_klsupport;

  for (coxtypes::CoxNbr y = 0; y < d_W->contextSize(); ++y) {
    if (!kls->isExtrAllocated(y))
      continue;

    bits::Permutation a(0);
    bits::sortI(d_W->extrList(y),a);

    kls->applyIPermutation(y,a);

    /* apply to the various klcontexts */

    if (d_W->d_kl) {
      d_W->d_kl->applyIPermutation(y,a);
    }
    if (d_W->d_invkl) {
      d_W->d_invkl->applyIPermutation(y,a);
    }
    if (d_W->d_uneqkl) {
      d_W->d_uneqkl->applyIPermutation(y,a);
    }
  }

  return;
}


/*
  This function is an auxiliary to permute; it takes care of checking
  that the smaller one of the pairs (y,y_inverse) is allocated in the
  extrList and in all lists which depend on that.
*/
void CoxGroup::CoxHelper::checkInverses()
{
  klsupport::KLSupport& kls = d_W->d_klsupport;

  for (coxtypes::CoxNbr y = 0; y < d_W->contextSize(); ++y) {
    coxtypes::CoxNbr yi = d_W->inverse(y);
    if (yi <= y)
      continue;
    if (kls.isExtrAllocated(y))
      continue;
    /* if we get here, we should transfer lists from yi to y */
    kls.move_extr_list_from_inverse(y);
    if (d_W->d_kl)
      d_W->d_kl->applyInverse(y);
    if (d_W->d_invkl)
      d_W->d_invkl->applyInverse(y);
    if (d_W->d_uneqkl)
      d_W->d_uneqkl->applyInverse(y);
  }

  return;
}

};
