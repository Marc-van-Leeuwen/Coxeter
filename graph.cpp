/*
  This is graph.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "graph.h"

#include "sl_list.h"
#include "constants.h"
#include "directories.h"
#include "interactive.h"

namespace {
  using namespace graph;

  coxtypes::CoxSize dihedralOrder(CoxGraph& G, bits::Lflags I);
  coxtypes::ParSize extrQuotOrder(CoxGraph& G, bits::Lflags I, coxtypes::Generator s);
  void fillCoxAMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxBMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxDMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxEMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxFMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxGMatrix(CoxMatrix& m);
  void fillCoxHMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxIMatrix(CoxMatrix& m);
  void fillCoxaMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxbMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxcMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxdMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxeMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxfMatrix(CoxMatrix& m, coxtypes::Rank l);
  void fillCoxgMatrix(CoxMatrix& m);
  void fillCoxXMatrix(CoxMatrix& m, const coxtypes::Rank& l, const type::Type& t);
  void fillCoxYMatrix(CoxMatrix& m, coxtypes::Rank l);
  coxtypes::CoxSize finiteOrder(const type::Type& type, const coxtypes::Rank& rank);
  Ulong gcd(Ulong a, Ulong b);
  const type::Type& irrType(CoxGraph& G, bits::Lflags I);
  coxtypes::Generator lastGenerator(CoxGraph& G, bits::Lflags I);
  coxtypes::ParSize lastQuotOrder(const type::Type& type, coxtypes::Rank rank);
  CoxMatrix Coxeter_matrix(const type::Type& x, const coxtypes::Rank& l);
  containers::vector<bits::Lflags> vertex_stars
    (const CoxMatrix& m, const coxtypes::Rank& l);
  containers::vector<bits::Lflags> finite_edge_list
    (const CoxMatrix& m, const coxtypes::Rank& l);
  CoxEntry maxCoefficient(CoxGraph& G, bits::Lflags I);
  CoxEntry minCoefficient(CoxGraph& G, bits::Lflags I);
  coxtypes::CoxSize A_order(coxtypes::Rank rank);
  coxtypes::CoxSize B_order(coxtypes::Rank rank);
  coxtypes::CoxSize D_order(coxtypes::Rank rank);
};

/****************************************************************************

        Chapter I -- The CoxGraph class.

  The CoxGraph class provides access to the various data contained in the
  Coxeter matrix : the matrix itself, and the underlying graph structure.

  The following functions are provided :

   - CoxGraph(x,l) : constructs a CoxGraph of type x and rank l;
   - ~CoxGraph() : not implemented yet;

   - component(I,s) : returns the connected component of s in I;
   - extremities(I) : returns the extremities of I;
   - nodes(I) : returns the nodes of I;

 ****************************************************************************/

namespace graph {


// Construct a Coxeter graph of type x and rank l.
CoxGraph::CoxGraph(const type::Type& x, const coxtypes::Rank& l)
  : d_type(x)
  , d_rank(l)
  , d_matrix(Coxeter_matrix(x,d_rank))
  , d_S()
  , d_star()
{
  if (error::ERRNO)
    return;

  /* the restriction on the rank should be removed eventually */

  if (l <= coxtypes::MEDRANK_MAX) // required to have |bits::Lflags| large enough
  {
    d_S = constants::leqmask[d_rank-1]; // all |d_rank| initial bits set
    d_star = vertex_stars(d_matrix,d_rank);
  }

  d_finite_edges = finite_edge_list(d_matrix,l);
}


// Automatic destruction does everything required
CoxGraph::~CoxGraph() {}


// Return the connected component of |s| in |I| as a |bits::Lflags|
bits::Lflags CoxGraph::component(bits::Lflags I, coxtypes::Generator s) const
{
  bits::Lflags nf = constants::lmask[s];
  bits::Lflags f = 0;

  while (nf)  /* there are new elements to be considered */
    {
      f |= nf;
      for (bits::Lflags f1 = nf; f1; f1 &= f1-1)
	nf |= (I & d_star[constants::firstBit(f1)]);
      nf &= ~f;
    }

  return f;
}


bits::Lflags CoxGraph::extremities(bits::Lflags I) const

/*
  This function returns a bitmap of the set of points in I which are extremal
  in the induced graph; this means that the valency of the point is one.
*/

{
  bits::Lflags f = 0;
  bits::Lflags f1 = I;

  while (f1)
    {
      coxtypes::Generator s = constants::firstBit(f1);
      if (bits::bitCount(d_star[s]&I) == 1)  /* s is an extremity */
	f |= constants::lmask[s];
      f1 &= f1-1;
    }

  return f;
}


bits::Lflags CoxGraph::nodes(bits::Lflags I) const

/*
  This function returns a bitmap of the set of points in I which are nodes
  for the induced graph; this means that the valency of the point is > 2.
*/

{
  bits::Lflags f,f1;
  coxtypes::Generator s;

  f = 0;
  f1 = I;

  while (f1)
    {
      s = constants::firstBit(f1);
      if (bits::bitCount(d_star[s]&I) > 2)  /* s is a node */
	f |= constants::lmask[s];
      f1 &= f1-1;
    }

  return f;
}

};


/****************************************************************************

        Chapter II -- Coxeter matrix functions.

  This section provides the functions which will fill in the Coxeter matrices
  in the predefined types, and those reading a matrix from a file, or reading
  it in interactively. These functions are private to graph.c.

  The following functions are provided :

  for filling the matrix of a finite Weyl group :

   - fillCoxAMatrix(m,l) : type A;
   - fillCoxBMatrix(m,l) : type B = C;
   - fillCoxDMatrix(m,l) : type D;
   - fillCoxEMatrix(m,l) : type E;
   - fillCoxFMatrix(m,l) : type F;
   - fillCoxGMatrix(m,l) : type G;
   - fillCoxHMatrix(m,l) : type H;
   - fillCoxIMatrix(m,l) : type I (the remaining dihedrals) (interactive);

   for filling the matrix of an affine Weyl group :

   - fillCoxaMatrix(m,l) : type a;
   - fillCoxbMatrix(m,l) : type b;
   - fillCoxcMatrix(m,l) : type c;
   - fillCoxdMatrix(m,l) : type d;
   - fillCoxeMatrix(m,l) : type e;
   - fillCoxfMatrix(m,l) : type f;
   - fillCoxgMatrix(m,l) : type g;

   for reading a Coxeter matrix from a file :

   - fillCoxXMatrix(m,l);

   for reading in a Coxeter matrix interactively :

   - fillCoxYMatrix(m,l);

  The functions filling in the actual data in the CoxGraph structure are :

   - m = Coxeter_matrix(x,l) : computes the Coxeter matrix;
   - star = vertex_stars(m,l) : computes the star array;
   - ops = finite_edge_list(m,l) : enumerate the finite edges

 ****************************************************************************/

namespace {

void fillCoxAMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  for (coxtypes::Rank j = 1; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxBMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  m[1] = 4;
  m[l] = 4;

  for (coxtypes::Rank j = 2; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxDMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  if (l == 2)
    return;

  m[2] = 3;
  m[l + 2] = 3;
  m[2*l] = 3;
  m[2*l + 1] = 3;

  for (coxtypes::Rank j = 3; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxEMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  m[2] = 3;
  m[2*l] = 3;

  if (l == 3)
    return;

  m[l + 3] = 3;
  m[2*l + 3] = 3;
  m[3*l + 1] = 3;
  m[3*l + 2] = 3;

  for (coxtypes::Rank j = 4; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxFMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  for (coxtypes::Rank j = 1; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[l + 2] = 4;
  m[2*l + 1] = 4;

  return;
}


void fillCoxGMatrix(CoxMatrix& m)

{
  m[1] = 6;
  m[2] = 6;

  return;
}

void fillCoxHMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  m[1] = 5;
  m[l] = 5;

  for (coxtypes::Rank j = 2; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxIMatrix(CoxMatrix& m)

{
  CoxEntry m_12;

  m_12 = interactive::getCoxEntry(1,2);

  if (error::ERRNO)
    return;

  m[1] = m_12;
  m[2] = m_12;

  return;
}


void fillCoxaMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  if (l == 2)
    {
      m[1] = 0;
      m[2] = 0;
      return;
    }

  for (coxtypes::Rank j = 1; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[l-1] = 3;
  m[(l-1)*l] = 3;

  return;
}

void fillCoxbMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  if (l == 3) { // go over to type c3
    fillCoxcMatrix(m,3);
    return;
  }

  m[1] = 4;
  m[l] = 4;

  for (coxtypes::Rank j = 2; j < l-1; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[(l-3)*l + (l-1)] = 3;
  m[(l-1)*l + (l-3)] = 3;

  return;
}

void fillCoxcMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  m[1] = 4;
  m[l] = 4;

  for (coxtypes::Rank j = 2; j < l-1; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[(l-2)*l + (l-1)] = 4;
  m[(l-1)*l + (l-2)] = 4;

  return;
}

void fillCoxdMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  m[2] = 3;
  m[2*l] = 3;

  for (coxtypes::Rank j = 2; j < l-1; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[(l-3)*l + (l-1)] = 3;
  m[(l-1)*l + (l-3)] = 3;

  return;
}

void fillCoxeMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  m[2] = 3;
  m[2*l] = 3;

  m[l + 3] = 3;
  m[2*l + 3] = 3;
  m[3*l + 1] = 3;
  m[3*l + 2] = 3;

  for (coxtypes::Rank j = 4; j < l-1; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  switch (l)
    {
    case 5:
      m[l-1] = 3;
      m[l + (l-1)] = 3;
      m[(l-1)*l] = 3;
      m[(l-1)*l + 1] = 3;
      break;
    case 6:
      m[2*l + (l-1)] = 3;
      m[(l-1)*l + 2] = 3;
      break;
    case 7:
      m[l + (l-1)] = 3;
      m[(l-1)*l + 1] = 3;
      break;
    case 8:
      m[l-1] = 3;
      m[(l-1)*l] = 3;
      break;
    case 9:
      m[(l-2)*l + (l-1)] = 3;
      m[(l-1)*l + (l-2)] = 3;
      break;
    }

  return;
}

void fillCoxfMatrix(CoxMatrix& m, coxtypes::Rank l)

{
  for (coxtypes::Rank j = 1; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[l + 2] = 4;
  m[2*l + 1] = 4;

  return;
}


void fillCoxgMatrix(CoxMatrix& m)

{
  m[1] = 6;
  m[3] = 6;
  m[5] = 3;
  m[7] = 3;

  return;
}

/*
  Recall that in type X the type is really a string, where the name
  of a valid input file follows X, lying under coxeter_matrices. The
  program reads the first l entries of the first l lines in the file;
  this allows us sometimes to define a whole family of groups with
  a single file.

  NOTE : this should perhaps be completed by a recognition test, in order
  to renumber the generators if necessary (and then modify the matrix
  accordingly). For the time being, we leave it as is.
*/
void fillCoxXMatrix(CoxMatrix& m, const coxtypes::Rank& l, const type::Type& t)
{
  static std::string buf;
  const std::string& name = t.name();

  buf = directories::coxmatrix_dir + "/" + name;
  FILE *inputfile = fopen(buf.c_str(),"r");

  for (coxtypes::Rank i = 0; i < l; i++) {
    for (coxtypes::Rank j = 0; j < l; j++) {

      /* check for EOL */

      if (interactive::endOfLine(inputfile)) {
	Error(error::BAD_LINE,name.c_str()+1,l,i,j);
	error::ERRNO = error::ABORT;
	return;
      }

      m[i*l + j] = interactive::readCoxEntry(i,j,inputfile);

      if (error::ERRNO) {
	error::Error(error::ERRNO,i,j);
	error::ERRNO = error::ABORT;
	return;
      }

      /* check for symmetry */

      if (j < i)
	if (m[i*l + j] != m[j*l + i]) {
	  error::Error(error::NOT_SYMMETRIC,name.c_str()+1,&m,l,i,j);
	  error::ERRNO = error::ABORT;
	  return;
	}
    }

    /* flush remaining line of the inputfile */

    char c;

    while((c = getc(inputfile)) != EOF)
      if (c == '\n')
	break;
  }

  fclose(inputfile);
}

void fillCoxYMatrix(CoxMatrix& m, coxtypes::Rank l)

/*
  This is the type for arbitrary input, where the coxeter matrix is gotten
  interactively from the user. He is prompted only for entries i,j for
  which i < j.

  NOTE : this should be completed by a recognition test, in order
  to renumber the generators if necessary (and then modify the matrix
  accordingly). For the time being, we leave it as is.
*/

{
  for (coxtypes::Rank i = 0; i < l; i++)
    for (coxtypes::Rank j = i+1; j < l; j++) {
      m[i*l + j] = interactive::getCoxEntry(i+1,j+1);
      if (error::ERRNO) {
	error::Error(error::ERRNO);
	error::ERRNO = error::ERROR_WARNING;
	return;
      }
      m[j*l + i] = m[i*l + j];
    }

  return;
}



/*
  Return the Coxeter matrix for |x| and |l|. In the case of type X,
  this includes checking the input file for correct values; in the
  case of type I, we need to get to get the only non-trivial entry
  of the matrix from the user.

  In the case of failure, it sets an error, and returns.
*/
CoxMatrix Coxeter_matrix(const type::Type& x, const coxtypes::Rank& l)
{
  CoxMatrix result(l*l,2); // generators commute "by default"

  // set diagonal to entries 1
  for (Ulong j = 0; j < l; ++j)
    result[j*(l+1)] = 1;

  switch (x[0])
    {
    case 'A':
      fillCoxAMatrix(result,l);
      break;
    case 'B':
      fillCoxBMatrix(result,l);
      break;
    case 'D':
      fillCoxDMatrix(result,l);
      break;
    case 'E':
      fillCoxEMatrix(result,l);
      break;
    case 'F':
      fillCoxFMatrix(result,l);
      break;
    case 'G':
      fillCoxGMatrix(result);
      break;
    case 'H':
      fillCoxHMatrix(result,l);
      break;
    case 'I':
      fillCoxIMatrix(result);
      if (error::ERRNO)
	return result; // must return something even on error
      break;
    case 'a':
      fillCoxaMatrix(result,l);
      break;
    case 'b':
      fillCoxbMatrix(result,l);
      break;
    case 'c':
      fillCoxcMatrix(result,l);
      break;
    case 'd':
      fillCoxdMatrix(result,l);
      break;
    case 'e':
      fillCoxeMatrix(result,l);
      break;
    case 'f':
      fillCoxfMatrix(result,l);
      break;
    case 'g':
      fillCoxgMatrix(result);
      break;
    case 'X':
      fillCoxXMatrix(result,l,x);
      break;
    case 'Y':
      fillCoxYMatrix(result,l);
      break;
    }

  return result;
}



/*
  Makes the star-array of the Coxeter graph. This is an array of l
  bits::Lflags, flagging the "stars" of each generator in the Coxeter diagram.
*/
containers::vector<bits::Lflags> vertex_stars
  (const CoxMatrix& m, const coxtypes::Rank& l)
{
  containers::vector<bits::Lflags> result;

  for(coxtypes::Generator s = 0; s < l; s++) {
    bits::Lflags star = 0;
    for (coxtypes::Generator t = 0; t < l; t++)
      if ((m[s*l + t] > 2) || (m[s*l + t] == 0))
	star |= constants::lmask[t];
    result.push_back(star);
  }

  return result;
}


/*
  Returns a list of the finite edges of the graph, computed from the matrix. For
  each finite edge in the graph (in some order) this list holds 2-element subset
  that definies the edge, represented as |bits::Lflags|
*/
containers::vector<bits::Lflags> finite_edge_list
  (const CoxMatrix& m, const coxtypes::Rank& l)
{
  /* count number of finite edges */

  containers::sl_list<bits::Lflags> result;

  for (coxtypes::Generator s = 0; s < l; ++s)
    for (coxtypes::Generator t = s+1; t < l; ++t)
      if ((m[s*l + t] > 2) && (m[s*l + t] != infty)) // finite edge condition
	result.push_back(constants::lmask[s] | constants::lmask[t]);

  return result.to_vector();
}

};

/****************************************************************************

      Chapter III -- Graph analysis.

  This section defines various functions for the analysis of subsets of the
  Coxeter graph. The following functions are provided :

   - isAffine(G,I);
   - isConnected(G,I);
   - isCrystallographic(G,I);
   - isFinite(G,I);
   - isLoop(G,I);
   - isSimplyLaced(G,I);
   - isTree(G,I);

  In fact, the real work for isFinite and isAffine is provided in :

   - irrType(G,I) : returns the type of an irreducible subset;

 ****************************************************************************/

namespace graph {

bool isAffine(CoxGraph& G, bits::Lflags I)

/*
  Returns true if the group generated by I is affine, false otherwise. Uses the
  classification of the graphs of affine Coxeter groups.

  It is assumed that I is already irreducible.
*/

{
  const type::Type& type = irrType(G,I);

  if (strchr("abcdefg",type[0]))  /* group is affine */
    return true;
  else
    return false;
}


bool isConnected(CoxGraph& G, bits::Lflags I)

/*
  Returns true if the graph induced on I is connected, false otherwise.
*/

{
  if (I == 0)
    return false;

  coxtypes::Generator s = constants::firstBit(I);

  if (G.component(I,s) == I)
    return true;
  else
    return false;
}


bool isCrystallographic(CoxGraph& G, bits::Lflags I)

/*
  Checks if the restriction of the Coxeter graph to I is crystallographic,
  i.e., if the entries in the Coxeter matrix are all 2,3,4,6 or infinity.
*/

{
  for (coxtypes::Generator s = 0; s < G.rank(); s++)
    for (coxtypes::Generator t = s+1; t < G.rank(); t++)
      {
	switch (G.M(s,t)) {
	case 0:
	case 2:
	case 3:
	case 4:
	case 6:
	  continue;
	default:
	  return false;
	};
      }

  return true;
}


bool isFinite(CoxGraph& G, bits::Lflags I)

/*
  Returns true if the group generated by I is finite, false otherwise. Uses the
  classification of the graphs of finite Coxeter groups.
*/

{
  while (I)
    {
      coxtypes::Generator s = constants::firstBit(I);
      bits::Lflags f = G.component(I,s);
      const type::Type& type = irrType(G,f);
      if (strchr("ABCDEFGHI",type[0]) == NULL)
	return false;
      I &= ~f;
    }

  return true;
}


bool isLoop(CoxGraph& G, bits::Lflags I)

/*
  Returns 1 if the graph induced on I is a loop, 0 otherwise. Uses the
  characterization that a loop is a connected graph for which all
  valencies are equal to 2.
*/

{
  if (!isConnected(G,I))
    return false;

  for (bits::Lflags f = I; f; f &= f-1)
    {
      coxtypes::Generator s = constants::firstBit(f);
      if (bits::bitCount(G.star(I,s)) != 2)
	return false;
    }

  return true;
}


/*
  Return whether the Coxeter graph restricted to I is simply laced (i.e., all
  edges have label 3).
*/
bool isSimplyLaced(CoxGraph& G, bits::Lflags I)
{
  for (bits::Lflags fs = I; fs; fs &= fs-1)
    {
      coxtypes::Generator s = constants::firstBit(fs);
      for (bits::Lflags ft = fs & (fs-1); ft; ft &= ft-1)
	{
	  coxtypes::Generator t = constants::firstBit(ft);
	  if ((G.M(s,t) == 0) || (G.M(s,t) > 3))
	    return false;
	}
    }

  return true;
}

bool isTree(CoxGraph& G, bits::Lflags I)

/*
  Returns 1 if the graph induced on I is a tree, 0 otherwise. Uses the
  characterization that a tree is a connected graph for which the
  number of edges is the number of vertices minus one.
*/

{
  if (!isConnected(G,I))
    return false;

  unsigned edgecount = 0;

  for (bits::Lflags f = I; f; f &= f-1)
    {
      coxtypes::Generator s = constants::firstBit(f);
      edgecount += bits::bitCount(G.star(I,s));
    }

  edgecount /= 2;  /* each edge was counted twice */

  if (edgecount == (bits::bitCount(I) - 1))
    return true;
  else
    return false;
}

};

namespace {

const type::Type& irrType(CoxGraph& G, bits::Lflags I)

/*
  Returns the type of the subgraph induced on I, if this subgraph is
  irreducible, finite or affine. Assumes that irreducibility has already
  been checked. Returns type "X" if the type is not defined.

  The result is returned as type, which is a safe place until the next call
  to irrType.
*/

{
  static type::Type type("X");

  if (bits::bitCount(I) == 1) {
    type[0] = 'A';
    return type;
  }

  if (bits::bitCount(I) == 2)
    {
      coxtypes::Generator s = constants::firstBit(I);
      coxtypes::Generator t = constants::firstBit(I & (I-1));
      CoxEntry m = G.M(s,t);

      switch (m)
	{
	case 0:
	  type[0] = 'a';
	  return type;
	case 3:
	  type[0] = 'A';
	  return type;
	case 4:
	  type[0] = 'B';
	  return type;
	case 5:
	  type[0] = 'H';
	  return type;
	case 6:
	  type[0] = 'G';
	  return type;
	default:
	  type[0] = 'I';
	  return type;
	};
    }

  /* from here on the rank is at least three */

  if (!isTree(G,I))  /* type must be a_n */
    {
      if (!isLoop(G,I))  /* unknown type */
	return type;
      if (!isSimplyLaced(G,I))  /* unknown type */
	return type;
      type[0] = 'a';
      return type;
    }

  /* from here on the graph is a tree */

  CoxEntry m = maxCoefficient(G,I);

  switch (m) {
  case 3: { /* simply laced : type is A, D, E, d, e if known */
    bits::Lflags fn = G.nodes(I);
    switch (bits::bitCount(fn))
      {
      case 0: /* type A */
	type[0] = 'A';
	return type;
      case 1: { /* type is D, E, e, or d5, if known */
	coxtypes::Generator n = constants::firstBit(G.nodes(I));
	switch (bits::bitCount(G.star(n)))
	  {
	  case 3: { /* type is D, E or e */
	    bits::Lflags f = G.extremities(I);
	    switch (bits::bitCount(f & G.star(n)))  /* short branches */
	      {
	      case 3:  /* type is D4 */
		type[0] = 'D';
		return type;
	      case 2:  /*  type is Dn, n >= 5 */
		type[0] = 'D';
		return type;
	      case 1: { /* type is E6, E7, E8, e8 or e9 */
		/* trim branches by one */
		bits::Lflags J = I & ~f;
		f = G.extremities(J);
		switch (bits::bitCount(f & G.star(n)))
		  {
		  case 0:  /* two branches of length > 2 */
		    if (bits::bitCount(I) == 8)  /* type e8 */
		      type[0] = 'e';
		    return type;
		  case 1:  /* one branch of length 2 */
		    switch (bits::bitCount(I))
		      {
		      case 7:  /* type E7 */
		      case 8:  /* type E8 */
			type[0] = 'E';
			return type;
		      case 9:  /* type e9 */
			type[0] = 'e';
			return type;
		      default:  /* unknown type */
			return type;
		      };
		  case 2:  /* two branches of length 2 */
		    if (bits::bitCount(I) == 6)  /* type E6 */
		      type[0] = 'E';
		    return type;
		  };
	      }
	      case 0:  /* type has to be e7 */
		if (bits::bitCount(I) == 7)
		  type[0] = 'e';
		return type;
	      };
	  }
	  case 4:  /* type is d5 */
	    if (bits::bitCount(I) == 5)
	      type[0] = 'd';
	    return type;
	  default:  /* unknown type */
	    return type;
	  };
      }
      case 2: {
	bits::Lflags f = G.extremities(I);
	if (bits::bitCount(f) > 4)  /* unknown type */
	  return type;
	/* from here on each node has three branches */
	bits::Lflags J = I & ~f;
	f = G.extremities(J);
	if (f == fn)  /* type d */
	  type[0] = 'd';
	return type;
      }
      default:  /* unknown type */
	return type;
      };
  }
  case 4: { /* type is B, F, b, c or f if known */
    switch (bits::bitCount(G.nodes(I)))
      {
      case 0: { /* graph is a string : type is B, F, c or f */
	bits::Lflags f = G.extremities(I);
	bits::Lflags J = I & ~f;
	switch (maxCoefficient(G,J))
	  {
	  case 1:
	  case 3: { /* type is B or c */
	    type[0] = 'B';
	    coxtypes::Generator s = constants::firstBit(f);
	    coxtypes::Generator t = constants::firstBit(G.star(s));
	    CoxEntry m1 = G.M(s,t);
	    if (m1 == 3)
	      return type;
	    f &= f-1;
	    s = constants::firstBit(f);
	    t = constants::firstBit(G.star(s));
	    m1 = G.M(s,t);
	    if (m1 == 4)
	      type[0] = 'c';
	    return type;
	  }
	  case 4:  /* type is F or f, if known */
	    switch (bits::bitCount(I))
	      {
	      case 4:  /* type F4 */
		type[0] = 'F';
		return type;
	      case 5: {
		CoxEntry m1 = minCoefficient(G,J);
		if (m1 == 3)  /* type f5 */
		  type[0] = 'f';
		return type;
	      }
	      default:  /* unknown type */
		return type;
	      };
	  default: /* unknown type */
	    return type;
	  };
      }
      case 1: { /* type is b if known */
	bits::Lflags f = G.extremities(I);
	if (bits::bitCount(f) > 3)  /* more than three branches */
	  return type;
	bits::Lflags J = I & ~f;
	if (!isSimplyLaced(G,J))  /* unknown type */
	  return type;
	coxtypes::Generator n = constants::firstBit(G.nodes(I));
	f &= G.star(n);
	switch (bits::bitCount(f))
	  {
	  case 2: /* exactly one long branch */
	    J = f | constants::lmask[n];
	    if (isSimplyLaced(G,J))  /* type is b */
	      type[0] = 'b';
	    return type;
	  case 3: /* type is b4 */
	    type[0] = 'b';
	    return type;
	  default: /* more than one long branch */
	    return type;
	  };
      }
      default:  /* unknown type */
	return type;
      };
  }
  case 5: { /* type must be H3 or H4 if known */
    switch (bits::bitCount(I))
      {
      case 3: {
	CoxEntry m1 = minCoefficient(G,I);
	if (m1 == 3)
	  type[0] = 'H';
	return type;
      }
      case 4: {
	if (G.nodes(I))  /* graph is not a string */
	  return type;
	bits::Lflags f = G.extremities(I);
	bits::Lflags J = I & ~f;
	if (!isSimplyLaced(G,J))  /* unknown type */
	  return type;
	J = 0;
	for (; f; f &= f-1)
	  {
	    coxtypes::Generator s = constants::firstBit(f);
	    J |= G.star(s);
	  }
	CoxEntry m1 = minCoefficient(G,J);
	if (m1 == 3)
	  type[0] = 'H';
	return type;
      }
      default:
	return type;
      };
    break;
  }
  case 6: { /* type must be g3 if known */
    switch(bits::bitCount(I))
      {
      case 3: {
	CoxEntry m1 = minCoefficient(G,I);
	if (m1 == 3)
	  type[0] = 'g';
	return type;
      }
      default:
	return type;
      };
  }
  default:  /* unknown type */
    return type;
  };

  return type; // unreachable
}

};

namespace graph {

/*
  Returns the type of the group generated by I as a |type::Type| containing one
  letter for each component of I, in increasing order : A-I if the component is
  finite, a-g if it is affine, X otherwise. So the length of the string is the
  number of components of the group. Returns the empty |type::Type| if I = 0.

  The result is returned as the static variable |type|, which reference is
  safe until the next call to this function (also called |type|).
*/
const type::Type& type(CoxGraph& G, bits::Lflags I)
{
  static type::Type type(0);
  type.name().resize(G.rank()); // long enough, and '\0' filled

  for (Ulong j = 0; I!=0; j++)  /* run through connected components */
    {
      coxtypes::Generator s = constants::firstBit(I);
      bits::Lflags f = G.component(I,s);
      type[j] = (irrType(G,f))[0];
      I &= ~f;
    }

  return type;
}

};

/****************************************************************************

      Chapter IV -- Order computations.


 This section regroups functions for computing the order of subgroups
 generated by subsets of S.

 The functions are :

 - A_order(rank),B_order(rank),D_order(rank) : order functions for the
   infinite families;
 - dihedralOrder(G,I) : returns the order for dihedral groups;
 - extrQuotOrder(G,I,s) : returns the order of the quotient of I by I\{s},
   where I is irreducible and s extremal;
 - finiteOrder(type,rank) : returns the order of the standard irreducible
   finite Coxeter groups;
 - order(G,I) : returns the order of the subgroup generated by I;
 - quotOrder(G,I,J) : returns |W_I/W_J|, assuming that J is included in I;
 - lastQuotOrder(type,rank) : returns the order of the privileged
   quotient in the given type and rank;

   NOTE : the return value 0 should probably be changed to overflow.

 ****************************************************************************/

namespace {

coxtypes::CoxSize A_order(coxtypes::Rank rank)

/*
  Returns the order of the Coxeter group of type A and rank given, if
  this fits into a CoxSize, 0 otherwise. Of course the answer is
  (rank+1)!
*/

{
  coxtypes::CoxSize a = 1;

  for (coxtypes::Rank j = 1; j <= rank; j++)
    {
      if (a > coxtypes::COXSIZE_MAX/(j+1))
	return 0;
      a *= j+1;
    }

  return a;
}

coxtypes::CoxSize B_order(coxtypes::Rank rank)

/*
  Returns the order of the Coxeter group of type B and rank given, if
  this fits into a CoxSize, 0 otherwise. The answer is 2^rank*(rank!)
*/

{
  coxtypes::CoxSize a = 2;

  for (coxtypes::Rank j = 2; j <= rank; j++)
    {
      if (a > coxtypes::COXSIZE_MAX/(2*j))
	return 0;
      a *= 2*j;
    }

  return a;
}


coxtypes::CoxSize D_order(coxtypes::Rank rank)

/*
  Returns the order of the Coxeter group of type D and rank given, if
  this fits into a CoxSize, 0 otherwise. The answer is 2^(rank-1)*(rank!).
*/

{
  coxtypes::CoxSize a = 24;

  for (coxtypes::Rank j = 4; j <= rank; j++)
    {
      if (a > coxtypes::COXSIZE_MAX/(2*j))
	return 0;
      a *= 2*j;
    }

  return a;
}

coxtypes::CoxSize dihedralOrder(CoxGraph& G, bits::Lflags I)

/*
  Assuming that |I| = 2, returns the order of the subgroup generated
  by I. The answer is 2*m, where m=m_{s,t} is the Coxeter coefficient
  determined by I.
*/

{
  coxtypes::CoxSize m;
  coxtypes::Generator s, t;

  s = constants::firstBit(I);
  I &= I-1;
  t = constants::firstBit(I);
  m = (coxtypes::CoxSize)(G.M(s,t));

  if (m > coxtypes::COXSIZE_MAX/2)
    return 0;

  return 2*m;
}


coxtypes::ParSize extrQuotOrder(CoxGraph& G, bits::Lflags I, coxtypes::Generator s)

/*
  Assuming I irreducible and s extremal, this function returns
  the order of the quotient of I by I\{s}.

  It is assumed that LPARNBR_MAX >= 2^31
*/

{
  coxtypes::Rank l;
  bits::Lflags I1;
  coxtypes::Generator s1;
  CoxEntry m;

  const type::Type& t = irrType(G,I);
  l = bits::bitCount(I);

  if (l == 1)
    return (coxtypes::ParSize)2;

  I1 = I & ~constants::lmask[s];
  const type::Type& t1 = irrType(G,I1);

  switch (t[0])
    {
    case 'A':
      return (coxtypes::ParSize)(l+1);
    case 'B':
      switch (t1[0])
	{
	case 'A':  /* return 2^l */
	  if (l == BITS(coxtypes::ParSize))
	    return 0;
	  else
	    return (coxtypes::ParSize)1 << l;
	case 'B':
	  return (coxtypes::ParSize)(2*l);
	};
    case 'D':
      switch (t1[0])
	{
	case 'A':  /* return 2^(l-1) */
	  return (coxtypes::ParSize)1 << (l-1);
	case 'D':
	  return (coxtypes::ParSize)(2*l);
	};
    case 'E':
      switch (l)
	{
	case 6:
	  switch (t1[0])
	    {
	    case 'A':
	      return (coxtypes::ParSize)72;
	    case 'D':
	      return (coxtypes::ParSize)27;
	    };
	case 7:
	  switch (t1[0])
	    {
	    case 'A':
	      return (coxtypes::ParSize)576;
	    case 'D':
	      return (coxtypes::ParSize)126;
	    case 'E':
	      return (coxtypes::ParSize)56;
	    };
	case 8:
	  switch (t1[0])
	    {
	    case 'A':
	      return (coxtypes::ParSize)17280;
	    case 'D':
	      return (coxtypes::ParSize)2160;
	    case 'E':
	      return (coxtypes::ParSize)240;
	    };
	};
    case 'F':
      return (coxtypes::ParSize)24;
    case 'G':
      return (coxtypes::ParSize)6;
    case 'H':
      switch (l)
	{
	case 2:
	  return (coxtypes::ParSize)5;
	case 3:
	  switch (t1[0])
	    {
	    case 'H':
	      return (coxtypes::ParSize)12;
	    case 'A':
	      return (coxtypes::ParSize)20;
	    };
	case 4:
	  switch (t1[0])
	    {
	    case 'H':
	      return (coxtypes::ParSize)120;
	    case 'A':
	      return (coxtypes::ParSize)600;
	    };
	};
    case 'I':
      I &= ~(constants::lmask[s]);
      s1 = constants::firstBit(I);
      m = G.M(s,s1);
      return (coxtypes::ParSize)m;
    default:  /* group is not finite */
      return 0;
    };
}


coxtypes::CoxSize finiteOrder(const type::Type& type, const coxtypes::Rank& rank)

/*
  This function returns the order of the group of the given type
  and rank, by dispatching it to the appropriate sub-function.
  It is assumed that type[0] is one of A-H. Type I is handled
  separately. It is assumed that the rank has been scanned so
  that it is >= 1 in type A, >= 2 in type B, >= 4 in type D,
  6,7,8 in type E, 4 in type F, 2 in type G, 2,3,4 in type H.

  The returnvalue is the order of the group if this fits into
  a CoxSize, 0 otherwise.
*/

{
  switch (type[0]) {
  case 'A':
    return A_order(rank);
  case 'B':
  case 'C':
    return B_order(rank);
  case 'D':
    return D_order(rank);
  case 'E':
    switch (rank) {
    case 6:
      return static_cast<coxtypes::CoxSize>(51840);
    case 7:
      return static_cast<coxtypes::CoxSize>(2903040);
    case 8:
      return static_cast<coxtypes::CoxSize>(696729600);
    };
  case 'F':
    return static_cast<coxtypes::CoxSize>(1152);
  case 'G':
    return static_cast<coxtypes::CoxSize>(12);
  case 'H':
    switch (rank) {
    case 2:
      return static_cast<coxtypes::CoxSize>(10);
    case 3:
      return static_cast<coxtypes::CoxSize>(120);
    case 4:
      return static_cast<coxtypes::CoxSize>(14400);
    };
  default: // unreachable
    return 0;
  };
}


coxtypes::ParSize lastQuotOrder(const type::Type& type, coxtypes::Rank rank)

/*
  Returns the order of the privileged quotient in the given type and
  rank. Ther cannot be any overflow. It is assumed that the type is
  one of A-H.
*/

{
  switch (type[0]) {
  case 'A':
    return static_cast<coxtypes::ParSize>(rank+1);
  case 'B':
  case 'C':
  case 'D':
    return static_cast<coxtypes::ParSize>(2*rank);
  case 'E':
    switch (rank) {
    case 6:
      return static_cast<coxtypes::ParSize>(27);
    case 7:
      return static_cast<coxtypes::ParSize>(56);
    case 8:
      return static_cast<coxtypes::ParSize>(240);
    };
  case 'F':
    return static_cast<coxtypes::ParSize>(24);
  case 'G':
    return static_cast<coxtypes::ParSize>(6);
  case 'H':
    switch (rank) {
    case 2:
      return static_cast<coxtypes::ParSize>(5);
    case 3:
      return static_cast<coxtypes::ParSize>(12);
    case 4:
      return static_cast<coxtypes::ParSize>(120);
    };
  default: // unreachable
    return 0;
  };
}

};


namespace graph {

coxtypes::CoxSize order(CoxGraph& G, bits::Lflags I)

/*
  Returns the order of the subgroup generated by I, if this fits
  into a CoxSize, 0 otherwise.
*/

{
  if (I == 0)
    return 1;

  coxtypes::Generator s = constants::firstBit(I);
  bits::Lflags J = G.component(I,s);

  if (J != I)  /* group is not irreducible */
    {
      coxtypes::CoxSize c1 = order(G,J);
      coxtypes::CoxSize c2 = order(G,I&~J);
      if (c1 & c2 & (c2 > coxtypes::COXSIZE_MAX/c1))  /* overflow */
	return 0;
      return c1*c2;
    }

  const type::Type& t = irrType(G,I);
  coxtypes::Rank l = bits::bitCount(I);

  if (t[0] == 'I')
    return dihedralOrder(G,I);
  else
    return finiteOrder(t,l);
}


coxtypes::ParSize quotOrder(CoxGraph& G, bits::Lflags I, bits::Lflags J)

/*
  Returns the number of elements of W_I/W_J, assuming that J is contained
  in I, that W_I is finite, and that the results does not exceed
  LPARNBR_MAX; returns 0 otherwise.
*/

{
  if (I == J)
    return 1;

  coxtypes::Generator s = constants::firstBit(I);
  bits::Lflags I1 = G.component(I,s);

  if (I1 != I)  /* argue by induction */
    {
      bits::Lflags J1 = J & I1;
      bits::Lflags I2 = I & ~I1;
      bits::Lflags J2 = J & ~J1;
      coxtypes::ParSize c1 = quotOrder(G,I1,J1);
      coxtypes::ParSize c2 = quotOrder(G,I2,J2);
      if (c1 & c2 & (c2 > coxtypes::LPARNBR_MAX/c1))  /* overflow */
	return 0;
      return c1*c2;
    }

  /* now I is irreducible */

  const type::Type& type = irrType(G,I);
  if (strchr("ABCDEFGHI",type[0]) == NULL)  /* group is infinite */
    return 0;

  coxtypes::Rank l = bits::bitCount(I);

  if (l == 2)  /* dihedral case */
    {
      coxtypes::Generator s = constants::firstBit(I);
      coxtypes::Generator t = constants::firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      if (m == 0)  /* group is infinite */
	return 0;

      switch(bits::bitCount(J))
	{
	case 0:
	  return static_cast<coxtypes::ParSize>(2*m);
	case 1:
	  return static_cast<coxtypes::ParSize>(m);
	};
    }

  s = lastGenerator(G,I);

  I1 = I & ~(constants::lmask[s]);
  bits::Lflags J1 = J & ~(constants::lmask[s]);

  coxtypes::ParSize c1 = lastQuotOrder(type,l);
  coxtypes::ParSize c2 = quotOrder(G,I1,J1);
  if (c2 == 0)
    return 0;

  if ((J & constants::lmask[s]) == 0)  /* s is not in J */
    goto exit;

  J = G.component(J,s);

  {
    coxtypes::ParSize q = extrQuotOrder(G,J,s);
    coxtypes::ParSize d = gcd((Ulong)c1,(Ulong)q);
    c1 /= d;
    q /= d;
    c2 /= q;  /* now c2 must be divisible by q */
  }

 exit:
  if (c2 > coxtypes::LPARNBR_MAX/c1)  /* overflow */
    return 0;

  return c1*c2;
}

};

/****************************************************************************

      Chapter V -- Utility functions.

  This section defines various utility functions :

   - gcd(a,b) : yes, the classic Euclidian algorithm;
   - getConjugacyClasses(cl) : puts in cl the conjugacy classes of generators;
   - lastGenerator(G,I) : returns the last generator in the standard
     enumeration of I;
   - maxCoefficient(G,I) : returns the largest coefficient in m|I;
   - minCoefficient(G,I) : returns the smallest coefficient in m|I;

 ****************************************************************************/

namespace {

Ulong gcd(Ulong a, Ulong b)

/*
  The classic Euclidian algorithm. It is assumed that a and b are
  non-zero.
*/

{
  if (a < b)  /* exchange a and b */
    return gcd(b,a);

  Ulong r = a%b;

  while (r != 0)
    {
      a = b;
      b = r;
      r = a%b;
    }

  return b;
}

};

namespace graph {


/*
  This function returns in cl the conjugacy classes of generators in W (i.e.
  the partition of S induced by the partition of W in conjugacy classes.) It
  is known that these are the connected components of the graph obtained by
  removing from the Coxeter graph all edges with even or infinite labels.
*/
containers::vector<bits::Lflags> conjugacy_classes(const CoxGraph& G)
{
  containers::vector<bits::Lflags> odd_star; // odd-edge neighbours of each node
  odd_star.reserve(G.rank());

  for(coxtypes::Generator s = 0; s < G.rank(); ++s) {
    odd_star.push_back(0);
    for (coxtypes::Generator t = 0; t < G.rank(); ++t)
      if ((G.M(s,t)%2) && (G.M(s,t) > 1))
	odd_star.back() |= constants::lmask[t];
  }

  containers::vector<bits::Lflags> classes;

  for (bits::Lflags fS = G.supp(); fS; /* |fS &= ~f| */) {
    bits::Lflags nf = constants::lmask[constants::firstBit(fS)];
    bits::Lflags f = 0;
    while (nf)  /* there are new elements to be considered */
      {
	f |= nf;
	for (bits::Lflags f1 = nf; f1; f1 &= f1-1)
	  nf |= (odd_star[constants::firstBit(f1)]);
	nf &= ~f;
      }
    classes.push_back(f);
    fS &= ~f;
  } // |for(fS)|

  return classes;
}

};

namespace {

coxtypes::Generator lastGenerator(CoxGraph& G, bits::Lflags I)

/*
  Assuming that I is irreducible, this function returns an element
  in I which can be taken as a last generator in the standard enumeration.
*/

{
  coxtypes::Rank l = bits::bitCount(I);
  if (l <= 2)
    return constants::firstBit(I);

  /* from now on the rank is at least three */

  const type::Type& x = irrType(G,I);
  bits::Lflags f = G.extremities(I);

  switch (x[0])
    {
    case 'A':
      return constants::firstBit(f);
    case 'B': {
      coxtypes::Generator s = constants::firstBit(f);
      coxtypes::Generator t = constants::firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      switch (m)
	{
	case 3:
	  return s;
	case 4:
	  f &= ~(constants::lmask[s]);
	  return constants::firstBit(f);
	};
    }
    case 'D': {
      coxtypes::Generator s = constants::firstBit(f);
      coxtypes::Generator n = constants::firstBit(G.nodes(I));
      f &= ~(G.star(n));
      if (f)
	return constants::firstBit(f);
      else  /* rank is 4 */
	return s;
    }
    case 'E': {
      coxtypes::Generator n = constants::firstBit(G.nodes(I));
      f &= ~(G.star(n));
      coxtypes::Generator s = constants::firstBit(f);
      switch (l)
	{
	case 6:
	  return s;
	case 7:
	case 8: {
	  coxtypes::Generator t = constants::firstBit(G.star(I,s));
	  if (constants::lmask[t] & G.star(n)) {
	    f &= ~(constants::lmask[s]);
	    return constants::firstBit(f);
	  }
	  else
	    return s;
	}
	};
    }
    case 'F':
      return constants::firstBit(f);
    case 'H': {
      coxtypes::Generator s = constants::firstBit(f);
      coxtypes::Generator t = constants::firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      switch (m)
	{
	case 3:
	  return s;
	case 5:
	  f &= ~(constants::lmask[s]);
	  return constants::firstBit(f);
	};
    }
    case 'a':
      return constants::firstBit(I);
    case 'b': {
      coxtypes::Generator s = constants::firstBit(f);
      coxtypes::Generator t = constants::firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      switch (m)
	{
	case 3:
	  return s;
	case 4:
	  f &= ~(constants::lmask[s]);
	  return constants::firstBit(f);
	};
    }
    case 'c':
      return constants::firstBit(f);
    case 'd':
      return constants::firstBit(f);
    case 'e': {
      switch (l)
	{
	case 7:
	  return constants::firstBit(f);
	case 8: {
	  coxtypes::Generator n = constants::firstBit(G.nodes(I));
	  f &= ~(G.star(n));
	  return constants::firstBit(f);
	}
	case 9: {
	  coxtypes::Generator n = constants::firstBit(G.nodes(I));
	  f &= ~(G.star(n));
	  coxtypes::Generator s = constants::firstBit(f);
	  coxtypes::Generator t = constants::firstBit(G.star(I,s));
	  if (constants::lmask[t] & G.star(n))
	    {
	      f &= ~(constants::lmask[s]);
	      return constants::firstBit(f);
	    }
	  else
	    return s;
	}
	};
    }
    case 'f': {
      coxtypes::Generator s = constants::firstBit(f);
      bits::Lflags I1 = I & ~(constants::lmask[s]);
      switch ((irrType(G,I1))[0])
	{
	case 'B':
	  f &= ~(constants::lmask[s]);
	  return constants::firstBit(f);
	case 'F':
	  return s;
	};
    }
    case 'g': {
      coxtypes::Generator s = constants::firstBit(f);
      coxtypes::Generator t = constants::firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      switch (m)
	{
	case 3:
	  return s;
	case 6:
	  f &= ~(constants::lmask[s]);
	  return constants::firstBit(f);
	};
    }
    default:
      return constants::lastBit(I);
    };
}


CoxEntry maxCoefficient(CoxGraph& G, bits::Lflags I)

/*
  Returns the maximal coefficient in the Coxeter matrix restricted to
  I (we assume that I is not empty).
*/

{
  if (bits::bitCount(I) == 1)
    return 1;

  CoxEntry m = 2;

  for (bits::Lflags fs = I; fs; fs &= fs-1)
    {
      coxtypes::Generator s = constants::firstBit(fs);
      for (bits::Lflags ft = fs&G.star(s); ft; ft &= ft-1)
	{
	  coxtypes::Generator t = constants::firstBit(ft);
	  if (G.M(s,t) == 0)
	    return 0;
	  if (G.M(s,t) > m)
	    m = G.M(s,t);
	}
    }

  return m;
}


CoxEntry minCoefficient(CoxGraph& G, bits::Lflags I)

/*
  Returns the minimal coefficient > 2 in the Coxeter matrix restricted
  to I, if there is such; otherwise, if |I| > 1, returns 2; otherwise,
  returns 1. It is assumed that I is not empty.
*/

{
  if (bits::bitCount(I) == 1)
    return 1;

  CoxEntry m = maxCoefficient(G,I);
  if (m == 2)
    return 2;

  for (coxtypes::Generator s = 0; s < G.rank(); s++)
    for (bits::Lflags f = I&G.star(s); f; f &= f-1)
      {
	coxtypes::Generator t = constants::firstBit(f);
	if ((G.M(s,t) != 0) && (G.M(s,t) < m))
	  m = G.M(s,t);
      }

  return m;
}

};
