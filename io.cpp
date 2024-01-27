/*
  This is io.cpp

  Coxeter version 3.0_demo  Copyright (C) 2001 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include <stdarg.h>
#include <ctype.h>

#include "io.h"

/****************************************************************************

  This file used to be one of the biggest in the program; now it contains
  just some elementary output functions that really pertain to the standard
  libraries; in the current stage of the transition, and also because of
  nagging worries about memory allocation performance, I have preferred
  to keep my own clumsy string types.

  The interaction with the user has been farmed out to interactive.c; the
  output functions for the Coxeter group have gone mostly to interface.c;
  and the output functions for the various types are in the corresponding
  files (with the exception of polynomials, because I got stuck in the
  template problems.)

 ****************************************************************************/

/*****************************************************************************

        Chapter I --- Output functions returning strings

  This section regroups various output functions returning strings. We have
  not strived for speed, but for safety and robustness. The main aim has been
  to procure output functions that are easy to use, and relieve the user
  from the worries of having to allocate memory, or watch out for overflow.
  So usually each function posesses its own "safe" place, and the result is
  returned as a pointer to there. The size of the safe area is grown as
  necessary.

  The functions provided here are :

  - append(l,c) : appends the character c to the varstring l;
  - append(l,s) : appends the (c)string s to the varstring l;
  - append(l1,l2) : appends the varstring l2 to the varstring l1;
  - append(l,n) : appends the Ulong (int, unsigned) n to the varstring l;
  - erase(l,n) : erases the last n characters from l;
  - pad(l,n) : pads l with spaces to length n;
  - reset(l) : resets a varstring to the empty line;
  - setString(l,s,first,r) : sets l to a substring of s;

******************************************************************************/

namespace io {



std::string& append(std::string& str, const Ulong& n)
{
  return str.append(std::to_string(n));
}


std::string& append(std::string& str, const long& m)
{
  return str.append(std::to_string(m));
}



std::string& append(std::string& str, const int& n)
{
  return str += std::to_string(n);
}



// Pads the string with white spaces to length n.
std::string& pad(std::string& l, const Ulong& n)
{
  if (n <= l.length()) /* do nothing */
    return l;

  return l.append(n-l.length(),' ');
}


};


/*****************************************************************************

        Chapter II --- Output to files

  This chapter regroups some basic input/output functions at the level
  of strings and such. Of course this should have come from the STL.

  The following functions are provided :

  - foldLine(file,str,ls,h,hyphens) : fold a long output line into file;
  - print(file,str) : prints the (var)string to the file;
  - print(file,v,n) : prints a list of integers;
  - printFile(file,name) : prints the contents of the file name;
  - printFile(file,name,dir_name) : prints the contents of the file
    dir_name/name;

******************************************************************************/

namespace io {


/*
  This function breaks up the string str into lines of at most |ls| characters
  to improve output of long lines. It outputs the extra lines with an
  indentation of |h| white spaces. It chooses a breakpoint just before one of
  the characters in the string |hyphens|, whenever possible; otherwise it
  breaks the line brutally where it must.
*/
void foldLine(FILE* file, const std::string& str, const Ulong& ls,
	      const Ulong& h, const char* hyphens)
{
  if (str.length() <= ls) { /* str fits on one line */
    io::print(file,str);
    return;
  }

  /* search for hyphenation point */

  Ulong bp = 0; // break point

  for (Ulong j = 0; j < ls; j += // distance to next |hyphen|, or '\0'
	 std::strcspn(str.c_str()+j,hyphens))
  {
    bp = j;
    ++j; // skip over |hyphen| charcater before searching again
  }

  if (bp == 0) // no |hyphen| found
    bp = ls; // then break brutally

  print(file,str.substr(0,bp)); // initial part of |str|

  /* print continuation lines */

  Ulong p = bp; // continue where we left off

  while (p < str.length()-ls+h)
  {
    bp = 0; // repeat previous computation, from position |p|, targeting |ls-h|
    for (Ulong j = 0; j < ls-h; j += strcspn(str.c_str()+p+j,hyphens)) {
      bp = j;
      j++;
    }
    if (bp == 0)
      bp = ls-h;
    fprintf(file,"\n%*s",static_cast<int>(h),""); // spaces
    print(file,str.substr(p,bp)); // line between two break points
    p += bp;
  }

  /* print last line */

  fprintf(file,"\n%*s",static_cast<int>(h),""); // spaces
  print(file,str.substr(p)); // remainder of |str|

  return;
}



/*
  Print to |file| the string representation of the first |n| elements pointed
  by |v|, as a comma-separated and square-bracket-enclosed list.
*/

void print(FILE *file, const int* const& v, const Ulong& n)
{
  fprintf(file,"[");

  for (Ulong j = 0; j < n; j++)
    {
      fprintf(file,"%d",v[j]);
      if (j+1 < n)  /* more to come */
	fprintf(file,",");
    }

  fprintf(file,"]");

  return;
}


// Print to |file| the contents of another file named |dir_name|/|name|
void printFile(FILE* file, const char *name, const char *dir_name)
{
  std::string buf = dir_name;
  buf.push_back('/');
  buf.append(name);

  FILE* inputfile;
  char c;

  inputfile = fopen(buf.c_str(),"r");

  if (inputfile == 0) {
    Error(FILE_NOT_FOUND,buf.c_str());
    return;
  }

  while((c = getc(inputfile)) != EOF)
    putc(c,file);

  fclose(inputfile);
}

// Print to |file| the contents of another file named |name|
void printFile(FILE* file, const char *name)
{
  FILE* inputfile;
  char c;

  inputfile = fopen(name,"r");

  if (inputfile == 0) {
    Error(FILE_NOT_FOUND,name);
    return;
  }

  while((c = getc(inputfile)) != EOF)
    putc(c,file);

  fclose(inputfile);
}

};


/*****************************************************************************

        Chapter III -- General input functions.

  This section contains general input functions.

  The following functions are defined :

  - getInput(file,buf) : the universal input function;

******************************************************************************/

namespace io {

/*
  Read from |inputfile| until either |EOF| or a newline is reached;
  append the result to |buf| starting at position |len|; resize |buf|
  as needed while writing. The newline is not written on the |buf|.
*/

const char* getInput(FILE *inputfile, std::string& buf, Ulong len)
{
  buf.resize(len);
  while(true)
  {
    int c = getc(inputfile);
    if ((c == EOF) || (c == '\n'))
      break;
    buf.push_back(c);
  }

  return buf.c_str();
}

};


/*****************************************************************************

        Chapter III --- Miscellaneous.

This chapter regroups various functions that did not seem to fall into any
of the other categories :

  - alphabeticDigits(c,b) : returns the number of digits in
    alphabetic representation

******************************************************************************/

namespace io {

int alphabeticDigits(Ulong c, Ulong b)

/*
  Returns the number of digits for the alphabetic representation of c
  with b letters. In this representation, the numbers n s.t.

     1 + b + ... + b^(d-1) <= n < 1 + b + ... + b^d

  are representable with d digits (and 0 is represented by the empty
  word.)
*/

{
  Ulong j = 0;

  for(; c; c = (c-1)/b)
    ++j;

  return j;
}


int digits(Ulong c, Ulong b)

/*
  Returns the number of digits in the representation of c in base b.
*/

{
  int j = 1;

  for (c /= b; c; c /= b)
    ++j;

  return j;
}


/*
  Skips from character position |p| over white space (characters recognized
  by |isspace|), return number of places advanced
*/
Ulong skipSpaces(const std::string& l, Ulong p)
{
  Ulong j = p;

  while (j<l.length() and isspace(l[j]))
    ++j;

  return j-p;
}

};
