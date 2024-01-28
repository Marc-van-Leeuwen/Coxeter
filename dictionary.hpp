/*
  This is dictionary.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

/****************************************************************************

  This module contains an implementation of dictionaries. For us, a dictionary
  is a search structure (in practice, a tree) recognizing strings. The
  operations which are allowed are inserting and deleting a string, and
  searching for a given string. The stored value is a shared pointer; which
  will be returned from search requests; it is shared because many nodes will
  share the same data, due to automatic interpretation of unique prefixes.

 ****************************************************************************/

namespace dictionary {

/****************************************************************************

        Chapter I -- The DictCell structure.

 ****************************************************************************/

/*
  This is a basic structure, holding data via a shared pointer to to our
  template type |T|, some character data pertaining to the place this node has
  in the search tree, and two pointers |left|, |right| to form a binary tree.

  The unique constructor given in the class definition just sets the field from
  the arguments; by default it set the pointerss to |nullptr|, but if different
  pointers are explicitly supplied, the constructor takes posession of them.

  The cell owns the entire subtree it gives access to, with the excpetion of the
  user data for which it has shared ownership. This destructor will therefore
  recursively remove the cells of that subtree, and the user data if no external
  sharing is left. The data pointer |ptr| being a smart pointer, the latter is
  taken care of automatically, so we only need to |delete| our descendants.
*/
template <class T> DictCell<T>::~DictCell()
{
  delete left; delete right; // this is implicitly doubly recursive
}


/****************************************************************************

        Chapter II -- The Dictionary class.

 ****************************************************************************/


/*
  Default an only constructor for the Dictionary class. A dictionary is always
  created with a root cell, corresponding to the empty string "". Note that
  the |letter| field is set to |'\0'| to signal a leaf node (iniitally).
*/

template <class T> Dictionary<T>::Dictionary()
  : d_root(new DictCell<T>('\0',nullptr,true,false)) {}


/*
  Destructor for the Dictionary class. The destructor has to run through
  the cells, and destroy each of them. The recursion is taken care of by
  |DictCell<T>::~DictCell| so we simply need to |delete| the root pointer.
*/
template <class T> Dictionary<T>::~Dictionary()
{
  delete d_root;
}


/*
  Search for a value in the dictionary. Return |nullptr| in case of failure.

  The dictionary has a curious layout, where the string represented by a cell
  has its stored character as final character, the previous one is the closest
  ancestor cell from whic it was readched by LEFT descent, and so forth; the
  root node has a null character that is excuded from this sequence of
  characters , and a right descendent that is forever |nullptr|, so every
  nonzero string corresponds to a cell in the left subtree of the root.

  In a search if we successfully found a part of |str| and are at its
  correspondin cell, but the is more to find, we always descend to the left
  cell, which is the first ode of effectively a linked list via the RIGHT
  descendents of cells for strings that have onemore letter at the end.

  NOTE : in fact, recognizes if |str| is a prefix of a word in the dictionary,
  as any such prefix has a corresponding cell with data stored. It is up to the
  client to decide what payload data tot store there, but the condition of
  whether this prefix is actually a stored word (|fullname|) and if as a prefix
  it is unique (|uniquePrefix|) are stored in the cell proper. It is not
  possible to restrict to only store leaf nodes.
*/
template <class T> DictCell<T>* Dictionary<T>::findCell(const std::string& str)
  const
{
  DictCell<T> *cell = d_root;

  for (Ulong j = 0; str[j]; ++j)
  { char c = str[j]; // we need to find this character
    if (cell->left == nullptr) /* leaf reached */
      return nullptr;
    cell = cell->left;
    while (cell->right!=nullptr and  c > cell->letter)
      cell = cell->right;
    if (c != cell->letter)
      return nullptr;
  }

  return cell;
}

template <class T>
  std::shared_ptr<T> Dictionary<T>::find(const std::string& str) const

{
  DictCell<T>* dc = findCell(str);

  if (dc)
    return dc->value();
  else
    return nullptr;
}


/*
  Insert a new word in the dictionary. The root of the tree always
  corresponds to the empty string "" (which may or may not be considered
  to be in the dictionary.) The nodes are classified in three types :
  dictionary words, unique prefixes (i.e., strings which are prefixes
  to a unique dictionary word) and non-unique prefixes.
*/
template <class T> void Dictionary<T>::insert(const std::string& str,
					      std::shared_ptr<T> value)
{
  DictCell<T>* cell = findCell(str);

  if (cell && cell->fullname) { /* word was already in the dictionary */
    cell->ptr = value; // overwrite any stored value
    return;
  }

  /* from now on we are sure that the word is not already in the
     dictionary */

  cell = d_root; // restart

  for (Ulong j = 0; j<str.length(); ++j)
    {
      const bool final = j+1 == str.length();
      if (cell->left == nullptr) { // end of stored prefix reached
	if (final)  // one more character to record: add a leaf |fullname| node
	  cell->left = new DictCell<T>(str[j],value,true,false);
	else // at least 2 more chars, add a non |fullname|
	  cell->left = new DictCell<T>(str[j],nullptr,false,true);
	cell = cell->left; // descend into newly created node
	continue; // w.r.t. |for(j)|; remainder is effectively |else| branch
      }

      // now |cell->left!=nullptr|, we are prefix of some previous word

      if (str[j] < cell->left->letter) { /* insert at beginning */
	if (final)  // one more character to record: add a leaf |fullname| node
	  cell->left = new DictCell<T>(str[j],value,true,false,
				       nullptr,cell->left);
	else // at least 2 more chars, add a non |fullname|
	  cell->left = new DictCell<T>(str[j],nullptr,false,true,
				       nullptr,cell->left);
	cell = cell->left;
	continue;
      }
      cell = cell->left;
      while (cell->right!=nullptr and cell->right->letter <= str[j])
	cell = cell->right;
      if (cell->letter < str[j]) { /* add new cell */
	if (final)  // one more character to record: add a leaf |fullname| node
	  cell->right = new DictCell<T>(str[j],value,true,false,
					nullptr,cell->right);
	else // at least 2 more chars, add a non |fullname|
	  cell->right = new DictCell<T>(str[j],0,false,true,
					nullptr,cell->right);
	cell = cell->right;
	continue;
      }

      /* if we reach this point cell->letter = str[j] */

      cell->uniquePrefix = false;
      if (final) { // word is complete
	cell->fullname = true;
	cell->ptr = value;
      }
    }
}

template <class T> void Dictionary<T>::remove(const std::string& str)

/*
  Not implemented.
*/

{}

};



/****** Auxiliary functions ************************************************/

namespace dictionary {


/*
  This function prints to |file| all the possible extensions of the string
  corresponding to |cell|. The string |name| should contain the strict prefix,
  which excludes |cell->letter|, it is used as a working variable but will be
  restored to its original value (but not necessesarily original capacity) upon
  returning from this recursive function. The output will start with |sep|
  unless |first|, which after so suppressing |sep| once is set to |false|.

  The actual behaviour reflects the weird structure of dictionaries: not only
  are the extensions for the current |cell| printed, but also the extensions of
  its strict prefix |name| whose final letter is alphabetically beyond
  |cell->letter|, since these are present in the right subtree. This makes the
  recurion easy, and makes no difference when called with |cell| equal to the
  root of the tree, which has an empty right subtree, and does the right thing
  when |cell| is |parent->left|, where |parent| is the node for |name|.
*/
template <class T>
  void printExtensions(FILE* file, DictCell<T>* cell, std::string& name,
		       bool &first, const char* sep)
{
  if (cell == nullptr)
    return;
  name.push_back(cell->letter);
  if (cell->fullname) { // print prefix for current cell
    if (first) /* first time a name is found */
      first = false;
    else
      fprintf(file,"%s",sep);
    fprintf(file,"%s",name.c_str());
  }
  printExtensions(file,cell->left,name,first,sep);
  name.pop_back(); // remove last character from |name|
  printExtensions(file,cell->right,name,first,sep);
}

};
