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

template <class T> Dictionary<T>::Dictionary(std::shared_ptr<T> v)
  : d_root(new DictCell<T>('\0',std::move(v))) {}


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
  client to decide what payload data to store there, but the condition of
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
    while (c > cell->letter and cell->right!=nullptr)
      cell = cell->right;
    if (c != cell->letter)
      return nullptr;
  }

  return cell;
}

template <class T>
T* Dictionary<T>::find(const std::string& str, bool& absent_action) const

{
  DictCell<T>* dc = findCell(str);

  if (dc == nullptr)
    return absent_action=false,nullptr;
  absent_action = dc->ptr==nullptr;
  return dc->ptr.get();
}


/*
  Insert a new word in the dictionary. The root of the tree always corresponds
  to the empty string "" (which may or may not be considered to be in the
  dictionary.) The nodes are classified in three types : dictionary words
  (|fullname| holds), unique prefixes (|not fullname and uniquePrefix|) and
  non-unique prefixes (|not fullname and not uniquePrefix|). These attributes
  are defined or modified on the cells on the path to our final cell.

  We use a "triple ref" trick, which is Algol68 lingo for a variable pointer to
  a pointer. Here that variable |p| points to the link that either points to an
  existing cell we want to modify, or after which we want to insert a fresh cell
  if none is present. The work of inserting can then be neatly split into
  advancing this to the correct link, and then doing our stuff depending on
  whether a node is already present and whether we are at the |final| letter of
  |str|. One nice aspect is that cases where we insert into a left or right
  descendant link require no distinction at all.
*/
template <class T> void Dictionary<T>::insert(const std::string& str,
					      std::shared_ptr<T> value)
{
  DictCell<T>** p = &d_root->left;
  for (const char& c: str)
  {
    bool final = &c==&str[str.length()-1];
    while(*p != nullptr and (*p)->letter<c)
      p = &(*p)->right; // skip over cells with |letter<c|
    if (*p != nullptr and (*p)->letter==c) // whether a proper cell is present
    { // so far we match an existing prefix
      if (final) // |str| is an existing prefix, replace any data
	(*p)->ptr = value;
    }
    else // not looking at letter |c|
    {
      if (final)
	*p = new DictCell<T>(c,value,nullptr,*p);
      // having the |uniquePrefix| argument be |false| seems wrong here, but
      // the convention appears to be |fullname| implies |not uniquePrefix|
      else
	*p = new DictCell<T>(c,nullptr,nullptr,*p);
    } // |if (letter==c)|
    p = &(*p)->left; // either way continue with left link for this cell

  } // |for(c:str)|
}


/*
  Not implemented.
*/
template <class T> void Dictionary<T>::remove(const std::string& str)
{}


// here is how one can define an embedded class method out-of-line

/* Advance from |p| pointing to a cell to making it point to the next cell.
   In order to be able to move back up, ancestors for which we are a left
   descendent are maintained on |stack|.
*/
template <typename T>
  typename Dictionary<T>::const_iterator
  Dictionary<T>::const_iterator::operator++()
{
  if (p->left!=nullptr)
  {
    stack.push(p);
    p=p->left;
    return *this;
  }

  // now we have no left descendent to advance to,
  // move to right descendent, possibly from an ancestor when lacking any heir
  while (p->right==nullptr and not stack.empty())
  {
    p=stack.top();
    stack.pop();
  }

  p = p->right; // move to next heir in line, or |nullptr| when there is none
  return *this;
}

}; // |namespace dictionary|



/****** Auxiliary functions ************************************************/

namespace dictionary {

template <class T>
  void DictCell<T>::make_complete()
{
  if (left!=nullptr)
    left->make_complete();
  if (ptr==nullptr) // do nothing for prefixes that already have an action
  {
     assert(left!=nullptr); // if not final, there is an extension
     if (left->right==nullptr) // then there is a unique extension letter
       ptr = left->ptr; // copy its action here (or nothing, if ambiguous)
  }
  if (right!=nullptr)
    right->make_complete();
}

template <class T>
  bool DictCell<T>::has_own_action() const
{ if (ptr==nullptr)
    return false;
  if (left==nullptr)
    return true;
  return (left->ptr!=ptr);
}

/*
  This method lists all the known extensions of the string |name|, which should
  by the one corresponding to our cell.
*/
template <class T>
  containers::sl_list<std::string> DictCell<T>::extensions
    (std::string name) const
{
  containers::sl_list<std::string> result;
  for (auto* cell = left;  // use only left subtree
       cell != nullptr;
       cell = cell->right) // walk down right-list of left subtree
  {
    name.push_back(cell->letter);
    if (cell->has_own_action())  // then append prefix for current cell
      result.push_back(name);
    result.append(cell->extensions(name));
    name.pop_back(); // restore original |name| for next iteration
  }
  return result;
}


};
