/*
  This is commands.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef COMMANDS_H  /* guard against multiple inclusions */
#define COMMANDS_H

#include <string>
#include "globals.h"
#include "dictionary.h"
#include "io.h"

#include "coxgroup.h"

/******** type declarations ************************************************/

namespace commands {
  struct CommandData;
  class CommandTree;
};

/******** constants ********************************************************/

/******** function declarations ********************************************/

namespace commands {

  // these are exported implicitly, since used as default function arguments
  extern void (*const default_help)(); // defined here, its value in help
  void default_error(const char* str); // report |COMMAND_NOT_FOUND|

  // the follwing are exported so that help.cpp can print command lists
  CommandTree* mainCommandTree(); // for ordinary computation commands
  CommandTree* uneqCommandTree(); // for commands for unequal-parameter groups
  CommandTree* interfaceCommandTree();
  namespace interf {
    // this namespace used to be called |interface|, but that makes it
    // hard to refer to the global namespace of the same name, so don't
    CommandTree* inCommandTree();
    CommandTree* outCommandTree();
  };
  void printCommands(FILE* file, const CommandTree& tree);

  void relax_f(); // no-op, to be used as action function
  void run(); // the function called by |main| to run the main command loop
};

/******** Type definitions *************************************************/


namespace commands {

struct CommandData {
  std::string name;
  std::string tag;
  void (*action)();
  void (*help)();
  bool autorepeat; // whether this command can be autorepeated
/* Constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(CommandData));}
  CommandData(const char* str, const char* t, void (*a)(),
	      void (*h)(), bool rep);
  ~CommandData();
};

class CommandTree:public dictionary::Dictionary<CommandData> {
 private:
  std::string d_prompt;
  CommandTree* d_help;
  void (*d_entry)();
  void (*d_error)(const char* str);
  void (*d_exit)();
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(CommandTree));}
  CommandTree(const char *str, void (*action)(),
	      void (*filler)(CommandTree& tree), // function filling the tree
	      void (*entry)() = &relax_f,
	      void (*error)(const char*) = &default_error,
	      void (*exit)() = &relax_f, void (*h)() = nullptr
    );
  ~CommandTree();
/* modifiers */
  void add(const char* name, const char* tag, void (*action)(),
	   void (*help)() = default_help, bool rep = true);
  void set_default_action(void (*a)()) { root_action()->action=a; }

/* accessors */
  void prompt() const; // print |d_prompt|
  void call_entry() const { d_entry(); }
  void call_error(const char *str) const { d_error(str); }
  void call_exit() const { d_exit(); }
  CommandTree* helpMode() const { return d_help; } // may return |nullptr|
};

};

#endif
