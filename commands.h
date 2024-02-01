/*
  This is commands.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef COMMANDS_H  /* guard against multiple inclusions */
#define COMMANDS_H

#include <string>
#include "globals.h"

namespace commands {
  using namespace globals;
};

#include "coxgroup.h"

/******** type declarations ************************************************/

namespace commands {
  struct CommandData;
  class CommandTree;
};

/******** constants ********************************************************/

namespace commands {
  extern void (*default_help)();
};

/******** function declarations ********************************************/

namespace commands {
  coxgroup::CoxGroup* currentGroup(); // get group fixed upon main mode entry
  void default_error(const char* str); // report |COMMAND_NOT_FOUND|
  CommandTree* mainCommandTree(); // for ordinary computation commands
  CommandTree* uneqCommandTree(); // for commands for unequal-parameter groups
  CommandTree* interfaceCommandTree();
  namespace interface {
    CommandTree* inCommandTree();
    CommandTree* outCommandTree();
  };
  void printCommands(FILE* file, const CommandTree& tree);
  void relax_f(); // no-op, to be used as action function
  void run();
};

/******** Type definitions *************************************************/

#include "dictionary.h"
#include "io.h"

namespace commands {
  using namespace dictionary;
  using namespace io;
};

namespace commands {

struct CommandData {
  std::string name;
  std::string tag;
  void (*action)();
  void (*help)();
  bool autorepeat;
/* Constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(CommandData));}
  CommandData(const char* str, const char* t, void (*a)(),
	      void (*h)(), bool rep);
  ~CommandData();
};

class CommandTree:public Dictionary<CommandData> {
 private:
  std::string d_prompt;
  CommandTree* d_help;
  void (*d_entry)();
  void (*d_error)(const char* str);
  void (*d_exit)();
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(CommandTree));}
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
  void set_default_action(void (*a)()); // assign |a| to |d_root->action|
  void setRepeat(const char* str, bool b);
/* accessors */
  void prompt() const; // print |d_prompt|
  void call_entry() const { d_entry(); }
  void call_error(const char *str) const { d_error(str); }
  void call_exit() const { d_exit(); }
  CommandTree* helpMode() const { return d_help; } // may return |nullptr|
};

};

#endif
