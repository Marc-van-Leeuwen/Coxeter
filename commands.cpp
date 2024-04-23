/*
  This is commands.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "commands.h"

#include <memory> // for |std::unique_ptr|

#include "directories.h"
#include "sl_list.h"
#include "error.h"
#include "fcoxgroup.h"
#include "help.h"
#include "interactive.h"
#include "special.h"
#include "typeA.h"

namespace commands {

  using namespace error;

namespace {

  bool wgraph_warning = true; // a user preference about repeating warnings

  /* used in the definition of command trees */
  struct Empty_tag {};
  struct Interface_tag {};
  struct Main_tag {};
  struct Uneq_tag {};

  containers::stack<CommandTree*> mode_stack;

  coxgroup::CoxGroup* W = nullptr; // the current working Coxeter group
  // this is set depending on modes on the stack, by entry and exit functions

  void activate(CommandTree& tree); // push |tree| onto the command stack
  void report_ambiguous_command
    (const CommandTree& tree, const std::string& str);

  void empty_error(const char* str); // "error" function in empty mode
  CommandTree* emptyCommandTree(); // get pointer to static mode variable

  template<class C> // C is just a tag to select one of these
    void initCommandTree(CommandTree&); // command initialization functions

  // the following are functions whose pointers will go into command trees

  void startup(); // will be installed as action for an initial empty command

  void interface_entry();
  void interface_exit();

  void main_entry();
  void main_exit();

  void uneq_entry();
  void uneq_exit();

  void author_f();
  void betti_f();
  void coatoms_f();
  void compute_f();
  void descent_f();
  void duflo_f();
  void extremals_f();
  void fullcontext_f();
  void help_f();
  void ihbetti_f();
  void inorder_f();
  void interface_f();
  void interval_f();
  void invpol_f();
  void klbasis_f();
  void lcorder_f();
  void lcells_f();
  void lcwgraphs_f();
  void lrcorder_f();
  void lrcells_f();
  void lrcwgraphs_f();
  void lrwgraph_f();
  void lwgraph_f();
  void matrix_f();
  void mu_f();
  void pol_f();
  void q_f();
  void qq_f();
  void rank_f();
  void rcorder_f();
  void rcells_f();
  void rcwgraphs_f();
  void rwgraph_f();
  void schubert_f();
  void show_f();
  void showmu_f();
  void slocus_f();
  void sstratification_f();
  void type_f();
  void uneq_f();

  // the short help strings for the above commands
  const char* author_tag = "prints a message about the author";
  const char* betti_tag = "prints the ordinary betti numbers";
  const char* coatoms_tag = "prints out the coatoms of an element";
  const char* compute_tag = "prints out the normal form of an element";
  const char* descent_tag = "prints out the descent sets";
  const char* duflo_tag = "prints out the Duflo involutions";
  const char* extremals_tag =
    "prints out the Kazhdan-Lusztig polynomials for the extremal pairs";
  const char* fullcontext_tag = "sets the context to the full group";
  const char* help_tag = "enters help mode";
  const char* ihbetti_tag = "prints the IH betti numbers";
  const char* input_tag = "(in help mode only) explains the input conventions";
  const char* interface_tag = "changes the interface";
  const char* interval_tag = "prints an interval in the Bruhat ordering";
  const char* intro_tag =
    "(in help mode only) prints a message for first time users";
  const char* inorder_tag = "tells whether two elements are in Bruhat order";
  const char* invpol_tag = "prints a single inverse Kazhdan-Lusztig polynomial";
  const char* klbasis_tag = "prints an element of the Kazhdan-Lusztig basis";
  const char* lcorder_tag = "prints the left cell order";
  const char* lcells_tag = "prints out the left Kazhdan-Lusztig cells";
  const char* lcwgraphs_tag = "prints out the W-graphs of the left Kazhdan-Lusztig cells";
  const char* lrcorder_tag = "prints the two-sided cell order";
  const char* lrcells_tag = "prints out the tow-sided Kazhdan-Lusztig cells";
  const char* lrcwgraphs_tag =
    "prints out the W-graphs of the two-sided Kazhdan-Lusztig cells";
  const char* lrwgraph_tag = "prints out the two-sided W-graph";
  const char* lwgraph_tag = "prints out the left W-graph";
  const char* matrix_tag = "prints the current Coxeter matrix";
  const char* mu_tag = "prints a single mu-coefficient";
  const char* pol_tag = "prints a single Kazhdan-Lusztig polynomial";
  const char* q_tag = "exits the current mode";
  const char* qq_tag = "exits the program";
  const char* rank_tag = "resets the rank";
  const char* rcorder_tag = "prints the right cell order";
  const char* rcells_tag = "prints out the right Kazhdan-Lusztig cells";
  const char* rcwgraphs_tag = "prints out the W-graphs of the right Kazhdan-Lusztig cells";
  const char* rwgraph_tag = "prints out the right W-graph";
  const char* schubert_tag = "prints out the kl data for a schubert variety";
  const char* show_tag = "maps out the computation of a Kazhdan-Lusztig polynomial";
  const char* showmu_tag = "maps out the computation of a mu coefficient";
  const char* slocus_tag =
    "prints the rational singular locus of the Schubert variety";
  const char* sstratification_tag =
    "prints the rational singular stratification of the Schubert variety";
  const char* type_tag =
    "resets the type and rank (hence restarts the program)";
  const char* uneq_tag = "puts the program in unequal-parameter mode";

  namespace uneq {
    // these commands have a modified definition in 'uneq' mode

    void klbasis_f();
    void lcorder_f();
    void lrcorder_f();
    void lcells_f();
    void lrcells_f();
    void mu_f();
    void pol_f();
    void rcells_f();
    void rcorder_f();

    // their help tags could be different in this mode
    const char* lcorder_tag = "prints the left cell order";
    const char* lrcorder_tag = "prints the two-sided cell order";
    const char* lcells_tag = "prints out the left Kazhdan-Lusztig cells";
    const char* lrcells_tag = "prints out the two-sided Kazhdan-Lusztig cells";
    const char* mu_tag = "prints out a mu-coefficient";
    const char* pol_tag = "prints out a single Kazhdan-Lusztig polynomial";
    const char* rcells_tag = "prints out the right Kazhdan-Lusztig cells";
    const char* rcorder_tag = "prints the right cell order";
    }; // |namespace uneq|
  }; // |namespace|

  // define the function pointer constant |default_help|, set to |default_h|
  void (*const default_help)() = &help::default_h;

  namespace interf {
    // functions to be installed in the 'interface' mode

    std::unique_ptr<interface::GroupEltInterface> in_buf; // mode variable

    struct In_tag {};
    struct Out_tag {};

    void in_entry();
    void in_exit();
    void out_entry();
    void out_exit();

    void abort_f();
    void alphabetic_f();
    void bourbaki_f();
    void default_f();
    void decimal_f();
    void gap_f();
    void hexadecimal_f();
    void in_f();
    void out_f();
    void permutation_f();
    void symbol_f();
    void terse_f();
    void ordering_f();

    const char* abort_tag = "leaves without modifying the interface";
    const char* alphabetic_tag = "sets alphabetic generator symbols";
    const char* bourbaki_tag = "sets Bourbaki conventions for i/o";
    const char* decimal_tag = "sets decimal generator symbols";
    const char* default_tag = "sets i/o to default mode";
    const char* gap_tag = "sets i/o to GAP mode";
    const char* hexadecimal_tag = "sets hexadecimal generator symbols";
    const char* in_tag = "enters reset-input mode";
    const char* out_tag = "enters reset-output mode";
    const char* permutation_tag =
      "sets permutation notation for i/o (in type A only)";
    const char* ordering_tag = "modifies the ordering of the generators";
    const char* terse_tag = "sets i/o to terse mode";

    namespace in {
      void alphabetic_f();
      void bourbaki_f();
      void decimal_f();
      void default_f();
      void gap_f();
      void hexadecimal_f();
      void permutation_f();
      void postfix_f();
      void prefix_f();
      void separator_f();
      void terse_f();
      const char* alphabetic_tag =
        "sets alphabetic generator symbols for input";
      const char* bourbaki_tag = "sets Bourbaki conventions for input";
      const char* decimal_tag = "sets decimal generator symbols for input";
      const char* default_tag = "sets default conventions for input";
      const char* gap_tag = "sets GAP conventions for input";
      const char* hexadecimal_tag =
        "sets hexadecimal generator symbols for input";
      const char* permutation_tag =
        "sets permutation notation for input (in type A only)";
      const char* postfix_tag = "resets the input postfix";
      const char* prefix_tag = "resets the input prefix";
      const char* separator_tag = "resets the input separator";
      const char* symbol_tag = "resets an input symbol";
      const char* terse_tag = "sets terse conventions for input";
    }; // |namespace in|

    namespace out {
      void alphabetic_f();
      void bourbaki_f();
      void decimal_f();
      void default_f();
      void gap_f();
      void hexadecimal_f();
      void permutation_f();
      void postfix_f();
      void prefix_f();
      void separator_f();
      void terse_f();
      const char* alphabetic_tag =
        "sets alphabetic generator symbols for output";
      const char* bourbaki_tag = "sets Bourbaki conventions for output";
      const char* decimal_tag = "sets decimal generator symbols for output";
      const char* default_tag = "sets default conventions for output";
      const char* gap_tag = "sets GAP conventions for output";
      const char* hexadecimal_tag =
        "sets hexadecimal generator symbols for output";
      const char* permutation_tag =
        "sets permutation notation for output (in type A only)";
      const char* postfix_tag = "resets the output postfix";
      const char* prefix_tag = "resets the output prefix";
      const char* separator_tag = "resets the output separator";
      const char* symbol_tag = "resets an output symbol";
      const char* terse_tag = "sets terse conventions for output";
    }; // |namespace out|
  }; // |namespace interf|


/*****************************************************************************

  This module contains the code for the command interface. Although overall
  I'm happy with the way it works, it suffers from a certain amount of
  clumsiness.

  The idea is that at each point in time, there is a certain active CommandTree
  object. This is basically a dictionary of recognized command names, together
  with the functions that will executed for them; in other words, something that
  should be a map in STL parlance. Actually, the active command tree is the top
  of the command mode stack |mode_stack|; exiting the current mode means popping
  the stack; entering a new mode means pushing it onto the stack.

  Each mode has an associated entry and exit function, which take care of
  initialization and clean-up duties. Actually, there is mostly one main mode;
  the entry function for this is the one which gets type and rank for the user;
  the exit function destroys the current group. Redefining type or rank means
  exiting and re-entering the main mode. In addition, there is the "empty" mode,
  active on startup only, where nothing is defined yet, and some auxiliary modes
  which temporarily hide the main mode in order to perform certain duties :
  interface mode to set the i/o preferences of the user, help mode for help, and
  also unequal-parameter mode which sets unequal parameters for the
  Kazhdan-Lusztig functions; this is in fact a sort of duplicate main mode.

  Command completion is implemented to the extent that incomplete commands are
  recognized when non-ambiguous. This is realized by augmenting the command
  trees with commands for each unique prefix, sharing the pointer to the
  |CommondCell| so that the prefix becomes an alias of the complete command.

 *****************************************************************************/

/*****************************************************************************

        Chapter I -- Running the program

  This section contains the following functions :

  - report_ambiguous_command(str) : what to do with an ambiguous command;
  - mainCommandTree() : returns a pointer to the initial command tree (and
    builds it on first call);
  - relax_f() : does nothing;
  - run() : runs an interactive session;

 *****************************************************************************/

void relax_f() {}

/*
  This function runs an interactive session of the program.
*/
void run()
{
  std::string name;

  activate(*emptyCommandTree());

  if (ERRNO) {     // if already something went wrong
    Error (ERRNO); // report it
    return;        // and quit program
  }

  while (true) { // to exit from this loop either use "qq" command, or EOF
    CommandTree& mode = *mode_stack.top();
    mode.prompt();
    if (not io::getInput(stdin,name))
      break;

    const auto* cell = mode.find_cell(name);
    if (cell == nullptr) {
      mode.call_error(name.c_str());
      continue;
    }
    const auto cd = cell->action();
    if (cd==nullptr) {
      report_ambiguous_command(mode,name);
      continue;
    }

    if (name!=cd->name)
      fprintf(stdout,"%s\n",cd->name.c_str()); // show completed command

    cd->action();

    if (cd!=mode.root()->action()) // whether command differs from that of ""
      mode.set_default_action(cd->autorepeat ? cd : nullptr);
  }
} // |run|


// Default response to an unknown command.
void default_error(const char* str)
{
  Error(COMMAND_NOT_FOUND,str);
}

namespace {


/*
  Puts the tree on top of mode_stack, and executes the initialization function.

  If an error occurs, report just that and return with |MODECHANGE_FAIL|

  Curiously |tree| is already on the stack before |call_entry| is invoked,
  but it probably makes little difference since the main command loot that
  uses the stack is not visited before |call_entry| finishes.
*/
void activate(CommandTree& mode)
{
  mode_stack.push(&mode);
  mode.call_entry();

  if (ERRNO) { /* an error occured during initialization */
    Error(ERRNO);
    mode_stack.pop();
    ERRNO = MODECHANGE_FAIL;
  }
}


/*
  Response to ambiguous commands. Prints a warning and the list of possible
  completions in the current tree on stderr.
*/
void report_ambiguous_command(const CommandTree& mode, const std::string& str)
{
  fprintf(stderr,"%s : ambiguous (",str.c_str());

  bool first = true;
  for (const auto& ext :  mode.find_cell(str)->extensions(str))
    fprintf(stderr,first ? first=false, "%s" : ",%s",ext.c_str());

  fprintf(stderr,")\n");
}

// the error function for the empty mode pushes the main command mode!
void empty_error(const char* str)
{
  static auto* const type_node=mainCommandTree()->find("type");
  static auto* const rank_node=mainCommandTree()->find("rank");

  assert(type_node!=nullptr and rank_node!=nullptr);

  CommandTree& mode = *mainCommandTree();

  const auto* cell = mode.find_cell(str);
  if (cell == nullptr) {
    default_error(str); // report "command not found" after all
    return;
  }

  const auto cd = cell->action();
  if (cd==nullptr) {
    report_ambiguous_command(mode,str);
    return;
  }

  // now that an existing command was unambiguously identified, enter main mode
  activate(mode);
  if (ERRNO) { /* something went wrong during initialization */
    Error(ERRNO);
    return;
  }

  if (str!=cd->name)
    fprintf(stdout,"%s\n",cd->name.c_str()); // show completed command

  // type and rank are already set during |activate(mode)|; skip their action
  if ((cd.get() != type_node) && (cd.get() != rank_node))
    cd->action();

  if (cd!=mode.root()->action()) // whether command differs from that of ""
    mode.set_default_action(cd->autorepeat ? cd : nullptr);
}


/*
  The response to the first carriage return. Sets the response to "help"
  to a less verbose version, and starts up the program.
*/
void startup()
{
  activate(*mainCommandTree());

  if (ERRNO)
    Error(ERRNO);
}

}; // |namespace|

/*****************************************************************************

        Chapter II -- The CommandTree class.

  The purpose of a CommandTree is to get a command-name from the user (or
  perhaps from a file), and execute the corresponding command. For this,
  it maintains a tree of CommandCell s, one for each initial subword of
  each recognized command. Each CommandCell knows which command it should
  execute.

  Recognizing initial subwords allows for command completion : when the
  completion is unique, the command is executed as if the full name were
  typed. When the completion is not unique, the function for ambiguous
  commands is executed: as defined here, it prints the list of all possible
  completions in the current tree, and prompts the user again.

  The case of the empty command is special : either it does nothing, or,
  for most commands, it repeats the previous command.

  This setup supports the concept of mode : at all times, there is a current
  command tree, and in some situations this will change : new commands can
  become available, commands can change behaviour or can become unavailable.
  Help mode is an example of this. Another example is the interface command,
  which loads the interface mode tree, so that the user can set the various
  i/o parameters.

  NOTE : even though I like the actual behaviour of the setup, it is all
  rather clumsy and should be re-done. The command tree could be replaced
  with some associative container like map. The current behaviour could
  be pretty much kept as is, until we include the functionalities of
  readline.

  The following functions are defined :

   - constructors and destructors :

     - CommandTree(prompt,a,entry,error,exit,h) : builds a command tree with
       the given prompt, action a, help function h, given entry and exit
       functions, and error function (called when a command is not found);
     - ~CommandTree();

   - accessors :

     - prompt : prints the prompt;

   - manipulators :

     - add(name,tag,a,h,rep) : adds a command with the given name, tag (used
       in the online help), action a, help-action h and repetition flag rep;
     - setAction(str,a) : resets the action of the command for str;
     - setRepeat(str,b) : resets the repetition flag of the command for b;

******************************************************************************/

void help_filler(CommandTree& tree)
{
  tree.add("q",q_tag,&q_f,0,false);
}

/*
  Initialize a command tree with the given prompt and action |a| for the
  empty command, and call |filler| to construct the remainder of the tree.

  The reason that |filler| is passed as an argument, rather than letting it be
  called by our caller after completing the constructor, is that |CommandTree|
  variables will be |static| inside a function |f|, so that their constructor is
  called only once, but all following code every time that |f| is called.
*/
CommandTree::CommandTree(const char* prompt,
			 void (*a)(),
			 void (*filler)(CommandTree&),
			 void (*entry)(),
			 void (*error)(const char*),
			 void (*exit)(),
			 void (*h)())
  : dictionary::Dictionary<CommandData>(nullptr)
  , d_prompt(prompt), d_entry(entry), d_error(error), d_exit(exit)
{
  if (h!=nullptr) { /* add help functionality */
    d_help = new CommandTree("help",&help::cr_h,&help_filler,h);
    add("help",help_tag,&help_f,&help::help_h,false);
  }
  filler(*this); // add all mode-specific nodes
  install_command_completion();
  if (auto* p=helpMode())
    p->install_command_completion();
}


/*
  The memory allocated by a CommandTree object is hidden in the dictionary
  and in the d_help pointer.
*/

CommandTree::~CommandTree()
{
  delete d_help;
}

/******** accessors *********************************************************/

void CommandTree::prompt() const

/*
  Prints the prompt for the command tree.
*/

{
  printf("%s : ",d_prompt.c_str());
}

/******** manipulators ******************************************************/


/*
  This function adds a new command to the tree, adding new cells as
  necessary.
*/

void CommandTree::add(const char* name, const char* tag, void (*a)(),
		      void (*h)(), bool rep)
{
  auto cd = std::make_shared<CommandData>(name,tag,a,h,rep);

  insert(std::string(name),cd);
  if (d_help!=nullptr && h!=nullptr) { /* add help functionality */
    d_help->add(name,tag,h,0,false);
  }
}


/*****************************************************************************

        Chapter III -- The CommandData class.

  The CommandData structure collects the data associated to a given command
  name. The function a defines the action associated with the command; the
  function h defines the action associated with the command in help mode.
  The flag b is set if the command should be repeated on a carriage return.

 *****************************************************************************/

CommandData::CommandData(const char* str, const char* t,
			 void (*a)(), void (*h)(), bool rep)
  :name(str), tag(t), action(a), help(h), autorepeat(rep)
{
  assert(action!=nullptr);
}


/*
  No memory is allocated directly
*/
CommandData::~CommandData()
{}


/*****************************************************************************

        Chapter IV -- Building the command tree

  This section contains the functions used for the construction of the primary
  command tree, i.e., the initialization of the command module.

  The following functions are defined :

  - commandCompletion(tree) : finishes off the command tree;
  - initCommandTree<Empty_tag> : builds the empty command tree;
  - initCommandTree<Interface_tag> : builds the interface command tree;
  - initCommandTree<Main_tag> : builds the main command tree;
  - initCommandTree<Uneq_tag> : builds the command tree for unequal parameters;
  - emptyCommandTree() : returns a pointer to the initial command tree;
  - interfaceCommandTree() : returns a pointer to the interface command tree;
  - mainCommandTree() : returns a pointer to the main command tree;
  - uneqCommandTree() : returns a pointer to the unequal-parameter command
    tree;

 *****************************************************************************/

namespace {

/*
  This function builds the initial (empty mode) command tree of the program.
  The idea is that all commands on the main command tree will be considered
  entry commands, and so will do the necessary initialization. This is achieved
  through the |empty_error| function, since almost eny command will not be
  recognized in empty mode, so typing it will call the mode's error function,
  which will then install the main command mode.

  Like all |initCommandTree| instances, this function should be called only
  once, to construct the tree held in a |static| variable. This is achieved by
  passing a pointer to it to the constructor for that variable, which construtor
  then calls the function with the object under construction as |tree| argument.
*/
template<> void initCommandTree<Empty_tag>(CommandTree& tree)
{
  tree.add("author","author_tag",&author_f,&relax_f,false);
  tree.add("qq",qq_tag,&qq_f,&help::qq_h,false);

  tree.helpMode()->add("intro",intro_tag,&help::intro_h,0,false);
}


/*
  Return a pointer to the initial command tree of the program, building it on
  the first call.
*/
CommandTree* emptyCommandTree()
{
  static CommandTree empty_tree
    ("empty",&startup,&initCommandTree<Empty_tag>,
     &relax_f,&empty_error,&relax_f,&help::intro_h);
  return &empty_tree;
}


/*
  This function builds the interface command tree; this makes available the
  various little commands that are needed to reset the interface, and that
  have no reason to be clogging up the main command tree.
*/

template<> void initCommandTree<Interface_tag>(CommandTree& tree)
{
  tree.add("alphabetic",commands::interf::alphabetic_tag,
	   &commands::interf::alphabetic_f,&help::interface::alphabetic_h);
  tree.add("bourbaki",commands::interf::bourbaki_tag,
	   &commands::interf::bourbaki_f,&help::interface::bourbaki_h);
  tree.add("decimal",commands::interf::decimal_tag,
	   &commands::interf::decimal_f,&help::interface::decimal_h);
  tree.add("default",commands::interf::default_tag,
	   &commands::interf::default_f,&help::interface::default_h);
  tree.add("gap",commands::interf::out::gap_tag,
	   &commands::interf::out::gap_f, &help::interface::gap_h);
  tree.add("hexadecimal",commands::interf::hexadecimal_tag,
	   &commands::interf::hexadecimal_f,
	   &help::interface::hexadecimal_h);
  tree.add("in",commands::interf::in_tag,&commands::interf::in_f,
	   help::interface::in_h,false);
  tree.add("ordering",commands::interf::ordering_tag,
	   &commands::interf::ordering_f,help::interface::ordering_h,false);
  tree.add("out",commands::interf::out_tag,&commands::interf::out_f,
	   help::interface::out_h,false);
  tree.add("permutation",commands::interf::permutation_tag,
	   &commands::interf::permutation_f,
	   &help::interface::permutation_h);
  tree.add("q",q_tag,&q_f,0,false);
  tree.add("terse",commands::interf::out::terse_tag,
	   &commands::interf::out::terse_f, &help::interface::out::terse_h);
}


/*
  This function builds the command tree for the input-modification mode.
*/
template<> void initCommandTree<commands::interf::In_tag>(CommandTree& tree)
{
  using namespace commands::interf;

  tree.add("q",q_tag,&q_f,0,false);

  tree.add("abort",abort_tag,&abort_f,&help::interface::abort_h);
  tree.add("alphabetic",in::alphabetic_tag,&in::alphabetic_f,
	   &help::interface::in::alphabetic_h,false);
  tree.add("bourbaki",in::bourbaki_tag,&in::bourbaki_f,
	   &help::interface::in::bourbaki_h);
  tree.add("decimal",in::decimal_tag,&in::decimal_f,
	   &help::interface::in::decimal_h,false);
  tree.add("default",in::default_tag,&in::default_f,
	   &help::interface::in::default_h);
  tree.add("gap",in::gap_tag,&in::gap_f,&help::interface::in::gap_h);
  tree.add("hexadecimal",in::hexadecimal_tag,&in::hexadecimal_f,
	   &help::interface::in::hexadecimal_h,false);
  tree.add("permutation",in::permutation_tag,&in::permutation_f,
	   &help::interface::in::permutation_h,false);
  tree.add("postfix",in::postfix_tag,&in::postfix_f,
	   &help::interface::in::postfix_h);
  tree.add("prefix",in::prefix_tag,&in::prefix_f,
	   &help::interface::in::prefix_h);
  tree.add("separator",in::separator_tag,
	   &in::separator_f,&help::interface::in::separator_h);
  tree.add("symbol",in::symbol_tag,&symbol_f,
	   &help::interface::in::symbol_h);
  tree.add("terse",in::terse_tag,&in::terse_f,&help::interface::in::terse_h);
}


/*
  This function builds the command tree for the output-modification mode.
*/

template<> void initCommandTree<commands::interf::Out_tag>(CommandTree& tree)
{
  using namespace commands::interf;

  tree.add("q",q_tag,&q_f,0,false);

  tree.add("alphabetic",out::alphabetic_tag,&out::alphabetic_f,
	   &help::interface::out::alphabetic_h,false);
  tree.add("bourbaki",out::bourbaki_tag,&out::bourbaki_f,
	   &help::interface::out::bourbaki_h);
  tree.add("decimal",out::decimal_tag,&out::decimal_f,
	   &help::interface::out::decimal_h,false);
  tree.add("default",out::default_tag,&out::default_f,
	   &help::interface::out::default_h);
  tree.add("gap",out::gap_tag,&out::gap_f,
	   &help::interface::out::gap_h);
  tree.add("hexadecimal",out::hexadecimal_tag,&out::hexadecimal_f,
	   &help::interface::out::hexadecimal_h,false);
  tree.add("permutation",out::permutation_tag,&out::permutation_f,
	   &help::interface::out::permutation_h,false);
  tree.add("postfix",out::postfix_tag,&out::postfix_f,
	   &help::interface::out::postfix_h);
  tree.add("prefix",out::prefix_tag,&out::prefix_f,
	   &help::interface::out::prefix_h);
  tree.add("separator",out::separator_tag,
	   &out::separator_f,&help::interface::out::separator_h);
  tree.add("symbol",out::symbol_tag,&symbol_f,
	   &help::interface::out::symbol_h);
  tree.add("terse",out::terse_tag,&out::terse_f,
	   &help::interface::out::terse_h);
}

}; // |namespace|

/*
  Return a pointer to the interface command tree, building it on the first
  call.
*/

CommandTree* interfaceCommandTree()
{
  static CommandTree interface_tree
    ("interface",&relax_f,&initCommandTree<Interface_tag>,
     &interface_entry,&default_error,&interface_exit,&help::interface_help);
  return &interface_tree;
}

CommandTree* interf::inCommandTree()
{
  static CommandTree in_tree
    ("in",&relax_f,&initCommandTree<In_tag>,
     &in_entry,&default_error,&in_exit,&help::interface::in_help);
  return &in_tree;
}

CommandTree* interf::outCommandTree()
{
  static CommandTree out_tree
    ("out",&relax_f,&initCommandTree<Out_tag>,
     &out_entry,&default_error,
     &out_exit,&help::interface::out_help);
  return &out_tree;
}


namespace {


/*
  This function builds the main command tree, the one that is being run on
  startup. Auxiliary trees may be grafted onto this one (thru the pushdown
  stack mode_stack) by some functions needing to be in special modes.
*/
template<> void initCommandTree<Main_tag>(CommandTree& tree)
{
  using namespace help;
  tree.add("author",author_tag,&author_f,&relax_f,false);
  tree.add("betti",betti_tag,&betti_f,&betti_h,false);
  tree.add("coatoms",coatoms_tag,&coatoms_f,&coatoms_h);
  tree.add("compute",compute_tag,&compute_f,&compute_h);
  tree.add("descent",descent_tag,&descent_f,&descent_h);
  tree.add("duflo",duflo_tag,&duflo_f,&duflo_h);
  tree.add("extremals",extremals_tag,&extremals_f,&extremals_h);
  tree.add("fullcontext",fullcontext_tag,&fullcontext_f,&fullcontext_h);
  tree.add("ihbetti",ihbetti_tag,&ihbetti_f,&ihbetti_h,false);
  tree.add("interface",interface_tag,&interface_f,&interface_h,false);
  tree.add("interval",interval_tag,&interval_f,&interval_h,false);
  tree.add("inorder",inorder_tag,&inorder_f,&inorder_h);
  tree.add("invpol",invpol_tag,&invpol_f,&invpol_h);
  tree.add("lcorder",lcorder_tag,&lcorder_f,&lcorder_h,false);
  tree.add("lcells",lcells_tag,&lcells_f,&lcells_h,false);
  tree.add("lcwgraphs",lcwgraphs_tag,&lcwgraphs_f,&lcwgraphs_h,false);
  tree.add("lrcorder",lrcorder_tag,&lrcorder_f,&lrcorder_h,false);
  tree.add("lrcells",lrcells_tag,&lrcells_f,&lrcells_h,false);
  tree.add("lrcwgraphs",lrcwgraphs_tag,&lrcwgraphs_f,&lrcwgraphs_h,false);
  tree.add("lrwgraph",lrwgraph_tag,&lrwgraph_f,&lrwgraph_h,false);
  tree.add("lwgraph",lwgraph_tag,&lwgraph_f,&lwgraph_h,false);
  tree.add("klbasis",klbasis_tag,&klbasis_f,&klbasis_h,true);
  tree.add("matrix",matrix_tag,&matrix_f,&matrix_h);
  tree.add("mu",mu_tag,&mu_f,&mu_h);
  tree.add("pol",pol_tag,&pol_f,&pol_h);
  tree.add("q",q_tag,&q_f,0,false);
  tree.add("qq",qq_tag,&qq_f,&qq_h,false);
  tree.add("rank",rank_tag,&rank_f,&rank_h,false);
  tree.add("rcorder",rcorder_tag,&rcorder_f,&rcorder_h,false);
  tree.add("rcells",rcells_tag,&rcells_f,&rcells_h,false);
  tree.add("rcwgraphs",rcwgraphs_tag,&rcwgraphs_f,&rcwgraphs_h,false);
  tree.add("rwgraph",rwgraph_tag,&rwgraph_f,&rwgraph_h,false);
  tree.add("schubert",schubert_tag,&schubert_f,&schubert_h);
  tree.add("show",show_tag,&show_f,&show_h);
  tree.add("showmu",showmu_tag,&showmu_f,&showmu_h);
  tree.add("slocus",slocus_tag,&slocus_f,&slocus_h);
  tree.add("sstratification",sstratification_tag,&sstratification_f,
	   &sstratification_h);
  tree.add("type",type_tag,&type_f,&type_h,false);
  tree.add("uneq",uneq_tag,&uneq_f,&uneq_h,false);

  special::addSpecialCommands(&tree);

  tree.helpMode()->add("intro",intro_tag,&intro_h,0,false);
  tree.helpMode()->add("input",input_tag,&input_h,0,false);
}

}; // |namespace|


/*
  Return a pointer to the main command tree of the program, building it on
  the first call.
*/

CommandTree* mainCommandTree()
{
  static CommandTree main_tree
    ("coxeter",&relax_f,&initCommandTree<Main_tag>,
     &main_entry,&default_error, &main_exit,&help::main_help);
  return &main_tree;
}

namespace {


/*
  This function builds the unequal-parameter command tree. It contains
  essentially the same functions as the main command tree, except that the
  unequal-parameter versions have been substituted for the Kazhdan-Lusztig
  functions.
*/
template<> void initCommandTree<Uneq_tag>(CommandTree& tree)
{
  using namespace help;
  tree.add("author",author_tag,&author_f,&relax_f,false);
  tree.add("coatoms",coatoms_tag,&coatoms_f,&coatoms_h);
  tree.add("compute",compute_tag,&compute_f,&compute_h);
  tree.add("descent",descent_tag,&descent_f,&descent_h);
  tree.add("fullcontext",fullcontext_tag,&fullcontext_f,&fullcontext_h);
  tree.add("interface",interface_tag,&interface_f,&interface_h,false);
  tree.add("klbasis",klbasis_tag,&uneq::klbasis_f,&help::uneq::klbasis_h,true);
  tree.add("lcorder",uneq::lcorder_tag,&uneq::lcorder_f,
	   &help::uneq::lcorder_h,false);
  tree.add("lrcorder",uneq::lrcorder_tag,&uneq::lrcorder_f,
	   &help::uneq::lrcorder_h,false);
  tree.add("lcells",uneq::lcells_tag,&uneq::lcells_f,&help::uneq::lcells_h,
	   false);
  tree.add("lrcells",uneq::lrcells_tag,&uneq::lrcells_f,&help::uneq::lrcells_h,
	   false);
  tree.add("matrix",matrix_tag,&matrix_f,&matrix_h);
  tree.add("mu",uneq::mu_tag,&uneq::mu_f,&help::uneq::mu_h);
  tree.add("pol",uneq::pol_tag,&uneq::pol_f,&help::uneq::pol_h);
  tree.add("rcells",uneq::rcells_tag,&uneq::rcells_f,&help::uneq::rcells_h,
	   false);
  tree.add("rcorder",uneq::rcorder_tag,&uneq::rcorder_f,
	   &help::uneq::rcorder_h,false);
  tree.add("q",q_tag,&q_f,0,false);
  tree.add("qq",qq_tag,&qq_f,&qq_h,false);
}

}; // |namespace|


/*
  Return a pointer to the uneq command tree, building it on the first call.
*/
CommandTree* uneqCommandTree()
{
  static CommandTree uneq_tree
    ("uneq",&relax_f,&initCommandTree<Uneq_tag>,
     &uneq_entry,&default_error,&uneq_exit,&help::uneq_help);
  return &uneq_tree;
}


/*****************************************************************************

        Chapter V -- Functions for the predefined commands.

  This section contains the functions defining the responses to the various
  commands which are provided by the program. The functions are placed in
  the unnamed namespace defined in this file.

  The following functions are defined :

  - author_f() : response to "author";
  - betti_f() : response to "betti";
  - coatoms_f() : response to "coatoms";
  - compute_f() : response to "compute";
  - descent_f() : response to "descent";
  - duflo_f() : response to "duflo";
  - extremals_f() : response to "extremals";
  - fullcontext_f() : response to "fullcontext";
  - ihbetti_f() : response to "ihbetti";
  - inorder_f() : response to "inorder";
  - interface_f() : response to "interface";
  - help_f() : repsonse to "help";
  - klbasis_f() : response to "klbasis";
  - lcorder_f() : response to "lcorder";
  - lcells_f() : response to "lcells";
  - lcwgraphs_f() : response to "lcwgraphs";
  - lrcorder_f() : response to "lrcorder";
  - lrcells_f() : response to "lrcells";
  - lrcwgraphs_f() : response to "lrcwgraphs";
  - lrwgraph_f() : response to "lrwgraph";
  - lwgraph_f() : response to "lwgraph";
  - matrix_f() : response to "matrix";
  - mu_f() : response to "mu";
  - not_implemented_f() : response for not (yet) implemented features;
  - q_f() : response to "q";
  - qq_f() : response to "qq";
  - rank_f() : response to "rank";
  - rcorder_f() : response to "rcorder";
  - rcells_f() : response to "rcells";
  - rcwgraphs_f() : response to "rcwgraphs";
  - rwgraph_f() : response to "rwgraph";
  - relax_f() : does nothing;
  - schubert_f() : response to "schubert";
  - show_f() : response to "show";
  - showmu_f() : response to "showmu";
  - slocus_f() : response to "slocus";
  - sstratification_f() : response to "sstratification";
  - type_f() : response to "type";
  - uneq_f() : response to "uneq";

  In uneq mode :

  - klbasis_f() : response to "klbasis";
  - lcorder_f() : response to "lcorder";
  - lcells_f() : response to "lcells";
  - lrcorder_f() : response to "lrcorder";
  - lrcells_f() : response to "lrcells";
  - mu_f() : response to "mu";
  - pol_f() : response to "pol";
  - rcells_f() : response to "rcells";
  - rcorder_f() : response to "rcorder";

  In interface mode :

  - abort_f() : aborts input interface modification;
  - alphabetic_f() : sets alphabetic generator symbols;
  - bourbaki_f() : sets bourbaki conventions;
  - decimal_f() : sets decimal generator symbols;
  - default_f() : sets default i/o;
  - gap_f() : sets GAP-style i/o;
  - hexadecimal_f() : sets hexadecimal generator symbols;
  - ordering_f() : changes the ordering of the generators;
  - symbol_f() : resets generator symbol;
  - terse_f() : sets terse style i/o;
  - in::alphabetic_f() : sets alphabetic generator symbols;
  - in::bourbaki_f() : sets bourbaki conventions;
  - in::decimal_f() : sets decimal generator symbols;
  - in::default_f() : sets default-style input;
  - in::gap_f() : sets GAP-style input;
  - in::hexadecimal_f() : sets hexadecimal generator symbols;
  - in::permutation_f() : sets permutation input;
  - in::postfix_f() : resets input postfix;
  - in::prefix_f() : resets input prefix;
  - in::separator_f() : resets input separator;
  - in::terse_f() : sets terse-style input;
  - out::alphabetic_f() : sets alphabetic generator symbols;
  - out::bourbaki_f() : sets bourbaki conventions;
  - out::decimal_f() : sets decimal generator symbols;
  - out::default_f() : sets default-style output;
  - out::gap_f() : sets GAP-style output;
  - out::hexadecimal_f() : sets hexadecimal generator symbols;
  - out::permutation_f() : sets permutation output;
  - out::postfix_f() : resets output postfix;
  - out::prefix_f() : resets output prefix;
  - out::separator_f() : resest output separator;
  - out::terse_f() : sets terse-style output;

 *****************************************************************************/

namespace {


// Print a message about the author.
void author_f()
{
  io::printFile(stderr,"author.mess",directories::MESSAGE_DIR);
}

void betti_f()

/*
  Prints out the ordinary betti numbers of [e,y].

  NOTE : could be *much* improved! In particular, we would want to have
  betti(x,y) for x <= y.
*/

{
  static coxtypes::CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  files::OutputTraits& traits = W->outputTraits();
  printBetti(stdout,y,W->schubert(),traits);
}

void coatoms_f()

/*
  Prints out the coatoms of a given element, computing them in elementary
  fashion.
*/

{
  static coxtypes::CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  list::List<coxtypes::CoxWord> c(0);
  W->coatoms(c,g);

  for (Ulong j = 0; j < c.size(); ++j) {
    W->print(stdout,c[j]);
    printf("\n");
  }
}


//  Get an element from the user, and prints out its normal form.
void compute_f()
{
  printf("enter your element (finish with a carriage return) :\n");
  coxtypes::Cox_word g = interactive::getCox_word(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  W->to_normal_form(g);
  W->print(stdout,g);
  if (auto* Ws = dynamic_cast<fcoxgroup::SmallCoxGroup*> (W)) {
    coxtypes::CoxNbr x = 0;
    Ws->prodD(x,g);
    printf(" (#%lu)",static_cast<Ulong>(x));
  }
  coxtypes::CoxNbr x = W->context_number(g);
  if (x != coxtypes::undef_coxnbr)
    printf(" (%s%lu)","%",static_cast<Ulong>(x));
  printf("\n");
}

void descent_f()

/*
  Prints the left and right descent sets.
*/

{
  static coxtypes::CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  GenSet f = W->ldescent(g);
  printf("L:");
  W->printFlags(stdout,f);
  printf("; R:");
  f = W->rdescent(g);
  W->printFlags(stdout,f);
  printf("\n");
}

void duflo_f()

/*
  Prints the Duflo involutions. Works for finite groups only.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"duflo.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::dufloH,traits);
  printDuflo(file.f(),Wf->duflo(),Wf->cell<'l'>(),Wf->kl(),
	     W->interface(),traits);
}

void extremals_f()

/*
  Prints out the list of extremal elements x <= y (this is part of the schubert
  command.
*/

{
  static coxtypes::CoxWord g(0);

  printf("Enter your element (finish with a carriage-return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),files::extremalsH,traits);
  printExtremals(file.f(),y,W->kl(),W->interface(),traits);
}

void fullcontext_f()

/*
  Response to the fullcontext command. This sets the context to the full
  group. Of course, it works only for finite groups.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"fullcontext.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
  }
}


// Response to the help command.
void help_f()
{
  activate(*mode_stack.top()->helpMode());
}

void ihbetti_f()

/*
  Prints out the IH betti numbers of [e,y].
*/

{
  static coxtypes::CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  files::OutputTraits& traits = W->outputTraits();
  printIHBetti(stdout,y,W->kl(),traits);
}

void interface_f()

/*
  Response to the interface command.
*/

{
  activate(*interfaceCommandTree());
}

void interval_f()
{

  fprintf(stdout,"first : ");
  coxtypes::Cox_word g = interactive::getCox_word(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  fprintf(stdout,"second : ");
  coxtypes::Cox_word h = interactive::getCox_word(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (not W->Bruhat_leq(g,h))
  {
    fprintf(stderr,"the two elements are not in order\n");
    return;
  }

  W->extend_context(h);

  coxtypes::CoxNbr x = W->context_number(g);
  coxtypes::CoxNbr y = W->context_number(h);

  interactive::OutputFile file;

  bitmap::BitMap b = W->closure(y);

  schubert::CoxNbrList res;

  for (size_t z = y+1; b.back_up(z); ) // |b| pruned in reverse iteration
    if (W->inOrder(x,z))
      res.push_back(z);
    else
      b.andnot(W->closure(z));

  schubert::NFCompare nfc(W->schubert(),W->ordering());
  bits::Permutation a = bits::inverse_standardization(res,nfc);

  for (size_t j = 0; j < res.size(); ++j)
  {
    W->print(file.f(),res[a[j]]);
    fprintf(file.f(),"\n");
  }
}

void inorder_f()

/*
  Response to the inorder command. This will tell whether two elements
  are comparable in Bruhat order, using only the elementary string operations
  (and hence not consuming any memory.)
*/

{
  coxtypes::CoxWord g(0);
  coxtypes::CoxWord h(0);
  list::List<coxtypes::Length> a(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  fprintf(stdout,"second : ");
  h = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (W->inOrder(a,g,h)) {
    fprintf(stdout,"true :   ");
    Ulong i = 0;
    for (Ulong j = 0; j < a.size(); ++j) {
      while (i < a[j]) {
	W->printSymbol(stdout,h[i]-1);
	++i;
      }
      fprintf(stdout,".");
      ++i;
    }
    while (i < h.length()) {
      W->printSymbol(stdout,h[i]-1);
      ++i;
    }
    fprintf(stdout,"\n");
  }
  else
    fprintf(stdout,"false\n");
}

void invpol_f()

/*
  Response to the invpol command. This prints out a single inverse
  Kazhdan-Lusztig polynomial, without details.
*/

{
  coxtypes::CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  const kl::KLPol& pol = W->invklPol(x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  print(stdout,pol,"q");
  printf("\n");
}


/*
  Print out one element in the Kazhdan-Lusztig basis of the group, in the
  format defined by the current output mode.
*/
void klbasis_f()
{
  coxtypes::CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  kl::HeckeElt h;

  W->cBasis(h,y);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),files::basisH,traits);
  files::printAsBasisElt(file.f(),h,W->schubert(),W->interface(),traits);
}

void lcorder_f()

/*
  Prints the left cell order of the current group in a file. Works only for
  finite groups.
*/

{
  if (!fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"lcorder.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lCOrderH,traits);
  printLCOrder(file.f(),Wf->kl(),Wf->interface(),traits);
}

void lcells_f()

/*
  This function prints out the left cells in the group.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"lcells.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lCellsH,traits);
  printLCells(file.f(),Wf->cell<'l'>(),Wf->kl(),Wf->interface(),traits);
}

void lcwgraphs_f()

/*
  This function prints out the W-graphs of the left cells in the group.
  It works only for finite groups currently.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"lcells.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lCellWGraphsH,traits);
  printLCellWGraphs(file.f(),Wf->cell<'l'>(),Wf->kl(),W->interface(),traits);
}

void lrcorder_f()

/*
  Prints the two-sided cell order of the current group in a file. Works only
  for finite groups.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"lrcorder.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lrCOrderH,traits);
  printLRCOrder(file.f(),Wf->kl(),Wf->interface(),traits);
}

void lrcells_f()

/*
  This function prints out the two-sided cells in the group, together with the
  corresponding W-graphs. Works only for finite groups.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"lrcells.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lrCellsH,traits);
  printLRCells(file.f(),Wf->cell<'b'>(),Wf->kl(),Wf->interface(),traits);
}

void lrcwgraphs_f()

/*
  This function prints out the W-graphs of the two-sided cells in the group.
  It works only for finite groups currently.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"lcells.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lrCellWGraphsH,traits);
  printLRCellWGraphs(file.f(),Wf->cell<'b'>(),Wf->kl(),W->interface(),traits);
}

void lrwgraph_f()

/*
  Prints out the two-sided wgraph of the current context.
*/

{
  if (!W->isFullContext() && wgraph_warning) {
    io::printFile(stderr,"wgraph.mess",directories::MESSAGE_DIR);
    printf("continue ? y/n\n");
    if (not interactive::yesNo())
      return;
    printf("print this message next time ? y/n\n");
    if (not interactive::yesNo())
      wgraph_warning = false;
  }

  W->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),files::lrWGraphH,traits);
  printLRWGraph(file.f(),W->kl(),W->interface(),traits);
}

void lwgraph_f()

/*
  Prints out the left wgraph of the current context.
*/

{
  if (!W->isFullContext() && wgraph_warning) {
    io::printFile(stderr,"wgraph.mess",directories::MESSAGE_DIR);
    printf("continue ? y/n\n");
    if (not interactive::yesNo())
      return;
    printf("print this message next time ? y/n\n");
    if (not interactive::yesNo())
      wgraph_warning = false;
  }

  W->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),files::lWGraphH,traits);
  printLWGraph(file.f(),W->kl(),W->interface(),traits);
}

void matrix_f()

/*
  Prints the Coxeter matrix.
*/

{
  interactive::printMatrix(stdout,W);
}

void mu_f()

/*
  Response to the mu command. This prints out a single mu-coefficient,
  without details.
*/

{
  static coxtypes::CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  coxtypes::CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  klsupport::KLCoeff mu = W->mu(x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  printf("%lu\n",static_cast<Ulong>(mu));
}

void pol_f()

/*
  Response to the pol command. This prints out a single polynomial, without
  details.
*/

{
  static coxtypes::CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  const kl::KLPol& pol = W->klPol(x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  print(stdout,pol,"q");
  printf("\n");
}

void q_f()

/*
  Exits the current mode. If there is a problem on exit, the exit function
  has the option of setting an error, thus preventing the exit.
*/

{
  const CommandTree& mode = *mode_stack.top();
  mode.call_exit();

  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  mode_stack.pop();
}


// Exit the program after calling exit function of all active command modes
void qq_f()
{
  while(not mode_stack.empty()) {
    const CommandTree& mode = *mode_stack.top();
    mode.call_exit();
    mode_stack.pop();
  }

  exit(0);
}

void rank_f()

/*
  Sets the rank.
*/

{
  coxgroup::CoxGroup* Wloc = interactive::allocCoxGroup(W->type());

  if (ERRNO) {
    Error(ERRNO);
  }
  else {
    W = Wloc;
  }
}

void rcorder_f()

/*
  Prints the right cell order of the current group in a file. Works only for
  finite groups.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"rcorder.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::rCOrderH,traits);
  printRCOrder(file.f(),Wf->kl(),Wf->interface(),traits);
}

void rcells_f()

/*
  This function prints out the right cells in the group, together with the
  corresponding W-graphs.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"rcells.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::rCellsH,traits);
  printRCells(file.f(),Wf->cell<'r'>(),Wf->kl(),Wf->interface(),traits);
}

void rcwgraphs_f()

/*
  This function prints out the W-graphs of the right cells in the group.
  It works only for finite groups currently.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"lcells.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::rCellWGraphsH,traits);
  printRCellWGraphs(file.f(),Wf->cell<'r'>(),Wf->kl(),W->interface(),traits);
}

void rwgraph_f()

/*
  Prints out the right wgraph of the current context.
*/

{
  if (!W->isFullContext() && wgraph_warning) {
    io::printFile(stderr,"wgraph.mess",directories::MESSAGE_DIR);
    printf("continue ? y/n\n");
    if (not interactive::yesNo())
      return;
    printf("print this message next time ? y/n\n");
    if (not interactive::yesNo())
      wgraph_warning = false;
  }

  W->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),files::rWGraphH,traits);
  printRWGraph(file.f(),W->kl(),W->interface(),traits);
}

void schubert_f()

/*
  Response to the schubert command. This will print out the information
  corresponding to one element in the Kazhdan-Lusztig basis, and the information
  on the singularities of the corresponding Schubert variety, in the format
  popularized by Goresky.
*/

{
  static coxtypes::CoxWord g(0);

  printf("Enter your element (finish with a carriage-return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),files::closureH,traits);
  printClosure(file.f(),y,W->kl(),W->interface(),traits);
}

void show_f()

/*
  Response to the show command. This maps out the computation of a Kazhdan-
  Lusztig polynomial, letting you choose the descent generator.

  If no generator is given, the default descent generator is used.
*/

{
  static coxtypes::CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  fprintf(stdout,"generator (carriage return for default) : ");
  Lflags f = W->descent(y);
  coxtypes::Generator s = interactive::getGenerator(W,f);
  if (ERRNO) {
    Error (ERRNO);
    return;
  }

  interactive::OutputFile file;
  showKLPol(file.f(),W->kl(),x,y,W->interface(),s);
}

void showmu_f()

/*
  Response to the showmu command. This maps out the computation of a
  mu-coefficient, letting you choose the descent generator.

  If no generator is given, the default descent generator is used.
*/

{
  static coxtypes::CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  interactive::OutputFile file;
  showMu(file.f(),W->kl(),x,y,W->interface());
}

void slocus_f ()

/*
  Response to the slocus command. Prints out the singular locus of the
  Schubert variety cl(X_y).
*/

{
  static coxtypes::CoxWord g(0);

  printf("Enter your element (finish with a carriage-return) :\n");

  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),files::slocusH,traits);
  printSingularLocus(file.f(),y,W->kl(),W->interface(),traits);
}


/*
  Response to the 'slocus' command. Prints out the singular locus of the
  Schubert variety cl(X_y).
*/
void sstratification_f ()
{
  static coxtypes::CoxWord g(0);

  printf("Enter your element (finish with a carriage-return) :\n");

  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),files::sstratificationH,traits);
  printSingularStratification(file.f(),y,W->kl(),W->interface(),traits);
}

void type_f()

/*
  This function sets the type of W, i.e., gets a type and rank from the
  user and sets W to a new group of that type and rank.
*/

{
  coxgroup::CoxGroup* Wloc = interactive::allocCoxGroup();

  if (ERRNO) {
    Error(ERRNO);
  }
  else {
    delete W;
    wgraph_warning = true;
    W = Wloc;
  }
}

void uneq_f()

/*
  Response to the uneq command.
*/

{
  activate(*uneqCommandTree());
}

namespace uneq {

void klbasis_f()

/*
  Prints out one element in the Kazhdan-Lusztig basis of the group, in the
  format defined by the current output mode.
*/

{
  coxtypes::CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  uneqkl::HeckeElt h(0);

  W->uneqcBasis(h,y);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),files::basisH,traits);
  files::printAsBasisElt(file.f(),h,W->schubert(),W->interface(),traits);
}


// Print out the left cells in the group.
void lcells_f()
{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"lcells.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lCellsH,traits);
  printLCells(file.f(),Wf->uneq_cell<'l'>(),Wf->uneqkl(),Wf->interface(),traits);
}

void lcorder_f()

/*
  Prints the left cell order of the closure of the current group in a file.
  Works only for finite groups.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"lcorder.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lCOrderH,traits);
  printLCOrder(file.f(),Wf->uneqkl(),Wf->interface(),traits);
}

void lrcorder_f()

/*
  Prints the two-sided cell order of the closure of the current group in a
  file. Works only for finite groups.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"uneq/lrcorder.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lrCOrderH,traits);
  printLRCOrder(file.f(),Wf->uneqkl(),Wf->interface(),traits);
}

void lrcells_f()

/*
  This function prints out the two-sided cells in the group. Works only for
  finite groups.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"uneq/lrcells.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::lrCellsH,traits);
  printLRCells(file.f(),Wf->uneq_cell<'b'>(),Wf->uneqkl(),
	       Wf->interface(),traits);
}

void mu_f()

/*
  Response to the "mu" command. Prints out a single mu-coefficient. When
  the generator s act on the left, we use the fact that mu(left_s,x,y) is
  equal to mu(right_s,x^{-1},y^{-1}) to go over to the right action. This
  is necessary because the mu-tables are kept only for right actions.
*/

{
  static coxtypes::CoxWord g(0);
  bool leftAction = false;

  fprintf(stdout,"generator : ");
  coxtypes::Generator s = interactive::getGenerator(W);

  if (s >= W->rank()) { // action is on the left
    s -= W->rank();
    leftAction = true;
  }

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (leftAction)
    W->inverse(g);
  if (!W->isDescent(g,s)) { // mu(s,x,y) is undefined
    fprintf(stderr,"xs is greater than x\n");
    return;
  }
  coxtypes::CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (leftAction)
    W->inverse(g);
  if (W->isDescent(g,s)) { // mu(s,x,y) is undefined
    fprintf(stderr,"ys is smaller than y\n");
    return;
  }
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (x == y) {
    fprintf(stderr,"the two elements are equal\n");
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  const uneqkl::MuPol mu = W->uneqmu(s,x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  print(stdout,mu,"v");
  printf("\n");
}

void pol_f()

/*
  Response to the "pol" command. Prints out a single polynomial.
*/

{
  static coxtypes::CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  coxtypes::CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  const uneqkl::KLPol& pol = W->uneqklPol(x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  print(stdout,pol,"q");
  printf("\n");
}

void rcells_f()

/*
  This function prints out the left cells in the group. Works only for finite
  groups.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"rcells.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::rCellsH,traits);
  printRCells(file.f(),Wf->uneq_cell<'r'>(),Wf->uneqkl(),Wf->interface(),traits);
}

void rcorder_f()

/*
  Prints the right cell order of the current group in a file. Works only for
  finite groups.
*/

{
  if (not fcoxgroup::isFiniteType(W)) {
    io::printFile(stderr,"rcorder.mess",directories::MESSAGE_DIR);
    return;
  }

  auto* Wf = dynamic_cast<fcoxgroup::FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  files::OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),files::rCOrderH,traits);
  printRCOrder(file.f(),Wf->uneqkl(),Wf->interface(),traits);
}

}; // |namespace uneq|

}; // |namespace|


/*
  Abort the interface modification. Bypasses the exit function,
  so that there is no further checking of the abandoned choices.
*/
void interf::abort_f()
{
  in_buf.reset();
  mode_stack.pop();
}


/*
  Sets i/o to the standard alphabetic conventions : the symbols are the
  alphabetic sequence, prefix and postfix are empty, and the separator is
  "." whenever the rank is > 26.
*/
void interf::alphabetic_f()
{
  in_buf.reset
    (new interface::GroupEltInterface(W->rank(),interface::Alphabetic()));
  W->interface().setIn(*in_buf);
  W->interface().setOut(*in_buf);
}


/*
  Set Bourbaki conventions. This means that the ordering is reversed in
  types B and D, and symbols as well.

  NOTE : currently not implemented for affine groups.
*/
void interf::bourbaki_f()
{
  in_buf.reset
    (new interface::GroupEltInterface(W->interface().inInterface()));
  in::bourbaki_f();
  W->interface().setIn(*in_buf);

  in_buf.reset
    (new interface::GroupEltInterface(W->interface().outInterface()));
  out::bourbaki_f();
  W->interface().setOut(*in_buf);
}


/*
  Set i/o to the standard decimal conventions : the symbols are the
  decimal sequence, prefix and postfix are empty, and the separator is
  "." whenever the rank is > 9.
*/
void interf::decimal_f()
{
  in_buf.reset
    (new interface::GroupEltInterface(W->rank(),interface::Decimal()));
  W->interface().setIn(*in_buf);
  W->interface().setOut(*in_buf);
}


/*
  Sets i/o settings to the default style. This means that we use decimal
  symbols, no prefix or postfix, and separator "." only when rank is >= 10.
  The ordering is the internal default ordering.
*/
void interf::default_f()
{
  in_buf.reset(new interface::GroupEltInterface(W->rank()));

  W->interface().setIn(*in_buf);
  W->interface().setOut(*in_buf);

  W->interface().setOrder(bits::Permutation(W->rank()));
  W->interface().setDescent(io::Default());
  W->setOutputTraits(io::Pretty());
}


/*
  Set i/o settings to GAP style. This means first of all that Bourbaki
  conventions are adopted; decimal symbols are used for i/o with prefix
  "[", separator "," and postfix "]". Furthermore, output to files is
  done in GAP style, which produces files that are directly legible by
  GAP3.
*/
void interf::gap_f()
{
  in_buf.reset(new interface::GroupEltInterface(W->rank(),io::GAP()));

  in::bourbaki_f();
  W->interface().setIn(*in_buf);

  out::bourbaki_f();
  W->interface().setOut(*in_buf);
  W->interface().setDescent(io::GAP());
  W->setOutputTraits(io::GAP());
}


/*
  Sets i/o to the standard hexadecimal conventions : the symbols are the
  hexadecimal sequence, prefix and postfix are empty, and the separator is
  "." iff the rank is > 15.
*/
void interf::hexadecimal_f()
{
  in_buf.reset
    (new interface::GroupEltInterface(W->rank(),interface::Hexadecimal()));
  W->interface().setIn(*in_buf);

  W->interface().setOut(*in_buf);
}

void interf::in_f()

{
  activate(*inCommandTree());
}

void interf::ordering_f()

/*
  Allows the user to change the generator ordering.
*/

{
  static bits::Permutation in_order(W->rank());

  interactive::changeOrdering(W,in_order);

  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  W->setOrdering(in_order);
}

void interf::out_f()

{
  activate(*outCommandTree());
}

void interf::permutation_f()

/*
  Activates permutation i/o.
*/

{
  using namespace typeA;

  if (!isTypeA(W->type())) {
    io::printFile(stderr,"permutation.mess",directories::MESSAGE_DIR);
    return;
  }

  TypeACoxGroup* WA = dynamic_cast<TypeACoxGroup*>(W);

  WA->setPermutationInput(true);
  WA->setPermutationOutput(true);

  W->interface().setOrder(bits::Permutation(W->rank()));
  W->interface().setDescent(io::Default());
  W->setOutputTraits(io::Pretty());
}


/*
  Reset a symbol in |in_buf| (this will become either an input or an output
  symbol).
*/
void interf::symbol_f()
{
  static std::string buf;

  const interface::Interface& I = W->interface();
  coxtypes::Generator s = coxtypes::undef_generator;
  buf.clear();

  do {
    if (ERRNO)
      Error(ERRNO);
    printf("enter the generator symbol you wish to change, ? to abort:\n");
    if (not io::getInput(stdin,buf,0) or buf[0] == '?')
      return;
    io::skipSpaces(buf,0);
    interface::Token tok;
    I.symbolTree().find(buf,0,tok);
    if (interface::tokenType(tok) != interface::generator_type)/* error */
      ERRNO = NOT_GENERATOR;
    else
      s = tok-1;
  } while (ERRNO);

  printf("enter the new symbol (finish with a carriage return):\n");
  if (io::getInput(stdin,buf,0))
    in_buf->setSymbol(s,buf);
}


/*
  Set i/o settings to terse style. This style is meant for outputting files
  that are easily parsed by computer. The ordering of the generators is left
  untouched, and can be set independently by the user. Decimal output symbols
  are chosen, prefix is set to "[", postfix to "]" and separator to ",".
*/
void interf::terse_f()
{
  in_buf.reset(new interface::GroupEltInterface(W->rank(),io::GAP()));
  W->interface().setIn(*in_buf);
  W->interface().setOut(*in_buf);

  W->interface().setDescent(io::Default());
  W->setOutputTraits(io::Terse());
}

void interf::in::alphabetic_f()

/*
  Sets the input symbols to the alphabetic sequence.
*/

{
  const std::string* alpha =
    interface::alphabeticSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = alpha[j];
  }
}

void interf::in::bourbaki_f()

/*
  Sets Bourbaki conventions for input. This means reverting the ordering
  of the input symbols in types B and D.
*/

{
  const type::Type& x = W->type();

  if (!isFiniteType(x))
    return;
  if (!(isTypeB(x) || isTypeD(x)))
    return;

  for (coxtypes::Generator s = 0; s < W->rank(); ++s) {
    in_buf->symbol[s] = W->interface().inSymbol(W->rank()-s-1);
  }
}

void interf::in::decimal_f()

/*
  Sets the input symbols to the decimal sequence.
*/

{
  const std::string* dec = interface::decimalSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = dec[j];
  }
}


/*
  Sets the input interface to the default style.
*/

void interf::in::default_f()
{
  in_buf.reset(new interface::GroupEltInterface(W->rank()));
}


/*
  Sets the input interface to GAP style, and enforces Bourbaki conventions.
*/
void interf::in::gap_f()
{
  in_buf.reset(new interface::GroupEltInterface(W->rank(),io::GAP()));
  in::bourbaki_f();
}

void interf::in::hexadecimal_f()

/*
  Sets the input symbols to the hexadecimal sequence.
*/

{
  const std::string* hex = interface::hexSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = hex[j];
  }
}

void interf::in::permutation_f()

/*
  Sets input to permutation mode (type A only.)
*/

{
  using namespace typeA;

  if (!isTypeA(W->type())) {
    io::printFile(stderr,"permutation.mess",directories::MESSAGE_DIR);
    return;
  }

  TypeACoxGroup* WA = dynamic_cast<TypeACoxGroup*>(W);
  WA->setPermutationInput(true);

  in_buf.reset();
}

void interf::in::postfix_f()

/*
  Resets the input postfix.
*/

{
  printf("Enter the new input postfix (finish with a carriage return):\n");
  std::string buf;
  if (io::getInput(stdin,buf,0))
    in_buf->setPostfix(buf);
}

void interf::in::prefix_f()

/*
  Resets the input prefix.
*/

{
  printf("Enter the new input prefix (finish with a carriage return):\n");
  std::string buf;
  if (io::getInput(stdin,buf,0))
    in_buf->setPrefix(buf);
}

void interf::in::separator_f()

/*
  Resets the input separator.
*/

{
  printf("Enter the new input separator (finish with a carriage return):\n");
  std::string buf;
  if (io::getInput(stdin,buf,0))
    in_buf->setSeparator(buf);
}


/*
  Sets the input interface to terse style (the same as GAP style).
*/

void interf::in::terse_f()
{
  in_buf.reset(new interface::GroupEltInterface(W->rank(),io::GAP()));
}

void interf::out::alphabetic_f()

/*
  Changes the output symbols to hexadecimal.
*/

{
  const std::string* alpha =
    interface::alphabeticSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = alpha[j];
  }
}

void interf::out::bourbaki_f()

/*
  Sets Bourbaki conventions for input. This means reverting the ordering
  of the output symbols in types B and D, and setting the output ordering
  to the reverse of the default one.
*/

{
  const type::Type& x = W->type();

  if (!isFiniteType(x))
    return;
  if (!(isTypeB(x) || isTypeD(x))) {
    W->setOrdering(bits::Permutation(W->rank()));
    return;
  }

  for (coxtypes::Generator s = 0; s < W->rank(); ++s) {
    in_buf->symbol[s] = W->interface().outSymbol(W->rank()-s-1);
  }

  bits::Permutation a(W->rank());

  for (coxtypes::Generator s = 0; s < W->rank(); ++s) {
    a[s] = W->rank()-1-s;
  }

  W->setOrdering(a);
}

void interf::out::decimal_f()

/*
  Changes the output symbols to decimal.
*/

{
  const std::string* dec = interface::decimalSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = dec[j];
  }
}


/*
  Sets output styles to the default style.
*/
void interf::out::default_f()
{
  in_buf.reset(new interface::GroupEltInterface(W->rank()));
  W->setOrdering(bits::Permutation(W->rank()));

  W->setOutputTraits(io::Pretty());
}


/*
  Sets output styles to GAP style, and enforces Bourbaki conventions.
*/
void interf::out::gap_f()
{
  in_buf.reset(new interface::GroupEltInterface(W->rank(),io::GAP()));
  W->setOrdering(bits::Permutation(W->rank()));
  out::bourbaki_f();

  W->interface().setDescent(io::GAP());
  W->interface().setOut(*in_buf); // has to be done here so that output traits
                                  // will be correct.
  W->setOutputTraits(io::GAP());
}

void interf::out::hexadecimal_f()

/*
  Changes the output symbols to hexadecimal.
*/

{
  const std::string* hex = interface::hexSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = hex[j];
  }
}


/*
  Sets output to permutation mode (type A only.)
*/
void interf::out::permutation_f()
{
  using namespace typeA;

  if (!isTypeA(W->type())) {
    io::printFile(stderr,"permutation.mess",directories::MESSAGE_DIR);
    return;
  }

  TypeACoxGroup* WA = dynamic_cast<TypeACoxGroup*>(W);

  WA->setPermutationOutput(true);

  W->interface().setOrder(bits::Permutation(W->rank()));
  W->interface().setDescent(io::Default());
  W->setOutputTraits(io::Pretty());

  in_buf.reset();
}

void interf::out::postfix_f()

/*
  Resets the output postfix.
*/

{
  printf("enter the new output postfix (finish with a carriage return):\n");
  std::string buf;
  if (io::getInput(stdin,buf,0))
    in_buf->setPostfix(buf);
}

void interf::out::prefix_f()

/*
  Resets the output prefix.
*/

{
  printf("Enter the new output prefix (finish with a carriage return):\n");
  std::string buf;
  if (io::getInput(stdin,buf,0))
    in_buf->setPrefix(buf);
}

void interf::out::separator_f()

/*
  Resets the output separator.
*/

{
  printf("Enter the new output separator (finish with a carriage return):\n");
  std::string buf;
  if (io::getInput(stdin,buf,0))
    in_buf->setSeparator(buf);
}


/*
  Sets output styles to terse style (the same as GAP style).
*/
void interf::out::terse_f()
{
  in_buf.reset(new interface::GroupEltInterface(W->rank(),io::GAP()));

  W->interface().setDescent(io::Default());
  W->interface().setOut(*in_buf); // has to be done here so that output
                                  // traits will be correct
  W->setOutputTraits(io::Terse());
}

/*****************************************************************************

        Chapter VI -- Miscellaneous.

  This section contains some auxiliary functions :

  - printCommands(file,tree) : prints info about the various commands
    on the tree;
  - interface_entry() : entry function for the interface mode;
  - interface_exit() : exit function for the interface mode;
  - interf::in_entry() : entry function for the interf::in mode;
  - interf::in_exit() : exit function for the interf::in mode;
  - interf::out_entry() : entry function for the interf::out mode;
  - interf::out_exit() : exit function for the interf::out mode;
  - main_entry() : entry function for the main mode;
  - main_exit() : exit function for the main mode;
  - uneq_entry() : entry function for the uneq mode;
  - uneq_exit() : exit function for the uneq mode;

 *****************************************************************************/


/*
  Print one line for each command on the tree (sorted in alphabetical order)
  with the name of the command and the information contained in the tag field.
  The tree root is skipped, even if it should have an own action.
*/
void printCommands(FILE* file, const CommandTree& tree)
{
  for (auto it = std::next(tree.begin()); it != tree.end(); ++it)
    if (it->has_own_action())
    {
      const auto cd = it->action();
      fprintf(file,"  - %s : %s;\n",cd->name.c_str(),cd->tag.c_str());
    }
}


namespace {

void interface_entry()
{
  commands::interf::in_buf.reset(new interface::GroupEltInterface(W->rank()));
}

void interface_exit()
{
  commands::interf::in_buf.reset();
}


/*
  Set |W| after getting the type and rank.
  This is used as entry function for the main mode,
  and as the function that restarts the program when we change the type.

  NOTE : error handling should be done by the calling function.

  NOTE : something should be done about deallocating |W| before reallocating!
*/
void main_entry()
{
  W = interactive::allocCoxGroup();

  /* an error may be set here */
}

};


/*
  Entry function to the input interface modification mode. The global variable
  in_buf is originally set to value for the current group.
*/
void interf::in_entry()
{
  bits::Permutation a = W->interface().order().inverse();

  printf("current input symbols are the following :\n\n");
  interactive::printInterface(stdout,W->interface().inInterface(),a);
  printf("\n");

  in_buf.reset
    (new interface::GroupEltInterface(W->interface().inInterface()));
}


/*
  Exit function from the input modification mode. It checks if the
  modifications made by the user are consistent with the ones that will
  remain from the old interface; if yes, it confirms the modifications
  and exits peacefully; if not, it prints out the problems and keeps
  the user in the mode.
*/
void interf::in_exit()
{
  if (in_buf == nullptr) // hack to prevent execution in special cases
    return;

  bits::Permutation a = W->interface().order().inverse();

  /* at this point in_buf holds the full putative new interface; we
   need to check for reserved or repeated non-empty symbols */

  const std::string* str = checkLeadingWhite(*in_buf);

  if (str) {
    Error(LEADING_WHITESPACE,in_buf.get(),&W->interface().inInterface(),&a,str);
    goto error_exit;
  }

  str = checkReserved(*in_buf,W->interface());

  if (str) {
    Error(RESERVED_SYMBOL,in_buf.get(),&W->interface().inInterface(),&a,str);
    goto error_exit;
  }

  if (!checkRepeated(*in_buf)) {
    Error(REPEATED_SYMBOL,in_buf.get(),&W->interface().inInterface(),&a);
    goto error_exit;
  }

  /* if we reach this point, the new interface is ok */

  printf("new input symbols:\n\n");
  interactive::printInterface(stdout,*in_buf,a);
  printf("\n");

  W->interface().setIn(*in_buf);

  return;

 error_exit:
  ERRNO = ERROR_WARNING;
}


/*
  Entry function to the output interface modification mode. The global variable
  in_buf is originally set to value for the current group.
*/
void interf::out_entry()
{
  in_buf.reset
    (new interface::GroupEltInterface(W->interface().outInterface()));

  bits::Permutation a = W->interface().order().inverse();

  printf("current output symbols are the following :\n\n");
  interactive::printInterface(stdout,*in_buf,W->interface().inInterface(),a);
  printf("\n");
}


/*
  Exit function for the output modification mode. No checking is necessary
  here.
*/
void interf::out_exit()
{
  if (in_buf == nullptr) // hack to prevent execution in special cases
    return;

  bits::Permutation a = W->interface().order().inverse();

  printf("new output symbols:\n\n");
  interactive::printInterface(stdout,*in_buf,W->interface().inInterface(),a);
  printf("\n");

  W->interface().setOut(*in_buf);
}


namespace {

void main_exit()

/*
  Symmetric function to main_entry. Should undo what main_entry did.
*/

{
  delete W;
  wgraph_warning = true;
}

void uneq_entry()

{
  W->activateUEKL();
}


/*
  We keep the unequal-parameter context, because one might want to go back
  and forth between the unequal and the ordinary case. So exit is a no-op
*/
void uneq_exit()
{}

}; // |namespace|

}; // |namespace commands|
