/*
  This is main.c

  Coxeter version 3.0 Copyright (c) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "constants.h"
#include "commands.h"
#include "version.h"

namespace {
  using namespace version;
  void printVersion();
};

int main()

/*
  In this version, the program can only run in interactive mode, and
  does not take any arguments.
*/

{
  constants::initConstants();
  printVersion();
  commands::run();

  exit(0);
}


namespace {

void printVersion()

/*
  Prints an opening message and the version number.
*/

{

  printf("This is %s version %s.\nEnter help if you need assistance,\
 carriage return to start the program.\n\n",NAME,VERSION);

  return;
}

};
