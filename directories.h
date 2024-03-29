/*
  This is directories.h

  Coxeter version 3.0  Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef DIRECTORIES_H  /* guard against multiple inclusions */
#define DIRECTORIES_H

#include <string>

/*
  This file tells where the directories can be found which contain some
  auxiliary files used by the program. The following directories are defined :

    - COXMATRIX_DIR : contains the files for predefined Coxeter matrices;
      these are the files that are loaded through the "X" group type.

    - HEADER_DIR : contains headers for the output to files done by some
      of the commands.

    - MESSAGE_DIR : contains the text of various error and warning messages.
      This is used mostly by the help facility, and also in some of the
      error handling.
*/

namespace directories {
  static const std::string root(".");
  static const std::string coxmatrix_dir = root + "/coxeter_matrices";
  static const std::string header_dir = root + "/headers";
  static const std::string message_dir = root + "/messages";
  static const char* const COXMATRIX_DIR = coxmatrix_dir.c_str();
  static const char* const HEADER_DIR = header_dir.c_str();
  static const char* const MESSAGE_DIR = message_dir.c_str();
};

#endif
