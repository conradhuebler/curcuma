#pragma once
#include "src/global_config.h"

#include <string>

const std::string git_commit_hash = "f4eb2d4"; /*! Current hash of the git branch */
const std::string git_branch = "rmsd_mtd";      /*! Current name of the branch as printable string    */
const std::string git_tag = "0.0.146-7-gf4eb2d4";      /*! Current name of the branch as printable string    */

#ifdef _UNIX
const  std::string git_date = "Commit Date: Wed, 31 Jul 2024 17:38:09 +0200"; /*! Current date of the git branch branch as printable string*/
#else
const  std::string git_date = ""; /*! Current date of the git branch branch as printable string*/
#endif

const  std::string branch = " Git Branch: rmsd_mtd -- ";      /*! Current name of the branch as printable string    */
const  std::string commit_hash = "Commit: f4eb2d4 --"; /*! Current hash of the git branch branch as printable string*/

#ifdef _UNIX
const  std::string date = "Commit Date: Wed, 31 Jul 2024 17:38:09 +0200"; /*! Current date of the git branch branch as printable string*/
#else
const  std::string date = "Forever young."; /*! Current date of the git branch branch as printable string*/
#endif

#ifdef _DEBUG
const  std::string conf_mode = "DEBUG Mode"; /*! Last compilation mode used */
#else
const  std::string conf_mode = "RELEASE Mode"; /*! Last compilation mode used */
#endif

// const  std::string version = std::string("Curcuma 0.1 pre-Alpha \n\t%1 %2 %3 -  %4").arg(branch).arg(commit_hash).arg(date).arg(conf_mode); /*! Version name */


const double qint_version = 0.1;
