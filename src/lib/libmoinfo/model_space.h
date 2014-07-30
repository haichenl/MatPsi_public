/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef _psi_src_lib_libmoinfo_model_space_h_
#define _psi_src_lib_libmoinfo_model_space_h_

/*! \file    model_space.h
    \ingroup LIBMOINFO
    \brief   This class stores all the basic info regarding the model space
*/

#include "slater_determinant.h"

namespace psi {

class MOInfo;

class ModelSpace{
public:
  ModelSpace(MOInfo* moinfo_obj_);
  ~ModelSpace();
  void print();
private:
  void startup();
  void cleanup();
  void build();
  void classify();

  int wfn_sym;
  std::vector<SlaterDeterminant> determinants;
  std::vector<int>  closed_to_all;  // closed-shell determinants
  std::vector<int>  opensh_to_all;  // open-shell   determinants
  std::vector<int>  unique_to_all;  // spin-unique  determinants
  MOInfo* moinfo_obj;
};

    extern ModelSpace  *model_space;

}

#endif // _psi_src_lib_libmoinfo_model_space_h_