/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "util.h"

#define XC_GGA_K_OL1          512 /* Ou-Yang and Levy v.1 */

#include "maple2c/gga_k_ol1.c"

#define func maple2c_func
#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const xc_func_info_type xc_func_info_gga_k_ol1 = {
  XC_GGA_K_OL1,
  XC_KINETIC,
  "Ou-Yang and Levy v.1",
  XC_FAMILY_GGA,
  {&xc_ref_OuYang1991_379, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  5e-26,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, work_gga_k, NULL
};
