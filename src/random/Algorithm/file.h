/* rng/file.c
 * 
 * Copyright (C) 2003 Olaf Lenz
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


#include <stdio.h>
#include <gsl_errno.h>
#include <RandGenerator.h>

using gslcpp::CBase_rgn::gsl_rng_type;

int
gsl_rng_fread (FILE * stream, gsl_rng * r)
{
  size_t n = r->type->size ;

  char * state = r->state;

  size_t items = fread (state, 1, n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
  return GSL_SUCCESS;
}

int
gsl_rng_fwrite (FILE * stream, const gsl_rng * r)
{
  size_t n = r->type->size ;

  char * state = r->state;
  
  size_t items = fwrite (state, 1, n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return GSL_SUCCESS;
}
