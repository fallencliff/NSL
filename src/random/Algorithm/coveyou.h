/* rng/coveyou.c
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

/*
 * This generator is taken from
 *
 * Donald E. Knuth
 * The Art of Computer Programming
 * Volume 2
 * Third Edition
 * Addison-Wesley
 * Section 3.2.2
 *
 * This implementation copyright (C) 2001 Carlo Perassi
 * and (C) 2003 Heiko Bauke.
 * Carlo Perassi reorganized the code to use the rng framework of GSL.
 */


#include <stdlib.h>
#include <RandGenerator.h>

using gslcpp::gsl_rng_type;

namespace coveyou
{
	#define COVEYOU_MM 0xffffffffUL         /* 2 ^ 32 - 1 */

static inline unsigned long int cove_rand_get (void *vstate);
static double ran_get_double (void *vstate);
static void ran_set (void *state, unsigned long int s);

typedef struct
{
  unsigned long int x;
}
ran_state_t_cove;

static inline unsigned long int
cove_rand_get (void *vstate)
{
  ran_state_t_cove *state = (ran_state_t_cove *) vstate;

  state->x = (state->x * (state->x + 1)) & COVEYOU_MM;

  return state->x;
}

static double
ran_get_double (void *vstate)
{
  ran_state_t_cove *state = (ran_state_t_cove *) vstate;

  return cove_rand_get (state) / 4294967296.0;
}

static void
ran_set (void *vstate, unsigned long int s)
{
  ran_state_t_cove *state = (ran_state_t_cove *) vstate;

  unsigned long int diff = ((s % 4UL) - 2UL) % COVEYOU_MM;

  if (diff)
    state->x = (s - diff) & COVEYOU_MM;
  else
    state->x = s & COVEYOU_MM;

  return;
}

static const gsl_rng_type cove_ran_type = {
  "coveyou",                    /* name */
  COVEYOU_MM-1,                         /* RAND_MAX */
  2,                            /* RAND_MIN */
  sizeof (ran_state_t_cove),
  &ran_set,
  &cove_rand_get,
  &ran_get_double
};
}


const gsl_rng_type *gsl_rng_coveyou = &coveyou::cove_ran_type;
