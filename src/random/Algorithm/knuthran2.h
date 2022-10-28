/* rng/knuthran2.c
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
 * Page 108
 *
 * This implementation  copyright (C) 2001 Carlo Perassi
 * and (C) 2003 Heiko Bauke.
 */

#include <stdlib.h>
#include <RandGenerator.h>

using gslcpp::gsl_rng_type;

namespace knuthran2
{
#define AA1      271828183UL
#define AA2     1833324378UL    /* = -314159269 mod (2 ^ 31 -1) */
#define KNUTHRAN2_MM      0x7fffffffUL    /* 2 ^ 31 - 1 */
#define CEIL_SQRT_MM 46341UL    /* sqrt(2 ^ 31 - 1) */

static inline unsigned long int ran_get (void *vstate);
static double ran_get_double (void *vstate);
static void ran_set (void *state, unsigned long int s);

typedef struct
{
  unsigned long int x0;
  unsigned long int x1;
}
ran_state_t;

static inline unsigned long int
schrage (unsigned long int a, unsigned long int b, unsigned long int m)
{
  /* This is a modified version of Schrage's method. It ensures that no
   * overflow or underflow occurs even if a=ceil(sqrt(m)). Usual 
   * Schrage's method works only until a=floor(sqrt(m)).
   */
  unsigned long int q, t;
  if (a == 0UL)
    return 0UL;
  q = m / a;
  t = 2 * m - (m % a) * (b / q);
  if (t >= m)
    t -= m;
  t += a * (b % q);
  return (t >= m) ? (t - m) : t;
}

static inline unsigned long int
schrage_mult (unsigned long int a, unsigned long int b,
              unsigned long int m,
              unsigned long int sqrtm)
{
  /* To multiply a and b use Schrage's method 3 times.
   * represent a in base ceil(sqrt(m))  a = a1*ceil(sqrt(m)) + a0  
   * a*b = (a1*ceil(sqrt(m)) + a0)*b = a1*ceil(sqrt(m))*b + a0*b   
   */
  unsigned long int t0 = schrage (sqrtm, b, m);
  unsigned long int t1 = schrage (a / sqrtm, t0, m);
  unsigned long int t2 = schrage (a % sqrtm, b, m);
  unsigned long int t = t1 + t2;
  return (t >= m) ? (t - m) : t;
}

static inline unsigned long int
ran_get (void *vstate)
{
  ran_state_t *state = (ran_state_t *) vstate;

  const unsigned long int xtmp = state->x1;
  state->x1 = schrage_mult (AA1, state->x1, KNUTHRAN2_MM, CEIL_SQRT_MM)
    + schrage_mult (AA2, state->x0, KNUTHRAN2_MM, CEIL_SQRT_MM);

  if (state->x1 >= KNUTHRAN2_MM)
    state->x1 -= KNUTHRAN2_MM;

  state->x0 = xtmp;

  return state->x1;
}

static double
ran_get_double (void *vstate)
{
  ran_state_t *state = (ran_state_t *) vstate;

  return ran_get (state) / 2147483647.0;
}

static void
ran_set (void *vstate, unsigned long int s)
{
  ran_state_t *state = (ran_state_t *) vstate;

  if ((s % KNUTHRAN2_MM) == 0)
    s = 1;                      /* default seed is 1 */

  state->x0 = s % KNUTHRAN2_MM;
  state->x1 = s % KNUTHRAN2_MM;

  return;
}

static const gsl_rng_type ran_type = {
  "knuthran2",                  /* name */
  KNUTHRAN2_MM - 1L,                      /* RAND_MAX */
  0,                            /* RAND_MIN */
  sizeof (ran_state_t),
  &ran_set,
  &ran_get,
  &ran_get_double
};
}



const gsl_rng_type *gsl_rng_knuthran2 = &knuthran2::ran_type;
