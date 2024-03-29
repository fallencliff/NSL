/* ode-initval/cscal.c
 * 
 * Copyright (C) 2002 Brian Gough
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


#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl_errno.h>
#include <MathConst.h>
#include <odeiv_struct.h>

using namespace gslcpp;

namespace scale_control
{

typedef struct
{
  double eps_abs;
  double eps_rel;
  double a_y;
  double a_dydt;
  double * scale_abs;
}
sc_control_state_t;

static void *
sc_control_alloc (void)
{
  sc_control_state_t * s = 
    (sc_control_state_t *) malloc (sizeof(sc_control_state_t));
  
  if (s == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for sc_control_state", 
                      GSL_ENOMEM);
    }

  return s;
}

static int
sc_control_init (void * vstate, 
                  double eps_abs, double eps_rel, double a_y, double a_dydt)
{
  sc_control_state_t * s = (sc_control_state_t *) vstate;
  
  if (eps_abs < 0)
    {
      GSL_ERROR ("eps_abs is negative", GSL_EINVAL);
    }
  else if (eps_rel < 0)
    {
      GSL_ERROR ("eps_rel is negative", GSL_EINVAL);
    }
  else if (a_y < 0)
    {
      GSL_ERROR ("a_y is negative", GSL_EINVAL);
    }
  else if (a_dydt < 0)
    {
      GSL_ERROR ("a_dydt is negative", GSL_EINVAL);
    }
  
  s->eps_rel = eps_rel;
  s->eps_abs = eps_abs;
  s->a_y = a_y;
  s->a_dydt = a_dydt;

  return GSL_SUCCESS;
}

static int
sc_control_hadjust(void * vstate, size_t dim, unsigned int ord, const double y[], const double yerr[], const double yp[], double * h)
{
  sc_control_state_t *state = (sc_control_state_t *) vstate;

  const double eps_abs = state->eps_abs;
  const double eps_rel = state->eps_rel;
  const double a_y     = state->a_y;
  const double a_dydt  = state->a_dydt;
  const double * scale_abs = state->scale_abs;

  const double S = 0.9;
  const double h_old = *h;

  double rmax = DOUBLE_MIN;
  size_t i;

  for(i=0; i<dim; i++) {
    const double D0 = 
      eps_rel * (a_y * fabs(y[i]) + a_dydt * fabs(h_old * yp[i])) 
      + eps_abs * scale_abs[i];
    const double r  = fabs(yerr[i]) / fabs(D0);
    rmax = CMath::max(r, rmax);
  }

  if(rmax > 1.1) {
    /* decrease step, no more than factor of 5, but a fraction S more
       than scaling suggests (for better accuracy) */
    double r =  S / pow(rmax, 1.0/ord);
    
    if (r < 0.2)
      r = 0.2;

    *h = r * h_old;

    return GSL_ODEIV_HADJ_DEC;
  }
  else if(rmax < 0.5) {
    /* increase step, no more than factor of 5 */
    double r = S / pow(rmax, 1.0/(ord+1.0));

    if (r > 5.0)
      r = 5.0;

    if (r < 1.0)  /* don't allow any decrease caused by S<1 */
      r = 1.0;
        
    *h = r * h_old;

    return GSL_ODEIV_HADJ_INC;
  }
  else {
    /* no change */
    return GSL_ODEIV_HADJ_NIL;
  }
}


static void
sc_control_free (void * vstate)
{
  sc_control_state_t *state = (sc_control_state_t *) vstate;
  free (state->scale_abs);
  free (state);
}

static const gsl_odeiv_control_type sc_control_type =
{"scaled",                      /* name */
 &sc_control_alloc,
 &sc_control_init,
 &sc_control_hadjust,
 &sc_control_free};

}
const gsl_odeiv_control_type *gsl_odeiv_control_scaled = &scale_control::sc_control_type;

