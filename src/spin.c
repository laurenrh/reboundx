/**
 * @file    tides_constant_time_lag.c
 * @brief   Add constant time lag tides raised on primary, orbiting bodies, or both
 * @author  Stanley A. Baronett <stanley.a.baronett@gmail.com>, Dan Tamayo <tamayo.daniel@gmail.com>
 *
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Tides$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 Stanley A. Baronett, D. Tamayo, Noah Ferich
 * Implementation Paper    `Baronett et al., 2021 (in review) <https://arxiv.org/abs/2101.12277>`_.
 * Based on                `Hut 1981 <https://ui.adsabs.harvard.edu/#abs/1981A&A....99..126H/abstract>`_, `Bolmont et al., 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...583A.116B/abstract>`_.
 * C Example               :ref:`c_example_tides_constant_time_lag`.
 * Python Example          `TidesConstantTimeLag.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesConstantTimeLag.ipynb>`_.
 * ======================= ===============================================
 *
 * This adds constant time lag tidal interactions between orbiting bodies in the simulation and the primary, both from tides raised on the primary and on the other bodies.
 * In all cases, we need to set masses for all the particles that will feel these tidal forces. After that, we can choose to include tides raised on the primary, on the "planets", or both, by setting the respective bodies' physical radius particles[i].r, k2 (potential Love number of degree 2), constant time lag tau, and rotation rate Omega. See Baronett et al. (2021), Hut (1981), and Bolmont et al. 2015 above.
 *
 * If tau is not set, it will default to zero and yield the conservative piece of the tidal potential.
 *
 * **Effect Parameters**
 *
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * particles[i].r (float)       Yes         Physical radius (required for contribution from tides raised on the body).
 * tctl_k2 (float)              Yes         Potential Love number of degree 2.
 * tctl_tau (float)             No          Constant time lag. If not set will default to 0 and give conservative tidal potential.
 * Omega (float)                No          Rotation rate. If not set will default to 0.
 * ============================ =========== ==================================================================
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "reboundx.h"

struct reb_vec3d rebx_calculate_spin_orbit_accelerations(struct reb_particle* source, struct reb_particle* target, const double G, const double k2, const double sigma, const double sx, const double sy, const double sz){

  const double ms = source->m;
  const double mt = target->m;
  const double Rt = target->r;
  const double mratio = ms / mt; // have already checked for 0 and inf
  const double mu_ij = ms * mt / (ms + mt);

  // distance vector FROM j TO i
  const double dx = source->x - target->x;
  const double dy = source->y - target->y;
  const double dz = source->z - target->z;
  const double d2 = dx * dx + dy * dy + dz * dz;
  const double dr = sqrt(d2);

  struct reb_vec3d tot_force = {0};

  if (k2 != 0){
    const double quad_prefactor = ms * k2 * (Rt * Rt * Rt * Rt * Rt) / mu_ij;
    const double omega_dot_d = sx * dx + sy * dy + sz * dz;
    const double omega_squared = sx * sx + sy * sy + sz * sz;

    const double t1 = 5. * omega_dot_d * omega_dot_d / (2 * (dr * dr * dr * dr * dr * dr * dr));
    const double t2 = omega_squared / (2 * (dr * dr * dr * dr * dr));
    const double t3 = omega_dot_d / (dr * dr * dr * dr * dr);
    const double t4 = 6. * G * ms / (dr * dr * dr * dr * dr * dr * dr * dr);

    tot_force.x += quad_prefactor * ((t1 - t2 - t4) * dx - t3 * sx);
    tot_force.y += quad_prefactor * ((t1 - t2 - t4) * dy - t3 * sy);
    tot_force.z += quad_prefactor * ((t1 - t2 - t4) * dz - t3 * sz);

  }

  if (sigma != 0){
    // EKH FRAMEWORK
    const double a2 = k2 * k2 * Rt * Rt * Rt * Rt * Rt * Rt * Rt * Rt * Rt * Rt;

    // Velocity vector: j to i
    const double dvx = source->vx - target->vx;
    const double dvy = source->vy - target->vy;
    const double dvz = source->vz - target->vz;
    const double vel2 = dvx*dvx + dvy*dvy + dvz*dvz;

    // h vector - EKH
    const double hx = dy * dvz - dz * dvy;
    const double hy = dz * dvx - dx * dvz;
    const double hz = dx * dvy - dy * dvx;

    // h - r^2 Omega
    const double vec1_x = hx - d2 * sx;
    const double vec1_y = hy - d2 * sy;
    const double vec1_z = hz - d2 * sz;

    // second vector
    const double vec2_x = vec1_y * dz - vec1_z * dy;
    const double vec2_y = vec1_z * dx - vec1_x * dz;
    const double vec2_z = vec1_x * dy - vec1_y * dx;

    const double prefactor = (-9 * sigma * ms * ms * a2) / (2 * mu_ij * (d2 * d2 * d2 * d2 * d2));

    tot_force.x += prefactor * (vec1_x + vec2_x);
    tot_force.y += prefactor * (vec1_y + vec2_y);
    tot_force.z += prefactor * (vec1_z + vec2_z);
  }

  return tot_force;
}

static void rebx_spin_orbit_accelerations(struct reb_particle* source, struct reb_particle* target, const double G, const double k2, const double sigma, const double sx, const double sy, const double sz){

    // Input params all associated with target
    struct reb_vec3d tot_force = rebx_calculate_spin_orbit_accelerations(source, target, G, k2, sigma, sx, sy, sz);

    const double ms = source->m;
    const double mt = target->m;
    const double mtot = ms + mt;

    target->ax -= ((ms / mtot) * tot_force.x);
    target->ay -= ((ms / mtot) * tot_force.y);
    target->az -= ((ms / mtot) * tot_force.z);

    source->ax += ((mt / mtot) * tot_force.x);
    source->ay += ((mt / mtot) * tot_force.y);
    source->az += ((mt / mtot) * tot_force.z);
}

static void rebx_spin_derivatives(struct reb_ode* const ode, double* const yDot, const double* const y, const double t){
    struct reb_simulation* sim = ode->ref;
    struct rebx_extras* const rebx = sim->extras;
    unsigned int Nspins = 0;
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* primary = &sim->particles[0];
    for (int i=0; i<N_real; i++){
        struct reb_particle* pi = &sim->particles[i]; // target particle
        double* k2 = rebx_get_param(rebx, pi->ap, "k2");
        double* sigma = rebx_get_param(rebx, pi->ap, "sigma");
        double* moi = rebx_get_param(rebx, pi->ap, "moi");

	// Set initial spin accelerations to 0
        yDot[3*Nspins] = 0;
        yDot[3*Nspins+1] = 0;
        yDot[3*Nspins+2] = 0;

        const double sx = y[3 * Nspins];
        const double sy = y[3 * Nspins + 1];
        const double sz = y[3 * Nspins + 2];

        for (int j=0; j<N_real; j++){
          if (k2 != NULL && sigma != NULL && i != j){ // Look into this conditional more
              struct reb_particle* pj = &sim->particles[j];

              double dx = pj->x - pi->x;
              double dy = pj->y - pi->y;
              double dz = pj->z - pi->z;

              double mi = pi->m;
              double mj = pj->m;
              double mu_ij = -(mi * mj) / ((mi + mj));

              struct reb_vec3d tf = rebx_calculate_spin_orbit_accelerations(pj, pi, sim->G, *k2, *sigma, sx, sy, sz);
              yDot[3*Nspins] += ((dy * tf.z - dz * tf.y) * (mu_ij / *moi));
              yDot[3*Nspins+1] += ((dz * tf.x - dx * tf.z) * (mu_ij / *moi));
              yDot[3*Nspins+2] += ((dx * tf.y - dy * tf.x) * (mu_ij / *moi));
          }
        }
        Nspins += 1;
    }
    if (ode->length != Nspins*3){
        reb_error(sim, "rebx_spin ODE is not of the expected length.\n");
    }
}

static void rebx_spin_sync_pre(struct reb_ode* const ode, const double* const y0){
    struct reb_simulation* sim = ode->ref;
    struct rebx_extras* const rebx = sim->extras;
    unsigned int Nspins = 0;
    const int N_real = sim->N - sim->N_var;
    for (int i=0; i<N_real; i++){
        struct reb_particle* p = &sim->particles[i];
        double* k2 = rebx_get_param(rebx, p->ap, "k2");
        if (k2 != NULL){
            double* sx = rebx_get_param(rebx, p->ap, "spin_sx");
            double* sy = rebx_get_param(rebx, p->ap, "spin_sy");
            double* sz = rebx_get_param(rebx, p->ap, "spin_sz");
            ode->y[3*Nspins] = *sx;
            ode->y[3*Nspins+1] = *sy;
            ode->y[3*Nspins+2] = *sz;
            Nspins += 1;
        }
    }

    if (ode->length != Nspins*3){
        reb_error(sim, "rebx_spin ODE is not of the expected length.\n");
    }
}

static void rebx_spin_sync_post(struct reb_ode* const ode, const double* const y0){
    struct reb_simulation* sim = ode->ref;
    struct rebx_extras* const rebx = sim->extras;
    unsigned int Nspins = 0;
    const int N_real = sim->N - sim->N_var;
    for (int i=0; i<N_real; i++){
        struct reb_particle* p = &sim->particles[i];
        double* k2 = rebx_get_param(rebx, p->ap, "k2");
        if (k2 != NULL){
            rebx_set_param_double(rebx, &p->ap, "spin_sx", y0[3*Nspins]);
            rebx_set_param_double(rebx, &p->ap, "spin_sy", y0[3*Nspins+1]);
            rebx_set_param_double(rebx, &p->ap, "spin_sz", y0[3*Nspins+2]);
            Nspins += 1;
        }
    }
    if (ode->length != Nspins*3){
        reb_error(sim, "rebx_spin ODE is not of the expected length.\n");
    }
}


void rebx_spin_initialize_ode(struct reb_simulation* sim, struct rebx_force* const effect){
    struct rebx_extras* const rebx = sim->extras;
    unsigned int Nspins = 0;
    const int N_real = sim->N - sim->N_var;
    for (int i=0; i<N_real; i++){
        struct reb_particle* p = &sim->particles[i];
        // Only track spin if particle has moment of inertia and valid spin axisset
        double* moi = rebx_get_param(rebx, p->ap, "moi");
        double* sx = rebx_get_param(rebx, p->ap, "spin_sx");
        double* sy = rebx_get_param(rebx, p->ap, "spin_sy");
        double* sz = rebx_get_param(rebx, p->ap, "spin_sz");
        if (moi != NULL && sx != NULL && sy != NULL && sz != NULL){
            Nspins += 1;
        }
    }

    if (Nspins > 0){
        struct reb_ode* spin_ode = reb_create_ode(sim, Nspins*3);
        printf("Spin ODE length: %d\n", spin_ode->length);
        spin_ode->ref = sim;
        spin_ode->derivatives = rebx_spin_derivatives;
        spin_ode->pre_timestep = rebx_spin_sync_pre;
        spin_ode->post_timestep = rebx_spin_sync_post;
        rebx_set_param_pointer(rebx, &effect->ap, "ode", spin_ode);
    }
}

void rebx_spin(struct reb_simulation* const sim, struct rebx_force* const effect, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    const double G = sim->G;

    for (int i=0; i<N; i++){
        struct reb_particle* target = &particles[i];
        // Particle must have a k2 set, otherwise we treat this body as a point particle
        double* k2 = rebx_get_param(rebx, target->ap, "k2");
        double* sigma = rebx_get_param(rebx, target->ap, "sigma");
        double* sx = rebx_get_param(rebx, target->ap, "spin_sx");
        double* sy = rebx_get_param(rebx, target->ap, "spin_sy");
        double* sz = rebx_get_param(rebx, target->ap, "spin_sz");
        if (sx == 0 && sy == 0 && sz == 0){
            // Particle needs at least one component of spin to exert additional forces
            continue;
        }
        for (int j=0; j<N; j++){
            if (i==j){
                continue;
            }
            struct reb_particle* source = &particles[j]; // planet raising the tides on the star
            if (source->m == 0){
                continue;
            }
            rebx_spin_orbit_accelerations(source, target, G, *k2, *sigma, *sx, *sy, *sz);
        }
    }
}

// Calculate potential of conservative piece of tidal interaction
static double rebx_calculate_spin_potential(struct reb_particle* source, struct reb_particle* target, const double G, const double k2){
    const double ms = source->m;
    const double mt = target->m;
    const double Rt = target->r;

    const double mratio = ms/mt; // have already checked for 0 and inf
    const double fac = mratio*k2*Rt*Rt*Rt*Rt*Rt;

    const double dx = target->x - source->x;
    const double dy = target->y - source->y;
    const double dz = target->z - source->z;
    const double dr2 = dx*dx + dy*dy + dz*dz;

    return -1./2.*G*ms*mt/(dr2*dr2*dr2)*fac;
}

double rebx_spin_potential(struct rebx_extras* const rebx){
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return 0;
    }
    struct reb_simulation* const sim = rebx->sim;
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    const double G = sim->G;
    double H=0.;

    for (int i=0; i<N_real; i++){
        struct reb_particle* target = &particles[i];
        // Particle must have a k2 set, otherwise we treat this body as a point particle
        double* k2 = rebx_get_param(rebx, target->ap, "k2");
        if (k2 == NULL || target->r == 0 || target->m == 0){
            continue;
        }
        for (int j=0; j<N_real; j++){
            if (i==j){
                continue;
            }
            struct reb_particle* source = &particles[j]; // planet raising the tides on the star
            if (source->m == 0){
                continue;
            }
            H += rebx_calculate_spin_potential(source, target, G, *k2);
        }
    }

    return H;
}