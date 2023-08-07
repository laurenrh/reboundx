/** * @file    gr_potential.c
 * @brief   Post-newtonian general relativity corrections using a simple potential that gets the pericenter precession right.
 * @author  Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo
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
 * $General Relativity$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 H. Rein, D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
 * Based on                `Nobili and Roxburgh 1986 <http://labs.adsabs.harvard.edu/adsabs/abs/1986IAUS..114..105N/>`_.
 * C Example               :ref:`c_example_gr`
 * Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
 * ======================= ===============================================
 * 
 * This is the simplest potential you can use for general relativity.
 * It assumes that the masses are dominated by a single central body.
 * It gets the precession right, but gets the mean motion wrong by :math:`\mathcal{O}(GM/ac^2)`.  
 * It's the fastest option, and because it's not velocity-dependent, it automatically keeps WHFast symplectic.  
 * Nice if you have a single-star system, don't need to get GR exactly right, and want speed.
 * 
 * **Effect Parameters**
 * 
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * c (double)                   Yes         Speed of light, needs to be specified in the units used for the simulation.
 * ============================ =========== ==================================================================
 *
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"


static void rebx_calculate_nfw_potential(struct reb_particle* const particles, const int N, const double R_s, const double nfw_const, const double G){
    const struct reb_particle source = particles[0];
    
    for (int i=1; i<N; i++){
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        
        particles[i].ax += nfw_const*(p.x / ((R_s + pow(r2,.5)) * r2) 
                	- (p.x * log(1 + (pow(r2,.5) / R_s))
        		/ (pow(r2,1.5))) );
        particles[i].ay += nfw_const*(p.y / ((R_s + pow(r2,.5)) * r2) 
                	- (p.y * log(1 + (pow(r2,.5) / R_s))
        		/ (pow(r2,1.5))) );
        particles[i].az += nfw_const*(p.z / ((R_s + pow(r2,.5)) * r2) 
                	- (p.z * log(1 + (pow(r2,.5) / R_s))
        		/ (pow(r2,1.5))) );
        particles[0].ax -= p.m/source.m * particles[i].ax;
        particles[0].ay -= p.m/source.m * particles[i].ay;
        particles[0].az -= p.m/source.m * particles[i].az;
    }
}

 void rebx_nfw_potential(struct reb_simulation* const sim, struct rebx_force* const nfw_potential, struct reb_particle* const particles, const int N){
    double* rho_0 = rebx_get_param(sim->extras, nfw_potential->ap, "nfw_rho_0");
    double* Rs = rebx_get_param(sim->extras, nfw_potential->ap, "nfw_Rs"); //why not formetting right
    if (rho_0 == NULL){
        reb_error(sim, "REBOUNDx Error: Need to set ____ density.  See examples in documentation.\n");
    }
    else{
        if (Rs == NULL){
        reb_error(sim, "REBOUNDx Error: Need to set scale height?.  See examples in documentation.\n");
    	}
    	else{ 
		//const double Rs3 = (*Rs)*(*Rs)*(*Rs); //this was so not even necessry
		const double nfw_const = 4*M_PI*G*(*rho_0)*(*Rs)*(*Rs)*(*Rs);
        	rebx_calculate_nfw_potential(particles, N, *Rs, nfw_const, sim->G);
	}
	}
}
//^ is the else up here sloppy should i do like just if i forget whats allowed. look up

//change name maybe.irrelevant ish
static double rebx_calculate_nfw_potential_potential(struct reb_simulation* const sim, const double rho_0, const double R_s){
    const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;      //why. collisions? find out
	const double G = sim->G;
    const struct reb_particle source = particles[0];
	const double nfw_con = 4*M_PI*G*(*rho_0)*R_s*R_s*R_s;
    double H = 0.;

	for (int i=1;i<_N_real;i++){
		struct reb_particle pi = particles[i];
        double dx = pi.x - source.x;
        double dy = pi.y - source.y;
        double dz = pi.z - source.z;
        double r2 = dx*dx + dy*dy + dz*dz;
        H -= nfw_con*pi.m*log(1 + pow(r2,.5) / R_s) / pow(r2,.5);
    }		
	
    return H;
}
//oo modify later so u could also put in like virial wahetevrs and c and things? that will be a pain lol. not worth?
//fix spacing idk kinda ugly

double rebx_gr_potential_potential(struct rebx_extras* const rebx, const struct rebx_force* const gr_potential){
    double* rho_0 = rebx_get_param(rebx, nfw_potential->ap, "rho0");
    double* Rs = rebx_get_param(rebx, nfw_potential->ap, "Rs");
    if (rho_0 == NULL){
        rebx_error(rebx, "Need to set ____ density.  See examples in documentation.\n");
    	return 0;
    }
    if (Rs == NULL){
	reb_error(sim, "REBOUNDx Error: Need to set scale height?.  See examples in documentation.\n");
	return 0;
    }
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return 0;
    }
    return rebx_calculate_gr_potential_potential(rebx->sim, rho0, Rs);
}
