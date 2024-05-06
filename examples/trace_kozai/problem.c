/**
 *
 * Kozai Cycles with TRACE. Figure 8 in Lu, Hernandez & Rein (2024)
 *
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include <time.h>
 #include "rebound.h"

void heartbeat(struct reb_simulation* r);
double tmax = 3e6 * 2 * M_PI;
double e_init;

char title[100] = "energy.txt";
char title_remove[100] = "energy.txt";

double obl(struct reb_vec3d v1, struct reb_vec3d v2){
  return acos(reb_vec3d_dot(v1,v2) / (sqrt(reb_vec3d_length_squared(v1)) * sqrt(reb_vec3d_length_squared(v2))));
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    r->integrator        = REB_INTEGRATOR_TRACE;
    r->heartbeat = heartbeat;

    struct reb_particle star = {0};
    star.m  = 0.32;
    star.r = 0.1 * 0.00465;
    reb_simulation_add(r, star);

    double planet_m  = 5.15e-5; 
    double planet_a = 2.;
    double planet_e = 0.01;
    double planet_omega = 0.;
    reb_simulation_add_fmt(r, "m a e", planet_m, planet_a, planet_e);

    struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[1], r->particles[0]);
    r->dt = o.P / 30.12345;

    // The perturber
    double perturber_mass = 10. * 9.55e-4;
    double perturber_a  = 50.;
    double perturber_e = 0.52;
    double perturber_inc = 80. * M_PI/180.;
    double perturber_omega = 0.;

    reb_simulation_add_fmt(r, "m a e inc omega", perturber_mass, perturber_a, perturber_e, perturber_inc, perturber_omega);
    reb_simulation_move_to_com(r);

    system(title_remove);
    FILE* f = fopen(title,"w");
    fprintf(f, "t,E,a1,i1,e1,a2,i2,e2,mi,d\n");
    fclose(f);

    e_init = reb_simulation_energy(r);
    
    clock_t begin = clock();
    reb_simulation_integrate(r, tmax);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    
    printf("Energy Error: %e Time Spent: %e\n", fabs((reb_simulation_energy(r) - e_init) / e_init), time_spent);
    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* r){
   // Output spin and orbital information to file

   if (reb_simulation_output_check(r, 100.*2.*M_PI)){
        reb_simulation_output_timing(r, tmax);
   }

   if (reb_simulation_output_check(r, 1e4 * 2.*M_PI)){
     struct reb_particle* sun = &r->particles[0];
     struct reb_particle* p1 = &r->particles[1];
     struct reb_particle* pert = &r->particles[2];

     // orbits
     struct reb_orbit o1 = reb_orbit_from_particle(r->G, *p1, *sun);
     double a1 = o1.a;
     double e1 = o1.e;
     double i1 = o1.inc;
     double Om1 = o1.Omega;
     double pom1 = o1.pomega;
     struct reb_vec3d n1 = o1.hvec;

     struct reb_orbit o2 = reb_orbit_from_particle(r->G, *pert, *sun);
     double a2 = o2.a;
     double e2 = o2.e;
     double i2 = o2.inc;
     double Om2 = o2.Omega;
     double pom2 = o2.pomega;
     struct reb_vec3d n2 = o2.hvec;
     double mi = obl(n1,n2);

     double dx = p1->x - sun->x;
     double dy = p1->y - sun->y;
     double dz = p1->z - sun->z;
     double d = sqrt(dx*dx + dy*dy + dz*dz);

     FILE* f = fopen(title,"a");
     fprintf(f, "%f,%e,%f,%f,%f,%f,%f,%f,%f,%f\n", r->t,fabs((reb_simulation_energy(r) - e_init) / e_init),a1,i1,e1,a2,i2,e2,mi,d); // print spins and orbits
     fclose(f);
   }

}
