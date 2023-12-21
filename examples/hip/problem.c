
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <unistd.h>
 #include <math.h>
 #include "rebound.h"

void heartbeat(struct reb_simulation* r);
double tmax;
int ind;
int stable=1;
int ntest=3;
double mf;
double me;

double planet_as[10] = {0.1283,0.2061,0.88,1.06,1.37};
double planet_aerrs[10] = {1.5e-3, 2.4e-3, 0.01, 0.03, 0.02};

//char title[100] = "low_ob_";
char title_stats[100] = "2body_stability";
//char title_remove[100] = "rm -v low_ob_";

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_simulation_create();
    sim->integrator         = REB_INTEGRATOR_WHFAST;
    sim->heartbeat          = heartbeat;

    ind = 0;
    if (argc == 2){
      //strcat(title, argv[1]);
      //strcat(title_remove, argv[1]);
      ind = atoi(argv[1]);
    }

    sim->rand_seed = ind;

    // Initial conditions
    // Santerne et al 2019
    struct reb_particle star = {0};
    star.m  = 1.16;
    star.r = 1.273 * 0.00465;
    reb_simulation_add(sim, star);

    // Planets
    double mearth = 3e-6;
    double ri = 0.05 * M_PI/180.;

    // b
    double mb = 6.89 * mearth;
    double eb = 0.07;
    double ab = 0.1283;
    double ib = 0.0;//reb_random_uniform(sim, -1. * ri, ri);
    double Mb = reb_random_uniform(sim, 0, 2 * M_PI);

    double mc = 4.4 * mearth;
    double ec = 0.04;
    double ac = 0.2061;
    double ic = 0.0;//reb_random_uniform(sim, -1. * ri, ri);
    double Mc = reb_random_uniform(sim, 0, 2 * M_PI);

    double md = 4.6 * mearth;
    double ed = 0.06;
    double ad = 0.88;
    double id = 0.0;//reb_random_uniform(sim, -1. * ri, ri);
    double Md = reb_random_uniform(sim, 0, 2 * M_PI);

    me = reb_random_uniform(sim, 12. - 5., 12. + 5.) * mearth;
    double ee = 0.14;
    double ae = 1.06;//reb_random_uniform(sim, 1.06 - 0.02, 1.06 + 0.03);
    double ie = reb_random_uniform(sim, -1. * ri, ri);
    double Me = reb_random_uniform(sim, 0, 2 * M_PI);

    // This is the one we care abotu
    mf = reb_random_uniform(sim, 12. - 3., 12. + 3.) * mearth;
    double ef = 0.004;
    double af = 1.37;//reb_random_uniform(sim, 1.37 - 0.02, 1.37 + 0.02);
    double incf = reb_random_uniform(sim, -1. * ri, ri);
    double Mf = reb_random_uniform(sim, 0, 2 * M_PI);

    //double rhof = 1.0 * pow(1.496e13,3) / (1.989e33); // 1 g/cm^3 to solar masses/AU^3
    //rf =  pow((3 * mf / (4 * M_PI * rhof)), 1./3.);

    double planet_as[10] = {ab, ac, ad, ae, af};

    //reb_simulation_add_fmt(sim, "primary m a e inc M", star, mb, ab, eb, ib, Mb);
    //reb_simulation_add_fmt(sim, "primary m a e inc M", star, mc, ac, ec, ic, Mc);
    //reb_simulation_add_fmt(sim, "m a e inc M", md, ad, ed, id, Md);
    reb_simulation_add_fmt(sim, "m a e inc M", me, ae, ee, ie, Me);
    reb_simulation_add_fmt(sim, "m a e inc M", mf, af, ef, incf, Mf);

    //double rin = 1.5 * rf;
    //double rho_ring = 0.5;
    //double rout = 2.45 * pow(1./rho_ring, 1./3.) * rf;
    //for (unsigned int i = 0; i < ntest; i++){
    //reb_add_fmt(sim, "primary a inc Omega", sim->particles[2], rout, 30.*M_PI/180.);
    //}

    reb_simulation_move_to_com(sim);

    //FILE* of = fopen(title, "w");
    /*fprintf(of, "t,nx\n");
    fclose(of);

    system(title_remove);
    FILE* of = fopen(title, "w");
    fprintf(of, "t");
    //for (unsigned int i = 0; i < ntest; i++){
    fprintf(of, ",at,it");
    //}
    fprintf(of, "\n");
    fclose(of);
    */
    struct reb_orbit o = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    tmax = o.P * 1e8;
    sim->dt = o.P / 15.12345;
    reb_simulation_integrate(sim, tmax);
/*
    for (unsigned i = 0; i < 5; i++){
      struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[i + 1], sim->particles[0]);
      if (fabs(orb.a - planet_as[i]) / planet_as[i] > 0.1){
        //printf("Unstable %d %f\n", i+1,orb.a);
        stable=0;
        system(title_remove);
        break;
      }
    }
*/
    FILE* sf = fopen(title_stats, "a");
    fprintf(sf, "%d,%f,%f,%d,%d,%f\n",ind,me,mf,1,-1,sim->t/tmax);
    fclose(sf);

    //rebx_free(rebx);
    reb_simulation_free(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file

    if(reb_simulation_output_check(sim, 1. * 2 * M_PI)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;

      for (unsigned int i = 1; i < sim->N; i++){
        struct reb_orbit o = reb_orbit_from_particle(sim->G, sim->particles[i], sim->particles[0]); // planet orbit
        double a = o.a;

        if (a < (planet_as[i-1] - 0.1 * planet_as[i-1]) || a > (planet_as[i-1] + 0.1 * planet_as[i-1])){
          FILE* sf = fopen(title_stats, "a");
          fprintf(sf, "%d,%f,%f,%d,%d,%f\n",ind,me,mf,0,i,sim->t/tmax);
          fclose(sf);
          exit(1);
        }
      }
    }


    //if(reb_simulation_output_check(sim, 100.)){        // outputs to the screen
    //    reb_simulation_output_timing(sim, tmax);
    //}
}
