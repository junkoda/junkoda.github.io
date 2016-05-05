//
// Computes n(z) for given HOD parameter
//

#include <iostream>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>

#include <gsl/gsl_integration.h>

#include "msg.h"
#include "const.h"      // -> delta_c
#include "cosmology.h"
#include "growth.h"
#include "power.h"
#include "sigma.h"
#include "mf.h"
#include "hod.h"
#include "nbar.h"

using namespace std;
using namespace boost::program_options;

static double integrand_n_hod(double nu, void* params);


NbarIntegration::NbarIntegration()
{

}

NbarIntegration::NbarIntegration(NbarIntegration const * const ni)
{
  // use same setting as ni (for different redshift)
  w= ni->w;
  s= ni->s;
  mf= ni->mf;
  hod= ni->hod;
  rho_m= ni->rho_m;
  D= 0.0;
  z= 0.0;
}

NbarIntegration* nbar_integration_alloc(PowerSpectrum const * const ps,
					 Hod* const hod)
{
  NbarIntegration* const ni= new NbarIntegration();
  
  const double M_min= 1.0e10;
  const double M_max= 1.0e16;

  ni->w= gsl_integration_cquad_workspace_alloc(100);
  ni->s= new Sigma(ps, M_min, M_max);
  ni->mf= mf_alloc();
  ni->hod= hod;
  ni->rho_m= cosmology_rho_m();
  ni->D= 0.0;
  ni->z= 0.0;

  msg_printf(msg_debug, "Allocated an nbar_integration\n");

  return ni;
}

void nbar_integration_free(NbarIntegration* const ni)
{
  gsl_integration_cquad_workspace_free(ni->w);
  delete ni->s;
  mf_free(ni->mf);

  delete ni;
}

void nbar_read(const char filename[], const double z_min, const double z_max,
	       vector<Nbar>& v)
{
  cerr << "nbar file " << filename << endl;
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Error: unable to open nbar file " << filename << endl;
    throw filename;
  }

  char buf[256];
  while(fgets(buf, 255, fp)) {
    Nbar n; n.nbar= 0;
    if(buf[0] == '#') continue;

    int ret= sscanf(buf, "%le %le", &n.z, &n.nbar); assert(ret == 2);

    if(z_min <= n.z && n.z < z_max)
      v.push_back(n);
  }

  fclose(fp);

  // set error in the chi2 fitting
  for(vector<Nbar>::iterator p= v.begin(); p != v.end(); ++p)
    p->dnbar= p->nbar; // relative error
}

double integrand_n_hod(double nu, void* params)
{
  // Returns number density dn/dnu of HOD parameters
  // requirement: D computed for z
  //              tinker_set_redshift(z)
  //
  // dn= P(n|HOD)*n(M) dnu
  // nu = delta_c/sigma(M,z), sigma(M,z) = D(z) sigma0(M)
  // rho_m

  NbarIntegration* const ni= (NbarIntegration*) params;
  
  const double sigma0= delta_c/(ni->D*nu);
  const double M= ni->s->M(sigma0);

 return ni->hod->ncen(M)*(1.0 + ni->hod->nsat(M))*
        mf_f(ni->mf, nu)*ni->rho_m/M;
}

/*
void nbar_compute_nz(const double c[], vector<Nbar>* const v)
{
  // Input: parameters c[10]
  // Ouput: v

  Hod hod;
  // Computes nbar(z) from given parametrisation c[10]
  for(vector<Nbar>::iterator p= v->begin(); p != v->end(); ++p) {
    double z= p->z;
  
    hod_compute(c, z, &hod);           // c, z   --> hod(z; c)
    p->nbar= compute_n_hod(&hod, z);   // hod(z) --> nbar(z)
  }  
}
*/


double nbar_compute(NbarIntegration* const ni, const double z)
{
  // Computes number density of HOD galaxies from HOD parameters 'hod'
  // n = \int P(n|HOD)*n(M) dnu

  const double a= 1.0/(1.0 + z);

  ni->hod->compute_param_z(z);


  if(ni->D == 0.0 || ni->z != z) {
    mf_set_redshift(ni->mf, a);    
    ni->D= growth_D(a);
    ni->z= z;
    msg_printf(msg_verbose,
	       "setting nbar_integration at z=%.3f, D=%.4f, mf->alpha=%.4f\n",
	       ni->z, ni->D, ni->mf->alpha);
  }

  gsl_function F;
  F.function= &integrand_n_hod;
  F.params= (void*) ni;

  const double nu_min= delta_c*ni->s->sinv_min/ni->D;
  const double nu_max= delta_c*ni->s->sinv_max/ni->D;

  double result;
  gsl_integration_cquad(&F, nu_min, nu_max, 1.0e-5, 1.0e-5, ni->w,
			&result, 0, 0);

  return result;
}

