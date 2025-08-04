/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Define the components of the diagonal resistive tensor. 

  Use this function to supply the resistivity in the three directions
  \f$ \eta_{x1}\f$, \f$ \eta_{x2}\f$ and \f$ \eta_{x3}\f$.
  
  \authors A. Mignone (mignone@ph.unito.it)\n
  \date    Feb 26, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Resistive_eta(double *v, double x1, double x2, double x3,
                   double *J, double *eta)
/*!
 * Compute the resistive tensor components as function of the primitive
 * variables, coordinates and currents.
 *
 * \param [in]  v    array of primitive variables
 * \param [in]  x1   coordinate in the X1 direction
 * \param [in]  x2   coordinate in the X2 direction
 * \param [in]  x3   coordinate in the X3 direction
 * \param [in]  J    current components, J[IDIR], J[JDIR], J[KDIR]
 * \param [out] eta  an array containing the three components of
 *                   \f$ \tens{\eta}\f$.
 *
 *********************************************************************** */
{
	
	double r;
	
	double hTR   = g_inputParam[heta];
	double eta_0 = g_inputParam[eta0];
	double w_eta = g_inputParam[weta];
        double t = g_time;
        double t_res = g_inputParam[tres];
	
	r = sqrt(pow(x1,2)+pow(x2-hTR,2));

  if (t >= 0.0  && t < t_res){ 
  eta[IDIR] = eta_0*exp(-pow(r/w_eta,2)); 
  eta[JDIR] = eta_0*exp(-pow(r/w_eta,2));
  eta[KDIR] = eta_0*exp(-pow(r/w_eta,2));
   } else {
  eta[IDIR] = 0.0;
  eta[JDIR] = 0.0;
  eta[KDIR] = 0.0;
  }
}
