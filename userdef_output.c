#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{

   int i,j,k;
   double ***tmp,T,mu,***beta, ***p_mag,mu0;
   double ***p, ***rho, ***Bx, ***By;
   double ***CS, ***CA, ***p_fis, ***rho_fis;
   double ***Wave_energy_flux, ***vx, ***vy;
   double ***mag_Vel,***Qrad,***Kinetic_energy_density;
   
   int    klo, khi, kmid;
   static int ntab;
   double  Tmid, scrh, dT;
   static double *L_tab, *T_tab, E_cost;
   
   double nH;
     

   FILE *fcool;

   tmp = GetUserVar("tmp");
   //p_fis = GetUserVar("p_fis");
   //rho_fis = GetUserVar("rho_fis");
   //mag_Vel = GetUserVar("mag_Vel");
   
   #if PHYSICS == MHD 
   beta = GetUserVar("beta");
   //p_mag = GetUserVar("p_mag");
   //CS = GetUserVar("CS");
   CA = GetUserVar("CA");
   Wave_energy_flux = GetUserVar("Wave_energy_flux");
   Kinetic_energy_density = GetUserVar("Kinetic_energy_density");
   
   Bx  = d->Vc[BX1];
   By  = d->Vc[BX2];
   
   #endif
   
    #if COOLING == TABULATED
    Qrad =  GetUserVar("Qrad");
	#endif
    
   
   
   rho = d->Vc[RHO];
   p   = d->Vc[PRS];
   vx  = d->Vc[VX1];
   vy  = d->Vc[VX2];
   
   mu = 0.6;
   mu0 = 1.25664e-6; 
   
  
    #if COOLING == TABULATED
   
/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

  if (T_tab == NULL){
    FILE *fcool;
    printLog (" > Reading table from disk...\n");
    fcool = fopen("CHIANTI_coronal.dat","r");
    if (fcool == NULL){
      printLog ("! Radiat: cooltable.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, 
                                       L_tab + ntab)!=EOF) {
      ntab++;
    }
    E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
    fclose(fcool);
  }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  DOM_LOOP(k,j,i){
	  
  klo = 0;
  khi = ntab - 1;

  tmp[k][j][i] = (p[k][j][i]/rho[k][j][i])*KELVIN*mu;
   
  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Tmid = T_tab[kmid];
    if (tmp[k][j][i] <= Tmid){
      khi = kmid;
    }else if (tmp[k][j][i] > Tmid){
      klo = kmid;
    }
  }

   if(tmp[k][j][i]<1.e4){
	   Qrad[k][j][i] = 0.0; 
	   }
	   
	 else { 
	 
	  nH = UNIT_DENSITY/CONST_amu*H_MASS_FRAC/CONST_AH*rho[k][j][i]; 	 
	  dT        = T_tab[khi] - T_tab[klo];
	  scrh      = L_tab[klo]*(T_tab[khi]-tmp[k][j][i])/dT + L_tab[khi]*(tmp[k][j][i]-T_tab[klo])/dT;	 
	  Qrad[k][j][i] = -nH*nH*scrh*E_cost;
		}
	}
	#endif


    DOM_LOOP(k,j,i){
		
    tmp[k][j][i] = (p[k][j][i]/rho[k][j][i])*KELVIN*mu;
    
	//p_fis[k][j][i] = p[k][j][i];
	
	//rho_fis[k][j][i] = 1.e-12*rho[k][j][i];
	
	//mag_Vel[k][j][i] = 1e6*sqrt(pow(vx[k][j][i],2)+pow(vy[k][j][i],2));
	
	#if PHYSICS == MHD 
	
    beta[k][j][i] = (2.0*p[k][j][i])/(pow(Bx[k][j][i],2)+pow(By[k][j][i],2));

    //p_mag[k][j][i] = 0.5*(pow(Bx[k][j][i],2)+pow(By[k][j][i],2));
	
	//CS[k][j][i] = sqrt((g_gamma*p[k][j][i])/(1.e-12*rho[k][j][i]));
	
	CA[k][j][i] = (11.21e-4*sqrt(pow(Bx[k][j][i],2)+pow(By[k][j][i],2)))/(sqrt(mu0*1.e-12*rho[k][j][i]));
	
	Wave_energy_flux[k][j][i] = 0.5*(1.e-12*rho[k][j][i])*CA[k][j][i]*1e6*sqrt(pow(vx[k][j][i],2)+pow(vy[k][j][i],2));
	
	Kinetic_energy_density[k][j][i] = 0.5*(1.e-12*rho[k][j][i])*1e6*sqrt(pow(vx[k][j][i],2)+pow(vy[k][j][i],2));
	#endif

    }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





