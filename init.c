/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  2D Macrospicule formation with thermal condiuction and radiative thin cooling from a function of CHIANTI coronal abundances
  \author J.J. Gonzalez-Aviles (jjgonzalez@igeofisica.unam.mx)
  \date   July 15, 2020
  \last update January 15, 2021
  \references   Singh, B., Sharma, K., and Srivastava, A. K., 2019, Ann. Geophys., 37, 891
	        	Muraswski, K., Srivastava, A. K., and Zaqarashvilli, T. V., 2011, A&A, 535, A58

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
 
    double B0,a,b,Av,x0,y0,sigma;
    
    
	Av = 0.1;
	x0 = 5.0;
	y0 = 1.75;
	sigma = 0.15;
    

    g_gamma = 1.66;
    
    
    v[VX1] = 0.0; 
    v[VX2] = Av*exp((-pow(x1-x0,2)-pow(x2-y0,2))/pow(sigma,2));
    
	
	#if PHYSICS == MHD
    B0 = 3000.0;
    a = 0.0;
    b = -40.0;
	
    v[BX1] = (-2.0*B0*(x1-a)*(x2-b))/(pow(pow(a-x1,2)+pow(b-x2,2),2));
    v[BX2] = (B0*(pow(a-x1,2)-pow(b-x2,2)))/(pow(pow(a-x1,2)+pow(b-x2,2),2));
	#endif
	

	#if COOLING == TABULATED
	g_minCoolingTemp = 10000.0; 
    #endif    
    
    
    
    
  
   
}

#if PHYSICS == MHD
 //Define magnetic field functions to use in boundary conditions
/* ********************************************************************* */
void mag_field (double x1, double x2, double x3, double mag[2])
/*
 *
 *
 *********************************************************************** */
{
    
    double a,b,B0;
    
    B0 = 3000.0;
    a = 0.0;
    b = -40.0;
    
    
    mag[0] = (-2.0*B0*(x1-a)*(x2-b))/(pow(pow(a-x1,2)+pow(b-x2,2),2));
    mag[1] = (B0*(pow(a-x1,2)-pow(b-x2,2)))/(pow(pow(a-x1,2)+pow(b-x2,2),2));
    
    
}
#endif

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
    FILE *fpTemp,*fprho,*fpprs;
    char ch;
    int i,j,k,before_z1,ii,m=0,mm=0,ll=0,id;
    int num_datos,num_datos_rho,num_datos_prs;
    float zC7[801],TempC7[801];
    float zC7_rho[801],rhoC7[801];
    float zC7_prs[801],prsC7[801];
    double Te1,Te2,x1C7,x2C7,Te,Te_adi;
    double rho1,rho2,rhoe,rhoe_adi;
    double prs1,prs2,prse;
    double K_B,mu,m_p,time;
    double Ap,x0,y0,sigma;
    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
    
    num_datos = 801;
    num_datos_rho = 801;
    num_datos_prs = 801;
    
    /* Parameter of the pressure pulse */
    Ap = 12.0;
    x0 = 0.0;
    y0 = 1.8;
    sigma = 0.2;
    
    m_p = 1.672e-27;
    K_B = 1.380e-23;
    mu = 0.6;
	
	
       fprho = fopen("/home/javiles/2D_Macrospicules_thermal_conduction_rad_thin_cooling_CHIANTI_PLUTO/radiative_cooling_case/PLUTO4.4_case/rho_C7_40Mm.dat","r");
	   fpprs = fopen("/home/javiles/2D_Macrospicules_thermal_conduction_rad_thin_cooling_CHIANTI_PLUTO/radiative_cooling_case/PLUTO4.4_case/p_C7_40Mm.dat","r");


		 while(fscanf(fprho, "%f %f", zC7_rho+mm, rhoC7+mm) == 2)
		 {
			 mm++;
		 }
    
		 while(fscanf(fpprs, "%f %f", zC7_prs+ll, prsC7+ll) == 2)
		 {
			 ll++;
		 }

    
    /*Linear interpolation of data to PLUTO grid*/
		 
		 TOT_LOOP(k,j,i){
   
			for(ii=0;ii<num_datos_rho;ii++){
				if(zC7_rho[ii]>x2[j]){
				before_z1 = ii-1;
				break;
						}
				}
        
				x1C7 = zC7_rho[before_z1];
				x2C7 = zC7_rho[before_z1+1];
				rho1 = rhoC7[before_z1];
				rho2 = rhoC7[before_z1+1];
        
        /*Density profile interpolated to PLUTO grid*/
        
		rhoe = ((rho2-rho1)/(x2C7-x1C7))*(x2[j]-x1C7) + rho1;
        
		rhoe_adi = (rhoe)/(1.e-15);
        
        
		for(ii=0;ii<num_datos_prs;ii++){
			if(zC7_prs[ii]>x2[j]){
				before_z1 = ii-1;
				break;
			}
		}
           
		x1C7 = zC7_prs[before_z1];
		x2C7 = zC7_prs[before_z1+1];
		prs1 = prsC7[before_z1];
		prs2 = prsC7[before_z1+1];
           
        /*Pressure profile interpolated to PLUTO grid*/
           
		prse = ((prs2-prs1)/(x2C7-x1C7))*(x2[j]-x1C7) + prs1;
		   
        
        /*Calculating the initial conditions for density and pressure */
        
        d->Vc[RHO][k][j][i] = rhoe_adi;
        
        d->Vc[PRS][k][j][i] = prse;
        
        //d->Vc[PRS][k][j][i] = prse*(1.0 + Ap*exp((-pow(x1[i]-x0,2)-pow(x2[j]-y0,2))/pow(sigma,2)));
		
	
	}
    
    fclose(fprho);
    fclose(fpprs);
	
	 
    
    
}


/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*!
 *    Fixed boundary conditions
 *********************************************************************** */
{
    int i,j,k,id;
    double den,prs,xb_Bx1,xb_Bx2,xb_Bx3;
    double xb_vx1,xb_vx2,xb_vx3,xb_Bx1s;
    FILE *fp_rho,*fp_prs;
    double xx,px1,px2,px3;
    double mag[3],v[256];
	double Ap,x0,y0,sigma;
    double V0[8][36][16][506];
    double *x1l = grid->xl[IDIR];
    double *x1r = grid->xr[IDIR];
    double *x2l = grid->xl[JDIR];
    double *x2r = grid->xr[JDIR];
    double *x3l = grid->xl[KDIR];
    double *x3r = grid->xr[KDIR];
    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
	
    int before_z1,ii,m=0,mm=0,ll=0;
    int num_datos_rho,num_datos_prs,nv;
    float zC7_rho[401],rhoC7[401];
    float zC7_prs[401],prsC7[401];
    double rho1,rho2,rhoe,rhoe_adi;
    double prs1,prs2,prse,x1C7,x2C7;
	
     num_datos_rho = 401;
     num_datos_prs = 401;
	 
 		Ap = 0.0;
 		x0 = 5.0;
 		y0 = 1.75;
 		sigma = 0.2;
	 
	 
    
		 
		if(g_time == 0.0){ 
		 
 	 	id = InputDataOpen("/home/javiles/initial_conditions_C7/rho_C7_40Mm_1200x1200.dbl","/home/javiles/2D_Macrospicules_thermal_conduction_rad_thin_cooling_CHIANTI_PLUTO/radiative_cooling_case/PLUTO4.4_case/grid.out"," ",0,CENTER);
        
 	        TOT_LOOP(k,j,i){
               
 		d->Vc[RHO][k][j][i] = InputDataInterpolate (id,x1[i],x2[j],x3[k]);
		
 		}
		
 		InputDataClose(id);


 		id = InputDataOpen("/home/javiles/initial_conditions_C7/prs_C7_40Mm_1200x1200.dbl","/home/javiles/2D_Macrospicules_thermal_conduction_rad_thin_cooling_CHIANTI_PLUTO/radiative_cooling_case/PLUTO4.4_case/grid.out"," ",0,CENTER);
       
 		TOT_LOOP(k,j,i){
               
		d->Vc[PRS][k][j][i] = InputDataInterpolate (id,x1[i],x2[j],x3[k]);
		
		//d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i]*(1.0+Ap*exp((-pow(x1[i]-x0,2)-pow(x2[j]-y0,2))/pow(sigma,2)));

 		}

 		InputDataClose(id);
		 
	     TOT_LOOP(k,j,i){
			 
			#if PHYSICS == MHD 
			mag_field(x1[i],x2[j],x3[k],mag);
			
			d->Vc[BX1][k][j][i] = mag[0];
			d->Vc[BX2][k][j][i] = mag[1];
			
			V0[BX1][k][j][i] = d->Vc[BX1][k][j][i];
			V0[BX2][k][j][i] = d->Vc[BX2][k][j][i];
			#endif
			
			#ifdef GLM_MHD
			V0[PSI_GLM][k][j][i] = d->Vc[PSI_GLM][k][j][i];
			#endif 
				
			V0[RHO][k][j][i] = d->Vc[RHO][k][j][i];
			V0[PRS][k][j][i] = d->Vc[PRS][k][j][i];
			
			
		 
	     }
		 	    
		 
	 }

	
   
    
        if (side == X1_BEG){ // X1_BEG
			
		if (box->vpos == CENTER) {

		 BOX_LOOP(box,k,j,i){
			 
		     #if PHYSICS == MHD 
 			 V0[BX1][k][j][i] = d->Vc[BX1][k][j][i];
 			 V0[BX2][k][j][i] = d->Vc[BX2][k][j][i];
			 #endif
			 
			 #ifdef GLM_MHD
			 V0[PSI_GLM][k][j][i] = d->Vc[PSI_GLM][k][j][i];
			 #endif 

 			 V0[RHO][k][j][i] = d->Vc[RHO][k][j][i];
 			 V0[PRS][k][j][i] = d->Vc[PRS][k][j][i];
			 
			 //d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
			 //d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
			 
			 //d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IBEG];
			 //d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IBEG];

			 d->Vc[VX1][k][j][i] = 0.0;
			 d->Vc[VX2][k][j][i] = 0.0;
			
		 
		 		}
			}
			
		}
            
       
    
        if (side == X1_END){ // X1_END
            
        if (box->vpos == CENTER) {
        
   		 BOX_LOOP(box,k,j,i){
			
			    #if PHYSICS == MHD 
 			 	V0[BX1][k][j][i] = d->Vc[BX1][k][j][i];
 				V0[BX2][k][j][i] = d->Vc[BX2][k][j][i];
				#endif
				
				#ifdef GLM_MHD
				V0[PSI_GLM][k][j][i] = d->Vc[PSI_GLM][k][j][i];
				#endif 
				
   			 	d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IEND];
   			 	d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IEND];
				
   			 	d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND];
   			 	d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IEND];
				
    			//V0[RHO][k][j][i] = d->Vc[RHO][k][j][i];
    			//V0[PRS][k][j][i] = d->Vc[PRS][k][j][i];
				
    			//d->Vc[VX1][k][j][i] = 0.0;
    			//d->Vc[VX2][k][j][i] = 0.0;
			

        	}
		
		}
		
	}
        
    
    
        if (side == X2_BEG){ // X2_BEG
			
        if (box->vpos == CENTER) {
			
   		 BOX_LOOP(box,k,j,i){
			 
			 #if PHYSICS == MHD 
			 V0[BX1][k][j][i] = d->Vc[BX1][k][j][i];
			 V0[BX2][k][j][i] = d->Vc[BX2][k][j][i];
			 #endif
			 
			 #ifdef GLM_MHD
 			 V0[PSI_GLM][k][j][i] = d->Vc[PSI_GLM][k][j][i];
			 #endif 

			 V0[RHO][k][j][i] = d->Vc[RHO][k][j][i];
			 V0[PRS][k][j][i] = d->Vc[PRS][k][j][i];

			 d->Vc[VX1][k][j][i] = 0.0;
			 d->Vc[VX2][k][j][i] = 0.0;
			

				}
            }
		}
    
    
        if (side == X2_END){ // X2_END
            
        if (box->vpos == CENTER) {
       
   		 BOX_LOOP(box,k,j,i){
			 
		     #if PHYSICS == MHD 
			 V0[BX1][k][j][i] = d->Vc[BX1][k][j][i];
			 V0[BX2][k][j][i] = d->Vc[BX2][k][j][i];
			 #endif
			 
			 #ifdef GLM_MHD
 			 V0[PSI_GLM][k][j][i] = d->Vc[PSI_GLM][k][j][i];
			 #endif 

			 V0[RHO][k][j][i] = d->Vc[RHO][k][j][i];
			 V0[PRS][k][j][i] = d->Vc[PRS][k][j][i];

			 d->Vc[VX1][k][j][i] = 0.0;
			 d->Vc[VX2][k][j][i] = 0.0;
			
    			
	
    				}
        
    			}
          	}
        
  

}

#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  
   
    g[IDIR] = 0.0;
    g[JDIR] = -2.74e-4; //(-27400.0*UNIT_LENGTH)/(pow(UNIT_VELOCITY,2));
    g[KDIR] = 0.0;
}
#endif
