/*
 * ===========================================================================
 *
 *       Filename:  spatcontrol.c
 *
 *    Description:  C functions for spatcontrol.R
 *    		should be compiled by:
 *    		R CMD SHLIB spatcontrol.c
 *
 *        Version:  1.0
 *        Created:  06/10/2011 02:21:15 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  C. Barbu (CB), corentin.barbu@gmail.com
 *        Company:  CCEB, University of Pennsylvania
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


/*************** floating code start conversions to utm *********/
#include <math.h>
#include <float.h>

#define EPSILON         (DBL_EPSILON)
#define DEG_TO_RAD      (M_PI / 180.0)      /* degrees to radians */
#define RAD_TO_DEG      (180.0 / M_PI)      /* radians to degrees */

/* for checking the equality of floating point numbers */
#define FABS(n)                                                          \
                (((n) < 0) ? -(n) : (n))

#define DBL_EQ(n1, n2)                                                   \
                (((n1) == 0 && (n2) == 0)                                \
                 || (((n1) != 0)                                         \
                     && FABS((n1) - (n2))/FABS((n1)) <= DBL_EPSILON)     \
                 || FABS((n1) - (n2)) <= DBL_EPSILON)
#define DBL_LTEQ(n1, n2)                                                 \
                (((n1) < (n2)) || DBL_EQ((n1), (n2)))
#define DBL_GTEQ(n1, n2)                                                 \
                (((n1) > (n2)) || DBL_EQ((n1), (n2)))
#define DBL_LT(n1, n2)                                                   \
                (((n1) < (n2)) && !DBL_EQ((n1), (n2)))
#define DBL_GT(n1, n2)                                                   \
                (((n1) > (n2)) && !DBL_EQ((n1), (n2)))

/* the return types from UTM<-->LL conversion routines */
struct pair
{
  double x;
  double y;
};


/*-----------------------------------------------------------------------------
  lonlat_to_utm:
    Convert a point from Lat/Lon to UTM. All angles are in radians.
    Memory should be allocated for the structure before the function
    is called. 

  Author:
    Chris Grandin

  Algorithm Source:
    National Mapping Agency of Great Britain Ordnance Survey 
    <http://www.ordsvy.gov.uk>
  ---------------------------------------------------------------------------*/
void lonlat_to_utm(double lon, double lat, double *eastingNorthing_x, 
                   double *eastingNorthing_y, int utmZone);

/*-----------------------------------------------------------------------------
  utm_to_lonlat:
    Convert a point from UTM to Lat/Lon. All angles are in radians.
    Memory should be allocated for the structure before the function
    is called.

  Author:
    Chris Grandin

  Algorithm Source:
    National Mapping Agency of Great Britain Ordnance Survey 
    <http://www.ordsvy.gov.uk>
  ---------------------------------------------------------------------------*/
void utm_to_lonlat(double easting, double northing, double *lonlat_x,
                   double *lonlat_y, int utmZone,double correctif_south);


/* constants of the earth (Ellipsoidal) and conversion factors */
#define A_PRE   6378137.0                 /* major axis of ellipsoid (WGS 84)*/
#define B_PRE   6356752.3142              /* minor axis of ellipsoid (WGS 84)*/
#define SF      0.9996                    /* scale factor */
#define AXE_MAJOR       (A_PRE * SF)              /* major after scaled */
#define AXE_MINOR       (B_PRE * SF)              /* minor after scaled */
#define AXES_RATIO       ((AXE_MAJOR - AXE_MINOR) / (AXE_MAJOR + AXE_MINOR))       /* ratio used in calculations */
#define E       sqrt(((AXE_MAJOR*AXE_MAJOR)-(AXE_MINOR*AXE_MINOR))/(AXE_MAJOR*AXE_MAJOR)) /* eccentricity */
#define E0      500000.0                  /* grid eastings of true origin */

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

/*-----------------------------------------------------------------------------
  get_nu:
   
  Description:
  - get a radius of curvature at latitude 'phi' perpendicular to a meridian
  ---------------------------------------------------------------------------*/
static double get_nu(double phi)
{
  double sin1, sin2;

  sin1 = sin(phi);
  sin2 = sin1 * sin1;

  return (AXE_MAJOR / (sqrt(1.0 - (E * E * sin2))));
}

/*-----------------------------------------------------------------------------
  get_rho:
   
  Description:
  - get a radius of curvature of a meridian at latitude 'phi'
  ---------------------------------------------------------------------------*/
static double get_rho(double nu, double phi)
{
  double sin1, sin2;

  sin1 = sin(phi);
  sin2 = sin1 * sin1;

  return ((nu * (1.0 - E * E)) / (1.0 - (E * E * sin2)));
}

/*-----------------------------------------------------------------------------
  get_eta:
  ---------------------------------------------------------------------------*/
static double get_eta(double nu, double rho)
{
  return (sqrt((nu / rho) - 1.0));
}

/*-----------------------------------------------------------------------------
  get_M:

  Description:
  - get developed arc of a meridian from 'phi' to zero
  ---------------------------------------------------------------------------*/
static double get_M(double phi)
{
  double c1, c2, c3, c4, n2, n3, sin1, cos1, M;

  c1 = 5.0/4.0;
  c2 = 21.0/8.0;
  c3 = 15.0/8.0;
  c4 = 35.0/24.0;
  n2 = AXES_RATIO*AXES_RATIO;
  n3 = n2*AXES_RATIO;
  sin1 = sin(phi);
  cos1 = cos(phi);
  M =     ((1.0 + AXES_RATIO + (c1*n2) + (c1*n3)) * phi);
  M = M - (((3.0*AXES_RATIO) + (3.0*n2) + (c2*n3)) * sin(phi) * cos(phi));
  M = M + (((c3*n2) + (c3*n3)) * sin(2.0*phi) * cos(2.0*phi));
  M = M - ((c4*n3) * sin(3.0*phi) * cos(3.0*phi));

  return (AXE_MINOR*M);
}

/*-----------------------------------------------------------------------------
  lonlat_to_utm:
  ---------------------------------------------------------------------------*/
void lonlat_to_utm(double lon, double lat, double *eastingNorthing_x, 
                   double *eastingNorthing_y, int utmZone)
{
	// adapatation by CB to get same as usual utm
	lon = DEG_TO_RAD * lon;
	lat = DEG_TO_RAD * lat;
	int correctif_south = 0; 
	if ( lat <0 ) {
		correctif_south = 10000000.;
	}
	// printf("lon RAD: %f ; lat RAD: %f\n",lon,lat);

	// official code
  double sin1, cos1, cos3, cos5, tan1, tan2, tan4, nu, rho, eta, eta2;
  double P, P2, P3, P4, P5, P6, I, II, III, IIIA, IV, V, VI;
  int lonOrig = -177 + ((utmZone - 1) * 6);

  sin1 = sin(lat);
  cos1 = cos(lat);
  cos3 = cos1*cos1*cos1;
  cos5 = cos3*cos1*cos1;
  tan1 = tan(lat);
  tan2 = tan1*tan1;
  tan4 = tan2*tan2;
  nu = get_nu(lat);
  rho = get_rho(nu,lat);
  eta = get_eta(nu,rho);
  eta2 = eta*eta;
  P = lon - (lonOrig * DEG_TO_RAD);
  P2 = P*P;
  P3 = P2*P;
  P4 = P3*P;
  P5 = P4*P;
  P6 = P5*P;
  I = get_M(lat);
  II = (nu/2.0)*sin1*cos1;
  III = (nu/24.0)*sin1*cos3*(5.0-tan2+9.0*eta2);
  IIIA = (nu/720.0)*sin1*cos5*(61.0-58.0*tan2+tan4);
  IV = nu*cos1;
  V = (nu/6.0)*cos3*((nu/rho)-tan2);
  VI = (nu/120.0)*cos5*(5.0-18.0*tan2+tan4+14.0*eta2-58.0*tan2*eta2);

  *eastingNorthing_y = I+P2*II+P4*III+P6*IIIA + correctif_south;
  // *eastingNorthing_y = I+P2*II+P4*III+P6*IIIA;
  *eastingNorthing_x = E0+P*IV+P3*V+P5*VI;
}

/*-----------------------------------------------------------------------------
  utm_to_lonlat:
  ---------------------------------------------------------------------------*/

void utm_to_lonlat(double easting, double northing, double *lonlat_x,
                   double *lonlat_y, int utmZone, double correctif_south)
{
  double E1, E2, E3, E4, E5, E6, E7, lat1, M, nu, rho, eta, eta2;
  double tan1, tan2, tan4, tan6, sec1, nu3, nu5, nu7, VII, VIII;
  double IX, X, XI, XII, XIIA;
  int lonOrig = -177 + ((utmZone - 1) * 6);
  int i;

  E1 = easting - E0;
  E2 = E1*E1;
  E3 = E2*E1;
  E4 = E3*E1;
  E5 = E4*E1;
  E6 = E5*E1;
  E7 = E6*E1;
  lat1 = (northing)/(AXE_MAJOR);
  M = get_M(lat1);

  for(i=0;i<10;i++){
    lat1 += (northing-M)/AXE_MAJOR;
    M = get_M(lat1);
  }

  nu = get_nu(lat1);
  rho = get_rho(nu,lat1);
  eta = get_eta(nu,rho);
  eta2 = eta*eta;
  tan1 = tan(lat1);
  tan2 = tan1*tan1;
  tan4 = tan2*tan2;
  tan6 = tan4*tan4;
  sec1 = 1.0/(cos(lat1));
  nu3 = nu*nu*nu;
  nu5 = nu3*nu*nu;
  nu7 = nu5*nu*nu;
  VII = tan1/(2.0*rho*nu);
  VIII = tan1/(24.0*rho*nu3)*(5.0+3.0*tan2+eta2-9.0*tan2*eta2);
  IX = tan1/(720.0*rho*nu5)*(61.0+90.0*tan2+45.0*tan4);
  X = sec1/nu;
  XI = sec1/(6.0*nu3)*(nu/rho+2.0*tan2);
  XII = sec1/(120.0*nu5)*(5.0+28.0*tan2+24.0*tan4);
  XIIA = sec1/(5040.0*nu7)*(61.0+662.0*tan2+1320.0*tan4+720.0*tan6);

  /* 11 Jul 2003 [Nicholas Boers]
     Swapped the x and y, because values being returned were backwards. */
  *lonlat_x = (E1*X-E3*XI+E5*XII-E7*XIIA + (lonOrig * DEG_TO_RAD))*RAD_TO_DEG;
  *lonlat_y = (lat1 - E2*VII+E4*VIII-E6*IX)*RAD_TO_DEG;
  // printf("lon: %f; lat: %f\n",*lonlat_x,*lonlat_y);
}

/********************End of floating code for conversions to utm************************/
void multiple_ll_to_utm(double *lon, double *lat,int *ncoord, double *X, double *Y, int *utmZone){
	for (int i=0;i< *ncoord;i++){
		lonlat_to_utm(lon[i],lat[i],X+i,Y+i,*utmZone);
	}
}

void multiple_utm_to_ll(double *X, double *Y,int *ncoord, double *lon, double *lat, int *utmZone,double *correctif_south){
	for (int i=0;i< *ncoord;i++){
		utm_to_lonlat(X[i],Y[i],lon+i,lat+i,*utmZone,*correctif_south);
	}
}

// change all entries bigger than limit to 0
void filter_spam(double *entries, int *pnb_ent,double *limit){
	int i=0;
	for ( i = 0; i < *pnb_ent; i += 1 ) { 
		// printf("i %i ent %f limit %f \n",i,entries[i],*limit);
		if(entries[i]>*limit){
			entries[i]=0.;
		}
	}
}

// multiply by f the values of dm that corresponds to an entry in AS
void specific_multiply(int *pdim,double *entries_dm,int *pnb_ent_dm,int *rpoint_dm,int *colind_dm,double *entries_AS,int *pnb_ent_AS,int *rpoint_AS,int *colind_AS,double *f){
	int i=0;
	int k=0;
	int col=0;
	int row=1;
	double new_dist=0.;
	// printf("f:%f",*f);
	for ( i = 0; i < *pnb_ent_AS; i += 1 ) { 
		// printf("i:%i AS:%f colAS:%i \n",i,entries_AS[i],colind_AS[i]);
		 // find corresponding column and row
		 col=colind_AS[i];
		 while(rpoint_AS[row]<=i+1 && row<*pdim){
		 	row++;
		 }
		// find corresponding entry in dm
		k = rpoint_dm[row-1]-1;
		while(k<rpoint_dm[row]){
			if(colind_dm[k]==col){
				break;
			}
			k++;
		}
		// if(i>=1690 && i<1707){
		// 	printf("i %i col %i row %i k %i \n",i,col,row,k);
		// }
		
		// change the entry
		if(k==rpoint_dm[row]){
			printf("WARNING: dm doesn't have the entry [%i,%i], keep 0\n",row,col);
		}else{
			// printf("i:%i row:%i col:%i AS:%f dm:%f\n",i,row,col,entries_AS[i],entries_dm[k]);
			entries_dm[k]=entries_dm[k]* *f;
		}
	}
}

// for a given Delta recalculate Dmat
void DmatFromDelta(int *pdim,double *Delta,double *tr, double *entries_dm,int *pnb_ent_dm,int *rpoint_dm,int *colind_dm,double *entries_AS,int *pnb_ent_AS,int *rpoint_AS,int *colind_AS){	
	int i=0;
	int k=0;
	int col=0;
	int row=1;
	double new_dist=0.;
	for ( i = 0; i < *pnb_ent_AS; i += 1 ) { 
		// printf("i:%i AS:%f colAS:%i \n",i,entries_AS[i],colind_AS[i]);
		 // find corresponding column and row
		 col=colind_AS[i];
		 while(rpoint_AS[row]<=i+1 && row<*pdim){
		 	row++;
		 }
		// find corresponding entry in dm
		k = rpoint_dm[row-1]-1;
		while(k<rpoint_dm[row]){
			if(colind_dm[k]==col){
				break;
			}
			k++;
		}
		// change the entry
		if(k==rpoint_dm[row]){
			printf("WARNING: dm doesn't have the entry [%i,%i], keep 0\n",row,col);
		}else{
			// printf("i:%i row:%i col:%i AS:%f dm:%f\n",i,row,col,entries_AS[i],entries_dm[k]);
			new_dist=*Delta+entries_dm[k];
			if(new_dist>*tr){
				entries_dm[k]=0;
			}else{
				entries_dm[k]=new_dist;
			}
		}
	}
}

void resample_old(double *d,int *pd_l){
	// d vector of double 
	// d_l size of d
	int from=0;
	int to=0;
	double inter=0.;
	int d_l=*pd_l;
	
	// proceede by permutations
	int nb_perm=(d_l/2)+1; // minimum to be able to move everybody
	// printf("nb_perm: %i \n",nb_perm);
	int i=0;
	for ( i = 0; i < nb_perm; i++ ) { 
		 from=rand()%d_l;
		 to=rand()%d_l;
		 // printf("from %i to %i \n",from,to);
		 double inter=d[to];
		 d[to]=d[from];
		 d[from]=inter;
	}
}

void resample(double *d,int *pd_l){
	// d vector of double 
	// d_l size of d
	int d_ind[*pd_l];
	double d_new[*pd_l];
	int loc=0;
	int i=0;
	int k=0;
	int num_loc=0;
	
	for ( i = 0; i < *pd_l; i += 1 ) { 
		 d_ind[i]=1; 
	}
	for ( i = *pd_l; i > 0; i -- ) { 
		 num_loc=rand()%i;
		 loc=-1;
		 k=0;
		 while(k<num_loc+1 && k<*pd_l){
			 loc++;
			 k=k+d_ind[loc];
		 }
		 // printf("i: %i d_l: %i num_loc: %i k: %i loc:%i\n",i,*pd_l,num_loc,k,loc);
		 d_ind[loc]=0;
		 d_new[i-1]=d[loc];
	}
	for ( i = 0; i < *pd_l; i += 1 ) { 
		d[i]=d_new[i];
	}
}

void resample_spam_entries_by_row(double *entries, int *pnb_ent,int *rpoint,int*pnb_rows){
	int r=0;
	int nb_ent_row=0;
	for(r=0;r<*pnb_rows;r++){
		nb_ent_row = rpoint[r+1]-rpoint[r];
		//printf("nb_ent_row %i\n",nb_ent_row);
		resample(&entries[rpoint[r]-1],&nb_ent_row);
	}
}

void symmetric_resample_spam_entries_by_row(double *entries,int *coli, int *pnb_ent,int *rpoint,int*pnb_rows){
	int error=0;
	int r=0;
	int nb_ent_row=0;
	// order of raws to be resampled
	double rows_reordered[*pnb_rows];
	for(r=0;r<*pnb_rows;r++){
		rows_reordered[r]=r; // r codes rows not like R it begins at 0
	}
	resample(rows_reordered,pnb_rows);
	// resample raws
	double Free[*pnb_ent]; // 1 if entry hasn't been changed yet
	int i=0;
	for ( i = 0; i < *pnb_ent; i += 1 ) { 
		Free[i]=1;
	}
	// printf("prep ok\n");
	int draw=0;
	for(draw=0;draw<*pnb_rows;draw++){
		r=(int)rows_reordered[draw];
		nb_ent_row = rpoint[r+1]-rpoint[r];
		//printf("nb_ent_row %i\n",nb_ent_row);
		if(nb_ent_row>1){
			// printf("r %i ;",r);
			// limit the change to cells that haven't been changed
			int select_free[nb_ent_row];
			int size_select_free=0;
			int nument=0;
			// printf("ent:");
			for ( nument = rpoint[r]-1; nument < rpoint[r+1]-1; nument += 1 ) { 
				if(Free[nument]==1){
					size_select_free++;
					select_free[size_select_free-1]=nument;
					Free[nument]=0;
					// printf("%i ",nument);
				}
			}
			// resample the values
			double resampled_values[size_select_free];
			// printf("initial values:");
			for ( i = 0; i < size_select_free; i += 1 ) { 
				// printf(" %.1f",entries[select_free[i]]);
				resampled_values[i]=entries[select_free[i]];
			}
			resample(resampled_values,&size_select_free);
			// put back resampled values in the corresponding row and col
			for ( i = 0; i < size_select_free; i += 1 ) { 
				// printf("\nval %.1f",resampled_values[i]);
				int init_ent=select_free[i];
				entries[init_ent]=resampled_values[i]; // raw
				// printf("raw ok, new raw %i",coli[init_ent]-1);
				int rawFirstEnt=rpoint[coli[init_ent]-1]-1;
				int rawLastEnt=rpoint[coli[init_ent]]-1;
				int newEnt=rawFirstEnt;
				// printf("fe:%i le:%i ",rawFirstEnt,rawLastEnt);
				while(newEnt<=rawLastEnt && (coli[newEnt]-1)!=r){
					newEnt++;
				}
				// printf("newEnt: %i\n",newEnt);
				if((coli[newEnt]-1)==r){
					entries[newEnt]=resampled_values[i]; // col
					Free[newEnt]=0; 
				}else{
					printf("ERROR: matrix is not symmetric\n");
					error=1;
					return ;
				}
			}
		}
		// printf("\n");
	}
}

//formerly used for testing
/*
int main(){
 	int stime;
 	long ltime;

  	//get the current calendar time
 	ltime = time(NULL);
	stime = (unsigned) ltime/2;
	srand(stime);
	double vect[]={1.,2.,3.};
	int vect_l=3;
	resample(vect,&vect_l);
	int i;
	for(i=0;i<3;i++){
		printf("%f ",vect[i]);
	}
	printf("\n");
	return(1);
}

*/	


