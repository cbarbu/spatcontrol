/*
 * =====================================================================================
 *
 *       Filename:  useC2.c
 *
 *    Description:  complements to spam package
 *    		should be compiled by:
 *    		R CMD SHLIB useC2.c
 *    		Then can execute useC2.r
 *
 *        Version:  1.0
 *        Created:  06/02/2011 11:13:13 AM
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
#include <time.h>

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
int main(){
	int stime;
	long ltime;

	/* get the current calendar time */
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
	
	
	

