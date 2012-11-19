/*
 * ===========================================================================
 *
 *       Filename:  filter_spam.c
 *
 *    Description:  set entries of spam over x to be 0
 *    		should be compiled by:
 *    		R CMD SHLIB filter_spam.c
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


