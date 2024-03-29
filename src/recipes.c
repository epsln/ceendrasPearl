#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "include/recipes.h"

void maskitRecipe(double complex ta, double complex* gens){
	//See pp. 259
		
	printf("mu:  %lf + %lf\n", creal(ta), cimag(ta));

	gens[(0 * 2 + 0) * 2 + 0] = ta;
	gens[(0 * 2 + 0) * 2 + 1] = -I; 
	gens[(0 * 2 + 1) * 2 + 0] = -I;
	gens[(0 * 2 + 1) * 2 + 1] = 0; 

	gens[(1 * 2 + 0) * 2 + 0] = 1;
	gens[(1 * 2 + 0) * 2 + 1] = 2; 
	gens[(1 * 2 + 1) * 2 + 0] = 0;
	gens[(1 * 2 + 1) * 2 + 1] = 1; 

	gens[(2 * 2 + 0) * 2 + 0] =  gens[(0 * 2 + 1) * 2 + 1];
	gens[(2 * 2 + 0) * 2 + 1] = -gens[(0 * 2 + 0) * 2 + 1];
	gens[(2 * 2 + 1) * 2 + 0] = -gens[(0 * 2 + 1) * 2 + 0];
	gens[(2 * 2 + 1) * 2 + 1] =  gens[(0 * 2 + 0) * 2 + 0];
                                                             
	gens[(3 * 2 + 0) * 2 + 0] =  gens[(1 * 2 + 1) * 2 + 1];
	gens[(3 * 2 + 0) * 2 + 1] = -gens[(1 * 2 + 0) * 2 + 1];
	gens[(3 * 2 + 1) * 2 + 0] = -gens[(1 * 2 + 1) * 2 + 0];
	gens[(3 * 2 + 1) * 2 + 1] =  gens[(1 * 2 + 0) * 2 + 0];
}



void grandmaRecipe(double complex ta, double complex tb, double complex* gens){
	if (ta == 0 && tb == 0){
		printf("Error ! ta and tb cannot be both == 0 !\n Exiting...\n");
		exit(-3);
	}
	printf("ta:  %lf + %lf\n", creal(ta), cimag(ta));
	printf("tb:  %lf + %lf\n", creal(tb), cimag(tb));

	double complex a = 1;
	double complex b = (-ta * tb);
	double complex c = ta * ta + tb * tb;
       	double complex delta = b*b - 4 * a * c; 
	double complex tab = (- b - csqrt(delta))/(2 * a); 
	double complex z0 = ((tab - 2) * tb)/(tb * tab - 2 * ta + 2 * I * tab);

	gens[(0 * 2 + 0) * 2 + 0] = ta/2;
	gens[(0 * 2 + 1) * 2 + 0] =  (ta*tab - 2 * tb + 4 * I)/(z0*(2 * tab + 4)); 
	gens[(0 * 2 + 0) * 2 + 1] = ((ta * tab - 2 * tb - 4 * I)*z0)/(2* tab - 4);
	gens[(0 * 2 + 1) * 2 + 1] = ta/2; 
                                
	gens[(1 * 2 + 0) * 2 + 0] = (tb - 2 * I)/2;
	gens[(1 * 2 + 1) * 2 + 0] = tb/2; 
	gens[(1 * 2 + 0) * 2 + 1] = tb/2;
	gens[(1 * 2 + 1) * 2 + 1] = (tb + 2 * I)/2; 
                                
	gens[(2 * 2 + 0) * 2 + 0] =  gens[(0 * 2 + 1) * 2 + 1 ];
	gens[(2 * 2 + 1) * 2 + 0] = -gens[(0 * 2 + 1) * 2 + 0 ];
	gens[(2 * 2 + 0) * 2 + 1] = -gens[(0 * 2 + 0) * 2 + 1 ];
	gens[(2 * 2 + 1) * 2 + 1] =  gens[(0 * 2 + 0) * 2 + 0 ];
                                                             
	gens[(3 * 2 + 0) * 2 + 0] =  gens[(1 * 2 + 1) * 2 + 1 ];
	gens[(3 * 2 + 1) * 2 + 0] = -gens[(1 * 2 + 1) * 2 + 0 ];
	gens[(3 * 2 + 0) * 2 + 1] = -gens[(1 * 2 + 0) * 2 + 1 ];
	gens[(3 * 2 + 1) * 2 + 1] =  gens[(1 * 2 + 0) * 2 + 0 ];
}

void grandmaSpecialRecipe(double complex ta, double complex tb, double complex tab, double complex* gens){
	//Grandma's Special four-alarm two generator group recipe 
	//See pp. 261
//	printf("ta:  %lf + %lf\n", creal(ta), cimag(ta));
//	printf("tb:  %lf + %lf\n", creal(tb), cimag(tb));
//	printf("tab: %lf + %lf\n", creal(tab), cimag(tab));
	double complex tc = ta * ta + tb * tb + tab * tab - ta * tb * tab - 2;
	double complex Q  = csqrt(2 - tc);	
	double complex R = 0;
	
	if (cabs(tab) == 2){ 
		printf("Error ! taB cannot be == +/- 2 !\n Exiting...\n");
		exit(-3);
	}
	if (cabs(tc + I * Q * csqrt(tc + 2)) >= 2)
		R = csqrt(tc + 2);
	else
		R = -csqrt(tc + 2);
	
	double complex z0 = ((tab - 2) * (tb + R))/(tb * tab - 2 * ta + I * Q * tab);
	
	gens[(0 * 2 + 0) * 2 + 0] = ta/2;
	gens[(0 * 2 + 1) * 2 + 0] =  (ta * tab - 2 * tb + 2 * I * Q)/(z0 * (2 * tab + 4)); 
	gens[(0 * 2 + 0) * 2 + 1] = ((ta * tab - 2 * tb - 2 * I * Q) * z0)/(2 * tab - 4);
	gens[(0 * 2 + 1) * 2 + 1] = ta/2; 
                                
	gens[(1 * 2 + 0) * 2 + 0] = (tb - I * Q)/2;
	gens[(1 * 2 + 1) * 2 + 0] = (tb * tab - 2 * ta - I * Q * tab)/(z0*(2 * tab + 4));
	gens[(1 * 2 + 0) * 2 + 1] = ((tb * tab - 2 * ta + I * Q * tab)*z0)/(2 * tab - 4);
	gens[(1 * 2 + 1) * 2 + 1] = (tb + I * Q)/2; 

	gens[(2 * 2 + 0) * 2 + 0] =  gens[(0 * 2 + 1) * 2 + 1];
	gens[(2 * 2 + 1) * 2 + 0] = -gens[(0 * 2 + 1) * 2 + 0];
	gens[(2 * 2 + 0) * 2 + 1] = -gens[(0 * 2 + 0) * 2 + 1];
	gens[(2 * 2 + 1) * 2 + 1] =  gens[(0 * 2 + 0) * 2 + 0];
                                          
	gens[(3 * 2 + 0) * 2 + 0] =  gens[(1 * 2 + 1) * 2 + 1];
	gens[(3 * 2 + 1) * 2 + 0] = -gens[(1 * 2 + 1) * 2 + 0];
	gens[(3 * 2 + 0) * 2 + 1] = -gens[(1 * 2 + 0) * 2 + 1];
	gens[(3 * 2 + 1) * 2 + 1] =  gens[(1 * 2 + 0) * 2 + 0];

}

void jorgensen(double complex ta, double complex tb, double complex* gens){
	printf("ta:  %lf + %lf\n", creal(ta), cimag(ta));
	printf("tb:  %lf + %lf\n", creal(tb), cimag(tb));
	double complex z   = 0.5 * csqrt(ta * ta * tb * tb - 4 * ta * ta - 4 * tb * tb);
	double complex tab = 0.5 * (ta * tb) - z;

	gens[(0 * 2 + 0) * 2 + 0] = ta - tb / tab;
	gens[(0 * 2 + 1) * 2 + 0] = ta / (tab * tab);
	gens[(0 * 2 + 0) * 2 + 1] = ta;
	gens[(0 * 2 + 1) * 2 + 1] = tb / (tab);
                                
	gens[(1 * 2 + 0) * 2 + 0] = tb - ta / tab;
	gens[(1 * 2 + 1) * 2 + 0] = -tb / (tab * tab); 
	gens[(1 * 2 + 0) * 2 + 1] = -tb; 
	gens[(1 * 2 + 1) * 2 + 1] = ta / tab; 

	gens[(2 * 2 + 0) * 2 + 0] =  gens[(0 * 2 + 1) * 2 + 1];
	gens[(2 * 2 + 1) * 2 + 0] = -gens[(0 * 2 + 1) * 2 + 0];
	gens[(2 * 2 + 0) * 2 + 1] = -gens[(0 * 2 + 0) * 2 + 1];
	gens[(2 * 2 + 1) * 2 + 1] =  gens[(0 * 2 + 0) * 2 + 0];
                                          
	gens[(3 * 2 + 0) * 2 + 0] =  gens[(1 * 2 + 1) * 2 + 1];
	gens[(3 * 2 + 1) * 2 + 0] = -gens[(1 * 2 + 1) * 2 + 0];
	gens[(3 * 2 + 0) * 2 + 1] = -gens[(1 * 2 + 0) * 2 + 1];
	gens[(3 * 2 + 1) * 2 + 1] =  gens[(1 * 2 + 0) * 2 + 0];



}
