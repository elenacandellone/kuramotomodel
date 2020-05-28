// Kuramoto model (Synchronization of N oscillators) 
// headers from Numerical Recipes (3rd edition)
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <complex>
#include "nr3.h"			
#include "stepper.h"
#include "stepperdopr5.h"   // fifth-order Dormand-Prince method
#include "odeint.h"


using namespace std;

const int n = 4;		// number of oscillators
double omega_array[n];  	// omega's array
int time_int = 1200; 		// integration time
double K = 0.2;			// coupling parameter
int counts = 0;

// set of n differential equations
struct derivates{
	derivates(){};
	void operator()(const double& t, const VecDoub &yvector, VecDoub &dydx){
		int i = 0;
		int k = 0;
		counts++;																	
		for(i = 0; i < n; i++){	
			dydx[i] = omega_array[i];
			for(k = 0; k < n; k++){
				dydx[i] += K/n * sin(yvector[k] - yvector[i]);
			}
		}		// 2 for loops in order to realize Kuramoto model
	}
};
Doub zero_2_PI(Doub x){		// phases between 0 and 2pi
	while(x <= 0){
		x += 2.0*M_PI;
	}
	while(x > 2*M_PI){
		x -= 2.0*M_PI;
	}return x;
}

int main(){
	const int N = n;					
	Doub x1 = 0.0, x2 = 0.0, dx = 0.003;			// integration boundaries, stepsize
	Doub atol = 1e-12, rtol = 1e-12, stepmin = 1e-18;	// abs e rel tolerances, min stepsize 
	Doub x1max = time_int * dx;				// stepsize * integration time
	VecDoub ystart(N);					// initial states array
	int i = 0;
	ofstream omega("omega.csv");				// write omegas on omega.csv
	ofstream phase("phases.csv");				// write phases on phases.csv
	std:: default_random_engine gen;			// random normal distribution of frequencies
	std:: normal_distribution<double> norm(2.5,1.0);
	for(i=0; i<n;i++){										
		omega_array[i]=norm(gen);
		omega << omega_array[i] <<" "<<endl;			
	}
	std::default_random_engine gen2;			// uniform distribution of phases between 0 and 2pi
	std::uniform_real_distribution<double> unif(0,2*M_PI);
	for(i=0; i<n;i++){									
		ystart[i]=unif(gen2);
		phase<< ystart[i] <<" "<<endl;					
	}
	omega . close();
	phase . close();
	counts = 0;										
	ofstream trajectory("trajectory.csv");			
	ofstream order_1("order.csv");
	ofstream order_2("Rorder.csv");	
	for( x1 = 0; x1 < x1max; x1 += dx){			
		x2 = x1 + dx;				// int. boundaries are x1 and x1+dx=x2
		cout << "time: " << x1 << " " << x2 << endl; 	
		derivates d;
		Output out(10);				// Output is an object that takes intermediate values
		Odeint <StepperDopr5<derivates> > ode(ystart, x1, x2, atol ,rtol ,dx ,stepmin , out, d);
		ode.integrate();			// integration of differential eqs.
		for(i = 0; i < n; i++){
			ystart[i] = zero_2_PI(ystart[i]);	// phases between 0 and 2pi
		}
		trajectory << x1;			// lower bound printed on trajectory.csv
		complex<double> im = -1;
		im = sqrt(im);				// imaginary unit
		for(i = 0; i < n; i++){
			trajectory << "," << ystart[i] ;	// new phases printed on trajectory 
		}
		trajectory << endl;
		complex<double> order = exp(im * ystart[0]);	// order parameter
		for(i = 1; i < n; i++){
			order += exp(im * ystart[i]);
		}
		double realorder = real(order) / n;					// real part
		double imagorder = imag(order) / n;					// imaginary part
		double order1 = sqrt(realorder * realorder + imagorder * imagorder);  	// module
		order_1 << x1 << ","<< realorder << ","<< imagorder <<endl;
		order_2 << x1 <<","<< order1 << endl;					
	}
	trajectory.close();
	order_1.close();
	order_2.close();
	return 0;
}





