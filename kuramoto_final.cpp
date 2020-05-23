// Kuramoto model (sincronizzazione di n oscillatori) 
// utilizzo gli header tratti dal libro Numerical Recipes (3rd edition)
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

const int n = 4;  				// numero di oscillatori
double omega_array[n];  			// array delle omega 
int time_int = 1200; 				// tempo d'integrazione
double K = 0.2;					// coupling const.
int counts = 0;

// all'interno della struttura derivates è presente il set di n equazioni differenziali
struct derivates{
	derivates(){};
	void operator()(const double& t, const VecDoub &yvector, VecDoub &dydx){
		int i = 0;
		int k = 0;
		counts++;																	
		for(i = 0; i < n; i++){	// due cicli for che realizzano il modello di Kuramoto
			dydx[i] = omega_array[i];
			for(k = 0; k < n; k++){
				dydx[i] += K/n * sin(yvector[k] - yvector[i]);
			}
		}
	}
};
Doub zero_2_PI(Doub x){			// funzione che permette di avere angoli compresi tra 0 e 2pi
	while(x <= 0){
		x += 2.0*M_PI;
	}
	while(x > 2*M_PI){
		x -= 2.0*M_PI;
	}return x;
}

int main(){
	const int N = n;										// inizo il main con una carrellata di declarations
	Doub x1 = 0.0, x2 = 0.0, dx = 0.003;					// estremi d'integrazione, stepsize
	Doub atol = 1e-12, rtol = 1e-12, stepmin = 1e-18;			// abs e rel tolerances, stepsize minima
	Doub x1max = time_int * dx;							// stepsize * tempo d'integrazione
	VecDoub ystart(N);									// vettore di stati iniziali
	int i = 0;
	ofstream omega("omega.csv");							// scrivo le omega su un file esterno
	ofstream phase("phases.csv");							// scrivo le fasi su un file esterno
	std:: default_random_engine gen;							// genero una distribuzione normale di pulsazioni
	std:: normal_distribution<double> norm(2.5,1.0);
	for(i=0; i<n;i++){										// riempio l'array delle pulsazioni
		omega_array[i]=norm(gen);
		omega << omega_array[i] <<" "<<endl;			// stampo su outfile_w l'array delle pulsazioni
	}
	std::default_random_engine gen2;							// genero una distribuzione uniforme di angoli tra 0 e 2pi
	std::uniform_real_distribution<double> unif(0,2*M_PI);
	for(i=0; i<n;i++){										// riempio l'array delle fasi
		ystart[i]=unif(gen2);
		phase<< ystart[i] <<" "<<endl;					// stampo su outfile_p l'array delle fasi
	}
	omega . close();
	phase . close();
	counts = 0;										
	ofstream trajectory("trajectory.csv");						// traiettoria, parametro d'ordine e il suo modulo
	ofstream order_1("order.csv");
	ofstream order_2("Rorder.csv");	
	for( x1 = 0; x1 < x1max; x1 += dx){					// inizia un ciclo for che aggiunge la stepsize ad x1 fino ad arrivare a runtime*stepsize
		x2 = x1 + dx;									// integro sempre tra x1 e x1+dx=x2
		cout << "time: " << x1 << " " << x2 << endl; 	
		derivates d;
		Output out(10);									// Output è l'oggetto che mi fa prendere valori intermedi
		Odeint <StepperDopr5<derivates> > ode(ystart, x1, x2, atol ,rtol ,dx ,stepmin , out, d);
		ode.integrate();									// integrazione del set di equazioni differenziali
		for(i = 0; i < n; i++){
			ystart[i] = zero_2_PI(ystart[i]);					// scrivo le fasi tra 0 e 2pi	
		}
		trajectory << x1;									// stampo le x1 in trajectory (tempo iniziale)
		complex<double> im = -1;
		im = sqrt(im);									// introduco l'unità immaginaria
		for(i = 0; i < n; i++){
			trajectory << "," << ystart[i] ;					// mi stampa ystart su trajectory 
		}
		trajectory << endl;
		complex<double> order = exp(im * ystart[0]);			// calcolo il parametro d'ordine
		for(i = 1; i < n; i++){
			order += exp(im * ystart[i]);
		}
		double realorder = real(order) / n;					// parte reale
		double imagorder = imag(order) / n;					// parte immaginaria
		double order1 = sqrt(realorder * realorder + imagorder * imagorder);  	// modulo
		order_1 << x1 << ","<< realorder << ","<< imagorder <<endl;
		order_2 << x1 <<","<< order1 << endl;							// stampo parte reale e imm. su order.dat, il modulo su Roder.dat
	}
	trajectory.close();
	order_1.close();
	order_2.close();
	return 0;
}





