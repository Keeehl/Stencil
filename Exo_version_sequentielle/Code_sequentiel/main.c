#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define SIMULATION_TIME 500 //seconds, can be tweaked
#define TIME_STEPS 300 //Can be tweaked
#define DT (SIMULATION_TIME+.0)/TIME_STEPS

#define LENGTH 20.0 //length of the area of simulation in meters, can be tweaked
#define SQUARES_PER_ROW 25 // Can be tweaked
#define DX (LENGTH+.0)/SQUARES_PER_ROW

#define FACTOR 0.5 // Can be tweaked (serves to tweak diffusion coefficient) but has to be in [0;1] range for stability
#define DIFFUSION_COEFFICIENT FACTOR*(DX*DX)/(2*DT)
#define T0 297.0 //Can be tweaked
#define TEMPERATURE_AMPLITUDE 30 //Can be tweaked
#define PI 3.14


double Update_temperature(double current_square, double i_m1, double j_m1, double i_p1, double j_p1, double D, double dx, double dt);

int main() {
	//For displaying results
	FILE* script_file=fopen("plotting_script.txt","w+"); 
	fprintf(script_file,"unset key\nset pm3d at s\nsplot 'temperatures_file.txt' with pm3d");
	fclose(script_file);
	
	double**temperature_map=malloc(sizeof(double*)*SQUARES_PER_ROW); //temperature field at time step t
	double**new_temperature_map=malloc(sizeof(double*)*SQUARES_PER_ROW); //temperature field at time step t+1
	
	for (int i=0;i<SQUARES_PER_ROW;i++) { //Allocation the squares
		new_temperature_map[i]=malloc(sizeof(double)*SQUARES_PER_ROW);
		temperature_map[i]=malloc(sizeof(double)*SQUARES_PER_ROW);
	}
	
	//Boundary conditions and Initial conditions
	for (int i=0;i<SQUARES_PER_ROW;i++) { 
		for (int j=0;j<SQUARES_PER_ROW;j++) {
			if ((j==0)||(i==0)||(i==SQUARES_PER_ROW-1)||(j==SQUARES_PER_ROW-1)) {
				temperature_map[i][j]=T0; // Boundary conditions
				//You may change what goes into temperature_map[i][j], it will give different results
			}
			else {
				temperature_map[i][j]=T0+TEMPERATURE_AMPLITUDE*sin(2*PI*2*i*DX/LENGTH)*sin(2*PI*2*j*DX/LENGTH);
				//Initial conditions. You may change this line as well to observe how it affects the results
			}
			new_temperature_map[i][j]=temperature_map[i][j];
		}
	}
	
	//Simulating the evolution of the system
	for (int i=0;i<TIME_STEPS;i++) { 
		FILE* temperature_display=fopen("temperatures_file.txt","w+");
		for (int i=1;i<SQUARES_PER_ROW-1;i++) {
			for (int j=1;j<SQUARES_PER_ROW-1;j++) {
				new_temperature_map[i][j]=Update_temperature(temperature_map[i][j],temperature_map[i-1][j],temperature_map[i][j-1],temperature_map[i+1][j],temperature_map[i][j+1],DIFFUSION_COEFFICIENT,DX,DT);
				//Temperature is updated according to nearby temperatures values
			}
		}
		for (int i=0;i<SQUARES_PER_ROW;i++) { //Remplacing "old" temperatures values by new ones
			for (int j=0;j<SQUARES_PER_ROW;j++) {
				fprintf(temperature_display,"%f %f %f\n",i*DX,j*DX,temperature_map[i][j]);
				temperature_map[i][j]=new_temperature_map[i][j];
			}
			fprintf(temperature_display,"\n");
		}
		fclose(temperature_display);
	} 
	//End of simulation
	system("gnuplot -p 'plotting_script.txt'"); //Displaying results
	for (int i=0;i<SQUARES_PER_ROW;i++) {
		free(temperature_map[i]);
		free(new_temperature_map[i]);
	}
	free(temperature_map);
	free(new_temperature_map);
	return 0;
}


//This function takes temperatures values near a given square as argument and returns the temperature in that square for next time step
double Update_temperature(double current_square, double i_m1, double j_m1, double i_p1, double j_p1, double D, double dx, double dt) {
	//The arguments are, in order: temperature of square (i,j), of square (i-1,j),(i,j-1),(i+1,j),(i,j+1)
	// as well as the hyperparameters D,dx and dt (diffusion coefficient, space step and time step)
	double time_derivative=D*(i_p1+j_p1+i_m1+j_m1-4*current_square)/(dx*dx);
	double temperature=current_square+dt*time_derivative;
	return temperature;
}
