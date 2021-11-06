// Gravitation.cpp : Defines the entry point for the console application.
//
#include "stdio.h"
#include "math.h"

#include "util.h"

const int Nbodies = 3;
const double G = 6.674E-8;
const double Msun = 1.9885E33;
const double AU = 15000000000000;
const double Vfirst = sqrt(G*Msun/AU);
const int Niterations = 10000000;

void output(FILE* outFile, double& time, double** coordinate, double** velocity){
	fprintf(outFile, "%15.10g", time);
	for(int i = 0; i < Nbodies; ++i){
		for(int j = 0; j < 3; ++j){
			fprintf(outFile, " %15.10g", coordinate[i][j]);
		}
		for(int j = 0; j < 3; ++j){
			fprintf(outFile, " %15.10g", velocity[i][j]);
		}
	}
	fprintf(outFile, "\n");
}

void collectStat(double& E, double& px, double& py, double& pz, double& xc, double& yc, double& zc, double& Lx, double& Ly, double& Lz, double** coordinate, double** velocity, const double* mass){
	double distance[Nbodies][Nbodies];
	for(int i = 0; i < Nbodies; ++i){
		distance[i][i] = 0;
		for(int j = i+1; j < Nbodies; ++j){
			distance[i][j] = sqrt(sqr(coordinate[i][0] - coordinate[j][0]) + sqr(coordinate[i][1] - coordinate[j][1]) + sqr(coordinate[i][2] - coordinate[j][2]));
			distance[j][i] = distance[i][j];
		}
	}

	E = 0;
	double Ek = 0;
	double Ep = 0;
	px = 0;
	py = 0;
	pz = 0;
	xc = 0;
	yc = 0;
	zc = 0;
	Lx = 0;
	Ly = 0;
	Lz = 0;

	double M = 0;

	for(int i = 0; i < Nbodies; ++i){
		M += mass[i];

		Ek += 0.5*mass[i]*(velocity[i][0]*velocity[i][0] + velocity[i][1]*velocity[i][1] + velocity[i][2]*velocity[i][2]);
		for(int j = i+1; j < Nbodies; ++j){
			Ep -= G*mass[i]*mass[j]/distance[i][j];
		}

		px += mass[i]*velocity[i][0];
		py += mass[i]*velocity[i][1];
		pz += mass[i]*velocity[i][2];

		xc += mass[i]*coordinate[i][0];
		yc += mass[i]*coordinate[i][1];
		zc += mass[i]*coordinate[i][2];

		Lx += mass[i]*(coordinate[i][1]*velocity[i][2] - coordinate[i][2]*velocity[i][1]);
		Ly += mass[i]*(coordinate[i][2]*velocity[i][0] - coordinate[i][0]*velocity[i][2]);
		Lz += mass[i]*(coordinate[i][0]*velocity[i][1] - coordinate[i][1]*velocity[i][0]);
	}

	E = Ek + Ep;
	xc = xc/M;
	yc = yc/M;
	zc = zc/M;
}

void eulerIteration(double& dt, double** coordinate, double** velocity, const double* mass, double** tempCoordinate, double** tempVelocity){
	double distance[Nbodies][Nbodies];
	for(int i = 0; i < Nbodies; ++i){
		distance[i][i] = 0;
		for(int j = i+1; j < Nbodies; ++j){
			distance[i][j] = sqrt(sqr(coordinate[i][0] - coordinate[j][0]) + sqr(coordinate[i][1] - coordinate[j][1]) + sqr(coordinate[i][2] - coordinate[j][2]));
			distance[j][i] = distance[i][j];
		}
	}

	for(int i = 0; i < Nbodies; ++i){
		for(int j = 0; j < 3; ++j){
			tempCoordinate[i][j] = coordinate[i][j] + velocity[i][j]*dt;
			tempVelocity[i][j] = velocity[i][j];
		}

		double acceleration[3] = {0,0,0};

		for(int k = 0; k < Nbodies; ++k){
			if( i != k){
				for(int j = 0; j < 3; ++j){
					acceleration[j] += G*mass[k]*(coordinate[k][j] - coordinate[i][j])/cube(distance[i][k]);
				}
			}
		}

		for(int j = 0; j < 3; ++j){
			tempVelocity[i][j] += dt*acceleration[j];
		}
	}

	for(int i = 0; i < Nbodies; ++i){
		for(int j = 0; j < 3; ++j){
			coordinate[i][j] = tempCoordinate[i][j];
			velocity[i][j] = tempVelocity[i][j];
		}
	}
}

void evaluateRightPart(double** f, double** coordinate, double** velocity, const double* mass){
	double distance[Nbodies][Nbodies];
	for(int i = 0; i < Nbodies; ++i){
		distance[i][i] = 0;
		for(int j = i+1; j < Nbodies; ++j){
			distance[i][j] = sqrt(sqr(coordinate[i][0] - coordinate[j][0]) + sqr(coordinate[i][1] - coordinate[j][1]) + sqr(coordinate[i][2] - coordinate[j][2]));
			distance[j][i] = distance[i][j];
		}
	}

	for(int i = 0; i < Nbodies; ++i){
		f[i][0] = velocity[i][0];
		f[i][1] = velocity[i][1];
		f[i][2] = velocity[i][2];

		f[i][3] = 0;
		f[i][4] = 0;
		f[i][5] = 0;

		for(int k = 0; k < Nbodies; ++k){
			if( i != k){
				for(int j = 0; j < 3; ++j){
					f[i][3+j] += G*mass[k]*(coordinate[k][j] - coordinate[i][j])/cube(distance[i][k]);
				}
			}
		}
	}
}

void explicitRungeKuttaIteration(double& dt, double** coordinate, double** velocity, const double* mass, double** tempCoordinate, double** tempVelocity, double** rightPart){
	double k1[Nbodies][6];
	double k2[Nbodies][6];
	double k3[Nbodies][6];
	double k4[Nbodies][6];

	evaluateRightPart(rightPart, coordinate, velocity, mass);

	for(int i = 0; i < Nbodies; ++i){
		for(int j = 0; j < 6; ++j){
			k1[i][j] = rightPart[i][j];
		}

		tempCoordinate[i][0] = coordinate[i][0] + 0.5*dt*k1[i][0];
		tempCoordinate[i][1] = coordinate[i][1] + 0.5*dt*k1[i][1];
		tempCoordinate[i][2] = coordinate[i][2] + 0.5*dt*k1[i][2];

		tempVelocity[i][0] = velocity[i][0] + 0.5*dt*k1[i][3];
		tempVelocity[i][1] = velocity[i][1] + 0.5*dt*k1[i][4];
		tempVelocity[i][2] = velocity[i][2] + 0.5*dt*k1[i][5];
	}

	evaluateRightPart(rightPart, tempCoordinate, tempVelocity, mass);

	for(int i = 0; i < Nbodies; ++i){
		for(int j = 0; j < 6; ++j){
			k2[i][j] = rightPart[i][j];
		}

		tempCoordinate[i][0] = coordinate[i][0] + 0.5*dt*k2[i][0];
		tempCoordinate[i][1] = coordinate[i][1] + 0.5*dt*k2[i][1];
		tempCoordinate[i][2] = coordinate[i][2] + 0.5*dt*k2[i][2];

		tempVelocity[i][0] = velocity[i][0] + 0.5*dt*k2[i][3];
		tempVelocity[i][1] = velocity[i][1] + 0.5*dt*k2[i][4];
		tempVelocity[i][2] = velocity[i][2] + 0.5*dt*k2[i][5];
	}

	evaluateRightPart(rightPart, tempCoordinate, tempVelocity, mass);

	for(int i = 0; i < Nbodies; ++i){
		for(int j = 0; j < 6; ++j){
			k3[i][j] = rightPart[i][j];
		}

		tempCoordinate[i][0] = coordinate[i][0] + dt*k3[i][0];
		tempCoordinate[i][1] = coordinate[i][1] + dt*k3[i][1];
		tempCoordinate[i][2] = coordinate[i][2] + dt*k3[i][2];

		tempVelocity[i][0] = velocity[i][0] + dt*k3[i][3];
		tempVelocity[i][1] = velocity[i][1] + dt*k3[i][4];
		tempVelocity[i][2] = velocity[i][2] + dt*k3[i][5];
	}

	evaluateRightPart(rightPart, tempCoordinate, tempVelocity, mass);

	for(int i = 0; i < Nbodies; ++i){
		for(int j = 0; j < 6; ++j){
			k4[i][j] = rightPart[i][j];
		}

		coordinate[i][0] = coordinate[i][0] + dt*(k1[i][0] + 2*k2[i][0] + 2*k3[i][0] + k4[i][0])/6.0;
		coordinate[i][1] = coordinate[i][1] + dt*(k1[i][1] + 2*k2[i][1] + 2*k3[i][1] + k4[i][1])/6.0;
		coordinate[i][2] = coordinate[i][2] + dt*(k1[i][2] + 2*k2[i][2] + 2*k3[i][2] + k4[i][2])/6.0;

		velocity[i][0] = velocity[i][0] + dt*(k1[i][3] + 2*k2[i][3] + 2*k3[i][3] + k4[i][3])/6.0;
		velocity[i][1] = velocity[i][1] + dt*(k1[i][4] + 2*k2[i][4] + 2*k3[i][4] + k4[i][4])/6.0;
		velocity[i][2] = velocity[i][2] + dt*(k1[i][5] + 2*k2[i][5] + 2*k3[i][5] + k4[i][5])/6.0;
	}
}

int main()
{
	double dt = 0.0001*AU/Vfirst;
	double time = 0;

	double** coordinate;
	double** velocity;
	double** prevCoordinate;
	double** prevVelocity;
	double** tempCoordinate;
	double** tempVelocity;
	double** rightPart;
	double mass[Nbodies];
	srand(10);

	coordinate = new double*[Nbodies];
	velocity = new double*[Nbodies];
	prevCoordinate = new double*[Nbodies];
	prevVelocity = new double*[Nbodies];
	tempCoordinate = new double*[Nbodies];
	tempVelocity = new double*[Nbodies];
	rightPart = new double*[Nbodies];

	for(int i = 0; i < Nbodies; ++i){
		coordinate[i] = new double[3];
		velocity[i] = new double[3];
		prevCoordinate[i] = new double[3];
		prevVelocity[i] = new double[3];
		tempCoordinate[i] = new double[3];
		tempVelocity[i] = new double[3];
		rightPart[i] = new double[6];

		coordinate[i][0] = 2*AU*(uniformDistribution() - 0.5);
		coordinate[i][1] = 2*AU*(uniformDistribution() - 0.5);
		coordinate[i][2] = 2*AU*(uniformDistribution() - 0.5);

		velocity[i][0] = 1.0*Vfirst*(uniformDistribution() - 0.5);
		velocity[i][1] = 1.0*Vfirst*(uniformDistribution() - 0.5);
		velocity[i][2] = 1.0*Vfirst*(uniformDistribution() - 0.5);

		for(int j = 0; j < 3; ++j){
			prevCoordinate[i][j] = coordinate[i][j];
			prevVelocity[i][j] = velocity[i][j];
			tempCoordinate[i][j] = coordinate[i][j];
			tempVelocity[i][j] = velocity[i][j];
		}

		mass[i] = Msun*uniformDistribution();
	}

	double E, px, py, pz, xc, yc, zc, Lx, Ly, Lz;
	collectStat(E, px, py, pz, xc, yc, zc, Lx, Ly, Lz, coordinate, velocity, mass);

	double M = 0;
	for(int i = 0; i < Nbodies; ++i){
		M += mass[i];
	}

	double vx = px/M;
	double vy = py/M;
	double vz = pz/M;

	for(int i = 0; i < Nbodies; ++i){
		velocity[i][0] -= vx;
		velocity[i][1] -= vy;
		velocity[i][2] -= vz;
	}
	collectStat(E, px, py, pz, xc, yc, zc, Lx, Ly, Lz, coordinate, velocity, mass);

	FILE* general = fopen("general.dat","w");
	fprintf(general, "%g %g %g %g %g %g %g %g %g %g\n", E, px, py, pz, xc, yc, zc, Lx, Ly, Lz);
	fclose(general);

	FILE* outFile = fopen("output.dat", "w");
	output(outFile, time, coordinate, velocity);
	fclose(outFile);

	for(int k = 0; k < Niterations; ++k){
		printf("iteration numer %d\n", k);
		time += dt;

		//eulerIteration(dt, coordinate, velocity, mass, tempCoordinate, tempVelocity);
		explicitRungeKuttaIteration(dt, coordinate, velocity, mass, tempCoordinate, tempVelocity, rightPart);

		if(k%100 == 0){
			outFile = fopen("output.dat", "a");
			output(outFile, time, coordinate, velocity);
			fclose(outFile);

			general = fopen("general.dat","a");
			collectStat(E, px, py, pz, xc, yc, zc, Lx, Ly, Lz, coordinate, velocity, mass);
			fprintf(general, "%g %g %g %g %g %g %g %g %g %g\n", E, px, py, pz, xc, yc, zc, Lx, Ly, Lz);
			fclose(general);
		}
	}

	for(int i = 0; i < Nbodies; ++i){
		delete[] coordinate[i];
		delete[] velocity[i];
		delete[] prevCoordinate[i];
		delete[] prevVelocity[i];
		delete[] tempCoordinate[i];
		delete[] tempVelocity[i];
		delete[] rightPart[i];
	}
	delete[] coordinate;
	delete[] velocity;
	delete[] prevCoordinate;
	delete[] prevVelocity;
	delete[] tempCoordinate;
	delete[] tempVelocity;
	delete[] rightPart;

	return 0;
}

