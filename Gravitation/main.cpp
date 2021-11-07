// Gravitation.cpp : Defines the entry point for the console application.
//
#include "stdio.h"
#include "math.h"

#include "util.h"

const int Nbodies = 10;
const double G = 6.674E-8;
const double Msun = 1.9885E33;
const double AU = 14959787070000;
const double Vfirst = sqrt(G*Msun/AU);
const int Niterations = 60*876000;

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

void initSolarSystem(double** coordinate, double** velocity, double* const mass){
	if(Nbodies != 10){
		printf("Nbodies must be 10 in Solar System\n");
		exit(0);
	}

	//double r[Nbodies][2] = {{0,0},{0.47,0.47},{0.72,0.72},{0.98, 0.98},{1.39,1.39},{4.96,4.97},{9.18,9.18},{19.92,19.92},{30.12,30.12},{30.22,30.22}};
	//double beta[Nbodies][2] = {{0,0},{-2.9,-3.2},{3.3,3.3},{0,0},{-1.4,-1.4},{-1.2,-1.2},{-2.3,-2.3},{-0.7,-0.7},{0.2,0.2},{11.2,11.2}};
	//double lambda[Nbodies][2] = {{0,0},{252.4,255.2},{181.8,183.4},{99.9,100.9},{359.1,359.8},{36.2,36.3},{45.7,45.7},{316.4,316.4},{303.9,303.9},{250.5,250.5}};

	mass[0] = Msun;
	coordinate[0][0] = 0;
	coordinate[0][1] = 0;
	coordinate[0][2] = 0;
	velocity[0][0] = 0;
	velocity[0][1] = 0;
	velocity[0][2] = 0;

	mass[1] = 3.33E26;
	mass[2] = 4.86E27;
	mass[3] = (1.0 + 1.0/81.0)*5.97E27;
	mass[4] = 6.42E26;
	mass[5] = 1.9E30;
	mass[6] = 5.68E29;
	mass[7] = 8.68E28;
	mass[8] = 1.02E29;
	mass[9] = 1.3E25 + 1.52E24;

	/*for(int i = 1; i < Nbodies; ++i){
		coordinate[i][0] = AU*r[i][0]*sin(pi*(90 - beta[i][0])/180)*cos(pi*lambda[i][0]/180);
		double newx = AU*r[i][1]*sin(pi*(90 - beta[i][1])/180)*cos(pi*lambda[i][1]/180);
		coordinate[i][1] = AU*r[i][0]*sin(pi*(90 - beta[i][0])/180)*sin(pi*lambda[i][0]/180);
		double newy = AU*r[i][1]*sin(pi*(90 - beta[i][1])/180)*sin(pi*lambda[i][1]/180);
		coordinate[i][2] = AU*r[i][0]*cos(pi*(90 - beta[i][0])/180);
		double newz = AU*r[i][1]*cos(pi*(90 - beta[i][1])/180);

		velocity[i][0] = (newx - coordinate[i][0])/(24*3600);
		velocity[i][1] = (newy - coordinate[i][1])/(24*3600);
		velocity[i][2] = (newz - coordinate[i][2])/(24*3600);
	}*/
	//at 01 01 200 12-00
	coordinate[1][0] = -1.946172635585372E12;
	coordinate[1][1] = -6.691327526352400E12;
	coordinate[1][2] = -3.679854343749542E11;
	velocity[1][0] = 3.699499185727919E6;
	velocity[1][1] = -1.116441592561850E6;
	velocity[1][2] = -4.307628118658220E4;

	coordinate[2][0] = -1.074564940521906E13;
	coordinate[2][1] = -4.885014975872536E11;
	coordinate[2][2] = 6.135634299718402E11;
	velocity[2][0] = 1.381906029263447E5;
	velocity[2][1] = -3.514029517644670E6;
	velocity[2][2] = -5.600423382820807E4;

	coordinate[3][0] = -2.649903367743050E12;
	coordinate[3][1] = 1.446972967925493E13;
	coordinate[3][2] = -6.111494259536266E7;
	velocity[3][0] = -2.979426007043741E6;
	velocity[3][1] = -5.469294939770602E5;
	velocity[3][2] = 1.817836785027449E1;

	coordinate[4][0] = 2.080481406481886E13;
	coordinate[4][1] = -2.007052475139215E11;
	coordinate[4][2] = -5.156288876678155E11;
	velocity[4][0] = 1.162672383641686E5;
	velocity[4][1] = 2.629606454658444E6;
	velocity[4][2] = 5.222970037729198E4;

	coordinate[5][0] = 5.985676246570644E13;
	coordinate[5][1] = 4.396046799481729E13;
	coordinate[5][2] = -1.522686167298746E12;
	velocity[5][0] = -7.909860292172008E5;
	velocity[5][1] = 1.115621752636729E6;
	velocity[5][2] = 1.308656815823666E4;

	coordinate[6][0] = 9.583853590115459E13;
	coordinate[6][1] = 9.828562833594368E13;
	coordinate[6][2] = -5.521297984971136E12;
	velocity[6][0] = -7.431212960083195E5;
	velocity[6][1] = 6.736755749509767E5;
	velocity[6][2] = 1.777382899472442E4;

	coordinate[7][0] = 2.158975038179710E14;
	coordinate[7][1] = -2.054625247626725E14;
	coordinate[7][2] = -3.562546377838492E12;
	velocity[7][0] = 4.637273797531964E5;
	velocity[7][1] = 4.627599147267434E5;
	velocity[7][2] = -4.291789035828519E3;

	coordinate[8][0] = 2.515046471487719E14;
	coordinate[8][1] = -3.738714539167761E14;
	coordinate[8][2] = 1.903221741384149E12;
	velocity[8][0] = 4.465275177950522E5;
	velocity[8][1] = 3.075980180890997E5;
	velocity[8][2] = -1.662486428694513E4;

	coordinate[9][0] = -1.477331056749140E14;
	coordinate[9][1] = -4.182575030598021E14;
	coordinate[9][2] = 8.752155613302606E13;
	velocity[9][0] = 5.259811447117987E5;
	velocity[9][1] = -2.656425857509986E5;
	velocity[9][2] = -1.250553068550307E4;
}

int main()
{
	double dt = 0.0001*AU/Vfirst;
	dt = 60;
	double time = 0;

	double** coordinate;
	double** velocity;
	double** prevCoordinate;
	double** prevVelocity;
	double** tempCoordinate;
	double** tempVelocity;
	double** rightPart;
	double mass[Nbodies];
	srand(12);

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

	initSolarSystem(coordinate, velocity, mass);

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

		coordinate[i][0] -= xc;
		coordinate[i][1] -= yc;
		coordinate[i][2] -= zc;
	}
	collectStat(E, px, py, pz, xc, yc, zc, Lx, Ly, Lz, coordinate, velocity, mass);

	FILE* general = fopen("general.dat","w");
	fprintf(general, "%g %g %g %g %g %g %g %g %g %g\n", E, px, py, pz, xc, yc, zc, Lx, Ly, Lz);
	fclose(general);

	FILE* outFile = fopen("output.dat", "w");
	output(outFile, time, coordinate, velocity);
	fclose(outFile);

	FILE* mercFile = fopen("mercury.dat","w");
	fclose(mercFile);

	double mercPrevR = 0;
	double mercPrevPrevR = -1;

	double mercX0 = coordinate[1][0] - coordinate[0][0];
	double mercY0 = coordinate[1][1] - coordinate[0][1];
	double mercZ0 = coordinate[1][2] - coordinate[0][2];

	for(int k = 0; k < Niterations; ++k){
		printf("iteration numer %d\n", k);
		time += dt;

		//eulerIteration(dt, coordinate, velocity, mass, tempCoordinate, tempVelocity);
		double mercX = coordinate[1][0] - coordinate[0][0];
		double mercY = coordinate[1][1] - coordinate[0][1];
		double mercZ = coordinate[1][2] - coordinate[0][2];
		explicitRungeKuttaIteration(dt, coordinate, velocity, mass, tempCoordinate, tempVelocity, rightPart);

		double mercR = sqrt(sqr(coordinate[0][0] - coordinate[1][0]) + sqr(coordinate[0][1]- coordinate[1][1]) + sqr(coordinate[0][2] - coordinate[1][2]));
		if((mercPrevR < mercPrevPrevR) && (mercPrevR < mercR)){
			double theta = evaluateAngle(mercX0, mercY0, mercZ0, coordinate[1][0] - coordinate[0][0], coordinate[1][1] - coordinate[0][1], coordinate[1][2] - coordinate[0][2]);
			mercFile = fopen("mercury.dat","a");
			fprintf(mercFile,"%15.10g %15.10g %15.10g %15.10g %15.10g\n", time-dt, mercX, mercY, mercZ, theta);
			fclose(mercFile);
		}

		mercPrevPrevR = mercPrevR;
		mercPrevR = mercR;

		if(k%(12*60) == 0){
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

