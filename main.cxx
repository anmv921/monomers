#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h> 
#include <algorithm>
#include <map>
#include <stdexcept>

using namespace std;

#define NDIM  2
#define ALGORITHM "VERLET"

#include "helperMacros.h"
#include "prototypes.h"

#define DO_PART for (n = 0; n < numParticles; n ++)

typedef struct {
    VecR r, rTest, fa, fb, r0, vel, rPast;
    int nPart;
    real radius;
    real pressure;
    real f2, ddu; // duu is d^2[U]/dr^2
	real mass;
} Particle;

void pairForce(Particle& p1, Particle& p2, int in_stage);
void SRK_Step(int stage);
void InitVelocities();
void VerletStep();
void InitMasses();

// **************************
// Initalization of variables

vector<Particle> particles;
real dt, D, timeNow, uSum, T, r0, V0, lambda, rmax, L, msd, rcut, P, sigmaxy, ekSum, f2Med, dduMed, eTot;
int stepLimit, stepSample, stepCount, moreCycles;
int numParticles; // Must have an even square - ex. 25, 36, and so on. Best use powers of 2
unsigned seed;
ofstream ovito, propsDat;
VecR W, reg;

mt19937 generator;
normal_distribution<real> distribution(0.0, 1.0);

uniform_real_distribution<double> uniformDistribution(0.0, 1.0);


// ****************************
int main(int argc, char **argv)
{	
	auto start = chrono::steady_clock::now();
	cout << "Inicio do programa" << endl;
    ReadInput();
    
    cout << "Input lido" << endl;
    SetupJob();
    
    cout << "Setup finalizado. A executar o ciclo principal" << endl;
    PrintSummary();
    
    moreCycles = 1;
    PrintSummary();
    while (moreCycles) {
        SingleStep ();
        if (stepCount >= stepLimit) moreCycles = 0;
    }
    cout << "Fim do ciclo principal" << endl;
    
    // End process
    ovito.close();
    propsDat.close();
    
    cout << "Fim de execucao do programa" << endl;
    PrintElapsedTime(start);
} // End main
// ****************************

// ****************************
// Functions

void SingleStep ()
{
    ++stepCount;
    timeNow = stepCount * dt;
    
    if (ALGORITHM == "VERLET") {
		computeForces(1);
		VerletStep();
	} 
	else {
		computeForces(1);
		SRK_Step(1);
		computeForces(2);
		SRK_Step(2);
	}
	
    if (stepCount % stepSample == 0) {
        PrintSummary();
        writeProps();
    }
}

void SetupJob ()
{	
	cout << "A inicializar o workflow" << endl;
    AllocArrays();
    stepCount = 0;
    InitCoords();
    Init_F();
    VZero(W);
    ovito.open("ovito.xyz");
    timeNow = 0.0;
    uSum = 0.0;
    ekSum = 0.0;
    stepCount = 0;
    propsDat.open("props.txt");
    writePropsHeader();
    generator = mt19937(seed);
    f2Med = 0.0;
    
    if (ALGORITHM == "VERLET") InitVelocities();

	InitMasses();
    
    cout << "Workflow inicializado" << endl;
}

void AllocArrays ()
{	
	particles.resize(numParticles);
	cout << "Arrays allocados" << endl;
}



void InitVelocities() {
	int n;
	DO_PART {
		particles[n].vel.x = uniformDistribution(generator) - 0.5;
		particles[n].vel.y = uniformDistribution(generator) - 0.5;
	}
}

void InitMasses() {
	int n;
	DO_PART {
		particles[n].mass = 1.0;
	}
}


void CBD_Step ()
{
    int n;
	DO_PART {
		GenerateNoise( W, sqrt(2*dt*D) );
		//RStep( particles[n].r, particles[n].f, W);
	}
}

void VerletStep() {
	int n;
	ekSum = 0.0;
	eTot = 0.0;
	// uSum calcula-se na funcao computeForces
	DO_PART {
		real xNext = 2*particles[n].r.x - particles[n].rPast.x + particles[n].fa.x*dt*dt/particles[n].mass;
		real yNext = 2*particles[n].r.y - particles[n].rPast.y + particles[n].fa.y*dt*dt/particles[n].mass;
		
		particles[n].vel.x = (xNext - particles[n].rPast.x) / (2 * dt);
		particles[n].vel.y = (yNext - particles[n].rPast.y) / (2 * dt);

		real v2 = DotProd(particles[n].vel, particles[n].vel);

		ekSum += 0.5 * particles[n].mass * v2;

		particles[n].rPast.x = particles[n].r.x;
		particles[n].rPast.y = particles[n].r.y;

		particles[n].r.x = xNext;
		particles[n].r.y = yNext;
	}
	eTot = ekSum + uSum;
}

void SRK_Step (int stage)
{	
	int n;
	switch(stage) {
		// ++++++++
		case 1: {
			DO_PART {
				GenerateNoise( W, sqrt(2*dt*D) );
				particles[n].rTest.x = particles[n].r.x;
				particles[n].rTest.y = particles[n].r.y;
				//RStep( particles[n].R, particles[n].fa, W );
			}
			break;
		}
		// ++++++++
		case 2: {
			DO_PART {
				GenerateNoise( W, sqrt(2*dt*D) );
				
				VecR vec_f;
				vec_f.x = (particles[n].fa.x + particles[n].fb.x) / 2.0;
				vec_f.y = (particles[n].fa.y + particles[n].fb.y) / 2.0;
				
				particles[n].f2 = DotProd(vec_f, vec_f);
				
				RStep( particles[n].r, vec_f, W );
			}
			break;
		}
		// ++++++++
	}
}

void InitCoords() {
	// TODO Colocar isto a funcionar com qualquer num de particulas
	real delta = L / (sqrt(numParticles));
	int n = 0;
	for (int i=0; i<sqrt(numParticles); i++) {
		for (int j=0; j<sqrt(numParticles); j++) {
			particles[n].radius = r0, particles[n].nPart = n;
			particles[n].r.x = ((real)i+0.5)*delta;
			particles[n].r.y = ((real)j+0.5)*delta;
			particles[n].rPast.x = particles[n].r.x;
			particles[n].rPast.y = particles[n].r.y;
			n++;
		}
	}
	n = 0;
	DO_PART {
		particles[n].r0.x = particles[n].r.x;
		
		particles[n].r0.y = particles[n].r.y;
	}
	cout << "Fim da inicializacao das coordenadas" << endl;
}

void Init_F ()
{
    int n;
	DO_PART {
		VZero ( particles[n].fa );
		VZero ( particles[n].fb );
	};
}

void PrintSummary ()
{	
	ovito << numParticles << endl;
	ovito << "Lattice=\""<<L<<" 0.0 0.0 0.0 "<<L<<" 0.0 0.0 0.0 0.0\" Origin=\""<<0<<" ";
	ovito<<0<<" 0.0\"";
	ovito << " Properties=pos:R:2:Radius:R:1:Type:S:1:Identifier:I:1:Pressure:R:1";
	ovito << " t=" << timeNow << endl;
	
    int counter=0;
    for( counter=0; counter < numParticles; counter++ )
    {	
		Particle& part = particles[counter];
		
		
		real x_img = part.r.x;
		real y_img = part.r.y;

		// TODO escrever varias copias das partículas a volta das particulas
		// somando +-L, cores diferentes

		// Wrap periodic images
		//if (x_img > L ) x_img -= L;
		//if (x_img < 0 ) x_img += L;
		//if (y_img > L ) y_img -= L;
		//if (y_img < 0 ) y_img += L;

		ovito << x_img << "\t" <<  y_img << "\t" << part.radius << \
		 "\t" << counter << "\t" << part.nPart << "\t" << part.pressure  << "\t" << endl;
    }
}

void ReadInput() {
	cout << "A ler o input" << endl;
	ifstream infile("params.in");
	
	string strKey, strValue;
	
	map<string, string> config;
	
	while (infile >> strKey >> strValue) {
		config[strKey] = strValue;
	}
	// numCellsPerDim -> sqrt(numParticles)
	CNum(config["N"], numParticles);
	CNum(config["dt"], dt); // Time step
	CNum(config["D"], D); // Diffusion constant
	CNum(config["T"], T); // Temperature
	CNum(config["r0"], r0); // radius of particles
	CNum(config["stepLimit"], stepLimit);
	CNum(config["sampleStep"], stepSample);
	CNum(config["seed"], seed);
	CNum(config["lambda"], lambda); // Yukawa parameter
	CNum(config["V0"], V0); // Yukawa parameter
	CNum(config["L"], L); // Simulation domain lateral size
	CNum(config["rcut"], rcut); // Cutoff radius
	
	infile.close();
	config.clear();
	
	cout << "Input lido" << endl;
}

void writePropsHeader() {
	propsDat << "t" << "\t" << "u" << "\t" << "p" << "\t" << "msd" << \
	"\t" << "sigmaxy" << "\t" << "f2" << "\t" << "ek" << "\t" << "e" << endl;
}

void writeProps() {
	averageProps();
	propsDat << timeNow << "\t" << uSum / (real)numParticles << "\t" <<  P << "\t" << msd;
	propsDat << "\t" << sigmaxy << "\t" << f2Med << "\t";
	propsDat << ekSum / (real)numParticles << "\t" << eTot / (real)numParticles << endl;
}

void averageProps() {
	int n;
	P = 0.;
	msd = 0;
	real f2sum = 0;
	real du2sum = 0;
	f2Med = 0;
	DO_PART {
		P += (particles[n].pressure);
		
		VecR vec_d;
		SubtractVectors(vec_d, particles[n].r, particles[n].r0);
		
		msd += DotProd(vec_d, vec_d);
		
		f2sum += particles[n].f2;
		du2sum += particles[n].ddu;
	 }
	 P = T*numParticles/(L*L) + P / (L*L*NDIM*numParticles);
	 msd /= numParticles;
	 f2Med = (f2sum/numParticles) - T*du2sum/numParticles;
}

void PrintElapsedTime(chrono::steady_clock::time_point start) {
    auto end = chrono::steady_clock::now();
    
    if (chrono::duration_cast<chrono::seconds>(end - start).count() < 60) {
        cout << "Elapsed time : " \
            << chrono::duration_cast<chrono::seconds>(end - start).count() \
            << " sec";
    }
    else if (chrono::duration_cast<chrono::seconds>(end - start).count() > 60) {
        cout << "Elapsed time : " \
            << chrono::duration_cast<chrono::minutes>(end - start).count() << " min";
    }
    else if (chrono::duration_cast<chrono::minutes>(end - start).count() > 60) {
        cout << "Elapsed time : " \
            << chrono::duration_cast<chrono::hours>(end - start).count() << " hours";
    }
    cout << endl;
    
    ofstream tfile;
	tfile.open ("time.txt");
	tfile << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;
	tfile.close();
}

// ========================================================== //
void computeForces(int in_stage) {                            //
	int n, j, i;                                              //
	switch(in_stage) {                                        //
		// ++++++++++++++++++++++++++++++++++++++++++++++++++ //
		case 1:                                               //
		{                                                     //
			sigmaxy = 0.0;                                    //
			DO_PART {                                         //
				VZero ( particles[n].fa );                    //
				particles[n].pressure = 0;                    //
				particles[n].ddu = 0.0;                       //
				particles[n].f2 = 0.0;                        //
			}                                                 //
			uSum = 0.;                                        //
			                                                  //
			for (i = 0; i < numParticles; i++) {              //
				for (j = i + 1; j < numParticles-1; j++) {    //
					pairForce(particles[i], particles[j], 1); //
				}                                             //
			}                                                 //
			break;                                            //
		}                                                     //
		// ++++++++++++++++++++++++++++++++++++++++++++++++++ //
		case 2:                                               //
		{                                                     //
			DO_PART {                                         //
				VZero(particles[n].fb);                       //
				particles[n].pressure = 0;                    //
			}                                                 //
			uSum = 0.;                                        //
			                                                  //
			for (i = 0; i < numParticles; i++) {              //
				for (j = i + 1; j < numParticles-1; j++) {    //
					pairForce(particles[i], particles[j], 2); //
				}                                             //
			}                                                 //
		}                                                     //
		// ++++++++++++++++++++++++++++++++++++++++++++++++++ //
	}                                                         //
}                                                             //
// ========================================================== //

// ===============================================================================
void pairForce(Particle& par1, Particle& par2, int in_stage) {
	VecR vec_dr, vec_f12, vec_n;
	real F = 0.0;
	real V, dr;
	switch(in_stage) {
		// -----------------------------------------------------------------------
		case 1:
		{
			GetVecDr(vec_dr, par1.r, par2.r);
			ApplyPBC(vec_dr);
			dr = Norm(vec_dr.x, vec_dr.y);
			GetVecN(vec_n, vec_dr, dr);
			
			if (dr < rcut) {
				V = ( V0 * exp(-lambda*(dr-1)) / dr );
				F = V  * (1 + lambda*dr) / dr;
			}
			else {
				V = 0.0;
				F = 0.0;
			}
			uSum += V;
			GetVecF(vec_f12, F, vec_n);
			
			AddForce(par1.fa, vec_f12.x, vec_f12.y);
			AddForce(par2.fa, -vec_f12.x, -vec_f12.y);
			
			// Virial pressure
			// Minus signs in force and dr cancel
			// So that the pressure in both bodies is the same
			par1.pressure += DotProd(vec_f12, vec_dr); 
			par2.pressure += DotProd(vec_f12, vec_dr); 
			
			sigmaxy +=  vec_dr.x * vec_f12.y / 4.0;
			
			if (ALGORITHM != "VERLET") {
				throw runtime_error("TODO Verificar isto");
				par1.ddu += 0.5 * V * ( 2/(dr*dr) + 2*lambda/dr + lambda*lambda );
				par2.ddu += 0.5 * V * ( 2/(dr*dr) + 2*lambda/dr + lambda*lambda );
			}
			break;
		}
		// -----------------------------------------------------------------------
		case 2: {
			// GetVecDr(vec_dr, par1.R, par2.R);
			ApplyPBC(vec_dr);
			dr = Norm(vec_dr.x, vec_dr.y);
			GetVecN(vec_n, vec_dr, dr);
			
			if (dr < rcut) {
				V = ( V0 * exp(-lambda*(dr-1)) / dr );
				F = V * (1 + lambda * dr) / dr;
			}
			else {
				V = 0.0;
				F = 0.0;
			}
			uSum += V;
			GetVecF(vec_f12, F, vec_n);
			
			AddForce( par1.fb, vec_f12.x, vec_f12.y );
			AddForce( par2.fb, -vec_f12.x, -vec_f12.y );
			
			// Virial pressure
			// Or is it vec_n? I am almost sure it's vec_dr
			// Minus signs in force and dr cancel
			par1.pressure += DotProd( vec_f12, vec_dr );
			par2.pressure += DotProd( vec_f12, vec_dr );
			
			sigmaxy +=  vec_dr.x * vec_f12.y / 4.0;
			
			par1.ddu += 0.5 * V * ( 2/(dr*dr) + 2*lambda/dr + lambda*lambda );
			par2.ddu += 0.5 * V * ( 2/(dr*dr) + 2*lambda/dr + lambda*lambda );
			
			break;
		} // End case 2
		// -----------------------------------------------------------------------
	} // End switch
} // End function
// ===============================================================================
