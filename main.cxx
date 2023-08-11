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

using namespace std;

#define NDIM  2

typedef long double real;

#include "vectors.h"
#include "prototypes.h"

#define DO_PART for (n = 0; n < N; n ++)

typedef struct {
    VecR r, R, fa, fb, r0;
    int nPart;
    real rho; // radius
    real pa, pb; // pressure
    real f2, du2;
} Particle;

void pairForce(Particle& p1, Particle& p2, int stage);
void SRK_Step(int stage);

// *********************************************************************
// Initalization of variables

vector<Particle> particles;
vector<string> inVars;
real dt, D, timeNow, uSuma, uSumb, T, r0, V0, lambda, rmax, L, msd, rcut, P, sigmaxy, f2tdu2;
int stepLimit, printStep, stepCount, moreCycles;
int NParticles; // Must have an even square - ex. 25, 36...
unsigned seed;
ofstream ovito, propsDat;
VecR W, reg;

mt19937 generator;
normal_distribution<real> distribution(0.0, 1.0);

// *********************************************************************
// 

int main(int argc, char **argv)
{	
	auto start = chrono::steady_clock::now();
	cout << "Inicio do programa" << endl;
    ReadInput ();
    cout << "Input lido" << endl;
    SetupJob ();
    cout << "Setup finalizado. A executar o ciclo principal" << endl;
    
    moreCycles = 1;
    PrintSummary ();
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
}

// *********************************************************************
// Functions

void SingleStep ()
{
    ++stepCount;
    timeNow = stepCount * dt;
    
    computeForces(1);
    SRK_Step(1);
    computeForces(2);
    SRK_Step(2);
	
    if (stepCount % printStep == 0) {
        PrintSummary ();
        writeProps();
    }
}

void SetupJob ()
{	
	cout << "A inicializar o workflow" << endl;
    AllocArrays ();
    stepCount = 0;
    InitCoords();
    Init_F();
    VZero(W);
    ovito.open("ovito.xyz");
    timeNow = 0.0;
    uSuma = 0.0;
    uSumb = 0.0;
    stepCount = 0;
    propsDat.open("props.txt");
    writePropsHeader();
    generator = mt19937(seed);
    f2tdu2 = 0.0;
    cout << "Workflow inicializado" << endl;
}

void AllocArrays ()
{	
	cout << "A allocar os arrays" << endl;
	particles.resize(NParticles);
	cout << "Arrays allocados" << endl;
}

#define RStep(rinit, ff, ww) \
	(rinit).x += (D/T) * (ff).x * dt + (ww.x), \
	(rinit).y += (D/T) * (ff).y * dt + (ww.y)
	
#define GenerateNoise(ww, s) \
	(ww).x = distribution(generator) * s, \
	(ww).y = distribution(generator) * s

void CBD_Step ()
{
    int n;
	DO_PART {
		GenerateNoise( W, sqrt(2*dt*D) );
		//RStep( monomers[n].r, monomers[n].f, W);
	}
}

void SRK_Step (int stage)
{	
	int n;
	switch(stage) {
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		case 1: {
			DO_PART {
				GenerateNoise( W, sqrt(2*dt*D) );
				particles[n].R.x = particles[n].r.x;
				particles[n].R.y = particles[n].r.y;
				RStep( particles[n].R, particles[n].fa, W );
				//cout << monomers[n].R.x << " | " << monomers[n].R.y << endl;
			}
			break;
		}
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	}
}

#define GenRand01 \
	(double)rand() / RAND_MAX
	
void InitCoords() {
	real delta = L / (sqrt(NParticles));
	int n = 0;
	for (int i=0; i<sqrt(NMonomers); i++) {
		for (int j=0; j<sqrt(NMonomers); j++) {
			monomers[n].rho = r0, monomers[n].nMono = n;
			monomers[n].r.x = ((real)i+0.5)*delta;
			monomers[n].r.y = ((real)j+0.5)*delta;
			n++;
		}
	}
	
	n = 0;
	DO_MONO {
		monomers[n].r0.x = monomers[n].r.x;
		
		monomers[n].r0.y = monomers[n].r.y;
	}
	
	cout << "Fim da inicializacao das coordenadas" << endl;
}

void Init_F ()
{
    int n;
	DO_MONO {
		VZero (monomers[n].fa);
		VZero (monomers[n].fb);
	};
}

void PrintSummary ()
{	
	ovito << NMonomers << endl;
	ovito << "Lattice=\""<<L<<" 0.0 0.0 0.0 "<<L<<" 0.0 0.0 0.0 0.0\" Origin=\""<<0<<" ";
	ovito<<0<<" 0.0\"";
	ovito << " Properties=pos:R:2:Radius:R:1:Type:S:1:Identifier:I:1:Pressure:R:1";
	ovito << " t=" << timeNow << endl;
	
    int counter=0;
    for( counter=0; counter < NMonomers; counter++)
    {	
		Monomer& mon = monomers[counter];
		
		// Wrap periodic images
		real x_img = mon.r.x;
		//if (x_img > L ) x_img -= L;
		//if (x_img < 0 ) x_img += L;
		
		real y_img = mon.r.y;
		//if (y_img > L ) y_img -= L;
		//if (y_img < 0 ) y_img += L;

		ovito << x_img << "\t" <<  y_img << "\t" << mon.rho << \
		 "\t" << "3" << "\t" << mon.nMono << "\t" << (mon.pa + mon.pb)/2.0  << "\t" << endl;
    }
}

#define CNum(s, num) istringstream ( s ) >> num

void ReadInput() {
	cout << "A ler o input" << endl;
	ifstream infile("params.in");
	
	string strKey, strValue;
	
	map<string, string> config;
	
	while (infile >> strKey >> strValue) {
		config[strKey] = strValue;
	}
	
	CNum(config["NMonomers"], NMonomers);
	CNum(config["dt"], dt);
	CNum(config["D"], D);
	CNum(config["T"], T);
	CNum(config["r0"], r0);
	CNum(config["stepLimit"], stepLimit);
	CNum(config["sampleStep"], printStep);
	CNum(config["seed"], seed);
	CNum(config["lambda"], lambda);
	CNum(config["V0"], V0);
	CNum(config["L"], L);
	CNum(config["rcut"], rcut);
	
	cout << NMonomers << "|" << dt << "|" << D << "|" << \
	T << "|" << r0 << "|" << stepLimit << "|" << printStep << "|" << seed << "|" << \
	lambda << "|" << V0 << "|" << L << "|" << rcut << endl;
	
	infile.close();
	
	//config.~map();
	
	cout << "Input lido" << endl;
}

void writePropsHeader() {
	propsDat << "t" << "\t" << "U" << "\t" << "P" << "\t" << "msd" << \
	"\t" << "sigmaxy" << "\t" << "f2tdu2" << endl;
}

void writeProps() {
	averageProps();
	propsDat << timeNow << "\t" << (uSuma + uSumb) / (2*NMonomers) << "\t" <<  P << "\t" << msd << \
	"\t" << sigmaxy << "\t" << f2tdu2 << endl;
}

void averageProps() {
	int n;
	P = 0.;
	msd = 0;
	real f2sum = 0;
	real du2sum = 0;
	f2tdu2 = 0;
	DO_MONO {
		P += (monomers[n].pa + monomers[n].pb)/2;
		
		VecR vec_d;
		SubtractVectors(vec_d, monomers[n].r, monomers[n].r0);
		
		msd += DotProd(vec_d, vec_d);
		
		f2sum += monomers[n].f2;
		du2sum += monomers[n].du2;
	 }
	 P = T*NMonomers/(L*L) + P / (L*L*NDIM*NMonomers);
	 msd /= NMonomers;
	 f2tdu2 = (f2sum/NMonomers) - T*du2sum/NMonomers;
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


// =====================================================================
void computeForces(int stage) { 
	int n, j, i;
	switch(stage) {
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		case 1:
		{	
			sigmaxy = 0.0;
			DO_MONO {
				VZero (monomers[n].fa);
				monomers[n].pa = 0;
				monomers[n].du2 = 0.0;
				monomers[n].f2 = 0.0;
			}
			uSuma = 0.;
			
			for (i = 0; i < NMonomers; i++) {
				for (j = i + 1; j < NMonomers-1; j++) {
					pairForce(monomers[i], monomers[j], 1);
				}
			}
			break;
		}
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		case 2:
		{	
			DO_MONO {
				VZero (monomers[n].fb);
				monomers[n].pb = 0;
			}
			uSumb = 0.;
			
			for (i = 0; i < NMonomers; i++) {
				for (j = i + 1; j < NMonomers-1; j++) {
					pairForce(monomers[i], monomers[j], 2);
				}
			}
		}
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	}
}
// =====================================================================

#define ApplyPBC(rr) \
	(rr).x -= L * round(rr.x/L), \
	(rr).y -= L * round(rr.y/L)

#define GetVecDr(VecDr, RR1, RR2) \
	VecDr.x = RR1.x - RR2.x, \
	VecDr.y = RR1.y - RR2.y

#define GetVecN(VecN, VecDr, Dr) \
	VecN.x = VecDr.x/Dr, \
	VecN.y = VecDr.y/Dr
	
#define GetVecF(VecF, FF, VecN) \
	VecF.x = FF * VecN.x, \
	VecF.y = FF * VecN.y

// =====================================================================
void pairForce(Monomer& mon1, Monomer& mon2, int stage) {
	// TODO remove switch use args mon.r, mon.R
	VecR vec_dr, vec_f12, vec_n;
	real F = 0.0;
	real V, dr;
	switch(stage) {
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		case 1:
		{
			GetVecDr(vec_dr, mon1.r, mon2.r);
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
			uSuma += V;
			GetVecF(vec_f12, F, vec_n);
			
			AddForce(mon1.fa, vec_f12.x, vec_f12.y);
			AddForce(mon2.fa, -vec_f12.x, -vec_f12.y);
			
			// Virial pressure
			mon1.pa += DotProd(vec_f12, vec_dr); //
			mon2.pa += DotProd(vec_f12, vec_dr); // minus signs in force and dr cancel
			
			sigmaxy +=  vec_dr.x * vec_f12.y / 4.0;
			
			mon1.du2 += 0.5 * V * ( 2/(dr*dr) + 2*lambda/dr + lambda*lambda );
			
			break;
		}
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		case 2: {
			GetVecDr(vec_dr, mon1.R, mon2.R);
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
			uSumb += V;
			GetVecF(vec_f12, F, vec_n);
			
			AddForce(mon1.fb, vec_f12.x, vec_f12.y);
			AddForce(mon2.fb, -vec_f12.x, -vec_f12.y);
			
			// Virial pressure
			mon1.pb += DotProd(vec_f12, vec_dr); // or is it vec_n
			mon2.pb += DotProd(vec_f12, vec_dr); // minus signs in force and dr cancel
			
			sigmaxy +=  vec_dr.x * vec_f12.y / 4.0;
			
			
			mon1.du2 += 0.5 * V * ( 2/(dr*dr) + 2*lambda/dr + lambda*lambda );
			
			break;
		}
	} // end switch
}
// =====================================================================
