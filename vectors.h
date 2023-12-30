#ifndef V_DEFS
#define V_DEFS

typedef struct {
	real x, y;
} VecR;

#define VZero(v) \
	(v).x = 0., \
	(v).y = 0.
	
#define VProd(v, c) \
	v.x *= c, \
	v.y *= c
	
#define Dist(x, y) \
	sqrt(x*x + y*y)
	
#define Norm(x, y) \
	sqrt(x*x + y*y)
	
#define DotProd(v1, v2) \
	(v1).x * (v2).x + (v1).y * (v2).y

#define GetDr(DR, R1, R2) \
	DR.x = R1.x - R2.x, \
	DR.y = R1.y - R2.y

#define AddForce(f0, ffx, ffy) \
	(f0).x += ffx, \
	(f0).y += ffy
	
#define FLJ(dr, epsLJ, sigLJ) \
	(48 * epsLJ / (sigLJ*sigLJ)) * ( pow(sigLJ/dr, 14) - 0.5 * pow(sigLJ/dr, 8) );

#define ULJ(dr, epsLJ, sigLJ) \
	( 4 * epsLJ * ( pow(sigLJ/dr, 12) - pow(sigLJ/dr, 6) ) + epsLJ );

#define SubtractVectors(vv, vv1, vv2) \
	vv.x = vv1.x - vv2.x, \
	vv.y = vv1.y - vv2.y

#define GetVecDr(VecDr, RR1, RR2) \
	VecDr.x = RR1.x - RR2.x, \
	VecDr.y = RR1.y - RR2.y

#define ApplyPBC(rr) \
	(rr).x -= L * round(rr.x/L), \
	(rr).y -= L * round(rr.y/L)

#define GetVecN(VecN, VecDr, Dr) \
	VecN.x = VecDr.x/Dr, \
	VecN.y = VecDr.y/Dr
	
#define GetVecF(VecF, FF, VecN) \
	VecF.x = FF * VecN.x, \
	VecF.y = FF * VecN.y

#define CNum(s, num) istringstream ( s ) >> num

#define GenRand01 \
	(double)rand() / RAND_MAX

#define RStep(rinit, ff, ww) \
	(rinit).x += (D/T) * (ff).x * dt + (ww.x), \
	(rinit).y += (D/T) * (ff).y * dt + (ww.y)
	
#define GenerateNoise(ww, s) \
	(ww).x = distribution(generator) * s, \
	(ww).y = distribution(generator) * s

#endif

