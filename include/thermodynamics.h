/**@<thermodynamics.h>::*/

extern double
Rgas /* erg mole⁻¹ Kelvin⁻¹ */ ,
adbtc /* = 2 / 5 valid for diatomic gases only */ ,
AvogadrosNumber /* = 6.0221415e+23 particles per mole */ ,
//rho0 = 2 / 22.414e+03 /* g cm⁻³ for H2 */ ,
rho0 /* = 16 / 22.414e+03 g cm⁻³ for O2 */ ,
H2 /* = 2 / 6.0221415e+23*/, u0, Temp, O2/* = 32 / 6.0221415e+23 */;

extern double u[NMAX], udot[NMAX], P[NMAX], Paccel[NMAX][DIM];

extern void sph_quantities(int n);