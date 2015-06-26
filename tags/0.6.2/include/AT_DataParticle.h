#ifndef AT_DATAPARTICLE_H_
#define AT_DATAPARTICLE_H_

/**
 * @brief Particle properties
 */


/*
 *    AT_DataParticle.h
 *    ==============
 *
 *    Copyright 2006, 2010 The libamtrack team
 *
 *    This file is part of the AmTrack program (libamtrack.sourceforge.net).
 *
 *    AmTrack is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    AmTrack is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with AmTrack (file: copying.txt).
 *    If not, see <http://www.gnu.org/licenses/>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "AT_Error.h"
#include "AT_NumericalRoutines.h"

/**
 * Total number of particles in data table
 */
#define PARTICLE_DATA_N    96

/**
 * Particle no of proton
 */
#define PARTICLE_PROTON_NUMBER   1001

/**
 * Maximum length of particle name acronym
 */
#define PARTICLE_NAME_NCHAR 6


/**
 * TODO
 */
typedef struct {
  const long      n;
  const long      Z[PARTICLE_DATA_N];
  const double    atomic_weight[PARTICLE_DATA_N];
  const char*     element_name[PARTICLE_DATA_N];
  const char*     element_acronym[PARTICLE_DATA_N];
  const double    density_g_cm3[PARTICLE_DATA_N];
  const double    I_eV_per_Z[PARTICLE_DATA_N];
} AT_particle_data_struct;


/**
 * TODO
 */
static const AT_particle_data_struct AT_Particle_Data = {
    PARTICLE_DATA_N,
    {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
    61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
    71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
    81, 82, 83, 84, 86, 88, 89, 90, 91, 92,
    93, 94, 95, 96, 97, 98
    },
    {
    1.0079, 4.0026, 6.941, 9.0122, 10.811, 12.0107, 14.0067, 15.9994, 18.9984, 20.1797,
    22.9897, 24.305, 26.9815, 28.0855, 30.9738, 32.065, 35.453, 39.948, 39.0983, 40.078,
    44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845, 58.9332, 58.6934, 63.546, 65.39,
    69.723, 72.64, 74.9216, 78.96, 79.904, 83.8, 85.4678, 87.62, 88.9059, 91.224,
    92.9064, 95.94, 98, 101.07, 102.9055, 106.42, 107.8682, 112.411, 114.818, 118.71,
    121.76, 127.6, 126.9045, 131.293, 132.9055, 137.327, 138.9055, 140.116, 140.9077, 144.24,
    145, 150.36, 151.964, 157.25, 158.9253, 162.5, 164.9303, 167.259, 168.9342, 173.04,
    174.967, 178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.9665, 200.59,
    204.3833, 207.2, 208.9804, 209, 222, 226, 227, 232.0381, 231.0359, 238.0289,
    237, 244, 243, 247, 247, 251
    },
    {
    "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon",
    "Sodium", "Magnesium", "Aluminum", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium",
    "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc",
    "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium",
    "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium", "Indium", "Tin",
    "Antimony", "Tellurium", "Iodine", "Xenon", "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium",
    "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium", "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium",
    "Lutetium", "Hafnium", "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinum", "Gold", "Mercury",
    "Thallium", "Lead", "Bismuth", "Polonium", "Radon", "Radium", "Actinium", "Thorium", "Protactinium", "Uranium",
    "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium"
    },
    {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "Rn", "Ra", "Ac", "Th", "Pa", "U",
    "Np", "Pu", "Am", "Cm", "Bk", "Cf"
    },
    {
    0.09, 0.18, 0.53, 1.85, 2.34, 2.26, 1.25, 1.43, 1.7, 0.9,
    0.97, 1.74, 2.7, 2.33, 1.82, 2.07, 3.21, 1.78, 0.86, 1.55,
    2.99, 4.54, 6.11, 7.19, 7.43, 7.87, 8.9, 8.9, 8.96, 7.13,
    5.91, 5.32, 5.72, 4.79, 3.12, 3.75, 1.63, 2.54, 4.47, 6.51,
    8.57, 10.22, 11.5, 12.37, 12.41, 12.02, 10.5, 8.65, 7.31, 7.31,
    6.68, 6.24, 4.93, 5.9, 1.87, 3.59, 6.15, 6.77, 6.77, 7.01,
    7.3, 7.52, 5.24, 7.9, 8.23, 8.55, 8.8, 9.07, 9.32, 6.9,
    9.84, 13.31, 16.65, 19.35, 21.04, 22.6, 22.4, 21.45, 19.32, 13.55,
    11.85, 11.35, 9.75, 9.3, 9.73, 5.5, 10.07, 11.72, 15.4, 18.95,
    20.2, 19.84, 13.67, 13.5, 14.78, 15.1
    },
    {
    13.5984, 24.5874, 5.3917, 9.3227, 8.298, 11.2603, 14.5341, 13.6181, 17.4228, 21.5645,
    5.1391, 7.6462, 5.9858, 8.1517, 10.4867, 10.36, 12.9676, 15.7596, 4.3407, 6.1132,
    6.5615, 6.8281, 6.7462, 6.7665, 7.434, 7.9024, 7.881, 7.6398, 7.7264, 9.3942,
    5.9993, 7.8994, 9.7886, 9.7524, 11.8138, 13.9996, 4.1771, 5.6949, 6.2173, 6.6339,
    6.7589, 7.0924, 7.28, 7.3605, 7.4589, 8.3369, 7.5762, 8.9938, 5.7864, 7.3439,
    8.6084, 9.0096, 10.4513, 12.1298, 3.8939, 5.2117, 5.5769, 5.5387, 5.473, 5.525,
    5.582, 5.6437, 5.6704, 6.1501, 5.8638, 5.9389, 6.0215, 6.1077, 6.1843, 6.2542,
    5.4259, 6.8251, 7.5496, 7.864, 7.8335, 8.4382, 8.967, 8.9587, 9.2255, 10.4375,
    6.1082, 7.4167, 7.2856, 8.417, 10.7485, 5.2784, 5.17, 6.3067, 5.89, 6.1941,
    6.2657, 6.0262, 5.9738, 5.9915, 6.1979, 6.2817
    }
};


/**
 * TODO
 *
 * @param[in]  Z            atomic number
 * @param[in]  A            mass number
 * @return particle_no
 */
 long AT_particle_no_from_Z_and_A_single(  const long  Z,
    const long  A);


/**
 * Returns particle index number from given A and Z
 *
 * @param[in]  n            array sizes
 * @param[in]  Z            atomic numbers (array of size n)
 * @param[in]  A            mass number (array of size n)
 * @param[out] particle_no  corresponding particle index numbers (array of size n)
 * @return status code
 */
int AT_particle_no_from_Z_and_A( const long  n,
    const long  Z[],
    const long  A[],
    long  particle_no[]);


/**
 * Calculates mass number A for particle with given code number
 * @param[in] particle_no
 * @return A
 */
 long AT_A_from_particle_no_single(  const long  particle_no );


/**
 * Returns mass number for given particle number
 * @param[in]  n                length of arrays
 * @param[in]  particle_no      particle index number (array of size n)
 * @param[out] A                mass number (array of size n)
 * @return     return code
 */
int AT_A_from_particle_no(  const long  n,
    const long  particle_no[],
    long  A[]);


/**
 * Returns atomic weight for given Z
 * @param[in]  n                length of arrays
 * @param[in]  Z                atomic number (array of size n)
 * @param[out] atomic_weight    atomic weight (array of size n)
 * @return     return code
 */
int AT_atomic_weight_from_Z( const long  n,
    const long  Z[],
    double  atomic_weight[]);

/**
 * Calculates atomic number Z for particle with given code number
 * @param[in] particle_no
 * @return Z
 */
 long AT_Z_from_particle_no_single(  const long  particle_no );


/**
 * Returns atomic number for given particle number
 * @param[in]  n               length of arrays
 * @param[in]  particle_no     particle index number (array of size n)
 * @param[out] Z               atomic number (array of size n)
 * @return     return code
 */
int AT_Z_from_particle_no(  const long  n,
    const long  particle_no[],
    long  Z[]);


/**
 * Calculates atomic weight for particle with given code number
 * @param[in]  n              length of arrays
 * @param[in]  particle_no    particle index number (array of size n)
 * @param[out] atomic_weight  atomic weight (array of size n)
 * @return                    status code
 */
int AT_atomic_weight_from_particle_no(  const long  n,
    const long  particle_no[],
    double  atomic_weight[]);


/**
 * Return I value for given elements
 * @param[in]  n              number of elements
 * @param[in]  particle_no    particle index number (array of size n)
 * @param[out] I_eV           I value (array of size n)
 * @return                    status code
 */
int AT_I_eV_from_particle_no( const long  n,
    const long  particle_no[],
    double  I_eV[]);


/**
 * Returns nuclear spin from particle no
 * @param[in]  n              number of elements
 * @param[in]  particle_no    particle index number (array of size n)
 * @param[out] I              nuclear spin (array of size n)
 * @return                    status code
 */
int AT_nuclear_spin_from_particle_no_multi( const long  n,
    const long  particle_no[],
    double  I[]);


/**
 * Returns nuclear spin from particle no
 * @param[in]  particle_no    particle index number (array of size n)
 * @return     nuclear spin
 */
double AT_nuclear_spin_from_particle_no_single( const long  particle_no );


/**
 * Returns nuclear spin from Z and A
 * @param[in]  Z              atomic number
 * @param[in]  A              mass number
 * @return     nuclear spin
 */
double AT_nuclear_spin_from_Z_and_A( const long  Z,
    const long  A);


/**
 * Returns particle index name for given particle index number
 * @param[in]  particle_no    particle index number
 * @param[out] particle_name  corresponding particle name
 * @return status
 */
int AT_particle_name_from_particle_no_single( const long  particle_no,
    char* particle_name);

/**
 * Returns particle index number for given particle name
 * @param[in]  particle_name particle index name (array of size PARTICLE_NAME_NCHAR)
 * @return corresponding particle number
 */
long AT_particle_no_from_particle_name_single( const char particle_name[] );


/**
 * Returns particle index numbers for given particle names
 * @param[in]  n             TODO
 * @param[in]  particle_no   particle index numbers (array of size n)
 * @param[out] particle_name corresponding particle names
 * @return status
 */
int AT_particle_name_from_particle_no( const long  n,
    const long  particle_no[],
    char particle_name[][PARTICLE_NAME_NCHAR]);


/**
 * Returns particle names for given particle numbers
 * @param[in]  n              TODO
 * @param[in]  particle_name  particle names
 * @param[out] corresponding  particle numbers
 * @return status
 */
int AT_particle_no_from_particle_name( const long  n,
    char* particle_name[],
    long  particle_no[]);


/**
 * Returns Z for given elemental symbol
 * @param[in]  acronym        elemental symbols (array of size PARTICLE_NAME_NCHAR)
 * @return corresponding Z
 */

int AT_Z_from_element_acronym_single( const char acronym[] );



/**
 * Returns Z for given elemental symbols
 * @param[in]  n       number of elements
 * @param[in]  acronym elemental symbols (array of size n)
 * @param[out] Z       corresponding Z (array of size n)
 * @return status
 */

int AT_Z_from_element_acronym( const long n,
		char* acronym[],
		long Z[]);

/**
 * Returns elemental symbol for given Z
 * @param[in]  Z      atomic number
 * @param[out] acronym  corresponding elemental symbol (array of size PARTICLE_NAME_NCHAR)
 * return status
 */

int AT_element_acronym_from_Z_single( const long Z,
		char acronym[]);


/**
 * Returns elemental symbols for given Z
 * @param[in]  n       number of elements
 * @param[in]  Z       atomic numbers (array of size n)
 * @param[out] acronym corresponding elemental symbols (array of size n)
 * @return status
 */

int AT_element_acronym_from_Z( const long n,
		long Z[],
		char* acronym[]);




/**
 * Returns atomic weight for given elemental symbol
 * @param[in]  acronym   elemental symbols (array of size PARTICLE_NAME_NCHAR)
 * @return corresponding atomic weight
 */
double AT_atomic_weight_from_element_acronym_single( const char acronym[] );


/**
 * Returns atomic weight for given elemental symbols
 * @param[in]  n               number of elements in array
 * @param[in]  acronym         elemental symbols (array of size n)
 * @param[out] atomic_weight   corresponding A (array of size n)
 * @return status
 */
int AT_atomic_weight_from_element_acronym( const long n,
		char* acronym[],
		double atomic_weight[]);


/**
 * Returns the elemental density in g/cm3 for given elemental symbol
 * @param[in]  acronym   elemental symbols (array of size PARTICLE_NAME_NCHAR)
 * @return elemental density in g/cm3
 */
double AT_density_g_cm3_from_element_acronym_single( const char acronym[] );


/**
 * Returns the elemental densities in g/cm3 for given elemental symbols
 * @param[in]  n        number of elements in array
 * @param[in]  acronym  elemental symbols (array of size n)
 * @param[out] density  corresponding elemental densities in g/cm3 (array of size n)
 * @return status
 */
int AT_density_g_cm3_from_element_acronym( const long n,
		char* acronym[],
		double density[]);


/**
 * Returns the I value in eV for given elemental symbol
 * @param[in]  acronym  elemental symbol (array of size PARTICLE_NAME_NCHAR)
 * @return corresponding I value in eV
 */
double AT_I_eV_from_element_acronym_single( const char acronym[] );


/**
 * Returns the I values in eV for given elemental symbol
 * @param[in]  n        number of elements in array
 * @param[in]  acronym  elemental symbols (array of size n)
 * @param[out] I        corresponding I values in eV (array of size n)
 * @return status
 */
int AT_I_eV_from_element_acronym( const long n,
		char* acronym[],
		double I[]);


/**
 * Returns the electron density per cm3 for given elemental symbol
 * @param[in]   acronym  elemental symbol (array of size PARTICLE_NAME_NCHAR)
 * @return      corresponding electron density per cm3
 */
double AT_electron_density_cm3_from_element_acronym_single( const char acronym[] );


/**
 * Returns the electron densities per cm3 for given elemental symbols
 * @param[in]  n                       array length
 * @param[in]  acronym                 elemental symbols (array of size n)
 * @param[out] electron_density        corresponding electron densities per cm3 (array of size n)
 * @return status
 */
int AT_electron_density_cm3_from_element_acronym( const long n,
                                                 char* acronym[],
                                                 double electron_density[]);

#endif /* AT_DATAPARTICLE_H_ */
