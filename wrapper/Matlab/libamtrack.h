/* TEST HEADER FOR MATLAB
 * 28.01.2011
 */  

/**
 * Returns CSDA range (in m) from pstar tables for given energy.
 * In case of ions a simple scaling procedure (A/Z^2) will be used (even effective charge will be neglected)
 * @param[in]   E_MeV_u                  energy of particle
 * @param[in]   particle_no              type of the particle
 * @see          AT_DataParticle.h for definition
 * @param[in]   material_no              material index
 * @see          AT_DataMaterial.h for definition
 * @return      CSDA_range_m             CSDA range
 */
double AT_CSDA_range_m_single(  const double  E_MeV_u,
    const long    particle_no,
    const long    material_no);

/**
 * Returns CSDA range (in m) from pstar tables for given energy.
 * In case of ions a simple scaling procedure (A/Z^2) will be used (even effective charge will be neglected)
 * @param[in]   number_of_particles      number of particle types in the mixed particle field
 * @param[in]   E_MeV_u                  energy of particles in the mixed particle field (array of size number_of_particles)
 * @param[in]   particle_no              type of the particles in the mixed particle field (array of size number_of_particles)
 * @see          AT_DataParticle.h for definition
 * @param[in]   material_no              material index
 * @see          AT_DataMaterial.h for definition
 * @param[out]  CSDA_range_m            (array of size number_of_particles) to be allocated by the user which will be used to return the results
 */
void AT_CSDA_range_m(  const long  number_of_particles,
	    const double  E_MeV_u[],
	    const long    particle_no[],
	    const long    material_no,
	    double        CSDA_range_m[]);