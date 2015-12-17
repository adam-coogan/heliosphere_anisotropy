#ifndef DIFFUSIONTENSOR
#define DIFFUSIONTENSOR 

/*!
 * Diffusion tensor container.
 */
struct DiffusionTensor {
    //! Tensor elements 
    double rr, phph, rph, thth;

    //! Derivatives of tensor's elements 
    double rr_dr, rph_dr;

    // Zero in simple models.  Can implement in a subclass if necessary.
    //double phph_dr, thth_dr, rr_dth, phph_dth, rph_dth, thth_dth;
    //double rth_dr, thph_dr;
    //double rth_dth, thph_dth;
};

#endif


