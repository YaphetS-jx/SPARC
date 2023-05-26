/**
 * @file    exchangeCorrelation.h
 * @brief   This file declares the functions for calculating exchange-correlation components.
 *
 * @authors Qimen Xu <qimenxu@gatech.edu>
 *          Abhiraj Sharma <asharma424@gatech.edu>
 *          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
 * 
 * Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
 */
 

#ifndef EXCHANGECORRELATION_H
#define EXCHANGECORRELATION_H

#include "isddft.h"

// remove later
/**
* @brief Structure to hold the constants used in xc functional
**/
typedef struct _XCCST_OBJ {
double alpha_zeta2; 
double alpha_zeta;
double beta;
double fsec_inv;
double kappa_pbe;
double mu;
double rsfac;
double kappa;
double mu_divkappa_pbe;
double mu_divkappa;

double ec0_aa; double ec1_aa; double mac_aa;
double ec0_a1; double ec1_a1; double mac_a1;
double ec0_b1; double ec1_b1; double mac_b1;
double ec0_b2; double ec1_b2; double mac_b2;
double ec0_b3; double ec1_b3; double mac_b3;
double ec0_b4; double ec1_b4; double mac_b4;

double piinv;
double third;
double twom1_3;
double sixpi2_1_3;
double sixpi2m1_3;
double threefourth_divpi;
double gamma;
double gamma_inv;
double beta_gamma;
double factf_zeta;
double factfp_zeta;
double coeff_tt;
double sq_rsfac;
double sq_rsfac_inv;
} XCCST_OBJ;

void xc_constants_init(XCCST_OBJ *xc_cst, SPARC_OBJ *pSPARC);
void Calculate_Vxc_GGA_PBE(SPARC_OBJ *pSPARC, XCCST_OBJ *xc_cst, double *rho);
void Calculate_Vxc_GSGA_PBE(SPARC_OBJ *pSPARC, XCCST_OBJ *xc_cst, double *rho);
void Calculate_Exc_GGA_PBE(SPARC_OBJ *pSPARC, double *electronDens);
void Calculate_Exc_GSGA_PBE(SPARC_OBJ *pSPARC, double *electronDens);


////////////////////////////////////////////////////////////////////////////////////

/**
* @brief  Calculate exchange correlation potential
**/
void Calculate_Vxc(SPARC_OBJ *pSPARC);


/**
* @brief  Calculate exchange correlation energy
**/
void Calculate_Exc(SPARC_OBJ *pSPARC, double *electronDens);


/**
 * @brief   slater exchange
 */
void slater(int DMnd, double *rho, double *ex, double *vx);


/**
 * @brief   pw correaltion
 *          J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
 */
void pw(int DMnd, double *rho, double *ec, double *vc);


/**
 * @brief   pz correaltion
 *          J.P. Perdew and A. Zunger, PRB 23, 5048 (1981).
 */
void pz(int DMnd, double *rho, double *ec, double *vc);


/**
 * @brief   pbe exchange
 *
 * @param   iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
 * @param   iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
 * @param   iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
 * @param   iflag=4  Zhang-Yang Revised PBE: Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998)
 */
void pbex(int DMnd, double *rho, double *sigma, int iflag, double *ex, double *vx, double *v2x);


/**
 * @brief   pbe correlation
 *
 * @param   iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
 * @param   iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
 * @param   iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
 */
void pbec(int DMnd, double *rho, double *sigma, int iflag, double *ec, double *vc, double *v2c);


/**
 * @brief   slater exchange - spin polarized 
 */
void slater_spin(int DMnd, double *rho, double *ex, double *vx);


/**
 * @brief   pw correaltion - spin polarized 
 *          J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
 */
void pw_spin(int DMnd, double *rho, double *ec, double *vc);


/**
 * @brief   pz correaltion - spin polarized 
 *          J.P. Perdew and A. Zunger, PRB 23, 5048 (1981).
 */
void pz_spin(int DMnd, double *rho, double *ec, double *vc);


/**
 * @brief   pbe exchange - spin polarized
 *
 * @param   iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
 * @param   iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
 * @param   iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
 * @param   iflag=4  Zhang-Yang Revised PBE: Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998)
 */
void pbex_spin(int DMnd, double *rho, double *sigma, int iflag, double *ex, double *vx, double *v2x);


/**
 * @brief   pbe correlation - spin polarized
 *
 * @param   iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
 * @param   iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
 * @param   iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
 */
void pbec_spin(int DMnd, double *rho, double *sigma, int iflag, double *ec, double *vc, double *v2c);


/**
 * @brief   calculate square norm of gradient
 */
void calculate_square_norm_of_gradient(SPARC_OBJ *pSPARC, 
        double *rho, double *mag, int DMnd, int ncol, 
        double *sigma, double *Drho_x, double *Drho_y, double *Drho_z);

/**
 * @brief   calculate square norm of a set of vector
 */ 
void compute_norm_square(SPARC_OBJ *pSPARC, double *norm2, int DMnd, double *v1, double *v2, double *v3);


/**
 * @brief   add core electron density if needed
 */
void add_rho_core(SPARC_OBJ *pSPARC, double *rho_in, double *rho_out, int ncol);


/**
 * @brief   compute Drho times v2xc
 */
void Drho_times_v2xc(SPARC_OBJ *pSPARC, int DMnd, int ncol, double *Drho_x, double *Drho_y, double *Drho_z, double *v2xc);

/**
 * @brief   Calculate PBE short ranged exchange
 *          Taken from Quantum Espresson
 */
void pbexsr(double rho, double grho, double omega, double *e_xc_sr, double *XCPotential_sr, double *Dxcdgrho_sr);

/**
 * @brief   Calculate PBE short ranged enhancement factor
 *          Taken from Quantum Espresson
 */
void wpbe_analy_erfc_approx_grad(double rho, double s, double omega, double *Fx_wpbe, double *d1rfx, double *d1sfx);

/**
 * @brief  Calculate exchange correlation energy density
 **/
void Calculate_xc_energy_density(SPARC_OBJ *pSPARC, double *ExcRho);

#endif // EXCHANGECORRELATION_H


