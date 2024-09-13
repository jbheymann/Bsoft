/**
@file	mg_ctf_focal.h
@brief	Header file for CTF (contrast transfer function) functions
@author 	Bernard Heymann
@date	Created: 20000426
@date	Modified: 20240527
**/

#include "rwimg.h"
#include "ctf.h"

// Function prototypes
//vector<CTFparam>	ctf_series(CTFparam& cp, double& dstart, double& dend, double& dinc);
Bimage*		img_ctf_focal_series(CTFparam& cp, double def_min, double def_max, double def_inc,
				Vector3<long> size, Vector3<double> sam, double hires, double lores, int flag=0);
Bimage*		img_ctf_focal_series(CTFparam& cp, long nimg, Vector3<long> size, Vector3<double> sam,
				double hires, double lores, int flag);
Bimage*		img_ctf_focal_series(CTFparam& cp, vector<double>& dfocus,
				Vector3<long> size, Vector3<double> sam, double hires, double lores, int flag=0);
//Bimage*		img_ctf_focal_fit(Bimage* p, CTFparam& cp, vector<double>& dfocus, double hires, double lores,
//				double tmax, double Bfactor, long maxiter);
Bimage*		img_ctf_focal_fit(Bimage* p, CTFparam& cp, double hires, double lores,
				double tmax, double Bfactor, long maxiter);
int			img_fspace_weigh_sphere(Bimage* p, double volt);
Bimage*		img_fspace_extract_sphere(Bimage* p, CTFparam cp, vector<double>& dfocus);
Bimage*		img_ctf_apply_to_focal_series(Bimage* p, CTFparam cp, vector<double>& dfocus);
double		img_ctf_refine_focus_step(Bimage* p, double wl, double& focus_step,
				double focus_inc, double hires, double lores, double linear);
Bimage*		img_ctf_focal_reconstruct(Bimage* p, CTFparam cp, vector<double>& dfocus, double hires, int flag);
Bimage*		img_ctf_focal_reconstruct(Bimage* p, CTFparam cp, double hires, double lores=0, int flag=0);
Bimage*		img_ctf_focal_exit_wave_reconstruct(Bimage* p, CTFparam cp,
				double hires, long max_iter, double tol, int ctf_flag, int exit_flag);

