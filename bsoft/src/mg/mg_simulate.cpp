/**
@file	mg_simulate.cpp
@brief	Generates projects and simulate images 
@author	Bernard Heymann
@date	Created: 20030805
@date	Modified: 20250408
**/

#include "mg_processing.h"
#include "mg_multislice.h"
#include "mg_ctf.h"
#include "model_map.h"
#include "model_transform.h"
#include "model_util.h"
//#include "molecule_to_map.h"
#include "mol_transform.h"
#include "mol_edit.h"
#include "Complex.h"
#include "Vector3.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Applies a complex CTF function to a Fourier transform.
@param 	*p			complex Fourier transform (modified).
@param 	&cp			CTF parameters.
@return int 			0.


	The CTF is applied as a multiplication with a complex number.
	Both input and output are complex transforms.


**/
/*int			img_apply_complex_CTF(Bimage* p, CTFparam& cp)
{
	
	long			i, n, x, y, z;
	double			xx, yy, zz, sx2, sy2, sz2, s2, env;
	Complex<double>	ctf;
	Vector3<double>	freq_scale(1.0/p->real_size());
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Applying a CTF to a complex Fourier transform" << endl;

	if ( verbose & VERB_PROCESS ) {
		cp.show();
	}
	
	for ( i=n=0; n<p->images(); n++ ) {
		for ( z=0; z<p->sizeZ(); z++ ) {
			zz = z;
			if ( z > (p->sizeZ() - 1)/2 ) zz -= p->sizeZ();
			sz2 = freq_scale[2]*zz;
			sz2 *= sz2;
			for ( y=0; y<p->sizeY(); y++ ) {
				yy = y;
				if ( y > (p->sizeY() - 1)/2 ) yy -= p->sizeY();
				sy2 = freq_scale[1]*yy;
				sy2 *= sy2;
				for ( x=0; x<p->sizeX(); x++, i++ ) {
					xx = x;
					if ( x > (p->sizeX() - 1)/2 ) xx -= p->sizeX();
					sx2 = freq_scale[0]*xx;
					sx2 *= sx2;
					s2 = sx2 + sy2 + sz2;
					ctf = cp.calculate_complex(sqrt(s2), atan2(yy,xx));
					env = cp.partial_coherence_and_energy_spread(s2);
					p->set(i, ((*p)[i] * ctf) * env);
				}
			}
		}
	}
	
	return 0;
}*/

/**
@brief 	Generates a project for projection simulations.
@param 	nfield			number of fields-of-view.
@param 	nmg				number of micrographs per field-of-view.
@param 	npart			number of particles per micrograph.
@param 	pixel_size		micrograph pixel size.
@param 	img_origin		image origin within the simulation box.
@param 	cp				CTF parameters).
@param 	def_min			defocus minimum (angstrom).
@param 	def_max			defocus maximum (angstrom).
@param 	dose			electron dose (e/angstrom^2).
@param 	tsigma			translation standard deviation (pixels).
@param 	&fieldbase		field base name.
@param 	&mgbase			micrograph base name.
@param 	&partbase		particle image base name.
@param 	fieldnumber		field-of-view number.
@param 	mgnumber		micrograph number.
@param 	partnumber		particle number.
@return Bproject*			project structure.
**/
Bproject*	project_generate(int nfield, int nmg, int npart,
				Vector3<double> pixel_size, double img_origin,
				CTFparam& cp, double def_min, double def_max, double dose,
				double tsigma, Bstring& fieldbase, Bstring& mgbase, Bstring& partbase,
				int fieldnumber, int mgnumber, int partnumber)
{
	random_seed();
	
	double			def_range = def_max - def_min;
	if ( def_range < 0 ) def_range = 0;
	
	int				i, j, k, block = 0;
	double			irm = 1.0/get_rand_max(), irm2 = 2*irm;
	Vector3<float>	t;
	vector<View2<double>>	v(npart);
	Bproject*		project = new Bproject;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bstring			field_id, mg_id;
	
	if ( verbose ) {
		cout << "Generating a new project:" << endl;
		cout << "Number of fields:               " << nfield << endl;
		cout << "Number of micrographs:          " << nmg << endl;
		cout << "Number of particles:            " << npart << endl;
	}

	for ( i=0; i<nfield; i++ ) {
		field_id = fieldbase + Bstring(fieldnumber, "_%06d");
		if ( project->field ) field = field_add(&field, field_id);
		else field = field_add(&project->field, field_id);
		for ( k=0; k<npart; k++ ) {
			v[k] = View2<double>(random()*irm2 - 1, random()*irm2 - 1, random()*irm2 - 1, M_PI*(random()*irm2 - 1));
			v[k].normalize();
		}
		for ( j=0; j<nmg; j++ ) {
			mg_id = mgbase + Bstring(fieldnumber, "_%03d") + Bstring(mgnumber, "_%03d");
			if ( field->mg ) mg = micrograph_add(&mg, mg_id);
			else mg = micrograph_add(&field->mg, mg_id);
			mg->block = block;
			block++;
			mg->fpart = partbase + Bstring(fieldnumber, "_%03d") + Bstring(mgnumber, "_%03d.grd");
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->update(cp);
			mg->ctf->defocus_average(def_range*random()*irm + def_min);
			mg->pixel_size = pixel_size;
			mg->dose = dose;
			mg->box_size[0] = mg->box_size[1] = (int) (2*img_origin);
			mg->box_size[2] = 1;
			mg->ctf->zero(1);
			mg->ctf->baseline(0, 1);
			mg->ctf->envelope(0, 3);
			mg->ctf->envelope(1, -M_PI*M_PI*mg->ctf->alpha()*mg->ctf->alpha()*mg->ctf->defocus_average()*mg->ctf->defocus_average());
			for ( k=0; k<npart; k++ ) {
				if ( mg->part ) part = particle_add(&part, k+partnumber);
				else part = particle_add(&mg->part, k+partnumber);
				v[k] = View2<double>(random()*irm2 - 1, random()*irm2 - 1, random()*irm2 - 1, M_PI*(random()*irm2 - 1));
				v[k].normalize();
				if ( tsigma > 0 ) t = vector3_xy_random_gaussian(0.0, (double)tsigma);
				part->view2(v[k]);
				part->loc[0] = part->loc[1] = img_origin;
				part->ori = t + img_origin;
				part->ori[2] = 0;
			}
			mgnumber++;
		}
		fieldnumber++;
	}	
	
	return project;
}

/**
@brief 	Generates a project for projection simulations of an asymmetric unit.
@param 	&symmetry_string   symmetry designation.
@param 	pixel_size		micrograph pixel size.
@param 	img_origin		image origin within the simulation box.
@param 	theta_step		step size in theta (radians).
@param 	phi_step		step size in phi (radians).
@param 	cp				CTF parameters).
@param 	defocus			defocus minimum (angstrom).
@param 	dose			electron dose (e/angstrom^2).
@param 	&mgbase			micrograph base name.
@param 	&partbase		particle image base name.
@return Bproject*			project structure.
**/
Bproject*	project_generate_asu(string& symmetry_string,
				Vector3<double> pixel_size, double img_origin,
				double theta_step, double phi_step,
				CTFparam& cp, double defocus, double dose, Bstring& mgbase, Bstring& partbase)
{
	int				k, npart;
	Bproject*		project = new Bproject;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;

	Bsymmetry 		sym(symmetry_string);
	vector<View2<double>>	view = sym.asymmetric_unit_views(theta_step, phi_step, 0);
	npart = view.size();

	field = field_add(&project->field, mgbase);
	mg = micrograph_add(&field->mg, mgbase);
	mg->block = 0;
	mg->fpart = partbase + ".grd";
	if ( !mg->ctf ) mg->ctf = new CTFparam;
	mg->ctf->update(cp);
	mg->ctf->defocus_average(defocus);
	mg->pixel_size = pixel_size;
	mg->dose = dose;
	mg->box_size[0] = mg->box_size[1] = (int) (2*img_origin);
	mg->box_size[2] = 1;
	for ( k=0; k<npart; ++k ) {
		if ( mg->part ) part = particle_add(&part, k+1);
		else part = particle_add(&mg->part, k+1);
		part->view2(view[k]);
		part->loc = part->ori = Vector3<float>(img_origin, img_origin, 0);
	}
	
	return project;
}

/**
@brief 	Generates potential images using a multislice calculation.
@param 	*molgroup		molecule group structure.
@param 	*water			block of water as solvent.
@param 	*project		project structure with parameters.
@param 	&fieldname		selected field name (if "" do all).
@param 	&mgname			selected micrograph (if "" do all).
@param 	partselect		selected particle ( if <1 do all).
@param 	size			size of simulation block (angstrom).
@param 	thickness		thickness of slices for the multislice calculation (angstrom).
@param 	resolution		resolution for the multislice calculation (angstrom).
@param 	Bfactor			B-factor to apply to the multislice calculation (angstrom^2).
@param 	pottype			type of potential to calculate (???).
@param 	&paramfile		parameter file for scattering curves (???).
@return int				0.
**/
int			project_generate_potential(Bmolgroup* molgroup, Bmolgroup* water, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				Vector3<long> size, double thickness, double resolution,
				double Bfactor, int pottype, Bstring& paramfile)
{
	random_seed();
	
	Vector3<double>			sam, box;
	sam = project->field->mg->pixel_size;
	if ( resolution < 2*project->field->mg->pixel_size[0] )
		resolution = 2*project->field->mg->pixel_size[0];
	
	if ( size.volume() < 1 ) {
		size[0] = (long) (molgroup->box[0]/sam[0]);
		size[1] = (long) (molgroup->box[1]/sam[1]);
		size[2] = (long) (molgroup->box[2]/sam[2]);
	}
	if ( water ) {
		if ( molgroup->box.length2() < water->box.length2() ) {
			size[0] = (long) (water->box[0]/sam[0]);
			size[1] = (long) (water->box[1]/sam[1]);
			size[2] = (long) (water->box[2]/sam[2]);
		}
		box = sam * size;
		water->box = Vector3<float>(box[0], box[1], box[2]);
	}
	box = sam * size;
	molgroup->box = Vector3<float>(box[0], box[1], box[2]);
	
	Vector3<double>		sampling3(sam[0], sam[1], sam[2]);
	
	if ( verbose )
		cout << "Simulation box:                 " << molgroup->box << " = " << molgroup->box.volume() << " A^3" << endl;
	
	int				j;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bmolgroup*		mol_rot = NULL;
	Bmolgroup*		mol_sim = NULL;
	Bimage			*ppot;
	Bstring			potname;
	double			ranshift = molgroup->box[0];
	Vector3<float>	t, rot_ori;
	Vector3<double>	ori;
		
	for ( field = project->field; field; field = field->next ) 
		if ( fieldname.length() < 1 || field->id == fieldname ) {
		if ( verbose )
			cout << "Calculating potentials for field " << field->id << endl;
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mgname.length() < 1 || mg->id == mgname ) {
			if ( verbose )
				cout << "Calculating potentials for micrograph " << mg->id << endl;
			for ( part = mg->part; part; part = part->next )
				if ( partselect < 1 || part->id == partselect ) {
//				potname = number_filename(mg->fpart, part->id, 3);
				potname = mg->fpart;
				potname = potname.pre_rev('.') + Bstring(part->id, "_%03d.") + potname.post_rev('.');
				if ( access(potname.c_str(), F_OK) != 0 ) { // Skip if the file already exists
					if ( verbose )
						cout << "Calculating potentials for particle " << part->id << endl;
					t = (part->ori - part->loc) * sampling3;
//					cout << "origin=" << origin << " shift=" << t << endl;
					mol_rot = molgroup_rotate_from_view(molgroup, part->view2(), rot_ori, t);
					if ( water ) {
						mol_sim = molgroup_copy(water);
						molgroup_coor_shift_PBC(mol_sim, vector3_random(-ranshift, ranshift));
						molgroup_insert(mol_sim, mol_rot, 2);
						molgroup_kill(mol_rot);
					} else {
						mol_sim = mol_rot;
					}
					ori = {part->ori[0], part->ori[1], part->ori[2]};
					ppot = img_calc_potential(mol_sim, size, ori, sam, thickness,
						resolution, Bfactor, paramfile, pottype);
					ppot->complex_to_real();
					for ( j=0; j<ppot->images(); j++ ) {
						ppot->image[j].view(part->view);
						ppot->image[j].origin(part->ori);
					}
					part->fpart = potname;
					write_img(potname, ppot, 0);
					molgroup_kill(mol_sim);
					delete ppot;
				}
			}
		}
	}
	
	return 0;
}


/**
@brief 	Generates projection images for a project.
@param 	*molgroup		molecule group structure.
@param 	*water			block of water as solvent.
@param 	*project		project structure with parameters.
@param 	&fieldname		selected field name (if "" do all).
@param 	&mgname			selected micrograph (if "" do all).
@param 	partselect		selected particle ( if <1 do all).
@param 	size			size of simulation block (angstrom).
@param 	thickness		thickness of slices for the multislice calculation (angstrom).
@param 	resolution		resolution for the multislice calculation (angstrom).
@param 	Bfactor			B-factor to apply to the multislice calculation (angstrom^2).
@param 	type			type of projection to calculate: 0 = ewald; 1 = projection; 2 = multislice.
@param 	&paramfile		parameter file for scattering curves (???).
@return int				0.
**/
int			project_generate_projections(Bmolgroup* molgroup, Bmolgroup* water, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				Vector3<long> size, double thickness, double resolution,
				double Bfactor, int type, Bstring& paramfile)
{
	random_seed();
	
	Vector3<double>			sam, box;
	sam = project->field->mg->pixel_size;
	if ( resolution < 2*project->field->mg->pixel_size[0] )
		resolution = 2*project->field->mg->pixel_size[0];
	
	if ( size.volume() < 1 ) {
		size[0] = (long) (molgroup->box[0]/sam[0]);
		size[1] = (long) (molgroup->box[1]/sam[1]);
		size[2] = (long) (molgroup->box[2]/sam[2]);
	}
	if ( water ) {
		if ( molgroup->box.length2() < water->box.length2() ) {
			size[0] = (long) (water->box[0]/sam[0]);
			size[1] = (long) (water->box[1]/sam[1]);
			size[2] = (long) (water->box[2]/sam[2]);
		}
		box = sam * size;
		water->box = Vector3<float>(box[0], box[1], box[2]);
	}
	box = sam * size;
	molgroup->box = Vector3<float>(box[0], box[1], box[2]);
	
	Vector3<double>		sampling3(sam[0], sam[1], sam[2]);
	
	if ( verbose )
		cout << "Simulation box:                 " << molgroup->box << " = " << molgroup->box.volume() << " A^3" << endl;
	
	long			j, npart, pottype(0);
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bmolgroup*		mol_rot = NULL;
	Bmolgroup*		mol_sim = NULL;
	Bimage*			ppot = NULL;
	Bimage*			pone = NULL;
	Bimage*			p = NULL;
	Bstring			imgname;
	double			ranshift = molgroup->box[0];
	Vector3<float>	t, rot_ori;
	Vector3<double>	ori;
		
	for ( field = project->field; field; field = field->next ) 
		if ( fieldname.length() < 1 || field->id == fieldname ) {
		if ( verbose )
			cout << "Calculating potentials for field " << field->id << endl;
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mgname.length() < 1 || mg->id == mgname ) {
			if ( access(mg->fpart.c_str(), F_OK) != 0 ) { // Skip if the file already exists
				if ( verbose )
					cout << "Calculating projections for micrograph " << mg->id << endl;
				for ( npart=0, part = mg->part; part; part = part->next, npart++ ) ;
				p = new Bimage(Float, TSimple, size[0], size[1], 1, npart);
				p->sampling(mg->pixel_size);
				for ( j=0, part = mg->part; part; part = part->next, ++j ) {
					t = (part->ori - part->loc) * sampling3;
//					cout << "origin=" << origin << " shift=" << t << endl;
					mol_rot = molgroup_rotate_from_view(molgroup, part->view2(), rot_ori, t);
					if ( water ) {
						mol_sim = molgroup_copy(water);
						molgroup_coor_shift_PBC(mol_sim, vector3_random(-ranshift, ranshift));
						molgroup_insert(mol_sim, mol_rot, 2);
						molgroup_kill(mol_rot);
					} else {
						mol_sim = mol_rot;
					}
					ori = {part->ori[0], part->ori[1], part->ori[2]};
					if ( type < 1 ) {		// Ewald sphere projection
					
					} else if ( type == 1 ) {	// Projection approximation
					
					} else if ( type == 2 ) {
						ppot = img_calc_potential(mol_sim, size, ori, sam, thickness,
							resolution, Bfactor, paramfile, pottype);
						img_calc_phase_grating(ppot, mg->ctf->volt(), 0);
						pone = img_calc_multi_slice(ppot, thickness, mg->ctf->volt(), 0);
					} else {
					}
//					pone->fft();
//					img_ctf_apply_complex(pone, *mg->ctf, 0, 1, 0.1, 0, 0);		
					pone->fft(FFTW_BACKWARD, 1);
					pone->complex_to_intensities();
					p->replace(j, pone);
					p->image[j].view(part->view);
					p->image[j].origin(part->ori);
//					part->fpart = imgname;
					molgroup_kill(mol_sim);
					delete ppot;
					delete pone;
				}
				write_img(mg->fpart, p, 0);
				delete p;
			}
		}
	}
	
	return 0;
}

/**
@brief 	Generates projection images for a project from a model.
@param 	*model			molecule group structure.
@param 	*water			block of water as solvent.
@param 	*project		project structure with parameters.
@param 	&fieldname		selected field name (if "" do all).
@param 	&mgname			selected micrograph (if "" do all).
@param 	partselect		selected particle ( if <1 do all).
@param 	size			size of simulation block (angstrom).
@param 	thickness		thickness of slices for the multislice calculation (angstrom).
@param 	resolution		high resolution limit (angstrom).
@param 	Bfactor			B-factor to apply to the multislice calculation (angstrom^2).
@param 	type			type of projection: 0 = ewald; 1 = projection; 2 = multislice, 3 = conjugate ewald.
@param 	ab_flag			Flag to apply CTF:: 0 = not; 1 = apply; -1 = apply conjugate.
@param 	&paramfile		parameter file for scattering curves (???).
@return int				0.
**/
int			project_generate_projections(Bmodel* model, Bmodel* water, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				Vector3<long> size, double thickness, double resolution,
				double Bfactor, int type, int ab_flag, Bstring& paramfile)
{
	random_seed();
	
	Vector3<double>			sam;
	sam = project->field->mg->pixel_size;
	if ( resolution < 2*project->field->mg->pixel_size[0] )
		resolution = 2*project->field->mg->pixel_size[0];

	vector<Vector3<double>>	bounds = models_calculate_bounds(model);
	vector<Vector3<double>>	bounds2 = models_calculate_bounds(model);
	Vector3<double>			box = bounds[1] - bounds[0];
	Vector3<double>			box2 = bounds2[1] - bounds2[0];
	if ( box.length2() < box2.length2() ) box = box2;
	
	if ( size.volume() < 1 ) {
		size[0] = (long) (box[0]/sam[0]);
		size[1] = (long) (box[1]/sam[1]);
		size[2] = (long) (box[2]/sam[2]);
	}
	box = sam * size;
	
	if ( verbose )
		cout << "Simulation box:                 " << box << " = " << box.volume() << " A^3" << endl;
	
	long			j, npart;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bmodel*			mod_rot = NULL;
	Bmodel*			mod_sim = NULL;
	Bimage			*ppot, *pone, *p;
	Bstring			imgname;
	double			ranshift = box[0];
	Vector3<float>	t, rot_ori;
	Vector3<double>	ori;
		
	for ( field = project->field; field; field = field->next ) 
		if ( fieldname.length() < 1 || field->id == fieldname ) {
		if ( verbose )
			cout << "Calculating projections for field " << field->id << endl;
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mgname.length() < 1 || mg->id == mgname ) {
			if ( access(mg->fpart.c_str(), F_OK) != 0 ) { // Skip if the file already exists
				if ( verbose )
					cout << "Calculating projections for micrograph " << mg->id << endl;
				for ( npart=0, part = mg->part; part; part = part->next, npart++ ) ;
				p = new Bimage(Float, TSimple, size[0], size[1], 1, npart);
				p->sampling(mg->pixel_size);
				p->origin(p->size()/2);
				for ( j=0, part = mg->part; part; part = part->next, ++j ) {
					t = (part->ori - part->loc) * sam;
//					cout << "origin=" << origin << " shift=" << t << endl;
					Matrix3			mat = part->view2().matrix();
					mat = mat.transpose();
					mod_rot = model->copy();
					models_rotate(mod_rot, mat, rot_ori, t);
					if ( water ) {
						mod_sim = water->copy();
						mod_sim->shift(vector3_random(-ranshift, ranshift));
						model_insert(mod_sim, mod_rot, 2);
						delete mod_rot;
					} else {
						mod_sim = mod_rot;
					}
					ori = {part->ori[0], part->ori[1], part->ori[2]};
					if ( type < 2 ) {		// Projection
						pone = new Bimage(Float, TComplex, size[0], size[1], 1, 1);
						pone->sampling(mg->pixel_size);
						pone->origin(pone->size()/2);
						img_electron_scattering(mod_sim, 1, pone, *mg->ctf, 1, 0, paramfile.str(), 1-type, ab_flag);
						pone->combine_ewald();
						pone->phase_shift(pone->image->origin());
						pone->fft(FFTW_BACKWARD, 1);
						pone->complex_to_real();
					} else {				// Multislice
						long	nslices(size[2]/thickness);
						ppot = new Bimage(Float, TComplex, size[0], size[1], 1, nslices);
						ppot->sampling(sam[0], sam[1], thickness);
						ppot->origin(ppot->size()/2);
						ppot->fourier_type(Standard);
						img_potential_from_model_slices(mod_sim, ppot, paramfile.str(), resolution);
						ppot->phase_shift(ppot->image->origin());
						ppot->fft(FFTW_BACKWARD, 0);
						img_calc_phase_grating(ppot, mg->ctf->volt(), 0);
						pone = img_calc_multi_slice(ppot, thickness, mg->ctf->volt(), 0);
						if ( ab_flag )
							img_ctf_apply_complex(pone, *mg->ctf, 0, 1, 0.1, 0, 0);		
						pone->fft(FFTW_BACKWARD, 1);
						pone->complex_to_intensities();
						delete ppot;
					}
					p->replace(j, pone);
					p->image[j].view(part->view);
					p->image[j].origin(part->ori);
					delete mod_sim;
					delete pone;
				}
				write_img(mg->fpart, p, 0);
				delete p;
			}
		}
	}
	
	return 0;
}

/**
@brief 	Generates projection images for a project from a 3D map.
@param 	*map			3D map to project.
@param 	*project		project structure with parameters.
@param 	&fieldname		selected field name (if "" do all).
@param 	&mgname			selected micrograph (if "" do all).
@param 	partselect		selected particle ( if <1 do all).
@param 	size			size of simulation block (angstrom).
@param 	resolution		high resolution limit (angstrom).
@param 	Bfactor			B-factor to apply to the multislice calculation (angstrom^2).
@param 	ew_flag			type of ewald sphere application: 0 = central section, 1 = ewald; -1 = conjugate ewald.
@param 	ab_flag			Flag to apply CTF:: 0 = not; 1 = apply; -1 = apply conjugate.
@return int				0.
**/
int			project_generate_projections(Bimage* map, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				double resolution, int ew_flag, int ab_flag)
{
	random_seed();
	
	Vector3<double>			sam;
	sam = project->field->mg->pixel_size;
	if ( resolution < 2*project->field->mg->pixel_size[0] )
		resolution = 2*project->field->mg->pixel_size[0];

	Vector3<double>			size(map->sizeX(), map->sizeY(), 1);
	
	long				j;
	Bfield* 			field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;
	Bimage*				p;
	double				volt(0);
	
	FSI_Kernel* 	kernel = new FSI_Kernel(8, 2);
		
	for ( field = project->field; field; field = field->next ) 
		if ( fieldname.length() < 1 || field->id == fieldname ) {
		if ( verbose )
			cout << "Calculating projections for field " << field->id << endl;
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mgname.length() < 1 || mg->id == mgname ) {
			if ( access(mg->fpart.c_str(), F_OK) != 0 ) { // Skip if the file already exists
				if ( verbose )
					cout << "Calculating projections for micrograph " << mg->id << endl;
				if ( mg->ctf )
					volt = mg->ctf->volt();
				vector<View2<double>>	views;
				for ( part = mg->part; part; part = part->next )
					views.push_back(part->view2());
				p = map->project(views, resolution, kernel, volt, ew_flag, 0);
				if ( mg->ctf )
					img_ctf_apply_ewald(p, *mg->ctf, 0, resolution, ab_flag, 0);
				for ( j=0, part = mg->part; part; part = part->next, ++j )
					p->phase_shift(j, part->ori - part->loc);					
				p->fft(FFTW_BACKWARD);
				write_img(mg->fpart, p, 0);
				delete p;
			}
		}
	}

	if ( kernel ) delete kernel;
	
	return 0;
}

/**
@brief 	Generates final images from a multislice calculation.
@param 	*project		project structure with parameters.
@param 	thickness			thickness of slices for the multislice calculation (angstrom).
@param 	resolution		resolution for the multislice calculation (angstrom).
@return int						0.
**/
int			project_generate_image(Bproject* project, double thickness, double resolution)
{
	Bstring			potname = project->field->mg->fpart;
	potname = potname.pre_rev('.') + Bstring(project->field->mg->part->id, "_%03d.") + potname.post_rev('.');
	
	Bimage*			ppot = read_img(potname, 0, 0);
	Vector3<long>	size(ppot->sizeX(), ppot->sizeY(), 1);
	delete ppot;
	
	Vector3<double>	realsize(size);
	realsize *= project->field->mg->pixel_size;
	
	if ( resolution < 2*project->field->mg->pixel_size[0] )
		resolution = 2*project->field->mg->pixel_size[0];
	
	long			i, j, k, npart;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bimage			*p, *pone;
	int				spacegroup = 1;
		
	long	imgsize = 2*size[0]*size[1];
	UnitCell		unit_cell(realsize[0], realsize[1], realsize[2], M_PI_2, M_PI_2, M_PI_2);
	
	for ( field = project->field; field; field = field->next ) {
		if ( verbose )
			cout << "Calculating images for field " << field->id << endl;
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( verbose )
				cout << "Calculating images for micrograph " << mg->id << endl;
			for ( npart=0, part = mg->part; part; part = part->next, npart++ ) ;
			p = new Bimage(Float, TComplex, size[0], size[1], 1, npart);
			p->unit_cell(unit_cell);
			p->space_group(spacegroup);
			p->sampling(mg->pixel_size);
			for ( i=0, part = mg->part; part; part = part->next, i++ ) {
				if ( verbose )
					cout << "Calculating image for particle " << part->id << endl;
				potname = mg->fpart;
				potname = potname.pre_rev('.') + Bstring(part->id, "_%03d.") + potname.post_rev('.');
				ppot = read_img(potname, 1, -1);
				ppot->simple_to_complex();
				img_calc_phase_grating(ppot, mg->ctf->volt(), 0);
				pone = img_calc_multi_slice(ppot, thickness, mg->ctf->volt(), 0);
				for ( j=i*imgsize, k=0; k<imgsize; j++, k++ ) p->set(j, (*pone)[k]);
				p->image[i] = ppot->image[0];
				delete ppot;
				delete pone;
			}
			if ( mg->ctf->defocus_average() )
//				img_apply_complex_CTF(p, *mg->ctf);
				img_ctf_apply_complex(p, *mg->ctf, 0, 1, 0.1, 0, 0);		
			p->fft(FFTW_BACKWARD, 2, Real);
			p->complex_to_intensities();
			p->statistics();
			write_img(mg->fpart, p, 0);
			delete p;
		}
	}
	
	return 0;
}

/**
@brief 	Applies imaging distortions to the final images from a multislice calculation.
@param 	*project		project structure with parameters.
@param 	poisson			flag to add Poisson noise.
@param 	gauss			width of gaussian noise to add (0=no noise).
@param 	kmtf			mass transfer decay constant (0=no decay).
@return int				0.
**/
int			project_apply_distortions(Bproject* project, int poisson, double gauss, double kmtf)
{
	double			dose_per_pixel;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;

	Bimage*			p;
	Bstring			insert("_tf.");				// Transfer functions(s)
	if ( poisson || gauss ) insert = "_tfn.";	// Transfer functions(s) and noise
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( verbose )
				cout << "Applying distortions to micrograph " << mg->id << endl;
			p = read_img(mg->fpart, 1, -1);
			if ( p->compound_type() == TComplex ) {
				p->fft(FFTW_BACKWARD, 2, Real);
				p->complex_to_intensities();
			} else {
				p->change_type(Float);
			}
			dose_per_pixel = mg->dose*mg->pixel_size[0]*mg->pixel_size[1];
			if ( dose_per_pixel )
				p->rescale(dose_per_pixel/p->average(), 0);
			if ( poisson )
//				img_add_poisson_noise(p);
				p->noise_poisson(p->average());
			if ( gauss < 1000 )
//				img_add_gaussian_noise(p, gauss);
				p->noise_gaussian(0, p->standard_deviation()/sqrt(gauss));
			if ( kmtf > 0 )
				p->fspace_weigh_B_factor(4*kmtf);
			mg->fpart = mg->fpart.pre_rev('.') + insert + mg->fpart.post_rev('.');
			write_img(mg->fpart, p, 0);
			delete p;
		}
	}
	
	return 0;
}


