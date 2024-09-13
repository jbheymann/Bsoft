/**
@file	brecmod.cpp
@brief	A tool to merge data from micrograph reconstruction and model files.
@author Bernard Heymann
@date	Created: 20090828
@date 	Modified: 20110725
**/

#include "mg_processing.h"
#include "rwmg.h"
#include "rwimg.h"
#include "rwmodel.h"
#include "model_transform.h"
#include "marker.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bmodel*		project_model_generate(Bproject* project, Bmodel* temp, int flags);
int			project_model_consolidate(Bproject* project, Bmodel* model);

/* Usage assistance */
const char* use[] = {
" ",
"Usage: brecmod [options] in.star",
"--------------------------------",
"Merges data from micrograph and model files.",
" ",
"Actions:",
"-consolidate mod.star    Consolidate project and model parameters.",
"-images                  Convert particle images (use with -outmodel).",
"-submodels               Generate particle models (use with -outmodel).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling of images (A/pixel; a single value can be given).",
" ",
"Input:",
"-template model.star     Template to associate with each particle.",
" ",
"Output:",
"-output newmg.star       New micrograph parameter file.",
"-outmodel partmod.star   New model file from particles.",
"-outmark markmod.star    New model file from markers.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	Bstring			consolidate;			// Consolidate project and model parameters
	int				flags(0);				// Convert particle images to model maps
	Vector3<double>	sam;    				// Units for the three axes (A/pixel)
	Bstring			paramfile;				// Input parameter file name
	Bstring			tempfile;				// Input template model file name
	Bstring			outfile;				// Output micrograph parameter file name
	Bstring			partmodfile;			// Output particle model file name
	Bstring			markmodfile;			// Output marker model file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "consolidate" )
			consolidate = curropt->filename();
		if ( curropt->tag == "images" ) flags |= 1;
		if ( curropt->tag == "submodels" ) flags |= 2;
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "template" )
			tempfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "outmodel" )
			partmodfile = curropt->filename();
		if ( curropt->tag == "outmark" )
			markmodfile = curropt->filename();
    }
	option_kill(option);
	
	double			ti = timer_start();

	// Read all the micrograph parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);		
	string_kill(file_list);

	if ( !project ) {
		cerr << "Error: Input micrograph file not read!" << endl;
		bexit(-1);
	}
	
	if ( !project->rec ) {
		cerr << "Error: No reconstruction found!" << endl;
		bexit(-1);
	}

	
	// Read the model file
	Bmodel*			model = NULL;		
	Bmodel*			temp = NULL;
	
	if ( sam.volume() <= 0 )
		sam = project->rec->voxel_size;
	else
		project_set_mg_pixel_size(project, sam);
	
	if ( tempfile.length() ) {
		temp = read_model(tempfile.str(), paramfile.str());
		if ( !temp ) {
			cerr << "Error: Input template file " << tempfile << " not read!" << endl;
			bexit(-1);
		}
	}

	if ( consolidate.length() ) {
		model = read_model(consolidate.str(), paramfile.str());
		if ( !model ) {
			cerr << "Error: No model read!" << endl;
			bexit(-1);
		}
		project_model_consolidate(project, model);
	} else if ( partmodfile.length() ) {
		model = project_model_generate(project, temp, flags);
	} else if ( markmodfile.length() ) {
		model = model_from_markers(project->rec->mark, project->rec->origin, sam);
	}
	
	// Write an output parameter format file if a name is given
	if ( outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}

	// Write an output parameter format file if a name is given
    if ( partmodfile.length() ) {
		write_model(partmodfile.str(), model);
	}

	if ( markmodfile.length() ) {
		write_model(markmodfile.str(), model);
	}
	
	project_kill(project);
	model_kill(model);
	
	
		timer_report(ti);
	
	bexit(0);
}

/**
@brief 	Generates a model from each particle defined in a reconstruction.
@param 	*project	project parameter structure.
@param 	*temp		template model.
@param 	flags		conversion options: 1=images.
@return Bmodel*		new model referencing the particle map files.

	The model id is a combined reconstruction and particle id.
	The model type id is taken from the particle group id.
	The selection and FOM is taken from the particle properties.
	The model map file name is taken from the particle file name if it
	is defined, otherwise it is taken from the reconstruction structure
	and the image number is taken from the particle id.
	Flags:
		1	transform each particle images to standard orientation.
		2	write a model file for each particle.

**/
Bmodel*		project_model_generate(Bproject* project, Bmodel* temp, int flags)
{
	int					i;
	string				id("1");
	string				comptype("PRT");
	Breconstruction*	rec;
	Bparticle*			part;
	Bmodel*				model = NULL;
	Bmodel*				mp = NULL;
	Bmodel*				partmod = NULL;
	Bcomponent*			comp = NULL;
	Bcomptype*			ct = NULL;
	Bimage				*p;
	string				partfile, newpartfile, modfile;
	Vector3<double>		origin, shift;
	Matrix3				mat;
	
	if ( verbose ) 
		cout << "Generating models from particles:" << endl << "Reconstruction\tParticle\tModel" << endl;
	for ( i=1, rec = project->rec; rec; rec = rec->next, i++ ) {
		if ( mp ) mp = mp->add(rec->id.str());
		else model = mp = new Bmodel(rec->id.str());
		mp->FOM(rec->fom);
		mp->select(rec->select);
		mp->mapfile(rec->frec.str());
		mp->image_number(0);
//		mp->comment = project->comment;
		comp = NULL;
		for ( part = rec->part; part; part = part->next ) {
			if ( comp ) comp = comp->add(part->id);
			else mp->comp = comp = new Bcomponent(part->id);
			comp->location((part->loc - rec->origin) * rec->voxel_size);
			View2<double>	v(part->view2().backward());
			comp->view(View2<float>(v[0],v[1],v[2],v[3]));
			comp->FOM(part->fom[0]);
			comp->select(part->sel);
			ct = model->add_type(id);
			comp->type(ct);
			if ( part->fpart.length() ) {
				partfile = part->fpart.str();
//				modfile = partfile.pre_rev('.') + ".star";
//				newpartfile = partfile.pre_rev('.') + "_std." + partfile.post_rev('.');
				modfile = base(partfile) + ".star";
				newpartfile = base(partfile) + "_std." + extension(partfile);
			} else if ( rec->fpart.length() ) {
				partfile = rec->fpart.str();
//				modfile = partfile.pre_rev('.') + Bstring(part->id, "_%03d.star");
//				newpartfile = partfile.pre_rev('.') + Bstring(part->id, "_%03d.") + partfile.post_rev('.');
				modfile = base(partfile) + "_" + to_string(part->id) + ".star";
				newpartfile = base(partfile) + "_" + to_string(part->id) + "." + extension(partfile);
				ct->image_number(part->id - 1);
			} else {
//				modfile = rec->frec.pre_rev('.') + Bstring(part->id, "_%03d.star");
				modfile = rec->frec.str();
				modfile = base(modfile) + "_" + to_string(part->id) + ".star";
			}
			if ( partfile.length() ) {
				ct->file_name(newpartfile);
				if ( flags & 1 ) {
					if ( ( p = read_img(partfile, 1, ct->image_number()) ) == NULL ) {
						error_show("project_model_generate", __FILE__, __LINE__);
						return NULL;
					}
					ct->image_number(0);
					if ( p->background(long(0)) == 0 )
						p->calculate_background();
					origin = p->size()/2;
					shift = origin - part->ori;
					p->origin(part->ori);
					p->rotate(shift, part->view);
					part->view = View2<double>(0,0,1,0);
					part->ori = origin;
					write_img(ct->file_name(), p, 0);
					delete p;
				}
			}
			if ( flags & 2 ) {
				ct->file_name(modfile);
				ct->image_number(0);
				if ( temp ) {
//					partmod = model_copy(temp);
					partmod = temp->copy();
					partmod->identifier(comp->identifier());
//					origin = 0;
					mat = part->view.matrix();
					mat = mat.transpose();
					model_rotate(partmod, mat, origin, shift);
				} else partmod = new Bmodel(comp->identifier());
				partmod->FOM(part->fom[0]);
				partmod->select(part->sel);
				partmod->mapfile(partfile);
				if ( flags & 1 ) partmod->mapfile(newpartfile);
//				partmod->comment = mp->comment;
				write_model(modfile, partmod);
				model_kill(partmod);
			}
			if ( verbose )
				cout << rec->id << tab << part->id << tab << comp->identifier() << endl;
		}
	}
	
	return model;
}

/**
@brief 	Exchanges information between reconstruction and model structures.
@param 	*project	project parameter structure.
@param 	*model		model parameter structures.
@return int			0.

	This function assumes there is a 1-1 relationship between models and
	particles derived from reconstructions and that they are in corresponding
	sequence.
	The type id and hand from each model is converted into an integer and
	the particle group set to it.
	Each model id is reset to a combined reconstruction and particle id.

**/
int			project_model_consolidate(Bproject* project, Bmodel* model)
{
	Breconstruction*	rec;
	Bparticle*			part;
	Bmodel*				mp;
	Bstring				type_id;
	Vector3<double>		origin;
	
	if ( verbose )
		cout << "Consolidating particles with models:" << endl << "Reconstruction\tParticle\tModel" << endl;
	for ( rec = project->rec, mp = model; rec; rec = rec->next ) {
		for ( part = rec->part; part && mp; part = part->next, mp = mp->next ) {
			if ( verbose )
				cout << rec->id << tab << part->id << tab << mp->identifier() << endl;
//			type_id = mp->type_id.remove('_') + Bstring((mp->hand<0)?2:mp->hand, "%1d");
//			part->group = type_id.integer();
//			if ( part->group < 10 ) part->group = part->id;
//			mp->identifier() = rec->id + Bstring(part->id, "_%04d");
			if ( part->fpart.length() ) mp->mapfile(part->fpart.str());
//			model_shift(mp, (part->loc - rec->origin)*rec->voxel_size);
			model_rotate(mp, part->view2(), origin, (part->ori - rec->box_size/2)*rec->voxel_size);
		}
	}
	
	return 0;
}

