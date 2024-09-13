/**
@file	img_combine.cpp
@brief	Functions to combine two images in various ways
@author 	Bernard Heymann
@date	Created: 19990219
@date	Modified: 20240406
**/

#include "Bimage.h"
#include "img_combine.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Adds multiple images together with optional variance or standard deviation.
@param 	&file_list 	list of file names.
@param 	flags		flags to modify summation.
@return Bimage* 	resultant image (floating point).

	Images are read from a number files and added.
	All the images must be the same size, but could have different numbers of sub-images.
	The flags that can set are:
		1	calculate the average in stead of the sum.
		2	calculate the variance as FOM
		4	calculate the standard deviation as FOM (supercedes the variance)
	All images are converted to floating point.

**/
Bimage* 	img_add(vector<string>& file_list, int flags)
{
	if ( !file_list.size() ) return NULL;
	
	int				calcavg = flags & 1;
	int 			calcfom = (flags & 2) >> 1;
	long			nfiles(file_list.size());
	
	long			i, j, n(0), nf(0), nimg(0);
	double			v, va, w;
	Bimage*			p = NULL;

	Bimage*	 		psum = read_img(file_list[0], 1, 0);
	psum->clear();
	psum->change_type(Float);
	
	Bimage*			pfom = psum->next = psum->copy();
	
	long			imgsize = psum->sizeX()*psum->sizeY()*psum->sizeZ()*psum->channels();

	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		cout << endl << "Adding " << nfiles << " files together:" << endl;	
	
	for ( auto& filename: file_list ) {
		p = read_img(filename, 1, -1);
		if ( p != NULL ) {
			if ( verbose & VERB_LABEL )
				cout << "Adding file " << nf << " with images " << p->images() << endl;
			for ( n=j=0; n<p->images(); n++ ) {
				for ( i=n*imgsize; i<(n+1)*imgsize; i++, j++ ) {
					v = (*p)[j];
					psum->add(j, v);
					if ( calcfom ) pfom->add(j, v*v);
				}
			}
			nimg += p->images();
			delete p;
		}
		nf++;
	}

	w = 1.0/nimg;
	
	for ( i=n=0; n<psum->images(); n++ ) {
		for ( j=0; j<imgsize; j++, i++ ) {
			va = (*psum)[i] * w;
			if ( calcavg ) psum->set(i, va);
			else psum->set(i, (*psum)[i]);
			v = (*pfom)[i];
			if ( calcfom > 0 ) {
				v = v * w - va*va;
				if ( v < 0 ) v = 0;
			}
			if ( calcfom > 1 ) v = sqrt(v);
			pfom->set(i, v);
		}
	}
	
	psum->statistics();
	pfom->statistics();
	
	return psum;
}


/**
@brief 	Sets up a list of images for concatenation or summation.
@param 	&file_list		list of file names.
@param 	&nimg			number of concatenated images.
@param 	cat				flag to indicate concatenation.
@return Bimage*			new image into which to write data.

	The images can have different numbers of sub-images, sizes and data types.

**/
Bimage*	 	img_setup_combined(vector<string>& file_list, long& nimg, int cat)
{
	if ( !file_list.size() ) return NULL;
	
	long			nf(0), nn(0);
	Vector3<double>	sam(1,1,1);
	Bimage*			p;
	Bimage*			pc = new Bimage;
	
	nimg = 0;
	
	for ( auto& filename: file_list ) {
		p = read_img(filename, 0, -1);
		if ( p != NULL ) {
			if ( pc->channels() < p->channels() ) pc->channels(p->channels());
			if ( pc->data_type() < p->data_type() ) pc->data_type(p->data_type());
			if ( pc->compound_type() < p->compound_type() ) pc->compound_type(p->compound_type());
			pc->size(pc->size().max(p->size()));
			nimg += p->images();
			if ( nn < p->images() ) nn = p->images();
			if ( nf == 0 ) sam = p->image->sampling();
			delete p;
		}
		nf++;
	}

	if ( cat ) pc->images(nimg);
	else pc->images(nn);
	pc->origin(pc->default_origin());
	pc->sampling(sam);

	pc->data_alloc_and_clear();
	
	pc->information();
	
	return pc;
}



/**
@brief 	Catenates a list of images into a multi-image structure.
@param 	&file_list		list of file names.
@param 	&rawstring		format for re-interpretation of file.
@param 	nudatatype		new data type (default from first image).
@param 	nusize			new size (default from images).
@param 	setZslices		flag to create 2D images from slices.
@param 	fill_type		fill type for expanding images.
@param 	fill			fill value for expanding images.
@param 	newavg			new average to set each individual image.
@param 	newstd			new standard deviation to set each individual image.
@return Bimage*			catenated image.

	The images can have different numbers of sub-images, sizes and data types.

**/
/*Bimage*		img_catenate(vector<string>& file_list, string& rawstring, DataType nudatatype,
				Vector3<long> nusize, int setZslices, int fill_type, double fill,
				double newavg, double newstd)
{
	
	// Set up image structures
	Bimage*			p = NULL;
	Bimage*			pcat = NULL;
	int				err(0);
	long			i(0), j(0), k;

	if ( rawstring.length() > 0 && rawstring[0] != '#' ) rawstring = "#" + rawstring;
	
	// First read the file headers to identify problems
	for ( auto& name: file_list ) {
		string		filename = name + rawstring;
//		cout << filename << endl;
		p = read_img(filename, 0, -1);
		if ( p ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG img_catenate: " << p->file_name() << " " << 
						p->sizeX() << " " << p->sizeY() << " " << p->sizeZ() << " " << p->images() << " " << p->channels() << endl;
			if ( i == 0 ) {
				pcat = p;
				if ( fill_type == FILL_AVERAGE ) fill = p->average();
				if ( fill_type == FILL_BACKGROUND ) fill = p->background(long(0));
				if ( p->fourier_type() ) p->fourier_type(Standard);
				if ( setZslices ) i = p->sizeZ();
				else i = p->images();
			} else {
				if ( setZslices && p->images() > 1 ) {
					cerr << p->file_name() << " must be a single image file to pack slices" << endl;
					err += -1;
				}
//				err += img_compatibility(pcat, p);
				if ( !p->compatible(pcat) ) err += -1;
				if ( pcat->sizeX() < p->sizeX() ) pcat->sizeX(p->sizeX());
				if ( pcat->sizeY() < p->sizeY() ) pcat->sizeY(p->sizeY());
				if ( pcat->sizeZ() < p->sizeZ() ) pcat->sizeZ(p->sizeZ());
				if ( fill_type == FILL_AVERAGE ) fill += p->average();
				if ( fill_type == FILL_BACKGROUND ) fill += p->background(long(0));
				if ( setZslices ) i += p->sizeZ();
				else i += p->images();
				delete p;
			}
			j++;
		} else err += -1;
	}
	
	if ( err < 0 || j == 0 ) {
		cerr << "Files not concatenated! Number of errors: " << -err << endl;
		bexit(err);
	}
		
	// Set up for catenation
	if ( nusize.volume() > 0 ) pcat->size(nusize);
	nusize = pcat->size();
	if ( fill_type != FILL_USER ) fill /= j;
	if ( setZslices ) {
		pcat->sizeZ(i);
		pcat->images(1);
	} else {
		pcat->images(i);
	}		
	if ( nudatatype > Unknown_Type ) pcat->data_type(nudatatype);
	pcat->data_alloc();
	
	if ( verbose & VERB_FULL )
		cout << "Catenated size: n=" << pcat->images() << " size=" << pcat->size()
			<< " c=" << pcat->channels() << endl << endl;
	
	// Read the images including their data
	long			imgsize = (long) pcat->size().volume()*pcat->channels();
	Vector3<long>	translate;
	
	if ( setZslices ) imgsize = (long) pcat->sizeX()*pcat->sizeY()*pcat->channels();
	
	if ( verbose )
		cout << "File\tImages\tStart" << endl;
	i = 0;
	for ( auto& name: file_list ) {
		string		filename = name + rawstring;
		p = read_img(filename, 1, -1);
		if ( p != NULL ) {
			if ( verbose )
				cout << p->file_name() << tab << p->images()
					<< tab << i << endl;
			if ( setZslices ) {
				p->images_to_slices();
				nusize[2] = p->sizeZ();
			}
			if ( newstd ) p->rescale_to_avg_std(newavg, newstd);
			translate = (nusize - p->size())/2;
			p->resize(nusize, translate, fill_type, fill);
			p->change_type(nudatatype);
			if ( setZslices ) {
				for ( j=i*imgsize, k=0; k<imgsize*p->sizeZ(); j++, k++ ) pcat->set(j, (*p)[k]);
				i += p->sizeZ();
			} else {
				for ( j=i*imgsize, k=0; k<imgsize*p->images(); j++, k++ ) pcat->set(j, (*p)[k]);
				for ( j=0; j<p->images(); j++, i++ ) pcat->image[i] = p->image[j];
			}
			delete p;
		}
	}
	
	if ( verbose )
		cout << "Number of images:               " << pcat->images() << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_catenate: final_size=" << i*imgsize << endl;
	
	return pcat;
}
*/
Bimage*		img_catenate(vector<string>& file_list, string& rawstring, DataType nudatatype,
				Vector3<long> nusize, int setZslices, int fill_type, double fill,
				double newavg, double newstd)
{
	
	// Set up image structures
	Bimage*			p = NULL;
	Bimage*			pcat = NULL;
	int				err(0);
	long			i(0), j(0), k;

	if ( rawstring.length() > 0 && rawstring[0] != '#' ) rawstring = "#" + rawstring;
	
	// First read the file headers to identify problems
	for ( auto name: file_list ) {
		string		filename = name + rawstring;
//		cout << filename << endl;
		p = read_img(filename, 0, -1);
		if ( p ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG img_catenate: " << p->file_name() << " " <<
						p->sizeX() << " " << p->sizeY() << " " << p->sizeZ() << " " << p->images() << " " << p->channels() << endl;
			if ( i == 0 ) {
				pcat = p;
				if ( fill_type == FILL_AVERAGE ) fill = p->average();
				if ( fill_type == FILL_BACKGROUND ) fill = p->background(long(0));
				if ( p->fourier_type() ) p->fourier_type(Standard);
				if ( setZslices ) i = p->sizeZ();
				else i = p->images();
			} else {
				if ( setZslices && p->images() > 1 ) {
					cerr << p->file_name() << " must be a single image file to pack slices" << endl;
					err += -1;
				}
//				err += img_compatibility(pcat, p);
				if ( !p->compatible(pcat) ) err += -1;
				if ( pcat->sizeX() < p->sizeX() ) pcat->sizeX(p->sizeX());
				if ( pcat->sizeY() < p->sizeY() ) pcat->sizeY(p->sizeY());
				if ( pcat->sizeZ() < p->sizeZ() ) pcat->sizeZ(p->sizeZ());
				if ( fill_type == FILL_AVERAGE ) fill += p->average();
				if ( fill_type == FILL_BACKGROUND ) fill += p->background(long(0));
				if ( setZslices ) i += p->sizeZ();
				else i += p->images();
				delete p;
			}
			j++;
		} else err += -1;
	}
	
	if ( err < 0 || j == 0 ) {
		cerr << "Files not concatenated! Number of errors: " << -err << endl;
		bexit(err);
	}
		
	// Set up for catenation
	if ( nusize.volume() > 0 ) pcat->size(nusize);
	nusize = pcat->size();
	if ( fill_type != FILL_USER ) fill /= j;
	if ( setZslices ) {
		pcat->sizeZ(i);
		pcat->images(1);
	} else {
		pcat->images(i);
	}
	if ( nudatatype > Unknown_Type ) pcat->data_type(nudatatype);
	pcat->data_alloc();
	
	if ( verbose & VERB_FULL )
		cout << "Catenated size: n=" << pcat->images() << " size=" << pcat->size()
			<< " c=" << pcat->channels() << endl << endl;
	
	// Read the images including their data
	long			imgsize = (long) pcat->size().volume()*pcat->channels();
	Vector3<long>	translate;
	
	if ( setZslices ) imgsize = (long) pcat->sizeX()*pcat->sizeY()*pcat->channels();
	
	if ( verbose )
		cout << "File\tImages\tStart" << endl;
	i = 0;
	for ( auto name: file_list ) {
		string		filename = name + rawstring;
		p = read_img(filename, 1, -1);
		if ( p != NULL ) {
			if ( verbose )
				cout << p->file_name() << tab << p->images()
					<< tab << i << endl;
			if ( setZslices ) {
				p->images_to_slices();
				nusize[2] = p->sizeZ();
			}
			if ( newstd ) p->rescale_to_avg_std(newavg, newstd);
			translate = (nusize - p->size())/2;
			p->resize(nusize, translate, fill_type, fill);
			p->change_type(nudatatype);
			if ( setZslices ) {
				for ( j=i*imgsize, k=0; k<imgsize*p->sizeZ(); j++, k++ ) pcat->set(j, (*p)[k]);
				i += p->sizeZ();
			} else {
				for ( j=i*imgsize, k=0; k<imgsize*p->images(); j++, k++ ) pcat->set(j, (*p)[k]);
				for ( j=0; j<p->images(); j++, i++ ) pcat->image[i] = p->image[j];
			}
			delete p;
		}
	}
	
	if ( verbose )
		cout << "Number of images:               " << pcat->images() << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_catenate: final_size=" << i*imgsize << endl;
	
	return pcat;
}
/**
@brief 	Adds multiple images together with given weights.
@param 	&file_list 	list of file names.
@param 	weight		list of weights (can be NULL).
@param 	newavg		new average for rescaling.
@param 	newstd		new standard deviation for rescaling.
@param 	flags		flags to modify summation.
@return Bimage* 	resultant image (floating point).

	Images are read from a number files and added to each other, using
	the given weights to determine each contribution.
	The images are rescaled to a new average and standard deviation before 
	weighted addition. If the given standard deviation is zero or less,
	this step is omitted.
	The weighed average is calculated and returned as a new image.
	The flags that can set are:
		1	calculate the average in stead of the sum.
		2	calculate the variance as FOM
		4	calculate the standard deviation as FOM (supercedes the variance)
		8	center each image before summation
	All images are converted to floating point.

**/
Bimage* 	img_add_weighed(vector<string>& file_list, vector<double> weight,
					double newavg, double newstd, int flags)
{
	if ( !file_list.size() ) return NULL;
	
	int					calcavg = flags & 1;
	int 				calcfom = (flags & 2) >> 1;
	int					center = flags & 8;
	
	long				i, j, n(0), nf(0), nimg(0), c(0);
	double				v;
	Bimage*				p = NULL;

	Bimage*	 			psum = img_setup_combined(file_list, nimg, 0);
	long				nfiles(file_list.size());

	double				weightsum(0);
	if ( weight.size() < nfiles )
		for ( i=weight.size(); i<nfiles; i++ ) weight.push_back(1);

	for ( i=0; i<nfiles; i++ ) weightsum += weight[i];;

	double				degrees_of_freedom = weightsum*(1 - 1.0L/nfiles);
	
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Adding " << nimg << " images together:" << endl;
		cout << "New image size:                 " << psum->size() << endl;
		cout << "Degrees of freedom:             " << degrees_of_freedom << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << endl << "Adding " << nimg << " images together" << endl << endl;
	
	long			imgsize = psum->sizeX()*psum->sizeY()*psum->sizeZ();
	long			datasize = imgsize*psum->images();
	
	float*			fom = NULL;
	if ( calcfom ) {
		psum->next = new Bimage(Float, TSimple, psum->size(), psum->images());
		psum->next->sampling(psum->image->sampling());
		psum->next->origin(psum->image->origin());
		fom = (float *) psum->next->data_pointer();
	}

	vector<double>	imgweight(psum->images());
	
	for ( auto& filename: file_list ) {
		p = read_img(filename, 1, -1);
		if ( p != NULL ) {
			if ( newstd > 0 ) p->rescale_to_avg_std(newavg, newstd);
			if ( center ) p->center_wrap();
			if ( verbose & VERB_LABEL )
				cout << "Adding image " << nf << " with weight " << weight[nf] << endl << endl;
			for ( n=j=0; n<p->images(); n++ ) {
				imgweight[n] += weight[nf];
				for ( i=n*imgsize; i<(n+1)*imgsize; i++ ) {
					for ( c=0; c<psum->channels(); c++, j++ ) {
						v = (*p)[j];
						psum->add(j, weight[nf]*v);
						if ( calcfom ) fom[i] += weight[nf]*v*v;
					}
				}
			}
			delete p;
		}
		nf++;
	}
	
	if ( calcfom ) {
		for ( i=n=0; n<psum->images(); n++ )
			for ( j=0; j<imgsize; j++, i++ )
				fom[i] = (fom[i] - (*psum)[i]*(*psum)[i]/imgweight[n])/degrees_of_freedom;
		if ( calcfom > 1 ) {
			for ( i=0; i<datasize; i++ ) {
				if ( fom[i] > 0 ) fom[i] = sqrt(fom[i]);
				else fom[i] = 0;
			}
		}
	}
	
	if ( calcavg )
		for ( i=n=0; n<psum->images(); n++ )
			for ( j=0; j<imgsize; j++, i++ )
				psum->set(i, (*psum)[i] / imgweight[n]);

	psum->statistics();
	
	return psum;
}

int			img_pack_into_image(Bimage* p, long nn, Bimage* ppatch, Vector3<long> start, Vector3<long> overlap)
{
	long   			i, j, c, x, y, z, xx, yy, zz;
	Vector3<double>	w;
	
	overlap *= 0.5;
	
	if ( verbose & VERB_FULL )
		cout << "Placing tile " << nn << " of size " << p->size() << " at " << start << endl;
	
	for ( z=0; z<p->sizeZ(); z++ ) {
		zz = start[2] + z;
		w[2] = 0;
		if ( z >= overlap[2] && z < p->sizeZ() - overlap[2] ) w[2] = 1;
		if ( zz <= overlap[2] ) w[2] = 1;
		if ( zz >= ppatch->sizeZ() - overlap[2] ) w[2] = 1;
		for ( y=0; y<p->sizeY(); y++ ) {
			yy = start[1] + y;
			w[1] = 0;
			if ( y >= overlap[1] && y < p->sizeY() - overlap[1] ) w[1] = 1;
			if ( yy <= overlap[1] ) w[1] = 1;
			if ( yy >= ppatch->sizeY() - overlap[1] ) w[1] = 1;
			for ( x=0; x<p->sizeX(); x++ ) {
				xx = start[0] + x;
				w[0] = 0;
				if ( x >= overlap[0] && x < p->sizeX() - overlap[0] ) w[0] = 1;
				if ( xx <= overlap[0] ) w[0] = 1;
				if ( xx >= ppatch->sizeX() - overlap[0] ) w[0] = 1;
				if ( w.volume() ) {
//					i = ((z*p->sizeY() + y)*p->sizeX() + x)*p->channels();
//					j = ((zz*ppatch->sizeY() + yy)*ppatch->sizeX() + xx)*p->channels();
					i = p->index(0, x, y, z, nn);
					j = ppatch->index(0, xx, yy, zz, 0);
					for ( c=0; c<p->channels(); c++, i++, j++ ) ppatch->set(j, (*p)[i]);
				}
			}
		}
	}

	return(0);
}


/**
@brief 	Assembles overlapping tiles into a single image.
@param 	*file_list		list of files with tiles.
@param 	tile_file		text file specifying the tile size, overlap and locations.
@param 	nudatatype		new data type to convert to.
@param 	cutmin			minimum to truncate each tile before assembly.
@param 	cutmax			maximum to truncate each tile before assembly.
@param 	nuavg			new average to set each tile before assembly.
@param 	nustd			new standard deviation to set each tile before assembly.
@return Bimage*			new composite image.

	Each tile extents halfway into overlapped areas, representing a step
	function between tiles.

**/
Bimage*		img_patch(Bstring* file_list, Bstring tile_file, DataType nudatatype,
				double cutmin, double cutmax, double nuavg, double nustd)
{
	long				i, n(1);
	Bstring*			file_ptr = NULL;
	Vector3<long>		size;						// Overall size after patching
	Vector3<long>		overlap;					// Tile overlap
	
	// Get the number of tiles and find the last tile
	for ( file_ptr = file_list; file_ptr->next; file_ptr = file_ptr->next ) n++;

	// Read the tile size, overlap and locations
	ifstream		fd(tile_file.c_str());
	if ( fd.fail() ) {
		cerr << "Error: A text file specifying the size and tile origins must be given!" << endl;
		bexit(-1);
	}
	
	fd >> size[0] >> size[1] >> size[2];
	fd >> overlap[0] >> overlap[1] >> overlap[2];

//	Vector3<int>*	tile_origin = new Vector3<int>[n];
	vector<Vector3<long>>	tile_origin(n);
	
	for ( i=0; i<n; i++ )
		fd >> tile_origin[i][0] >> tile_origin[i][1] >> tile_origin[i][2];
		
	fd.close();

	// Get the overall size
	Bimage*			p = read_img(*file_list, 0, 0);
	
	Vector3<long>	tile_size(p->size());			// Tile size from first image
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG img_patch: tilesize=" << tile_size << " overlap=" << overlap << endl;
		cout << "DEBUG img_patch: tiles=" << n << endl;
	}
	
	delete p;
	
	// Get the overall size
	p = read_img(*file_ptr, 0, 0);
	
	if ( size != tile_origin[n-1] + p->size() ) {
		cerr << "Error: The expected size does not match the size given in file " << tile_file << endl;
		cerr << "	Expected size:  " << tile_origin[n-1] + p->size() << endl;
		cerr << "	Given size:     " << size << endl;
		cerr << "	Make sure all the tiles are present." << endl;
		bexit(-1);
	}
	
	// Set up the receiving image
	if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();
	Bimage*			ppatch = new Bimage(nudatatype, p->compound_type(), size, 1);
	ppatch->sampling(p->sampling(0));
	ppatch->image[0] = p->image[0];
	ppatch->origin(ppatch->default_origin());
	delete p;

	if ( verbose ) {
		cout << "Patching tiles:" << endl;
		cout << "Size:                           " << size << endl;
		cout << "Overlap:                        " << overlap << endl;
		cout << "Number of tiles:                " << n << endl << endl;
	}
	
	for ( i=0, file_ptr = file_list; i<n && file_ptr; i++, file_ptr = file_ptr->next ) {
		if ( ( p = read_img(*file_ptr, 1, 0) ) ) {
			if ( verbose & VERB_DEBUG )
				cout << p->file_name() << ": " << p->size() << " " << p->channels() << " " << p->images() << endl;
			if ( nustd ) p->rescale_to_avg_std(nuavg, nustd);
			if ( cutmin || cutmax )
				p->truncate_to_min_max(cutmin, cutmax);
			p->change_type(nudatatype);
//			ppatch->place_with_overlap(p, i, tile_size, n, tile_origin);
			img_pack_into_image(p, 0, ppatch, tile_origin[i], overlap);
			delete p;
		}
	}
	
	ppatch->calculate_background();

//	delete[] tile_origin;

	return ppatch;
}


