/**
@file	img_combine.h
@brief	Header file for combining two images in various ways
@author 	Bernard Heymann
@date	Created: 20000430
@date	Modified: 20240406
**/

#include "rwimg.h"

// Function prototypes
Bimage* 	img_add(vector<string>& file_list, int flags);
Bimage*	 	img_setup_combined(vector<string>& file_list, long& nimg, int cat=0);
Bimage*		img_catenate(vector<string>& file_list, string& rawstring, DataType newdatatype,
				Vector3<long> nusize, int setZslices=0, int fill_type=0, double fill=0,
				double newavg=0, double newstd=0);
Bimage* 	img_add_weighed(vector<string>& file_list, vector<double> weight,
					double newavg=0, double newstd=0, int flags=0);
int			img_pack_into_image(Bimage* p, long nn, Bimage* ppatch, Vector3<long> start, Vector3<long> overlap);
Bimage*		img_patch(Bstring* file_list, Bstring tile_file, DataType nudatatype,
				double cutmin, double cutmax, double nuavg, double nustd);
