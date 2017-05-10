/* Main functions of Panorama Image Stitching
*  Details: https://hypjudy.github.io/2017/05/10/panorama-image-stitching/
*  HYPJUDY, 2017.5.10
*/

#pragma once

/* I. Extract SIFT features from all n images */
/* Input: grayscale image
Return: features
map of descriptor and keypoints
*/
map<vector<float>, VlSiftKeypoint> extract_sift_features(
	const CImg<unsigned char> &gray_img) {
	const float scale_ratio = 0.8;
	map<vector<float>, VlSiftKeypoint> features;

	CImg<unsigned char> img(gray_img);
	img.resize(1.0 * img.width() / scale_ratio,  // enlarge
		1.0 * img.height() / scale_ratio,
		img.depth(), img.spectrum(), 3);
	const int w = img.width();
	const int h = img.height();
	const int dim = 128;

	/* Constructing a scale space */
	vl_sift_pix *img_data = new vl_sift_pix[w * h];
	cimg_forXY(img, x, y) { // x: width, y: height
		img_data[x + y * w] = img(x, y, 0); // one channel
	}

	// Create a new SIFT filter
	int noctaves = 4; // number of octaves. e.g. 4
	int nlevels = 5; // number of levels per octave. e.g. 5
	int o_min = 0; // first octave index
	VlSiftFilt* sift_filter = vl_sift_new(w, h, noctaves, nlevels, o_min);

	/* LoG Approximation */
	// starts processing a new image by computing its 
	// Gaussian scale space at the lower octave
	if (vl_sift_process_first_octave(sift_filter, img_data)\
		!= VL_ERR_EOF) { // has more octaves to process
		while (true) {
			/* Finding keypoints */
			/* Get rid of bad key points */
			// detect keypoints in the current octave 
			// filling the internal keypoint buffer
			vl_sift_detect(sift_filter);
			VlSiftKeypoint const *keypoints_list_pointer = \
				vl_sift_get_keypoints(sift_filter);

			/* Assigning an orientation to the keypoints */
			for (int i = 0; i < sift_filter->nkeys; ++i) {
				double angles[4];
				VlSiftKeypoint keypoint_pointer = *keypoints_list_pointer;
				keypoints_list_pointer++;
				int norientations = vl_sift_calc_keypoint_orientations(
					sift_filter, angles, &keypoint_pointer);

				/* Generate SIFT features */
				for (int j = 0; j < norientations; ++j) {
					vl_sift_pix descr[dim]; // SIFT descriptor (output)
											// computes the SIFT descriptor
					vl_sift_calc_keypoint_descriptor(
						sift_filter, descr, &keypoint_pointer, angles[j]);
					vector<float> descriptor(dim);
					for (int k = 0; k < dim; ++k) descriptor[k] = descr[k];

					keypoint_pointer.x *= scale_ratio;
					keypoint_pointer.y *= scale_ratio;
					keypoint_pointer.ix *= scale_ratio;
					keypoint_pointer.iy *= scale_ratio;

					features.insert(pair<vector<float>, VlSiftKeypoint>(
						descriptor, keypoint_pointer));
				}
			}
			if (vl_sift_process_next_octave(sift_filter) == VL_ERR_EOF) break;
		}
	}
	vl_sift_delete(sift_filter);
	delete img_data;
	img_data = NULL;
	return features;
}



/* II. Find 2 nearest-neighbours for each feature of query
image in searched image using a k-d tree.
If the best match much better than the next, accept. */
/* Input:
query_fea: all features of a query image
searched_fea: all features of a searched image
Output:
Matched features of the query image and searched image
*/
const int nneighbors = 2; // number of nearest neighbors to find
vector<KeypointPair> find_k_nearest_neighbours(
	const map<vector<float>, VlSiftKeypoint> &query_fea,
	const map<vector<float>, VlSiftKeypoint> &searched_fea) {

	vector<KeypointPair> res;
	const int len_s = searched_fea.size();

	const int dim = 128; // data dimension
	const int num_trees = 1;
	VlKDForest *kdforest = vl_kdforest_new(
		VL_TYPE_FLOAT, dim, num_trees, VlDistanceL1);

	// data: store each dimension of each feature descriptor
	//       in a searched image
	float* data = new float[dim * len_s];
	map<vector<float>, VlSiftKeypoint>::const_iterator it =
		searched_fea.begin();
	int i;
	for (it, i = 0; it != searched_fea.end(); ++it, ++i) {
		// each feature has dim(128) dimensions
		vector<float> descriptor = it->first;
		for (int j = 0; j < dim; ++j) {
			data[j + dim * i] = descriptor[j];
		}
	}

	// builds the KDTree for searched image by processing the data
	vl_kdforest_build(kdforest, len_s, data);

	VlKDForestSearcher* searcher = vl_kdforest_new_searcher(kdforest);

	// list of nearest neighbors found (output)
	VlKDForestNeighbor neighbors[nneighbors];
	it = query_fea.begin();
	for (it; it != query_fea.end(); ++it) {
		float query_point[128];
		for (int i = 0; i < dim; ++i) query_point[i] = (it->first)[i];
		// Query the forest
		vl_kdforestsearcher_query(searcher, neighbors, nneighbors, query_point);

		// VlKDForestNeighbor Struct Reference:
		// double distance: distance to the query point
		// vl_uindex index: index of the neighbor in the KDTree data
		const double ratio = 0.5;
		double NN_1 = neighbors[0].distance; // SSD of the closest match
		double NN_2 = neighbors[1].distance; // SSD of the second-closest match
											 //  Feature-space outlier rejection: SSD(patch1,patch2) < threshold
											 //  If the best match much better than the next, accept.
		if (NN_1 / NN_2 < ratio) {
			vector<float> descr(dim);
			for (int i = 0; i < dim; ++i)
				descr[i] = data[i + neighbors[0].index * dim];
			VlSiftKeypoint p1 = searched_fea.find(descr)->second;
			VlSiftKeypoint p2 = it->second;
			res.push_back(KeypointPair(p1, p2));
		}
	}
	vl_kdforestsearcher_delete(searcher);
	vl_kdforest_delete(kdforest);
	delete[]data;
	data = NULL;
	return res;
}


/* III. For each image: */
/* (i) Select m candidate matching images that have
the most feature matches to this image */



/* (ii) Find geometrically consistent feature matches
using RANSAC to solve for the homography between
pairs of images */
HomographyMatrix RANSAC_homography(const vector<KeypointPair> & pairs) {
	/* Calculate N for probability p of at least
	one sample with only inliers */
	int s = 4; // number of points to compute solution
	float p = 0.99; // chance of getting good sample
	float e = 0.5; // proportion outliers (outlier ratio)
				   /* P(sample set with all inliers) = (1 − e)^s
				   P(sample set will have at least one outlier) = (1-(1-e)^s)
				   P(all N samples have outlier) = (1-(1-e)^s)^N
				   We want P(all N samples have outlier) < (1-p)
				   So: (1-(1-e)^s)^N < (1-p) */
	int N = ceil(log(1 - p) / log(1 - pow(1 - e, s)));

	vector<int> inliers_idx, max_inliers_idx;
	/* RANSAC loop:*/
	while (N--) {
		srand(time(0));
		/* 1. Select four feature pairs (at random) */
		vector<KeypointPair> feature_pairs;
		set<int> idxs;
		while (feature_pairs.size() < s) {
			while (true) {
				int idx = rand() % pairs.size();
				if (idxs.find(idx) == idxs.end()) { // no duplicate
					idxs.insert(idx);
					feature_pairs.push_back(pairs[idx]);
					break;
				}
			}
		}

		/* 2. Compute homography H_k (exact): p1=Hp2 */
		HomographyMatrix H = get_homography_matrix(feature_pairs);

		/* 3. Compute inliers where SSD(p_i', H_k.p_i) < eps */
		float threshDist = 4.0;
		inliers_idx.clear();
		for (int i = 0; i < pairs.size(); ++i) {
			if (idxs.find(i) != idxs.end()) continue; // exclude original inliers

			float p1_x = pairs[i].p1.x;
			float p1_y = pairs[i].p1.y;
			float p2_x = pairs[i].p2.x;
			float p2_y = pairs[i].p2.y;
			float warped_x = get_warped_x(p1_x, p1_y, H);
			float warped_y = get_warped_y(p1_x, p1_y, H);

			float dist = sqrt((warped_x - p2_x)*(warped_x - p2_x) +
				(warped_y - p2_y)*(warped_y - p2_y));
			if (dist < threshDist) inliers_idx.push_back(i);
		}

		/* 4. Keep H_k, if C_k is the largest set of inliers */
		if (inliers_idx.size() > max_inliers_idx.size())
			max_inliers_idx = inliers_idx;

	} /* 5. For a while go to 1 */

	  /* 6. Re-compute least-squares H estimate on all of the C_k inliers */
	int inliers_num = max_inliers_idx.size();
	assert(inliers_num > 0);
	CImg<double> A(4, inliers_num, 1, 1, 0);
	CImg<double> bx(1, inliers_num, 1, 1, 0);
	CImg<double> by(1, inliers_num, 1, 1, 0);
	for (int i = 0; i < inliers_num; ++i) {
		int idx = max_inliers_idx[i];
		A(0, i) = pairs[idx].p1.x;
		A(1, i) = pairs[idx].p1.y;
		A(2, i) = pairs[idx].p1.x * pairs[idx].p1.y;
		A(3, i) = 1;

		bx(0, i) = pairs[idx].p2.x;
		by(0, i) = pairs[idx].p2.y;
	}
	CImg<double> x1 = bx.get_solve(A); // solve Ax=B
	CImg<double> x2 = by.get_solve(A);

	return HomographyMatrix(x1(0, 0), x1(0, 1), x1(0, 2),
		x1(0, 3), x2(0, 0), x2(0, 1), x2(0, 2), x2(0, 3));
}


/* (iii) Verify image matches using a probabilistic model */


/* IV. Find connected components of image matches */
void img_homography_warping(const CImg<unsigned char> &src,
	CImg<unsigned char> &dst, HomographyMatrix H,
	float offset_x, float offset_y) {
	// inverse warping
	cimg_forXY(dst, x, y) {
		float warped_x = get_warped_x(x + offset_x, y + offset_y, H);
		float warped_y = get_warped_y(x + offset_x, y + offset_y, H);
		if (warped_x >= 0 && warped_x < src.width() &&
			warped_y >= 0 && warped_y < src.height()) {
			cimg_forC(dst, c) {
				dst(x, y, c) = src(floor(warped_x), floor(warped_y), c);
			}
		}
	}
}

void img_shift(const CImg<unsigned char> &src,
	CImg<unsigned char> &dst, float offset_x, float offset_y) {
	cimg_forXY(dst, x, y) {
		int x0 = x + offset_x;
		int y0 = y + offset_y;
		if (x0 >= 0 && x0 < src.width() &&
			y0 >= 0 && y0 < src.height()) {
			cimg_forC(dst, c) {
				dst(x, y, c) = src(x0, y0, c);
			}
		}
	}
}

void feature_homography_warping(map<vector<float>, VlSiftKeypoint> &feature,
	HomographyMatrix H, float offset_x, float offset_y) {
	map<vector<float>, VlSiftKeypoint>::iterator it = feature.begin();
	for (it; it != feature.end(); ++it) {
		float px = it->second.x; // coordinate
		float py = it->second.y;
		it->second.x = get_warped_x(px, py, H) - offset_x;
		it->second.y = get_warped_y(px, py, H) - offset_y;
		it->second.ix = int(it->second.x); // unormolized coordinate
		it->second.iy = int(it->second.y);
	}
}

void feature_shift(map<vector<float>, VlSiftKeypoint> &feature,
	float offset_x, float offset_y) {
	map<vector<float>, VlSiftKeypoint>::iterator it = feature.begin();
	for (it; it != feature.end(); ++it) {
		it->second.x -= offset_x; // coordinate
		it->second.y -= offset_y;
		it->second.ix = int(it->second.x); // unormolized coordinate
		it->second.iy = int(it->second.y);
	}
}

/* V. For each connected component: */
/* (i) Perform bundle adjustment to solve for the rotation
theta1, theta2, theta3 and focal length f of all cameras */
// pass step (i)
/* (ii) Render panorama using multi-band blending */
CImg<unsigned char> multiband_blending(const CImg<unsigned char> &a, const CImg<unsigned char> &b) {
	/* Find overlapped area */
	int w = a.width(), h = a.height(); // a and b has the same size

	int sum_a_x = 0, sum_a_y = 0;
	int width_mid_a = 0; // width of image a's content (calculate in half of height)

	int sum_overlap_x = 0, sum_overlap_y = 0;
	int width_mid_overlap = 0;

	// only consider stitching image horizontally
	int mid_y = h / 2;
	int x = 0;
	while (a(x, mid_y) == 0) ++x; // avoid leading zero
	for (x; x < w; ++x) {
		if (a(x, mid_y) != 0) { // black
			sum_a_x += x;
			++width_mid_a;
			if (b(x, mid_y) != 0) {
				sum_overlap_x += x;
				++width_mid_overlap;
			}
		}
	}

	/* 0. Forming a Gaussian Pyramid */
	/* i. Start with the original image G0 */
	int max_len = w >= h ? w : h;
	int level_num = floor(log2(max_len));

	vector<CImg<float>> a_pyramid(level_num);
	vector<CImg<float>> b_pyramid(level_num);
	vector<CImg<float>> mask(level_num);

	a_pyramid[0] = a;
	b_pyramid[0] = b;
	mask[0] = CImg<float>(w, h, 1, 1, 0);
	assert(width_mid_a > 0);
	assert(width_mid_overlap > 0);
	float ratio = 1.0 * sum_a_x / width_mid_a;
	float overlap_ratio = 1.0 * sum_overlap_x / width_mid_overlap;
	// the x=overlap_ratio line should lie in the overlap area
	if (ratio < overlap_ratio) {
		for (int x = 0; x < overlap_ratio; ++x)
			for (int y = 0; y < h; ++y)
				mask[0](x, y) = 1;
	}
	else {
		for (int x = overlap_ratio + 1; x < w; ++x)
			for (int y = 0; y < h; ++y)
				mask[0](x, y) = 1;
	}


	/* ii. Perform a local Gaussian weighted averaging
	function in a neighborhood about each pixel,
	sampling so that the result is a reduced image
	of half the size in each dimension. */
	for (int i = 1; i < level_num; ++i) {
		int wp = a_pyramid[i - 1].width() / 2;
		int hp = a_pyramid[i - 1].height() / 2;
		int sp = a_pyramid[i - 1].spectrum();
		a_pyramid[i] = a_pyramid[i - 1].get_blur(2, true, true).
			get_resize(wp, hp, 1, sp, 3);
		b_pyramid[i] = b_pyramid[i - 1].get_blur(2, true, true).
			get_resize(wp, hp, 1, sp, 3);
		mask[i] = mask[i - 1].get_blur(2, true, true).
			get_resize(wp, hp, 1, sp, 3);
	} /* iii. Do this all the way up the pyramid Gl = REDUCE(Gl-1) */

	  /* iiii. Each level l node will represent a weighted
	  average of a subarray of level l. */



	  /* 1. Compute Laplacian pyramid of images and mask */
	  /* Making the Laplacians Li=Gi-expand(Gi+1)*/
	  /*  subtract each level of the pyramid from the next lower one
	  EXPAND:  interpolate new samples between those of
	  a given image to make it big enough to subtract*/
	for (int i = 0; i < level_num - 1; ++i) {
		int wp = a_pyramid[i].width();
		int hp = a_pyramid[i].height();
		int sp = a_pyramid[i].spectrum();
		a_pyramid[i] -= a_pyramid[i + 1].get_resize(wp, hp, 1, sp, 3);
		b_pyramid[i] -= b_pyramid[i + 1].get_resize(wp, hp, 1, sp, 3);
	}


	/* 2. Create blended image at each level of pyramid */
	/* Forming the New Pyramid
	A third Laplacian pyramid LS is constructed by copying
	nodes from the left half of LA to the corresponding
	nodes of LS and nodes from the right half of LB to the
	right half of LS.
	Nodes along the center line are set equal to
	the average of corresponding LA and LB nodes */
	vector<CImg<float>> blend_pyramid(level_num);
	for (int i = 0; i < level_num; ++i) {
		blend_pyramid[i] = CImg<float>(a_pyramid[i].width(),
			a_pyramid[i].height(), 1, a_pyramid[i].spectrum(), 0);
		cimg_forXYC(blend_pyramid[i], x, y, c) {
			blend_pyramid[i](x, y, c) =
				a_pyramid[i](x, y, c) * mask[i](x, y) +
				b_pyramid[i](x, y, c) * (1.0 - mask[i](x, y));
		}
	}

	/* 3. Reconstruct complete image */
	/* Using the new Laplacian Pyramid
	Use the new Laplacian pyramid with the reverse of how it
	was created to create a Gaussian pyramid. Gi=Li+expand(Gi+1)
	The lowest level of the new Gaussian pyramid gives the final
	result. */
	// float: cannot be unsigned char(invalid) type!
	CImg<float> expand = blend_pyramid[level_num - 1];
	for (int i = level_num - 2; i >= 0; --i) {
		expand.resize(blend_pyramid[i].width(),
			blend_pyramid[i].height(), 1, blend_pyramid[i].spectrum(), 3);
		cimg_forXYC(blend_pyramid[i], x, y, c) {
			expand(x, y, c) = blend_pyramid[i](x, y, c) + expand(x, y, c);
			if (expand(x, y, c) > 255) expand(x, y, c) = 255;
			else if (expand(x, y, c) < 0) expand(x, y, c) = 0;
		}
	}
	return expand;
}
