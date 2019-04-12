/* Main process of Panorama Image Stitching
*  Details: https://hypjudy.github.io/2017/05/10/panorama-image-stitching/
*  HYPJUDY, 2017.5.10
*/

/* Input: n unordered images
   Output: Panoramic image(s) */

#include "utils.h"
#include "stitching_functions.h"


int main() {
#ifdef __linux__
	string folder_path = "dataset3/";
#else
	string folder_path = "dataset3\\";
#endif
	const char *output_img_path = "output/flower.bmp";
	const int threshold = 30;
	CImg<unsigned char> output_img;
	vector<string> img_paths;
	get_files_in_folder(folder_path, img_paths);

	const int img_num = img_paths.size();
	vector<map<vector<float>, VlSiftKeypoint>> all_features(img_num);
	vector<CImg<unsigned char>> input_imgs(img_num);
	
	/* I. Extract SIFT features from all n images */
	for (int i = 0; i < img_num; ++i) {
		input_imgs[i] = CImg<unsigned char>(img_paths[i].c_str());
		input_imgs[i] = cylinderProjection(input_imgs[i]);
		//input_imgs[i].display();
		CImg<unsigned char> gray_img = rgb2gray(input_imgs[i]);
		
		all_features[i] = extract_sift_features(gray_img);
	}

	// indicate the adjacent relation
	vector<vector<int>> adjacent(img_num);
	vector<vector<bool>> adjacent_flag;
	// init
	for (int i = 0; i < img_num; ++i) {
		//input_imgs[i].display();
		vector<bool> v;
		for (int j = 0; j < img_num; ++j) {
			v.push_back(false);
		}
		adjacent_flag.push_back(v);
	}

	for (int i = 0; i < img_num; ++i) {
		for (int j = i + 1; j < img_num; ++j) {

			/* II. Find k nearest-neighbours for each feature
			       using a k-d tree */
			vector<KeypointPair> pairs = find_k_nearest_neighbours(
					all_features[i], all_features[j]);

			/* III. For each image: */
			/* (i) Select m candidate matching images that have
			       the most feature matches to this image */
			if (pairs.size() >= threshold) {
				adjacent[i].push_back(j);
				adjacent[j].push_back(i);
				adjacent_flag[i][j] = adjacent_flag[j][i] = true;
			}
		}
	}
	
	int last_idx = 0;
	queue<int> to_stitch_idxs;
	to_stitch_idxs.push(last_idx);
	CImg<unsigned char> stitched_img = input_imgs[last_idx];
	while (!to_stitch_idxs.empty()) {
		int to_stitch_idx = to_stitch_idxs.front();
		to_stitch_idxs.pop();

		for (int i = adjacent[to_stitch_idx].size() - 1; i >= 0; --i) {
			int adjacent_idx = adjacent[to_stitch_idx][i];
			if (!adjacent_flag[to_stitch_idx][adjacent_idx]) {
				continue;
			}
			else {
				adjacent_flag[to_stitch_idx][adjacent_idx] = false;
				adjacent_flag[adjacent_idx][to_stitch_idx] = false;
				to_stitch_idxs.push(adjacent_idx);
			}

			// feature will be shift when stitching, cannot be accessed 
			// from previous calculation (must calculate pairs again.
			vector<KeypointPair> feature_pairs = 
				find_k_nearest_neighbours(
					all_features[adjacent_idx], all_features[to_stitch_idx]);
			
			
			vector<KeypointPair> feature_pairs_inv;
			for (int j = 0; j < feature_pairs.size(); ++j) {
				feature_pairs_inv.push_back(
					KeypointPair(feature_pairs[j].p2, feature_pairs[j].p1));
			}

			/* (ii) Find geometrically consistent feature matches
			using RANSAC to solve for the homography between
			pairs of images */
			HomographyMatrix H = RANSAC_homography(feature_pairs);
			HomographyMatrix H_inv = RANSAC_homography(feature_pairs_inv);


			/* (iii) Verify image matches using a probabilistic model */
			// pass

			/* IV. Find connected components of image matches */
			CImg<unsigned char> adjacent_img(input_imgs[adjacent_idx]);
			float min_x = get_min_warped_x(adjacent_img, H_inv); // min <= 0
			float min_y = get_min_warped_y(adjacent_img, H_inv);
			float max_x = get_max_warped_x(adjacent_img, H_inv, stitched_img);
			float max_y = get_max_warped_y(adjacent_img, H_inv, stitched_img);
			int out_w = ceil(max_x - min_x);
			int out_h = ceil(max_y - min_y);
			CImg<unsigned char> last_stitch(out_w, out_h, 1, adjacent_img.spectrum(), 0);
			CImg<unsigned char> next_stitch(out_w, out_h, 1, adjacent_img.spectrum(), 0);

			img_shift(stitched_img, last_stitch, min_x, min_y);
			img_homography_warping(adjacent_img, next_stitch, H, min_x, min_y);
			feature_homography_warping(all_features[adjacent_idx], H_inv, min_x, min_y);
			feature_shift(all_features[last_idx], min_x, min_y);

			last_stitch.display();
			next_stitch.display();
			/* V. For each connected component: */
			/* (i) Perform bundle adjustment to solve for the rotation
			theta1, theta2, theta3 and focal length f of all cameras */
			// pass step (i)
			/* (ii) Render panorama using multi-band blending */
			stitched_img = multiband_blending(last_stitch, next_stitch);
			stitched_img.display();
			
		}
	}
	//stitched_img.display();
	stitched_img.save(output_img_path);

	return 0;
}
