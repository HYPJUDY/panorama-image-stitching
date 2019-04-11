/* Utils of Panorama Image Stitching
*  Details: https://hypjudy.github.io/2017/05/10/panorama-image-stitching/
*  HYPJUDY, 2017.5.10
*/

#pragma once

#include <iostream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <queue>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#ifdef __linux__
#include <dirent.h>
#include <dlfcn.h>
#else
#include <io.h>
#endif
#include <math.h>

#include "vl/generic.h"
#include "vl/sift.h"
#include "vl/kdtree.h"

#include "CImg.h"

using namespace cimg_library;
using namespace std;

struct KeypointPair {
	VlSiftKeypoint p1;
	VlSiftKeypoint p2;
	KeypointPair(VlSiftKeypoint _p1, VlSiftKeypoint _p2) :
		p1(_p1), p2(_p2) {}
};

struct HomographyMatrix {
	float a, b, c, d, e, f, g, h;
	HomographyMatrix(float _a, float _b, float _c,
		float _d, float _e, float _f, float _g, float _h) :
		a(_a), b(_b), c(_c), d(_d), e(_e), f(_f), g(_g), h(_h) {}
};

/* Get all files in a folder specified by path and store the file names in a vector */
void get_files_in_folder(string folder_path, vector<string> &file_paths) {
	string p;
#ifdef __linux__
    DIR * dir = opendir(folder_path.c_str());
    if (dir != NULL)
    {
        struct dirent * ent;
        while ((ent = readdir(dir)) != NULL)
        {
			/* Recuresively open files ? */
			if (ent->d_type == 4 && strcmp(ent->d_name, ".") != 0 && strcmp(ent->d_name, "..") != 0)
				get_files_in_folder(p.assign(folder_path).append("/").append(ent->d_name), file_paths);
			else if (ent->d_type == 8)
			{
				file_paths.push_back(p.assign(folder_path).append("/").append(ent->d_name));
				printf("adding %s\n", ent->d_name);
			}
        }
        closedir(dir);
    }
#else
	intptr_t hFile = 0;
	struct _finddata_t fileinfo;
	if ((hFile = _findfirst(p.assign(folder_path).append("\\*").c_str(), &fileinfo)) != -1) {
		do {
			if ((fileinfo.attrib & _A_SUBDIR)) {
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
				{
					file_paths.push_back(p.assign(folder_path).append("\\").append(fileinfo.name));
					get_files_in_folder(p.assign(folder_path).append("\\").append(fileinfo.name), file_paths);
				}
			}
			else {
				file_paths.push_back(p.assign(folder_path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
#endif
}

/* RGB to grayscale transformation */
CImg<unsigned char> rgb2gray(CImg<unsigned char> rgb_img) {
	CImg<unsigned char> gray_img(rgb_img._width, rgb_img._height, 1, 1, 0);
	cimg_forXY(rgb_img, x, y) {
		int r = rgb_img(x, y, 0);
		int g = rgb_img(x, y, 1);
		int b = rgb_img(x, y, 2);
		gray_img(x, y) = (unsigned char)(0.299 * r + 0.587 * g + 0.114 * b);
	}
	return gray_img;
}

/* Homography Estimation
   For homography we have: x'=Hx
   where p1=x'=[u,v] is a point in the reference image
   and p2=x=[x,y] is a point in the image that we want to warp
   H=[a,b,c;d,e,f;g,h,1];
*/
HomographyMatrix get_homography_matrix(const vector<KeypointPair>& pair) {
	assert(pair.size() == 4);

	float u0 = pair[0].p1.x, v0 = pair[0].p1.y;
	float u1 = pair[1].p1.x, v1 = pair[1].p1.y;
	float u2 = pair[2].p1.x, v2 = pair[2].p1.y;
	float u3 = pair[3].p1.x, v3 = pair[3].p1.y;

	float x0 = pair[0].p2.x, y0 = pair[0].p2.y;
	float x1 = pair[1].p2.x, y1 = pair[1].p2.y;
	float x2 = pair[2].p2.x, y2 = pair[2].p2.y;
	float x3 = pair[3].p2.x, y3 = pair[3].p2.y;

	float a, b, c, d, e, f, g, h;

	a = -(u0*v0*v1*x2 - u0*v0*v2*x1 - u0*v0*v1*x3 + u0*v0*v3*x1 - u1*v0*v1*x2 + u1*v1*v2*x0 + u0*v0*v2*x3 - u0*v0*v3*x2 + u1*v0*v1*x3 - u1*v1*v3*x0 + u2*v0*v2*x1 - u2*v1*v2*x0
		- u1*v1*v2*x3 + u1*v1*v3*x2 - u2*v0*v2*x3 + u2*v2*v3*x0 - u3*v0*v3*x1 + u3*v1*v3*x0 + u2*v1*v2*x3 - u2*v2*v3*x1 + u3*v0*v3*x2 - u3*v2*v3*x0 - u3*v1*v3*x2 + u3*v2*v3*x1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
			- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	b = (u0*u1*v0*x2 - u0*u2*v0*x1 - u0*u1*v0*x3 - u0*u1*v1*x2 + u0*u3*v0*x1 + u1*u2*v1*x0 + u0*u1*v1*x3 + u0*u2*v0*x3 + u0*u2*v2*x1 - u0*u3*v0*x2 - u1*u2*v2*x0 - u1*u3*v1*x0
		- u0*u2*v2*x3 - u0*u3*v3*x1 - u1*u2*v1*x3 + u1*u3*v1*x2 + u1*u3*v3*x0 + u2*u3*v2*x0 + u0*u3*v3*x2 + u1*u2*v2*x3 - u2*u3*v2*x1 - u2*u3*v3*x0 - u1*u3*v3*x2 + u2*u3*v3*x1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
			- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	c = (u0*v1*x2 - u0*v2*x1 - u1*v0*x2 + u1*v2*x0 + u2*v0*x1 - u2*v1*x0 - u0*v1*x3 + u0*v3*x1 + u1*v0*x3 - u1*v3*x0 - u3*v0*x1 + u3*v1*x0
		+ u0*v2*x3 - u0*v3*x2 - u2*v0*x3 + u2*v3*x0 + u3*v0*x2 - u3*v2*x0 - u1*v2*x3 + u1*v3*x2 + u2*v1*x3 - u2*v3*x1 - u3*v1*x2 + u3*v2*x1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
			- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	d = (u0*u1*v0*v2*x3 - u0*u1*v0*v3*x2 - u0*u2*v0*v1*x3 + u0*u2*v0*v3*x1 + u0*u3*v0*v1*x2 - u0*u3*v0*v2*x1 - u0*u1*v1*v2*x3 + u0*u1*v1*v3*x2 + u1*u2*v0*v1*x3 - u1*u2*v1*v3*x0 - u1*u3*v0*v1*x2 + u1*u3*v1*v2*x0
		+ u0*u2*v1*v2*x3 - u0*u2*v2*v3*x1 - u1*u2*v0*v2*x3 + u1*u2*v2*v3*x0 + u2*u3*v0*v2*x1 - u2*u3*v1*v2*x0 - u0*u3*v1*v3*x2 + u0*u3*v2*v3*x1 + u1*u3*v0*v3*x2 - u1*u3*v2*v3*x0 - u2*u3*v0*v3*x1 + u2*u3*v1*v3*x0)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
			- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	e = -(u0*v0*v1*y2 - u0*v0*v2*y1 - u0*v0*v1*y3 + u0*v0*v3*y1 - u1*v0*v1*y2 + u1*v1*v2*y0 + u0*v0*v2*y3 - u0*v0*v3*y2 + u1*v0*v1*y3 - u1*v1*v3*y0 + u2*v0*v2*y1 - u2*v1*v2*y0
		- u1*v1*v2*y3 + u1*v1*v3*y2 - u2*v0*v2*y3 + u2*v2*v3*y0 - u3*v0*v3*y1 + u3*v1*v3*y0 + u2*v1*v2*y3 - u2*v2*v3*y1 + u3*v0*v3*y2 - u3*v2*v3*y0 - u3*v1*v3*y2 + u3*v2*v3*y1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
			- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	f = (u0*u1*v0*y2 - u0*u2*v0*y1 - u0*u1*v0*y3 - u0*u1*v1*y2 + u0*u3*v0*y1 + u1*u2*v1*y0 + u0*u1*v1*y3 + u0*u2*v0*y3 + u0*u2*v2*y1 - u0*u3*v0*y2 - u1*u2*v2*y0 - u1*u3*v1*y0
		- u0*u2*v2*y3 - u0*u3*v3*y1 - u1*u2*v1*y3 + u1*u3*v1*y2 + u1*u3*v3*y0 + u2*u3*v2*y0 + u0*u3*v3*y2 + u1*u2*v2*y3 - u2*u3*v2*y1 - u2*u3*v3*y0 - u1*u3*v3*y2 + u2*u3*v3*y1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
			- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	g = (u0*v1*y2 - u0*v2*y1 - u1*v0*y2 + u1*v2*y0 + u2*v0*y1 - u2*v1*y0 - u0*v1*y3 + u0*v3*y1 + u1*v0*y3 - u1*v3*y0 - u3*v0*y1 + u3*v1*y0
		+ u0*v2*y3 - u0*v3*y2 - u2*v0*y3 + u2*v3*y0 + u3*v0*y2 - u3*v2*y0 - u1*v2*y3 + u1*v3*y2 + u2*v1*y3 - u2*v3*y1 - u3*v1*y2 + u3*v2*y1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
			- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	h = (u0*u1*v0*v2*y3 - u0*u1*v0*v3*y2 - u0*u2*v0*v1*y3 + u0*u2*v0*v3*y1 + u0*u3*v0*v1*y2 - u0*u3*v0*v2*y1 - u0*u1*v1*v2*y3 + u0*u1*v1*v3*y2 + u1*u2*v0*v1*y3 - u1*u2*v1*v3*y0 - u1*u3*v0*v1*y2 + u1*u3*v1*v2*y0
		+ u0*u2*v1*v2*y3 - u0*u2*v2*v3*y1 - u1*u2*v0*v2*y3 + u1*u2*v2*v3*y0 + u2*u3*v0*v2*y1 - u2*u3*v1*v2*y0 - u0*u3*v1*v3*y2 + u0*u3*v2*v3*y1 + u1*u3*v0*v3*y2 - u1*u3*v2*v3*y0 - u2*u3*v0*v3*y1 + u2*u3*v1*v3*y0)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
			- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	return HomographyMatrix(a, b, c, d, e, f, g, h);
}

float get_warped_x(float x, float y, HomographyMatrix H) {
	return H.a * x + H.b * y + H.c * x * y + H.d;
}

float get_warped_y(float x, float y, HomographyMatrix H) {
	return H.e * x + H.f * y + H.g * x * y + H.h;
}

float get_min_warped_x(const CImg<unsigned char> &input_img, HomographyMatrix H) {
	int w = input_img.width();
	int h = input_img.height();
	float min_x = cimg::min(get_warped_x(0, 0, H), get_warped_x(w - 1, 0, H),
		get_warped_x(0, h - 1, H), get_warped_x(w - 1, h - 1, H));
	return min_x < 0 ? min_x : 0;
}

float get_min_warped_y(const CImg<unsigned char> &input_img, HomographyMatrix H) {
	int w = input_img.width();
	int h = input_img.height();
	float min_y = cimg::min(get_warped_y(0, 0, H), get_warped_y(w - 1, 0, H),
		get_warped_y(0, h - 1, H), get_warped_y(w - 1, h - 1, H));
	return min_y < 0 ? min_y : 0;
}


float get_max_warped_x(const CImg<unsigned char> &input_img, HomographyMatrix H,
	const CImg<unsigned char> &stitched_img) {
	int w = input_img.width();
	int h = input_img.height();
	float max_x = cimg::max(get_warped_x(0, 0, H), get_warped_x(w - 1, 0, H),
		get_warped_x(0, h - 1, H), get_warped_x(w - 1, h - 1, H));
	return max_x >= stitched_img.width() ? max_x : stitched_img.width();
}

float get_max_warped_y(const CImg<unsigned char> &input_img, HomographyMatrix H,
	const CImg<unsigned char> &stitched_img) {
	int w = input_img.width();
	int h = input_img.height();
	float max_y = cimg::max(get_warped_y(0, 0, H), get_warped_y(w - 1, 0, H),
		get_warped_y(0, h - 1, H), get_warped_y(w - 1, h - 1, H));
	return max_y >= stitched_img.height() ? max_y : stitched_img.height();
}

template <class T>
T bilinear_interpolation(const CImg<T>& image, float x, float y, int channel) {
    /* This function comes from https://github.com/AmazingZhen/ImageStitching */
	assert(x >= 0 && x < image.width());
	assert(y >= 0 && y < image.height());
	assert(channel <= image.spectrum());

	int x_pos = floor(x);
	float x_u = x - x_pos;
	int xb = (x_pos < image.width() - 1) ? x_pos + 1 : x_pos;

	int y_pos = floor(y);
	float y_v = y - y_pos;
	int yb = (y_pos < image.height() - 1) ? y_pos + 1 : y_pos;

	float P1 = image(x_pos, y_pos, channel) * (1 - x_u) + image(xb, y_pos, channel) * x_u;
	float P2 = image(x_pos, yb, channel) * (1 - x_u) + image(xb, yb, channel) * x_u;

	return P1 * (1 - y_v) + P2 * y_v;
}

CImg<unsigned char> cylinderProjection(const CImg<unsigned char> &src) {
	/* This function comes from https://github.com/AmazingZhen/ImageStitching */
	int projection_width, projection_height;
	CImg<unsigned char> res(src.width(), src.height(), 1, src.spectrum(), 0);
	float r;
	float angle = 15.0;
	if (src.width() > src.height()) {
		projection_width = src.height();
		projection_height = src.width();

		r = (projection_width / 2.0) / tan(angle * cimg::PI / 180.0);

		for (int i = 0; i < res.width(); i++) {
			for (int j = 0; j < res.height(); j++) {
				float dst_x = j - projection_width / 2;
				float dst_y = i - projection_height / 2;

				float k = r / sqrt(r * r + dst_x * dst_x);
				float src_x = dst_x / k;
				float src_y = dst_y / k;

				if (src_x + projection_width / 2 >= 0 && src_x + projection_width / 2 < src.height()
					&& src_y + projection_height / 2 >= 0 && src_y + projection_height / 2 < src.width()) {
					for (int k = 0; k < res.spectrum(); k++) {
						res(i, j, k) = bilinear_interpolation(src, src_y + projection_height / 2, src_x + projection_width / 2, k);
					}
				}
			}
		}

	}
	else {
		projection_width = src.width();
		projection_height = src.height();

		r = (projection_width / 2.0) / tan(angle * cimg::PI / 180.0);

		for (int i = 0; i < res.width(); i++) {
			for (int j = 0; j < res.height(); j++) {
				float dst_x = i - projection_width / 2;
				float dst_y = j - projection_height / 2;

				float k = r / sqrt(r * r + dst_x * dst_x);
				float src_x = dst_x / k;
				float src_y = dst_y / k;

				if (src_x + projection_width / 2 >= 0 && src_x + projection_width / 2 < src.width()
					&& src_y + projection_height / 2 >= 0 && src_y + projection_height / 2 < src.height()) {
					for (int k = 0; k < res.spectrum(); k++) {
						res(i, j, k) = bilinear_interpolation(src, src_x + projection_width / 2, src_y + projection_height / 2, k);
					}
				}
			}
		}

	}

	return res;
}
