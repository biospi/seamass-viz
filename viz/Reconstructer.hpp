//
// Original authors: Andrew Dowsey <andrew.dowsey@manchester.ac.uk>
//
// Copyright 2014 CADET Bioinformatics Laboratory
//                University of Manchester, UK
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/gil/gil_all.hpp>
#include <Eigen/Core>

#include "PNGWriter.hpp"

using namespace std;
using namespace boost;
using namespace Eigen;

struct Coef
{
	int lx, ly, x, y;
	double v;
	
	bool operator()(Coef& a, Coef& b);
};


class Reconstructer
{
protected:
	filesystem::path out_path;

	Matrix<double,Dynamic,Dynamic> img; // accumulated image

	double mz_min;
	double mz_max;
	double rt_min;
	double rt_max;
	double max_counts;

	filesystem::path stream_path;
	ofstream         stream_ofs;
	int              stream_index;
	int              chunk_index;

	double read_start;
	double read_time;
	double viz_start;
	double viz_time;
	size_t chunk_count;
	PNGWriter writer;

	int nx,ny,nxy;

public:
	Reconstructer(const string& out_dir, int width, int height);

	void next_stream(double mz_min, double mz_max, double rt_min, double rt_max, double max_counts);
	void next_chunk(const vector<Coef>& cs);
};