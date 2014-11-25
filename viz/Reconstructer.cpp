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

#include "Reconstructer.hpp"

#include <omp.h>
#include <boost/filesystem.hpp>
#include <sstream>
#include <iomanip>
#include <iostream>

Reconstructer::
Reconstructer(const string& _out_dir) :
	out_dir(_out_dir),
	stream_i(0)
{
	boost::filesystem::path viz_path(out_dir);
	if (boost::filesystem::exists(viz_path))
	{
		for (boost::filesystem::directory_iterator end_dir_it, it(viz_path); it!=end_dir_it; ++it)
		{
			boost::filesystem::remove_all(it->path());
		}
	}
	else
	{
		boost::filesystem::create_directories(viz_path);
	}
}

void
Reconstructer::
next_stream()
{
	chunk_i = 0;

	boost::filesystem::path viz_path(out_dir);
	ostringstream oss; oss << stream_i;
	viz_path /= oss.str();
	boost::filesystem::create_directory(viz_path);

	viz_path /= "stream.csv";	
	ofs.open(viz_path.string());
	start = omp_get_wtime();

	stream_i++;
}

void
Reconstructer::
next_chunk(const vector<Coef>& cs)
{
	boost::filesystem::path viz_path(out_dir);
	ostringstream oss; oss << stream_i-1;
	viz_path /= oss.str();
	ostringstream oss2; oss2 << setfill('0') << setw(8) << chunk_i++ << ".txt";
	viz_path /= oss2.str();

	ofstream ofs2(viz_path.string());
	for (size_t i = 0; i < cs.size(); i++)
	{
		ofs2 << cs[i].v << ":[" << cs[i].lx << "," << cs[i].ly << "]:[" << cs[i].x << "," << cs[i].y << "]" << endl;
	}

	ofs << omp_get_wtime() - start << endl;
}
