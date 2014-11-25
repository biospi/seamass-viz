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

#include <iostream>
#include <boost/filesystem.hpp>

#include "SMVStreamer.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    cout << endl;
    cout << "seaMass Viz Client - Copyright (C) 2014 - biospI Laboratory" << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
    cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;  
	cout << endl;	
	if (argc != 8)
	{
		cout << "Usage" << endl;
		cout << "-----" << endl;
		cout << "viz_client <in_idx> <mz_min> <mz_max> <rt_min> <rt_max> <out_w> <out_h>" << endl;
		cout << endl;
		cout << "<in_file>:  Input idx file" << endl;
		cout << "<mz_min>:   Minimum m/z to display" << endl;
		cout << "<mz_max>:   Maxmimum m/z to display" << endl;
		cout << "<rt_min>:   Minimum retention time to display" << endl;
		cout << "<rt_max>:   Maximum retention time to display" << endl;
		cout << "<out_w>:    Output width in pixels" << endl;
		cout << "<out_h>:    Output height in pixels" << endl;
		return 0;
	}
	string in_file(argv[1]);
	double mz_min = atof(argv[2]);
	double mz_max = atof(argv[3]);
	double rt_min = atof(argv[4]);
	double rt_max = atof(argv[5]);
	int out_w  = atoi(argv[6]);
	int out_h  = atoi(argv[7]);

	boost::filesystem::path viz_path(in_file);
	boost::filesystem::path viz_dir = viz_path.parent_path();
	ostringstream oss; oss << viz_dir.stem().string() << ".out";
	Reconstructer* recon = new Reconstructer(oss.str());

	SMVStreamer* stream = new SMVStreamer(in_file);
	stream->stream_cs(mz_min, mz_max, rt_min, rt_max, 4096, recon);
	delete stream;

	delete recon;

    cout << endl;
	return 0;
}
