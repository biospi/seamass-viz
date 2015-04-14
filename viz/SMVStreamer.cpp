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

#include "SMVStreamer.hpp"

#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>
#include <math.h>
#include <limits>

using namespace SpatialIndex;


struct MyPair
{
	id_type id;
	Region* r;

	MyPair() {}

	MyPair(const id_type& _id, Region* _r) :
		id(_id),
		r(_r) {}

	void operator=(const MyPair& pair)
	{
		id = pair.id;
		r = pair.r;
	}
};


// yes, this is reversed
inline bool operator<(const MyPair& a, const MyPair& b)
{
	return a.r->getLow(2) > b.r->getLow(2);
}


// if all the nodes were internally sorted by intensity
// we could use a quicker sorting procedure than priority_queue
class MyQueryStrategy : public SpatialIndex::IQueryStrategy
{
private:
	const IShape& query;
	Reconstructer* recon;

	priority_queue<MyPair> node_q;
	priority_queue<MyPair> data_q;
	vector<MyPair> chunk;
	int chunk_index;
	int index;

	double start;

public:
	MyQueryStrategy(const IShape& _query,
		            Reconstructer* _recon,
					int _chunk_size) : 
	  query(_query),
	  recon(_recon),
	  chunk(_chunk_size),
	  index(0)
	{
	}

	void getNextEntry(const IEntry& entry, id_type& nextEntry, bool& hasNext)
	{
		if (!&entry)
		{
			hasNext = false;
			return;
		}

		const INode* n = dynamic_cast<const INode*>(&entry);

		IShape* ps;
		n->getShape(&ps);
		Region* pr = dynamic_cast<Region*>(ps);

		for (uint32_t cChild = 0; cChild < n->getChildrenCount(); cChild++)
		{
			IShape* ps;
			n->getChildShape(cChild, &ps);
			Region* r = dynamic_cast<Region*>(ps);
			if (query.intersectsShape(*r))
			{
				if (n->getLevel() == 0)
					data_q.push(MyPair(n->getChildIdentifier(cChild), r));
				else
					node_q.push(MyPair(n->getChildIdentifier(cChild), r));
			}
		}

		while (node_q.empty() && !data_q.empty() ||
			   !data_q.empty() && data_q.top().r->getLow(2) <= node_q.top().r->getLow(2))
		{
			MyPair d = data_q.top();
			data_q.pop();
			visitData(d, !(node_q.empty() && data_q.size() == 1));
		}

		if (!node_q.empty())
		{
			Region* r = node_q.top().r;
			nextEntry = node_q.top().id; node_q.pop();
			hasNext = true;
			delete r;
		}
		else
		{
			hasNext = false;
		}
	}

	void visitData(const MyPair& pair, bool hasNext)
	{
		chunk[index++] = pair;

		if (!hasNext)
		{
			chunk.resize(index);
		}
		if (index == chunk.size())		
		{
			vector<Coef> cs(chunk.size());
			for (size_t i = 0; i < chunk.size(); i++)
			{
				double w = 0.25 * (chunk[i].r->getHigh(0) - chunk[i].r->getLow(0));
				double h = 0.25 * (chunk[i].r->getHigh(1) - chunk[i].r->getLow(1));
				
				cs[i].x = (int) (chunk[i].r->getLow(0) / w) + 3;
				cs[i].y = (int) (chunk[i].r->getLow(1) / h) + 3;
				cs[i].lx = (int) (-log(w) / log(2.0));
				cs[i].ly = (int) (-log(h) / log(2.0));
				cs[i].v = (fp) -chunk[i].r->getLow(2) / (w * h);
			
				delete chunk[i].r;
			}

			// reconstruct
			recon->next_chunk(cs);

			cout << "." << flush;

			index = 0;
		}
	}
};


SMVStreamer::
SMVStreamer(const string& filename)
{
	int lastdot = filename.find_last_of("."); 
	basename = (lastdot == string::npos) ? basename : filename.substr(0, lastdot);

    cout << "Reading " << basename << ".txt" << endl;
	ostringstream oss; oss << basename << ".txt";
	ifstream ifs(oss.str().c_str());
	ifs >> mz_min >> mz_max >> rt_min >> rt_max; 

    cout << "Reading " << basename << ".idx" << endl;
	timer::cpu_timer time;
	time.start();

	diskfile = StorageManager::loadDiskStorageManager(basename);
	// this will try to locate and open an already existing storage manager.

	file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
	// applies a main memory random buffer on top of the persistent storage manager
	// (LRU buffer, etc can be created the same way).

	// If we need to open an existing tree stored in the storage manager, we only
	// have to specify the index identifier as follows
	tree = RTree::loadRTree(*file, 1);

    boost::chrono::duration<double> seconds = boost::chrono::nanoseconds(time.elapsed().user);
	cout << "Duration: " << seconds.count() << " seconds" << endl;
}


SMVStreamer::
~SMVStreamer()
{
	delete tree;
	delete file;
	delete diskfile;
}


// only supports cm with dimension 2 at present
void
SMVStreamer::
stream_cs(double _mz_min, double _mz_max, double _rt_min, double _rt_max,
          size_t chunk_size, Reconstructer* recon)
{
	_mz_min = mz_min > _mz_min ? mz_min : _mz_min;
	_mz_max = mz_max < _mz_max ? mz_max : _mz_max;
	_rt_min = rt_min > _rt_min ? rt_min : _rt_min;
	_rt_max = rt_max < _rt_max ? rt_max : _rt_max;

	cout << endl << "Reading " << basename << ".dat" << endl;
	cout << "Streaming " << _mz_min << " to " << _mz_max << " m/z, " << _rt_min << " to " << _rt_max << " mins" << endl;
	cout << "press Ctrl-C to stop" << endl << endl;
	recon->next_stream(_mz_min, _mz_max, _rt_min, _rt_max);

	double low[3] = {
		_mz_min * 60/1.0033548378,
		_rt_min,
		-numeric_limits<double>::max()
	};
	double high[3] = {
		_mz_max * 60/1.0033548378,
		_rt_max,
		0.0
	};
	Region r(low, high, 3);

	MyQueryStrategy qs(r, recon, chunk_size);
	tree->queryStrategy(qs);
}
