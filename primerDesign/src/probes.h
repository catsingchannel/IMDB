#pragma once

#include <algorithm>

#include "structs.h"
#include "utils.h"
#include "thread_pool.h"
#include "thal_wrapper.h"

#define probes singleton<std::vector<structs::probe>>::get_instance()

namespace probe {
	using namespace utils;
	using std::sort;
	using std::abs;
	using std::roundf;
	using std::cout;
	using std::endl;
	using std::min;
	using structs::probe;
	using structs::position;

	int filter_oligos(vector<probe>* p) {
		for (int i = 0; i < p->size();) {
			string oligo = (*p)[i].get_sequence();

			(*p)[i].tm = get_tm(get_tm_params(oligo), concProbe);
			if ((*p)[i].tm < probe_args.min_tm || (*p)[i].tm > probe_args.max_tm) {
				p->erase(p->begin() + i);
				continue;
			}

			(*p)[i].gc = gc_content(oligo);
			if ((*p)[i].gc < probe_args.min_gc || (*p)[i].gc > probe_args.max_gc) {
				p->erase(p->begin() + i);
				continue;
			}

			if (repeats(oligo) || three_end_runs(oligo) || five_end_G(oligo)) {
				p->erase(p->begin() + i);
				continue;
			}

			i++;
		}

		return 0;
	}

	int generate_oligos(const int& length, vector<probe>* result) {
		for (int i = 0; i < (tmp.size() - length + 1); i++) {
			result->push_back(probe(i, length, false));
			result->push_back(probe(i, length, true));
		}
		
		return filter_oligos(result);
	}

	void design_probes() {
		vector<vector<probe>> buckets = {};
		for (int i = probe_args.min_len; i < probe_args.max_len + 1; i++) {
			buckets.push_back(vector<probe>());
		}

		int l = probe_args.min_len;
		for (int i = 0; i < buckets.size(); i++) {
			thread_pool_.push(generate_oligos, l, &buckets[i]);
			l++;
		}
		thread_pool_.wait();

		for (int i = 0; i < buckets.size(); i++) {
			probes.insert(probes.end(), buckets[i].begin(), buckets[i].end());
		}

		sort(probes.begin(), probes.end(), [](const probe& p1, const probe& p2) { return p1.start > p2.start; });
		float avg_tm = (probe_args.max_tm + probe_args.min_tm) / 2;
		float avg_gc = (probe_args.max_gc + probe_args.min_gc) / 2;
		for (int i = 0; (i + 1) < probes.size();) {
			probe it = probes[i];
			probe it2 = probes[i + 1];
			if (it.reverse == it2.reverse) {
				if ((it.start >= it2.start && (it.start + it.length <= it2.start + it2.length))
					|| (it2.start >= it.start && (it2.start + it2.length <= it.start + it.length))) {

					float diff = roundf(abs(it.tm - avg_tm)), diff2 = roundf(abs(it2.tm - avg_tm));
					if (diff > diff2) {
						probes.erase(probes.begin() + i);
						continue;
					}
					else if (diff < diff2) {
						probes.erase(probes.begin() + (i + 1));
						i++;
						continue;
					}

					diff = roundf(abs(it.gc - avg_gc));
					diff2 = roundf(abs(it2.gc - avg_gc));
					if (diff > diff2) {
						probes.erase(probes.begin() + i);
						continue;
					}
					else if (diff < diff2) {
						probes.erase(probes.begin() + (i + 1));
						i++;
						continue;
					}
				}
			}

			i++;
		}
	}

	int pick_probe(const vector<int>& v, const position& p) {
		string p_fwd = p.get_fwd_seq();
		string p_rev = p.get_rev_seq();
		float fwd_tm = get_tm(get_tm_params(p_fwd), concPrimer);
		float rev_tm = get_tm(get_tm_params(p_rev), concPrimer);
		float tgt_tm_min = roundf((((fwd_tm + rev_tm) * 0.5) + 5.0) * 100) / 100;
		float tgt_tm_max = roundf((tgt_tm_min + 5.0) * 100) / 100;
		float lmin = probe_args.max_tm;
		vector<int> r = {};

		for (int i = 0; i < v.size(); i++) {
			float ptm = get_tm(get_tm_params(probes[v[i]].get_sequence()), concProbe);
			ptm = roundf(ptm * 100) / 100;

			if (ptm <= tgt_tm_max && ptm >= tgt_tm_min) {
				ptm = 0;
			}
			else if (ptm > tgt_tm_max) {
				ptm = roundf(abs(ptm - tgt_tm_max) * 100) / 100;
			}
			else if (ptm < tgt_tm_min) {
				ptm = roundf(abs(ptm - tgt_tm_min) * 100) / 100;
			}

			if (ptm > lmin) {
				continue;
			}
			else if (ptm == lmin) {
				r.push_back(v[i]);
			}
			else if (ptm < lmin) {
				lmin = ptm;
				r.clear();
				r.push_back(v[i]);
			}
		}

		if (r.size() == 1) {
			return r[0];
		}

		sort(r.begin(), r.end(), [](const int& r1, const int& r2) { return get_coverage({ probes[r1].start, (probes[r1].start + probes[r1].length - 1) }) < get_coverage({ probes[r2].start, (probes[r2].start + probes[r2].length - 1) }); });

		lmin = std::numeric_limits<float>::max();
		int min_id = 0;
		int limit = min(5, static_cast<int>(r.size()));
		for (int i = 0; i < limit; i++) {
			string seq = probes[r[i]].get_sequence();
			float total = self_dimer(seq);
			total += cross_dimer(seq, p_fwd);
			total += cross_dimer(seq, p_rev);
			total += hairpin(seq);
			total = roundf(total * 100) / 100;

			if (total >= lmin) {
				continue;
			}
			else {
				lmin = total;
				min_id = r[i];
			}
		}

		return min_id;
	}

	int get_probe(const position& p) {
		vector<int> avail_p = {};

		if (probes.empty()) {
			cout << "not available probes" << endl;
			return -1;
		}

		for (int i = 0; i < probes.size(); i++) {
			probe& temp = probes[i];
			if ((temp.start >= (p.fs + p.fl)) && (temp.start + temp.length) < (p.fs + p.pl - p.rl)) {
				avail_p.push_back(i);
			}
		}

		if (avail_p.empty()) {
			return -1;
		}
		else if (avail_p.size() == 1) {
			return avail_p[0];
		}
		else {
			return pick_probe(avail_p, p);
		}
	}
}