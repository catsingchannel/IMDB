#pragma once

#include <functional>

#include "utils.h"
#include "structs.h"
#include "thal_wrapper.h"

#define Cmp_ int(const primer_pair& p1, const primer_pair& p2)

namespace comparators {
	using namespace std::placeholders;
	using namespace utils;

	using std::abs;
	using std::roundf;
	using std::bind;
	using std::forward;
	using std::string;
	using std::function;
	using std::vector;
	using structs::primer_pair;
	using structs::position;
	using structs::primer_pair;

	template <typename T>
	int cmp(const T& a0, const T& a1) {
		if (a0 > a1) {
			return 0;
		}
		else if (a0 == a1) {
			return 1;
		}
		else {
			return 2;
		}
	}

	int fwd_tm(const primer_pair& p1, const primer_pair& p2) {
		float avg = (primer_args.max_tm + primer_args.min_tm) / 2;

		float p1_tm = roundf((abs(p1.fwd_tm - avg) * 100)) / 100;
		float p2_tm = roundf((abs(p2.fwd_tm - avg) * 100)) / 100;

		return cmp(p1_tm, p2_tm);
	}

	int rev_tm(const primer_pair& p1, const primer_pair& p2) {
		float avg = (primer_args.max_tm + primer_args.min_tm) / 2;

		float p1_tm = roundf((abs(p1.rev_tm - avg) * 100)) / 100;
		float p2_tm = roundf((abs(p2.rev_tm - avg) * 100)) / 100;

		return cmp(p1_tm, p2_tm);
	}
	
	int tm_diff(const primer_pair& p1, const primer_pair& p2) {
		float p1d = p1.tm_diff;
		float p2d = p2.tm_diff;

		if (p1.tm_diff <= primer_args.tm_diff) {
			p1d = 0;
		}
		else {
			p1d = roundf((abs(p1.tm_diff - primer_args.tm_diff) * 100)) / 100;
		}

		if (p2.tm_diff <= primer_args.tm_diff) {
			p2d = 0;
		}
		else {
			p2d = roundf((abs(p2.tm_diff - primer_args.tm_diff) * 100)) / 100;
		}

		return cmp(p1d, p2d);
	}

	int clamp_len(const primer_pair& p1, const primer_pair& p2) {
		return cmp(p1.clamp_len, p2.clamp_len);
	}

	int fwd_gc(const primer_pair& p1, const primer_pair& p2) {
		float avg = (primer_args.max_gc + primer_args.min_gc) / 2;

		float p1_gc = roundf((abs(p1.fwd_gc - avg) * 100)) / 100;
		float p2_gc = roundf((abs(p2.fwd_gc - avg) * 100)) / 100;

		return cmp(p1_gc, p2_gc);
	}

	int rev_gc(const primer_pair& p1, const primer_pair& p2) {
		float avg = (primer_args.max_gc + primer_args.min_gc) / 2;

		float p1_gc = roundf((abs(p1.fwd_gc - avg) * 100)) / 100;
		float p2_gc = roundf((abs(p2.fwd_gc - avg) * 100)) / 100;

		return cmp(p1_gc, p2_gc);
	}

	int dimer_check(const primer_pair& p1, const primer_pair& p2) {
		int c1 = 0, c2 = 0;

		if (p1.self_fwd > p2.self_fwd) {
			c1++;
		}
		else if (p1.self_fwd < p2.self_fwd) {
			c2++;
		}

		if (p1.self_rev > p2.self_rev) {
			c1++;
		}
		else if (p1.self_rev < p2.self_rev) {
			c2++;
		}

		if (p1.cross > p2.cross) {
			c1++;
		}
		else if (p1.cross < p2.cross) {
			c2++;
		}

		return cmp(c1, c2);
	}

	int hairpin_specificity(const primer_pair& p1, const primer_pair& p2) {
		int pc1 = 0, pc2 = 0;

		if (p1.hp_fwd > p2.hp_fwd) {
			pc1++;
		}
		else if (p1.hp_fwd < p2.hp_fwd) {
			pc2++;
		}

		if (p1.hp_rev > p2.hp_rev) {
			pc1++;
		}
		else if (p1.hp_rev < p2.hp_rev) {
			pc2++;
		}

		if (p1.spec_f > p2.spec_f) {
			pc1++;
		}
		else if (p1.spec_f < p2.spec_f) {
			pc2++;
		}

		if (p1.spec_r > p2.spec_r) {
			pc1++;
		}
		else if (p1.spec_r < p2.spec_r) {
			pc2++;
		}

		return cmp(pc1, pc2);
	}

	int coverage_compare(const primer_pair& p1, const primer_pair& p2) {
		float c1 = roundf(p1.coverage * 10000) / 10000;
		float c2 = roundf(p2.coverage * 10000) / 10000;

		return cmp(c1, c2);
	}

	int calculate_params(primer_pair& res) {
		string fwd = res.get_fwd_seq();
		string rev = res.get_rev_seq();

		res.fwd_tm = get_tm(get_tm_params(fwd), concPrimer);
		res.rev_tm = get_tm(get_tm_params(rev), concPrimer);
		res.tm_diff = abs(res.fwd_tm - res.rev_tm);
		res.fwd_gc = gc_content(fwd);
		res.rev_gc = gc_content(rev);
		res.gc_diff = abs(res.fwd_gc - res.rev_gc);

		if (!gc_clamp(fwd) && !gc_clamp(rev)) {
			res.clamp_len += 2;
		}
		else if (!gc_clamp(fwd) || !gc_clamp(rev)) {
			res.clamp_len++;
		}

		if (abs(res.fwd_length - res.rev_length) > primer_args.len_diff) {
			res.clamp_len += abs(abs(res.fwd_length - res.rev_length) - primer_args.len_diff);
		}

		res.self_fwd = self_dimer(fwd);
		res.self_rev = self_dimer(rev);
		res.cross = cross_dimer(fwd, rev);
		res.hp_fwd = hairpin(fwd);
		res.hp_rev = hairpin(rev);
		res.spec_f = KMP_search(tmp, fwd);
		string rrev = {};
		structs::complement(rev, rrev);
		res.spec_r = KMP_search(tmp, rrev);
		res.coverage = get_coverage({res.fwd_start, (res.fwd_length + res.fwd_start - 1), res.rev_start, (res.rev_start + res.rev_length - 1)});

		res.calc = 1;
		return 0;
	}

	vector<function<Cmp_>> get_all_cmps() {
		vector<function<Cmp_>> res = {};
		res.push_back(function<Cmp_>(bind(forward<Cmp_>(fwd_tm), _1, _2)));
		res.push_back(function<Cmp_>(bind(forward<Cmp_>(rev_tm), _1, _2)));
		res.push_back(function<Cmp_>(bind(forward<Cmp_>(tm_diff), _1, _2)));
		res.push_back(function<Cmp_>(bind(forward<Cmp_>(clamp_len), _1, _2)));
		res.push_back(function<Cmp_>(bind(forward<Cmp_>(fwd_gc), _1, _2)));
		res.push_back(function<Cmp_>(bind(forward<Cmp_>(rev_gc), _1, _2)));
		res.push_back(function<Cmp_>(bind(forward<Cmp_>(dimer_check), _1, _2)));
		res.push_back(function<Cmp_>(bind(forward<Cmp_>(hairpin_specificity), _1, _2)));
		res.push_back(function<Cmp_>(bind(forward<Cmp_>(coverage_compare), _1, _2)));

		return res;
	}

	int compare_all(primer_pair& p1, primer_pair& p2) {
		if (!(p1.calc)) {
			calculate_params(p1);
		}

		if (!(p2.calc)) {
			calculate_params(p2);
		}

		int ls = 0, eq = 0, gt = 0;
		vector<function<Cmp_>> c = get_all_cmps();
		for (auto it = c.begin(); it != c.end(); it++) {
			switch ((*it)(p1, p2))
			{
			case 0:
				ls++;
				continue;
			case 1:
				eq++;
				continue;
			case 2:
				gt++;
				continue;
			default:
				break;
			}
		}

		if (ls == 0 && gt >= 1) {
			return 2;
		}
		else if (gt == 0 && ls >= 1) {
			return 0;
		}
		else {
			return 1;
		}
	}

	auto compare = function<int(primer_pair&, primer_pair&)>(bind(forward<int(primer_pair&, primer_pair&)>(compare_all), _1, _2));
}

