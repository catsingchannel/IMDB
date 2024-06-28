#pragma once

#include <cstdlib>

#include "static_data.h"

namespace structs {
	using std::string;
	using std::vector;
	using std::pow;
	using std::sqrt;
	using std::roundf;
	using std::rand;
	struct primer_pair;

	int complement(const string& s, string& c) {
		for (auto rit = s.rbegin(); rit != s.rend(); rit++) {
			switch (*rit)
			{
			case 'A':
				c.push_back('T');
				break;
			case 'T':
				c.push_back('A');
				break;
			case 'C':
				c.push_back('G');
				break;
			case 'G':
				c.push_back('C');
				break;
			default:
				c.push_back('*');
				break;
			}
		}

		return 0;
	}

	struct position {
	public:
		int fs = 0;
		int fl = 0;
		int pl = 0;
		int rl = 0;

		position& operator=(const position&) = default;

		float get_distance(const position& p){
			return sqrt(pow((fs - p.fs), 2) + pow((fl - p.fl), 2) + pow((pl - p.pl), 2) + pow((rl - p.rl), 2));
		}

		string get_fwd_seq() const {
			return tmp.substr(fs, fl);
		}

		string get_rev_seq() const {
			string res = "";
			complement(tmp.substr((fs + pl - rl), rl), res);
			return res;
		}

		position& operator+=(const vector<float>& v) {
			fs = fs + roundf(v[0]); 
			fl = fl + roundf(v[1]);
			pl = pl + roundf(v[2]);
			rl = rl + roundf(v[3]);

			if (fl > primer_args.max_len || fl < primer_args.min_len) {
				fl = (rand() % (primer_args.max_len - primer_args.min_len)) + primer_args.min_len;
			}

			if (rl > primer_args.max_len || rl < primer_args.min_len) {
				rl = (rand() % (primer_args.max_len - primer_args.min_len)) + primer_args.min_len;
			}

			if (pl > primer_args.max_prolen || pl < primer_args.min_prolen) {
				pl = (rand() % (primer_args.max_prolen - primer_args.min_prolen)) + primer_args.min_prolen;
			}

			if (fs > tmp.size() - pl || fs < 0) {
				fs = rand() % (tmp.size() - pl);
			}

			return *this ;
		}

		position operator-(const position& p) const {
			position res = {};

			res.fs = fs - p.fs;
			res.fl = fl - p.fl;
			res.pl = pl - p.pl;
			res.rl = rl - p.rl;

			return res;
		}

		vector<float> operator*(const vector<float>& p) const {
			vector<float> res = {};

			res.push_back((float)fs * p[0]);
			res.push_back((float)fl * p[1]);
			res.push_back((float)pl * p[2]);
			res.push_back((float)rl * p[3]);

			return res;
		}

		bool operator<(const position& p) {
			return((fs < p.fs) && (fl < p.fl) && (pl < p.pl) && (rl < p.rl));
		}
	};

	struct primer_pair {
	public:
		int fwd_start = 0;
		int fwd_length = 0;
		int rev_start = 0;
		int rev_length = 0;
		int pro_length = 0;
		int clamp_len = 0;
		int spec_f = 0;
		int spec_r = 0;
		int calc = 0;
		float self_fwd = 0;
		float self_rev = 0;
		float cross = 0;
		float hp_fwd = 0;
		float hp_rev = 0;
		float fwd_tm = 0;
		float rev_tm = 0;
		float fwd_gc = 0;
		float rev_gc = 0;
		float tm_diff = 0;
		float gc_diff = 0;
		float coverage = 0;

		primer_pair(){}

		primer_pair(const position& p) {
			fwd_start = p.fs;
			fwd_length = p.fl;
			rev_start = p.fs + p.pl - p.rl;
			rev_length = p.rl;
			pro_length = p.pl;
		}

		~primer_pair(){}

		bool operator==(const primer_pair& p) {
			return ((fwd_start == p.fwd_start) && (fwd_length == p.fwd_length) && (pro_length == p.pro_length) && (rev_length == p.rev_length));
		}

		string get_fwd_seq() {
			return tmp.substr(fwd_start, fwd_length);
		}

		string get_rev_seq() {
			string res = "";
			complement(tmp.substr(rev_start, rev_length), res);
			return res;
		}

		position get_position() {
			position p;
			p.fs = fwd_start;
			p.fl = fwd_length;
			p.pl = pro_length;
			p.rl = rev_length;

			return p;
		}
	};

	struct probe {
	public:
		int start = 0;
		int length = 0;
		float tm = 0;
		float gc = 0;
		float coverage = 0;
		bool reverse = false;

		probe(){}

		probe(const int s, const int l, const bool r){
			start = s;
			length = l;
			reverse = r;
		}

		~probe() {}

		string get_sequence() {
			if (reverse) {
				string res = "";
				complement(tmp.substr(start, length), res);
				return res;
			}
			else {
				return tmp.substr(start, length);
			}
		}
	};

	struct assay {
		primer_pair pp;
		probe p;

		assay(const primer_pair& opp, const probe& op) {
			pp = opp;
			p = op;
		}
	};
}