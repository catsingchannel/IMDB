#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <utility>
#include <cmath>

#include "singleton.h"

namespace static_data {
	using std::unordered_map;
	using std::vector;
	using std::string;
	using std::make_pair;

	struct nn {
	public:
		unordered_map<string, float>* dH;
		unordered_map<string, float>* dS;

		nn() {
			dH = new unordered_map<string, float>;
			dS = new unordered_map<string, float>;

			vector<char*> nns = { "AA","TT","AT","TA","CA","TG","GT","AC","CT","AG","GA","TC","CG","GC","GG","CC","Initiation","AT_penalty","Symmetry_corr" };
			vector<float> dH_value = { -7.6, -7.6, -7.2, -7.2, -8.5, -8.5, -8.4, -8.4, -7.8, -7.8, -8.2, -8.2, -10.6, -9.8, -8.0, -8.0, 0.2, 2.2, 0 };
			vector<float> dS_value = { -21.3, -21.3, -20.4, -21.3, -22.7, -22.7, -22.4, -22.4, -21.0, -21.0, -22.2, -22.2, -27.2, -24.4, -19.9, -19.9, -5.7, 6.9, -1.4 };

			for (int i = 0; i < nns.size(); i++) {
				dH->insert(make_pair(nns[i], dH_value[i]));
				dS->insert(make_pair(nns[i], dS_value[i]));
			}
		};

		~nn() {
			delete dH;
			delete dS;
		};
	};

	struct probe_args {
	public:
		float max_tm = 65;
		float min_tm = 55;
		float max_gc = 0.55;
		float min_gc = 0.45;
		int max_len = 22;
		int min_len = 18;
	};

	struct primer_args {
	public:
		float max_tm = 70;
		float min_tm = 50;
		float tm_diff = 5;
		float max_gc = 0.55;
		float min_gc = 0.45;
		float gc_diff = 0.05;
		int max_len = 22;
		int min_len = 18;
		int len_diff = 5;
		int max_prolen = 150;
		int min_prolen = 90;
	};

}

static float concPrimer = 250;
static float concProbe = 250;
static float temperature = 37;
static float concNa = 0.05;

static int iter_max = 500;
static int particle_num = 16;
static float phi1 = 2.0;
static float phi2 = 2.0;
static int total_sample = 0;

std::string tmp = "";

#define nn singleton<static_data::nn>::get_instance()
#define probe_args singleton<static_data::probe_args>::get_instance()
#define primer_args singleton<static_data::primer_args>::get_instance()
#define profile singleton<std::vector<std::pair<int, std::vector<std::string>>>>::get_instance()
