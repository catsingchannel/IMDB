#pragma once

#include "static_data.h"

namespace utils {
	using std::pair;
	using std::string;
	using std::vector;
	using std::log;
	using std::unique;

	#define tmParams pair<float, float>

	tmParams get_tm_params(const string& s) {
		float dHsum = nn.dH->find("Initiation")->second;
		float dSsum = nn.dS->find("Initiation")->second;

		if (s[0] == 'A' || s[0] == 'T') {
			dHsum += nn.dH->find("AT_penalty")->second;
			dSsum += nn.dS->find("AT_penalty")->second;
		}

		if (s.back() == 'A' || s.back() == 'T') {
			dHsum += nn.dH->find("AT_penalty")->second;
			dSsum += nn.dS->find("AT_penalty")->second;
		}

		for (int i = 0; i < s.size() - 1; i++) {
			string dimer = s.substr(i, 2);
			dHsum += nn.dH->find(dimer)->second;
			dSsum += nn.dS->find(dimer)->second;
		}
		dSsum += 0.368 * (s.size() - 1) * log(concNa);

		return tmParams(dHsum, dSsum);
	}

	float get_tm(const tmParams& tm_params, const float& concOligo) {
		return (tm_params.first * 1000) / ((tm_params.second + (1.987 * log(concOligo * 1e-9)))) - 273.15;
	}

	float get_deltaG(const tmParams& tm_param) {
		return (tm_param.first * 1000 - (temperature + 273.15) * tm_param.second) / 1000;
	}

	float gc_content(const string& s) {
		float count = 0;
		for (int i = 0; i < s.size(); i++) {
			if (s[i] == 'G' || s[i] == 'C') {
				count++;
			}
		}

		return count / (float) s.size();
	}

	bool repeats(const string& s) {
		int sc = 1, dc1 = 1, dc2 = 1;
		char sr = s[0];
		string dm1 = s.substr(0, 2), dm2;
		
		for (int i = 1; i < s.size() - 1; i++) {
			if (s[i] == sr) {
				sc++;
				if (sc >= 5) {
					return true;
				}
			}
			else {
				sr = s[i];
				sc = 1;
			}

			if (i == 1) {
				dm2 = s.substr(1, 2);
				continue;
			}

			string dmer = s.substr(i, 2);
			if (i % 2 == 0) {
				if (dmer == dm1) {
					dc1++;
					if (dc1 >= 4) {
						return true;
					}
				}
				else {
					dm1 = dmer;
					dc1 = 1;
				}
			}
			else {
				if (dmer == dm2) {
					dc2++;
					if (dc2 >= 4) {
						return true;
					}
				}
				else {
					dm2 = dmer;
					dc2 = 1;
				}
			}
		}

		return false;
	}

	bool five_end_G(const string& s) {
		return(s[0] == 'G');
	}

	bool three_end_runs(const string& s){
		auto rit = s.rbegin();
		char ch = *rit;
		rit++;
		if (*rit != ch) {
			return false;
		}
		rit++;
		if (*rit != ch) {
			return false;
		}
		return true;
	}

	bool gc_clamp(const string& s) {
		int count = 0;
		auto rit = s.rbegin();
		for (int i = 0; i < 5; i++) {
			if (*rit == 'G' || *rit == 'C') {
				count++;
			}
			rit++;
		}

		if (count >= 2 && count <= 3) {
			return true;
		}
		else {
			return false;
		}
	}

	void get_LPS(const string& pat, const int& m, vector<int>& lps) {
		int len = 0;
		lps[0] = 0;

		for (int i = 1; i < m ;) {
			if (pat[i] == pat[len]) {
				len++;
				lps[i] = len;
				i++;
			}
			else {
				if (len != 0) {
					len = lps[len - 1];
				}
				else {
					lps[i] = len;
					i++;
				}
			}
		}
	}

	int KMP_search(const string& s, const string& pat) {
		int m = pat.size(), n = s.size();
		vector<int> LPS(m);

		get_LPS(pat, m, LPS);

		int res = 0;
		for (int i = 0, j = 0; i < n;) {
			if (pat[j] == s[i]) {
				i++;
				j++;
			}

			if (j == m) {
				j = LPS[j - 1];
				res++;
			}
			else if (i < n && pat[j] != s[i]) {
				if (j != 0) {
					j = LPS[j - 1];
				}
				else {
					i = i + 1;
				}
			}
		}

		return res;
	}

	float get_coverage(vector<int> p) {
		vector<string> all = {};

		if (p.size() == 6) {
			for (auto it = profile.begin(); it != profile.end(); it++) {
				if ((it->first >= p[0] && it->first <= p[1]) || (it->first >= p[2] && it->first <= p[3]) || (it->first >= p[4] && it->first <= p[5])) {
					all.insert(all.end(), it->second.begin(), it->second.end());
				}
				else if (it->first > p[5]) {
					break;
				}
			}
		}
		else if (p.size() == 4) {
			// for(auto it=p.begin(); it != p.end(); it++){
			// 	std::cout<<"param:"<<*it<<std::endl;
			// }

			for (auto it = profile.begin(); it != profile.end(); it++) {
				// std::cout<<"profile:"<<it->first<<"|"<<it->second.size()<<std::endl;
				// std::cout<<(it->first >= p[0])<<(it->first <= p[1])<<std::endl;
				// std::cout<<(it->first >= p[2])<<(it->first <= p[3])<<std::endl;

				if ((it->first >= p[0] && it->first <= p[1]) || (it->first >= p[2] && it->first <= p[3])) {
					// std::cout<<"hit:"<<it->first<<std::endl;
					all.insert(all.end(), it->second.begin(), it->second.end());
				}
				else if (it->first > p[3]) {
					break;
				}
			}
		}
		else {
			for (auto it = profile.begin(); it != profile.end(); it++) {
				if (it->first >= p[0] && it->first <= p[1]) {
					all.insert(all.end(), it->second.begin(), it->second.end());
				}
				else if (it->first > p[1]) {
					break;
				}
			}
		}

		all.erase(unique(all.begin(), all.end()), all.end());
		return (1.0 - (static_cast<float>(all.size()) / static_cast<float>(total_sample)));
	}
}