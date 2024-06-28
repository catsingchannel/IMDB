// primer design.cpp: 定义应用程序的入口点。
//

#include "primer_design.h"

void set_profile(const py::dict& p, const int& s) {
	profile.clear();
	vector<string> all = {};

	for (auto it = p.begin(); it != p.end(); it++) {
		py::list pl = p[it->first];
		vector<string> l = {};
		for (auto it = pl.begin(); it != pl.end(); it++) {
			string temp = it->cast<string>();
			l.push_back(temp);
			all.push_back(temp);
		}
		std::cout<<it->first.cast<int>();
		int cstemp = it->first.cast<int>();
		profile.push_back(make_pair((cstemp - s), l));
	}

	all.erase(unique(all.begin(), all.end()), all.end());
	total_sample = all.size();
	sort(profile.begin(), profile.end(), [](const pair<int, vector<string>>& p1, const pair<int, vector<string>>& p2) { return p1.first < p2.first; });
	// for(auto it : profile){
	// 	std::cout<<"C:"<<it.first<<std::endl;
	// }
}

void parser(const py::kwargs& kwargs) {
	assays.clear();

	set_profile(kwargs["profile"], kwargs["start"].cast<int>());

	if (kwargs) {
		for (auto it = kwargs.begin(); it != kwargs.end(); it++) {
			string varname = py::str(it->first);
			if (varname == "particles") {
				particle_num = it->second.cast<int>();
				continue;
			}

			if (varname == "iternum") {
				iter_max = it->second.cast<int>();
				continue;
			}

			if (varname == "phi1") {
				phi1 = it->second.cast<float>();
				continue;
			}

			if (varname == "phi2") {
				phi2 = it->second.cast<float>();
				continue;
			}

			if (varname == "primerMaxT") {
				primer_args.max_tm = it->second.cast<float>();
				continue;
			}

			if (varname == "primerMinT") {
				primer_args.min_tm = it->second.cast<float>();
				continue;
			}

			if (varname == "primerTDiff") {
				primer_args.tm_diff = it->second.cast<float>();
				continue;
			}

			if (varname == "primerMaxGC") {
				primer_args.max_gc = it->second.cast<float>();
				primer_args.max_gc = primer_args.max_gc * 0.01;
				continue;
			}

			if (varname == "primerMinGC") {
				primer_args.min_gc = it->second.cast<float>();
				primer_args.min_gc = primer_args.min_gc * 0.01;
				continue;
			}

			if (varname == "primerGCDiff") {
				primer_args.gc_diff = it->second.cast<float>();
				primer_args.gc_diff = primer_args.gc_diff * 0.01;
				continue;
			}

			if (varname == "primerMaxP") {
				primer_args.max_prolen = it->second.cast<int>();
				continue;
			}

			if (varname == "primerMinP") {
				primer_args.min_prolen = it->second.cast<int>();
				continue;
			}

			if (varname == "primerMaxL") {
				primer_args.max_len = it->second.cast<int>();
				continue;
			}

			if (varname == "primerMinL") {
				primer_args.min_len = it->second.cast<int>();
				continue;
			}

			if (varname == "primerLDiff") {
				primer_args.len_diff = it->second.cast<int>();
				continue;
			}

			if (varname == "temparature") {
				temperature = it->second.cast<float>();
				continue;
			}

			if (varname == "primerConc") {
				concPrimer = it->second.cast<int>();
				continue;
			}

			if (varname == "probeMaxT") {
				probe_args.max_tm = it->second.cast<float>();
				continue;
			}

			if (varname == "probeMinT") {
				probe_args.min_tm = it->second.cast<float>();
				continue;
			}

			if (varname == "probeMaxGC") {
				probe_args.max_gc = it->second.cast<float>();
				probe_args.max_gc = probe_args.max_gc * 0.01;
				continue;
			}

			if (varname == "probeMinGC") {
				probe_args.min_gc = it->second.cast<float>();
				probe_args.min_gc = probe_args.min_gc * 0.01;
				continue;
			}

			if (varname == "probeMaxL") {
				probe_args.max_len = it->second.cast<int>();
				continue;
			}

			if (varname == "probeMinL") {
				probe_args.min_len = it->second.cast<int>();
				continue;
			}

			if (varname == "probeConc") {
				concProbe = it->second.cast<float>();
				continue;
			}
		}
	}

	string nt = kwargs["seq"].cast<string>();
	if (nt == tmp) {
		if (kwargs["rerun"].cast<bool>()) {
			mopso.clear_particles();
			mopso.reset_iter();
		}
		else {
			mopso.reset();
		}
	}
	else {
		tmp.clear();
		tmp.insert(tmp.end(), nt.begin(), nt.end());
		mopso.reset();
	}
}

py::list get_result_() {
	py::list res;
	int count = 1;

	for (auto it = assays.begin(); it != assays.end(); it++) {
		py::dict temp;
		temp["id"] = count;
		count++;
		temp["start"] = it->pp.fwd_start;
		temp["length"] = it->pp.pro_length;
		temp["coverage"] = roundf(get_coverage({ it->pp.fwd_start, (it->pp.fwd_length + it->pp.fwd_start - 1), it->p.start, (it->p.start + it->p.length - 1), it->pp.rev_start, (it->pp.rev_start + it->pp.rev_length - 1) }) * 100) / 100;
		temp["forward_start"] = it->pp.fwd_start;
		temp["forward_length"] = it->pp.fwd_length;
		temp["forward_seq"] = py::str(it->pp.get_fwd_seq().c_str());
		temp["forward_tm"] = roundf(it->pp.fwd_tm * 100) / 100;
		temp["forward_gc"] = roundf(it->pp.fwd_gc * 100) / 100;
		temp["forward_deltaG"] = roundf(utils::get_deltaG(utils::get_tm_params(it->pp.get_fwd_seq())) * 100) / 100;
		temp["reverse_start"] = it->pp.rev_start;
		temp["reverse_length"] = it->pp.rev_length;
		temp["reverse_seq"] = py::str(it->pp.get_rev_seq().c_str());
		temp["reverse_tm"] = roundf(it->pp.rev_tm * 100) / 100;
		temp["reverse_gc"] = roundf(it->pp.rev_gc * 100) / 100;
		temp["reverse_deltaG"] = roundf(utils::get_deltaG(utils::get_tm_params(it->pp.get_rev_seq())) * 100) / 100;
		temp["tm_difference"] = roundf(it->pp.tm_diff * 100) / 100;
		temp["gc_difference"] = roundf(it->pp.gc_diff * 100) / 100;
		temp["forward_dimer"] = roundf(it->pp.self_fwd * 100) / 100;
		temp["forward_hairpin"] = roundf(it->pp.hp_fwd * 100) / 100;
		temp["reverse_dimer"] = roundf(it->pp.self_rev * 100) / 100;
		temp["reverse_hairpin"] = roundf(it->pp.hp_rev * 100) / 100;
		temp["cross_dimer"] = roundf(it->pp.cross * 100) / 100;
		temp["probe_start"] = it->p.start;
		temp["probe_length"] = it->p.length;
		temp["probe_seq"] = py::str(it->p.get_sequence().c_str());
		temp["probe_tm"] = roundf(it->p.tm * 100) / 100;
		temp["probe_gc"] = roundf(it->p.gc * 100) / 100;
		temp["probe_deltaG"] = roundf(utils::get_deltaG(utils::get_tm_params(it->p.get_sequence())) * 100) / 100;

		res.append(temp);
	}

	return res;
}

py::list primer_design_(const py::kwargs& kwargs) {
	if (!thal_wrapper) {
		thal_wrapper.init();
	}
	
	parser(kwargs);

	design_probes();
	mopso.step_to();

	vector<primer_pair>& res = mopso.get_result();
	for (auto it = res.begin(); it != res.end();) {
		primer_pair local = *it;

		if (local.tm_diff > primer_args.tm_diff || local.gc_diff > primer_args.gc_diff || abs(local.fwd_length - local.rev_length) > primer_args.len_diff) {
			it = res.erase(it);
			continue;
		}

		if (local.fwd_tm > primer_args.max_tm || local.rev_tm > primer_args.max_tm || local.fwd_tm < primer_args.min_tm || local.rev_tm < primer_args.min_tm) {
			it = res.erase(it);
			continue;
		}

		if (local.fwd_gc > primer_args.max_gc || local.rev_gc > primer_args.max_gc || local.fwd_gc < primer_args.min_gc || local.rev_gc < primer_args.min_gc) {
			it = res.erase(it);
			continue;
		}

		if (local.self_fwd >= primer_args.min_tm || local.self_rev >= primer_args.min_tm || local.hp_fwd >= primer_args.min_tm || local.hp_rev >= primer_args.min_tm || local.cross >= primer_args.min_tm) {
			it = res.erase(it);
			continue;
		}

		int t = get_probe(local.get_position());

		if (t < 0) {
			it = res.erase(it);
			continue;
		}
		assays.push_back(assay((*it), probes[t]));

		it++;
	}

	return get_result_();
}

//int main(int argc, char** argv)
//{
//	//std::cout << "Hello CMake." << std::endl;
//	
//	string ipt = "TCATGCTACTAGAGAAGCTGTTGGTACCAATTTACCTTTACAGCTAGGTTTTTCTACAGGTGTTAACCTAGTTGCTGTACCTACAGGTTATGTTGATACACCTAATAATACAGATTTTTCCAGAGTTAGTGCTAAACCACCGCCTGGAGATCAATTTAAACACCTCATACCACTTATGTACAAAGGACTTCCTTGGAATGTAGTGCGTATAAAGATTGTACAAATGTTAAGTGACACACTTAAAAATCTCTCTGACAGAGTCGTATTTGTCTTATGGGCACATGGCTTTGAGTTGACATCTATGAAGTATTTTGTGAAAATAGGACCTGAGCGCACCTGTTGTCTATGTGATAGACGTGCCACATGCTTTTCCACTGCTTCAGACACTTATGCCTGTTGGCATCATTCTATTGGATTTGATTACGTCTATAATCCGTTTATGATTGATGTTCAACAATGGGGTTTTACAGGTAACCTACAAAGCAACCATGATCTGTATTGTCAAGTCCATGGTAATGCACATGTAGCTAGTTGTGATGCAATCATGACTAGGTGTCTAGCTGTCCACGAGTGCTTTGTTAAGCGTGTTGACTGGACTATTGAATATCCTATAATTGGTGATGAACTGAAGATTAATGCGGCTTGTAGAAAGGTTCAACACATGGTTGTTAAAGCTGCATTATTAGCAGACAAATTCCCAGTTCTTCACGACATTGGTAACCCTAAAGCTATTAAGTGTGTACCTCAAGCTGATGTAGAATGGAAGTTCTATGATGCACAGCCTTGTAGTGACAAAGCTTATAAAATAGAAGAATTATTCTATTCTTATGCCACACATTCTGACAAATTCACAGATGGTGTATGCCTATTTTGGAATTGCAATGTCGATAGATATCCTGCTAATTCCATTGTTTGTAGATTTGACACTAGAGTGCTATCTAACCTTAACTTGCCTGGTTGTGATGGTGGCAGTTTGTATGTAAATAAACATGCATTCCACACACCAGCTTTTGATAAAAGTGCTTTTGTTAATTTAAAACAATTACCATTTTTCTATTACTCTGACAGTCCATGTGAGTCTCATGGAAAACAAGTAGTGTCAGATATAGATTATGTACCACTAAAGTCTGCTACGTGTATAACACGTTGCAATTTAGGTGGTGCTGTCTGTAGACATCATGCTAATGAGTACAGATTGTATCTCGATGCTTATAACATGATGATCTCAGCTGGCTTTAGCTTGTGGGTTTACAAACAATTTGATACTTATAACCTCTGGAACACTTTTACAAGACTTCAGAGTTTAGAAAATGTGGCTTTTAATGTTGTAAATAAGGGACACTTTGATGGACAACAGGGTGAAGTACCAGTTTCTATCATTAATAACACTGTTTACACAAAAGTTGATGGTGTTGATGTAGAATTGTTTGAAAATAAAAC";
//	thal_wrapper.init();
//
//	primer_design_(ipt);
//	vector<assay> result = assays;
//	return 0;
//}
