#pragma once

#include <algorithm>
//#include <chrono>

#include "pareto.h"
#include "comparators.h"
#include "thread_pool.h"

namespace mopso {
	using std::cerr;
	using std::endl;
	using std::make_pair;
	using std::sort;
	using std::roundf;
	using std::rand;
	using std::bind;
	using std::unique;
	using std::forward;
	using std::max;
	//using std::sqrtf;
	//using std::chrono::duration_cast;
	using std::pair;
	using std::vector;
	using std::unique_ptr;
	using std::function;
	//using std::chrono::high_resolution_clock;
	//using std::chrono::milliseconds;
	using structs::position;
	using structs::primer_pair;
	using pareto::pareto_set;
	using comparators::compare_all;

	int get_index(const int& size, const int& iter_n, const int& iter_max) {
		int d_idx = roundf(((float)(size - 1) * (float)(iter_max - iter_n) / (float)iter_max) + 1) - 1;
		return d_idx;
	}

	class particle {
	private:
		pareto_set<primer_pair, function<int(primer_pair&, primer_pair&)>> self_set = {};
		vector<float> v = { 0.0,0.0,0.0,0.0 };
		primer_pair p = {};
		position pbest = {};
		position gbest = {};

		int iter = 0;
		int iter_limit = 0;
	public:
		particle() {
			p.fwd_start = (rand() % (primer_args.max_len - primer_args.min_len)) + primer_args.min_len;
			p.rev_length = (rand() % (primer_args.max_len - primer_args.min_len)) + primer_args.min_len;
			p.pro_length = (rand() % (primer_args.max_prolen - primer_args.min_prolen)) + primer_args.min_prolen;
			p.fwd_length = rand() % (tmp.size() - p.pro_length);

			float phi = phi1 + phi2;
			float i = 2.0 / ((phi - 2.0) + sqrtf((phi * phi) - (4.0 * phi)));
			for (int i = 0; i < 4; i++) {
				v[i] = i;
			}
		};

		~particle() = default;
		particle(const particle&) = default;
		particle& operator=(const particle&) = default;

		pareto_set<primer_pair, function<int(primer_pair&, primer_pair&)>>& get_set() {
			return self_set;
		}

		primer_pair& get_pos() {
			return p;
		}

		int get_iter() const{
			return iter;
		}

		int set_iter(int i) {
			iter = 0;
			return iter;
		}

		int get_limit() const{
			return iter_limit;
		}

		int set_limit(int l) {
			iter_limit = l;
			return iter_limit;
		}

		int set_gbest(const vector<primer_pair>& global, const int& idx) {
			vector<pair<float, int>> dist;
			position lp = p.get_position();

			for (int i = 0; i < global.size(); i++) {
				primer_pair a = global[i];
				position temp = a.get_position();
				dist.push_back(make_pair(lp.get_distance(temp), i));
			}
			sort(dist.begin(), dist.end(), [](const pair<float, int>& p1, const pair<float, int>& p2) { return (p1.first > p2.first); });
			primer_pair pp = global[dist[idx].second];
			gbest = pp.get_position();

			return 0;
		}

		void update_pbest() {
			self_set.pareto_ops(p);
			const vector<primer_pair>& a = self_set.get_archive();
			vector<pair<float, int>> dist;
			position lp = p.get_position();

			for (int i = 0; i < a.size(); i++) {
				primer_pair tpp = a[i];
				position temp = tpp.get_position();
				dist.push_back(make_pair(lp.get_distance(temp), i));
			}
			sort(dist.begin(), dist.end(), [](const pair<float, int>& p1, const pair<float, int>& p2) { return (p1.first > p2.first); });
			primer_pair t = a[dist[get_index(dist.size(), iter, iter_limit)].second];
			pbest = t.get_position();
		}

		void step_once(const vector<primer_pair>& global, const int& idx) {
			if (iter >= iter_limit) {
				return ;
			}

			position lp = p.get_position();
			set_gbest(global, idx);

			float u1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / phi1));
			float u2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / phi2));
			vector<float> u1v(4, u1);
			vector<float> u2v(4, u2);

			vector<float> pv = (pbest - lp) * u1v;
			vector<float> gv = (gbest - lp) * u2v;

			for (int i = 0; i < 4; i++) {
				float ov = v[i];
				v[i] = ov + pv[i] + gv[i];
			}
			lp += v;
			p = primer_pair(lp);

			iter++;
		}
	};

	class mopso {
	private:
		pareto_set<primer_pair, function<int(primer_pair&, primer_pair&)>> global_archive_ = {};
		unique_ptr<particle[]> particles_ = nullptr;

		int iter = 0;
		int iter_limit = iter_max;
		int iter_to = iter_limit;
		int pnum = 0;

		void reset_particles() {
			particles_.reset();
			particles_ = std::unique_ptr<particle[]>(new particle[particle_num]);
			for (int i = 0; i < particle_num; i++) {
				particles_[i].set_limit(iter_max);
				particles_[i].get_set().add_comparator(comparators::compare);
			}
		}
	public:
		mopso(){
			global_archive_.add_comparator(comparators::compare);
			pnum = particle_num;
			reset_particles();
		}
		~mopso(){}

		int step_once(){
			//auto start = high_resolution_clock::now();

			if (iter > iter_max || iter > iter_to) {
				return 0;
			}

			for (int i = 0; i < particle_num; i++) {
				auto temp = bind(&particle::update_pbest, &particles_[i]);
				thread_pool_.push(temp);
			}
			thread_pool_.wait();

			/*std::cout<<" "<<global_archive_.get_archive().size()<<" ";*/

			if (iter == (iter_max - 1) || iter == (iter_to - 1)) {
				for (int i = 0; i < particle_num; i++) {
					global_archive_.pareto_ops(particles_[i].get_pos());
				}
			}
			else
			{
				vector<vector<int>> bucket = {};
				bucket.resize(particle_num);

				for (int i = 0; i < particle_num; i++) {
					auto temp = bind(&pareto_set<primer_pair, function<int(primer_pair&, primer_pair&)>>::para_ops, &global_archive_, particles_[i].get_pos(), &bucket[i]);
					thread_pool_.push(temp);
				}
				thread_pool_.wait();

				for (int i = 1; i < particle_num; i++) {
					bucket[0].insert(bucket[0].end(), bucket[i].begin(), bucket[i].end());
				}
				sort(bucket[0].begin(), bucket[0].end());
				bucket[0].erase(unique(bucket[0].begin(), bucket[0].end()), bucket[0].end());

				vector<primer_pair>& h = global_archive_.get_archive();
				for (int i = bucket[0].size() - 1, j = 0; i >= 0;) {
					j = bucket[0][i];
					h.erase(h.begin() + j);
					i--;
				}
			}

			int gbest_index = get_index(global_archive_.get_archive().size(), iter, iter_max);

			for (int i = 0; i < particle_num; i++) {
				auto temp = bind(&particle::step_once, &particles_[i], global_archive_.get_archive(), gbest_index);
				thread_pool_.push(temp);
			}
			thread_pool_.wait();

			iter++;
			return 0;
		}

		int step_to(const int& i){
			iter_to = i;
			int per = 0;

			for (; (iter <= iter_to) && (iter <= iter_limit);) {
				if (iter % max((iter_to / 100), 1) == 0) {
					std::cout << per << "%" << std::endl;
					per++;
				}

				step_once();
			}

			return 0;
		}

		int step_to() {
			return step_to(iter_limit);
		}

		int reset_iter() {
			iter = 0;
			iter_limit = iter_max;
			iter_to = iter_limit;

			if (particle_num != pnum) {
				reset_particles();
				pnum = particle_num;
			}
			else {
				for (int i = 0; i < particle_num; i++) {
					particles_[i].set_iter(0);
					particles_[i].set_limit(iter_limit);
				}
			}

			return iter;
		}

		void clear_global() {
			global_archive_.get_archive().clear();
		}

		void clear_particles() {
			for (int i = 0; i < particle_num; i++) {
				particles_[i].get_set().get_archive().clear();
			}
		}

		int clear() {
			clear_global();
			clear_particles();
			return 0;
		}

		int reset() {
			reset_iter();
			clear();
			return iter;
		}

		vector<primer_pair>& get_result() {
			return global_archive_.get_archive();
		}
	};
}

#define mopso singleton<mopso::mopso>::get_instance()
