#pragma once

#include <vector>
#include <mutex>
#include <functional>

namespace pareto {
	using std::vector;
	using std::mutex;
	using std::function;
	using std::bind;
	using std::forward;
	using std::cerr;
	using std::endl;

	template<typename T, typename Cmp>
	class pareto_set {
	private:
		vector<T> archive_ = {};
		vector<Cmp> comparators_ = {};
		mutex archive_lock_ = {};

	public:
		pareto_set() {};

		pareto_set(const pareto_set& p) {
			vector<T> cp = p.get_archive();
			archive_lock_.lock();
			archive_.insert(archive_.end(), cp.begin(), cp.end());
			archive_lock_.unlock();
			add_comparator(p.comparators_);
		};

		pareto_set& operator=(const pareto_set& p) {
			vector<T> cp = p.get_archive();
			archive_lock_.lock();
			archive_.insert(archive_.end(), cp.begin(), cp.end());
			archive_lock_.unlock();
			add_comparator(p.comparators_);
			return *this;
		};

		~pareto_set() {};

		int add_comparator(const Cmp& f) {
			comparators_.push_back(f);
			return 0;
		}

		int add_comparator(vector<Cmp>& v) {
			comparators_.insert(comparators_.end(), v.begin(), v.end());
			return 0;
		}

		int compare(T& t1, T& t2) {
			int less = 0, eq = 0, gt = 0;

			for (auto it = comparators_.begin(); it != comparators_.end(); it++) {
				switch ((*it)(t1, t2))
				{
				case 0:
					less++;
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

			if (less == 0 && gt >= 1) {
				return 2;
			}
			else if (gt == 0 && less >= 1) {
				return 0;
			}
			else {
				return 1;
			}
		}

		int pareto_ops(T& t) {
			if (comparators_.empty()) {
				cerr << "comparators list is empty" << endl;
				return -1;
			}

			bool flag = true;
			archive_lock_.lock();
			if (archive_.empty()) {
				archive_.push_back(t);
				archive_lock_.unlock();
				return 0;
			}

			for (auto it = archive_.begin(); it != archive_.end(); it++) {
				if ((*it) == t) {
					archive_lock_.unlock();
					return 0;
				}
			}

			for (auto it = archive_.begin(); it != archive_.end();) {
				switch (compare(t, (*it))) {
				case 0:
					it++;
					flag = false;
					continue;
				case 1:
					it++;
					continue;
				case 2:
					it = archive_.erase(it);
					continue;
				default:
					break;
				}
			}
			
			if (flag) {
				archive_.push_back(t);
			}
			archive_lock_.unlock();

			return 0;
		}

		void para_ops(T& t, vector<int>* r) {
			archive_lock_.lock();
			vector<T> local = archive_;
			archive_lock_.unlock();

			bool flag = true;

			for (int i = 0; i < local.size(); i++) {
				if (local[i] == t) {
					return;
				}
			}

			for (int i = 0; i < local.size(); i++) {
				switch (compare(t, local[i])) {
				case 0:
					flag = false;
					continue;
				case 1:
					continue;
				case 2:
					r->push_back(i);
					continue;
				default:
					break;
				}
			}

			if (flag) {
				archive_lock_.lock();
				archive_.push_back(t);
				archive_lock_.unlock();
			}
		}

		vector<T>& get_archive() {
			return archive_;
		}
	};
}
