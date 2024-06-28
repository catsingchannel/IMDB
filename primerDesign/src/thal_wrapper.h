#pragma once

#include <mutex>
#include <string>

#include "singleton.h"
#include "thal_parameters.h"
#include "thal.h"

namespace thal_wrapper {
	using std::mutex;
	using std::string;

	class thal_wrapper {
	private :
		mutex thal_global = {};

		thal_parameters param;
		thal_args a;
		thal_results r;
		thal_mode mode = THL_FAST;

		bool initial = false;
	public:
		thal_wrapper() = default;
		~thal_wrapper() = default;

		void init() {
			thal_set_null_parameters(&param);
			set_default_thal_parameters(&param);
			set_thal_default_args(&a);
			r.sec_struct = NULL;
			get_thermodynamic_values(&param, &r);
			initial = true;
		}

		explicit operator bool() const {
			return initial;
		}

		float self_dimer_f(const string& s) {
			const unsigned char* ls = reinterpret_cast<const unsigned char*>(s.c_str());
			thal_global.lock();
			thal(ls, ls, &a, mode, &r);
			r.sec_struct = NULL;
			float res = r.temp;
			thal_global.unlock();

			if (res <= 0) {
				return 0;
			}
			else {
				return res;
			}
		}

		float cross_dimer_f(const string& s1, const string& s2) {
			const unsigned char* ls1 = reinterpret_cast<const unsigned char*>(s1.c_str());
			const unsigned char* ls2 = reinterpret_cast<const unsigned char*>(s2.c_str());
			thal_global.lock();
			thal(ls1, ls2, &a, mode, &r);
			r.sec_struct = NULL;
			float res = r.temp;
			thal_global.unlock();

			if (res <= 0) {
				return 0;
			}
			else {
				return res;
			}
		}

		float hairpin_f(const string& s) {
			const unsigned char* ls = reinterpret_cast<const unsigned char*>(s.c_str());
			thal_global.lock();
			a.dimer = 0;
			a.type = thal_hairpin;
			thal(ls, ls, &a, mode, &r);
			float res = r.temp;
			a.dimer = 1;
			a.type = thal_any;
			r.sec_struct = NULL;
			thal_global.unlock();

			if (res <= 0) {
				return 0;
			}
			else {
				return res;
			}
		}
	};
}

#define thal_wrapper singleton<thal_wrapper::thal_wrapper>::get_instance()
#define self_dimer thal_wrapper.self_dimer_f
#define cross_dimer thal_wrapper.cross_dimer_f
#define hairpin thal_wrapper.hairpin_f