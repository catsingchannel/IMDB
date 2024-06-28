#pragma once

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <functional>
#include <vector>
#include <string>

#include "singleton.h"

namespace thread_pool {
	static size_t actually_thread_num = std::thread::hardware_concurrency();

	class ThreadPool {
	private:
		void worker(int i) {
			while (running) {
				std::function<void()> task;
				std::unique_lock<std::mutex> latch(q_latch_);
				if (waiting) {
					done_.notify_one();
				}
				cv_.wait(latch, [this] {return !tasks_.empty() || !running; });
				if (running) {
					task = std::move(tasks_.front());
					tasks_.pop();
					latch.unlock();
					running_tasks++;
					task();
					running_tasks--;
					if (waiting) {
						done_.notify_one();
					}
				}
			}
		}

		std::queue<std::function<void()>> tasks_ = {};
		std::mutex q_latch_ = {};
		std::condition_variable cv_ = {};
		std::condition_variable done_ = {};
		std::unique_ptr<std::thread[]> threads = nullptr;

		size_t thread_count_total = actually_thread_num;
		std::atomic<bool> running = {false};
		std::atomic<bool> waiting = {false};
		std::atomic<size_t> running_tasks = {0};
	public:
		ThreadPool() {
			threads = std::unique_ptr<std::thread[]>(new std::thread[thread_count_total]);
			running = true;
			for (int i = 0; i < thread_count_total; i++) {
				threads[i] = std::thread(&ThreadPool::worker, this, i);
			}
		}

		ThreadPool(const ThreadPool& t) = delete;
		ThreadPool& operator=(const ThreadPool& t) = delete;

		~ThreadPool() {
			running = false;
			cv_.notify_all();
			for (int i = 0; i < thread_count_total; i++) {
				threads[i].join();
			}
		};

		template<typename F, typename... A>
		void push(F&& f, A&&... args) {
			std::function<void()> task = std::bind(std::forward<F>(f), std::forward<A>(args)...);

			q_latch_.lock();
			//log_.push_back("main pushing");
			tasks_.push(task);
			//log_.push_back("main exit pushing");
			q_latch_.unlock();

			cv_.notify_one();
		}

		void wait() {
			waiting = true;
			std::unique_lock<std::mutex> latch(q_latch_);
			//log_.push_back("main waiting");
			done_.wait(latch, [this] {return running_tasks == 0 && tasks_.empty(); });
			waiting = false;
			//log_.push_back("main exit waiting");
			//for (auto it = worker_logs_.begin(); it != worker_logs_.end(); it++) {
			//	it->clear();
			//}
		}

		size_t get_threads_total() const {
			return thread_count_total;
		}
	};
}

#define thread_pool_ singleton<thread_pool::ThreadPool>::get_instance()
