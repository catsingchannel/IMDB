#pragma once

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <bitset>
#include <cstring>
#include <string>
#include <vector>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <functional>

class Row {
private:
	std::string refpos_,
		refvar_,
		qvar_,
		qpos_,
		qlength_,
		rlength_,
		rname_,
		sample_,
		protein_,
		variant_,
		varclass_,
		annotation_;

public:
	Row(const std::string& s){
		protein_ = std::string();
		variant_ = std::string();
		varclass_ = std::string();
		annotation_ = std::string();

		uint32_t tabs = 0;
		uint32_t f = 0;
		uint32_t r = 0;

		for (uint32_t i = 0; i < s.length(); i++) {
			if (s[i] == '\t') {
				switch (tabs){
				case 0:
					r = i;
					refpos_ = std::string(s, f, r);
					f = i + 1;
					break;
				case 1:
					r = i - f;
					refvar_ = std::string(s, f, r);
					f = i + 1;
					break;
				case 2:
					r = i - f;
					qvar_ = std::string(s, f, r);
					f = i + 1;
					break;
				case 3:
					r = i - f;;
					qpos_ = std::string(s, f, r);
					f = i + 1;
					break;
				case 7:
					f = i + 1;
					break;
				case 8:
					r = i - f;
					rlength_ = std::string(s, f, r);
					f = i + 1;
					break;
				case 9:
					r = i - f;
					qlength_ = std::string(s, f, r);
					f = i + 1;
					break;
				case 11:
					f = i + 1;
					break;
				case 12:
					r = i - f;
					rname_ = std::string(s, f, r);
					f = i + 1;
					sample_ = std::string(s, f);
					break;
				default:
					break;
				}

				tabs++;
			}
		}
	}

	std::string& operator[](int index) {
		switch (index)
		{
		case 0 :
			return refpos_;
		case 1:
			return refvar_;
		case 2:
			return qvar_;
		case 3:
			return qpos_;
		case 4:
			return rlength_;
		case 5:
			return qlength_;
		case 6:
			return rname_;
		case 7:
			return sample_;
		case 8:
			return protein_;
		case 9:
			return variant_;
		case 10:
			return varclass_;
		case 11:
			return annotation_;
		default:
			break;
		}

		return sample_;
	}

	std::ostream& operator<<(std::ostream& os) {
		uint8_t i = 0;
		std::string output;
		while (i < 7) {
			output.append((*this)[i]);
			output.append("\t");
			i++;
		}
		output.append((*this)[7]);

		os << output << std::endl;
		return os;
	}
};

struct Gff
{
	uint32_t UTR5_ = 0;
	uint32_t UTR3_ = 0;
	std::bitset<192> strand_;
	std::vector<uint32_t> starts_;
	std::vector<uint32_t> ends_;
	std::vector<std::string> varclass_;
	std::vector<std::string> annotation_;
};

struct Group {
	std::vector<Row> snp;
	std::vector<Row> ins;
	std::vector<Row> del;

	Group() {
		snp = std::vector<Row>();
		ins = std::vector<Row>();
		del = std::vector<Row>();
	}
};

uint32_t NumConvert(const std::string& snum) {
	uint32_t rt = 0;
	for (int i = 0; i < snum.size(); i++) {
		rt = snum[i] - '0' + (rt * 10);
	}
	return rt;
}

class ThreadPool {
private:
	void worker(){
		while (running) {
			std::function<void()> task;
			std::unique_lock<std::mutex> latch(q_latch_);
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

	size_t thread_count_total = std::min(static_cast<unsigned int>(16), std::thread::hardware_concurrency());
	std::atomic<bool> running{false};
	std::atomic<bool> waiting{false};
	std::atomic<size_t> running_tasks{0};

public:
	ThreadPool(){
		threads = std::unique_ptr<std::thread[]>(new std::thread[4]);
		std::cout << "threads : " << thread_count_total << std::endl;
		running = true;
		for (int i = 0; i < thread_count_total; i++) {
			threads[i] = std::thread(&ThreadPool::worker, this);
		}
	}

	ThreadPool(const ThreadPool& t) = delete;
	ThreadPool(ThreadPool&& t) = delete;

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
		tasks_.push(task);
		q_latch_.unlock();

		cv_.notify_one();
	}

	void wait(){
		waiting = true;
		std::unique_lock<std::mutex> latch(q_latch_);
		done_.wait(latch, [this] {return running_tasks == 0 && tasks_.empty(); });
		waiting = false;
	}

	size_t get_threads_total() const{
		return thread_count_total;
	}
};

class Process {
private:
	std::vector<Row>* merge(std::vector<Row>* nucmer, char type) {
		if (nucmer->empty()) {
			return nucmer;
		}

		auto prev_it = nucmer->begin();
		uint32_t prevqpos = NumConvert((*nucmer)[0][3]);
		uint8_t corrector = 0;

		for (auto it = nucmer->begin(); it != nucmer->end();) {
			std::string r = (*it)[1];
			std::string q = (*it)[2];
			if ((r != "A" && r != "T" && r != "C" && r != "G" && r != ".") || (q != "A" && q != "T" && q != "C" && q != "G" && q != ".")) {
				it = nucmer->erase(it);

				if (it == nucmer->end()) {
					break;
				}

				prev_it = it;
				prevqpos = NumConvert((*it)[3]);
				corrector = 0;
				continue;
			}

			if (it == nucmer->begin()) {
				it++;
				continue;
			}
			
			uint32_t qpos = NumConvert(it->operator[](3));
			if (type != 'D') {
				if (qpos != 1 && qpos == (prevqpos + 1 + corrector)) {
					if (type == 'I') {
						prev_it->operator[](2).append(it->operator[](2));
					}else{
						prev_it->operator[](1).append(it->operator[](1));
						prev_it->operator[](2).append(it->operator[](2));
					}

					corrector++;
					it = nucmer->erase(it);
					continue;
				} else {
					prev_it = it;
					prevqpos = NumConvert(it->operator[](3));
					corrector = 0;
					it++;
				}
			}
			else {
				if (qpos != 1 && qpos == prevqpos) {
					prev_it->operator[](1).append(it->operator[](1));
					it = nucmer->erase(it);
					continue;
				} else {
					prev_it = it;
					prevqpos = NumConvert(it->operator[](3));
					it++;
				}
			}
		}

		return nucmer;
	}

	std::string reverse(const std::string& s) {
		std::string reverse;
		for (int i = 0; i < s.size(); i++) {
			switch (s[i])
			{
			case 'A':
				reverse.push_back('T');
				break;
			case 'T':
				reverse.push_back('A');
				break;
			case 'C':
				reverse.push_back('G');
				break;
			case 'G':
				reverse.push_back('C');
				break;
			default:
				reverse.push_back('*');
				break;
			}
		}

		std::reverse(reverse.begin(), reverse.end());
		return reverse;
	}

	std::string translate(const std::string& s) {
		std::string aas;

		for (uint32_t i = 0; i < s.size(); i += 3) {
			std::string temp(s, i, 3);
			if (temp == "GGT" || temp == "GGC" || temp == "GGA" || temp == "GGG") {
				aas.push_back('G');
			}
			else if (temp == "GCT" || temp == "GCC" || temp == "GCA" || temp == "GCG") {
				aas.push_back('A');
			}
			else if (temp == "GTT" || temp == "GTC" || temp == "GTA" || temp == "GTG") {
				aas.push_back('V');
			}
			else if (temp == "CTT" || temp == "CTC" || temp == "CTA" || temp == "CTG" || temp == "TTA" || temp == "TTG") {
				aas.push_back('L');
			}
			else if (temp == "ATT" || temp == "ATC" || temp == "ATA") {
				aas.push_back('I');
			}
			else if (temp == "CCT" || temp == "CCA" || temp == "CCG" || temp == "CCC") {
				aas.push_back('P');
			}
			else if (temp == "TTT" || temp == "TTC") {
				aas.push_back('F');
			}
			else if (temp == "TAT" || temp == "TAC") {
				aas.push_back('Y');
			}
			else if (temp == "TGG") {
				aas.push_back('W');
			}
			else if (temp == "TCT" || temp == "TCA" || temp == "TCC" || temp == "TCG" || temp == "AGT" || temp == "AGC") {
				aas.push_back('S');
			}
			else if (temp == "ACT" || temp == "ACC" || temp == "ACG" || temp == "ACA") {
				aas.push_back('T');
			}
			else if (temp == "TGT" || temp == "TGC") {
				aas.push_back('C');
			}
			else if (temp == "ATG") {
				aas.push_back('M');
			}
			else if (temp == "AAT" || temp == "AAC") {
				aas.push_back('N');
			}
			else if (temp == "CAA" || temp == "CAG") {
				aas.push_back('Q');
			}
			else if (temp == "GAT" || temp == "GAC") {
				aas.push_back('D');
			}
			else if (temp == "GAA" || temp == "GAG") {
				aas.push_back('E');
			}
			else if (temp == "AAA" || temp == "AAG") {
				aas.push_back('K');
			}
			else if (temp == "CGT" || temp == "CGC" || temp == "CGG" || temp == "CGA" || temp == "AGA" || temp == "AGG") {
				aas.push_back('R');
			}
			else if (temp == "CAT" || temp == "CAC") {
				aas.push_back('H');
			}
			else {
				aas.push_back('*');
			}
		}
		return aas;
	}

	std::vector<Row>* annotation(std::vector<Row>* nucmer, const Gff& annot_, const std::string& refseq) {
		if (nucmer->size() == 0) {
			return nucmer;
		}

		for (auto it = nucmer->begin(); it != nucmer->end(); it++) {
			//uint32_t rpos_ = NumConvert((*it)[0]);
			float rpos = static_cast<float>(NumConvert((*it)[0]));
			if (rpos < annot_.UTR5_) {
				(*it)[8] = "5'UTR";
				(*it)[9] = (*it)[0];
				(*it)[10] = "extragenic";
				(*it)[11] = "extragenic";
				continue;
			}
			else if (rpos > annot_.UTR3_) {
				(*it)[8] = "3'UTR";
				(*it)[9] = (*it)[0];
				(*it)[10] = "extragenic";
				(*it)[11] = "extragenic";
				continue;
			}

			for (int j = 0; j < annot_.ends_.size(); j++) {
				if (rpos <= annot_.ends_[j]) {
					if (rpos < annot_.starts_[j]) {
						(*it)[8] = "intergenic";
						(*it)[9] = (*it)[0];
						(*it)[10] = "extragenic";
						(*it)[11] = "extragenic";
						break;
					}
					else{
						(*it)[8] = annot_.varclass_[j];
						(*it)[11] = annot_.annotation_[j];
						float start = annot_.starts_[j];
						float end = annot_.ends_[j];
						std::string ref(refseq, static_cast<size_t>(start - 1), static_cast<size_t>(end - start + 1));
						std::string refpep;
						if (annot_.strand_[j] == false) {
							refpep = translate(reverse(ref));
						}
						else {
							refpep = translate(ref);
						}


						if ((*it)[2] == ".") {
							if (((*it)[1].size()) % 3 != 0) {
								uint32_t mutpos = std::round((rpos - start + 1) / 3);
								if (rpos - start + 1 < 2) {
									mutpos = 1;
								}

								if (annot_.strand_[j] == false) {
									mutpos = ((end - start + 1) / 3) - mutpos + 1;
								}
								mutpos--;

								(*it)[9].push_back(refpep[mutpos]);
								(*it)[9].append(std::to_string(mutpos + 1));
								(*it)[10] = "deletion_frameshift";
								break;
							}
							else {
								std::string var = ref;
								var.erase(static_cast<size_t>(rpos - start), (*it)[1].size());
								if (annot_.strand_[j] == false) {
									var = reverse(var);
								}
								std::string varpep = translate(var);

								for (int k = 0; k < varpep.size(); k++) {
									if (varpep[k] != refpep[k]) {
										(*it)[9].push_back(refpep[k]);
										(*it)[9].append(std::to_string(k + 1));

										if (varpep[k] == '*') {
											(*it)[10] = "deletion_stop";
										}
										else {
											(*it)[10] = "deletion";
										}
										break;
									}
								}
							}
						}
						else if ((*it)[1] == ".") {
							if (((*it)[2].size()) % 3 != 0) {
								uint32_t mutpos = std::round((rpos - start + 1) / 3);
								if (rpos - start + 1 < 2) {
									mutpos = 1;
								}

								if (annot_.strand_[j] == false) {
									mutpos = ((end - start + 1) / 3) - mutpos + 1;
								}
								mutpos--;

								(*it)[9].push_back(refpep[mutpos]);
								(*it)[9].append(std::to_string(mutpos + 1));
								(*it)[10] = "insertion_frameshift";
								break;
							}
							else {
								std::string var = ref;
								var.insert(static_cast<size_t>(rpos - start + 1), (*it)[2]);
								if (annot_.strand_[j] == false) {
									var = reverse(var);
								}
								std::string varpep = translate(var);

								for (int k = 0; k < refpep.size(); k++) {
									if (varpep[k] != refpep[k]) {
										uint32_t aa_inserted_ = ((*it)[2].size()) / 3;
										std::string pep_inserted_(varpep, k, aa_inserted_);
										(*it)[9] = pep_inserted_;
										(*it)[9].append(std::to_string(k + 1));

										for (int l = 0; l < pep_inserted_.size(); l++) {
											if (pep_inserted_[l] == '*') {
												(*it)[10] = "insertion_stop";
												break;
											}
											else if (l + 1 == pep_inserted_.size()) {
												(*it)[10] = "insertion";
											}
										}
										break;
									}
								}
							}
						}
						else {
							std::string var = ref;
							if ((*it)[2].size() == 1) {
								var[static_cast<size_t>(rpos - start)] = (*it)[2][0];
								if (annot_.strand_[j] == false) {
									var = reverse(var);
								}
								std::string varpep = translate(var);

								if (varpep == refpep) {
									uint32_t mutpos = std::round((rpos - start + 1) / 3);
									if (rpos - start + 1 < 2) {
										mutpos = 1;
									}

									if (annot_.strand_[j] == false) {
										mutpos = ((end - start + 1) / 3) - mutpos + 1;
									}
									mutpos--;

									(*it)[9].push_back(refpep[mutpos]);
									(*it)[9].append(std::to_string(mutpos + 1));
									(*it)[9].push_back(varpep[mutpos]);
									(*it)[10] = "SNP_silent";
								}
								else{
									for (int k = 0; k < varpep.size(); k++) {
										if (varpep[k] == '*') {
											(*it)[9].push_back(refpep[k]); 
											(*it)[9].append(std::to_string(k + 1));
											(*it)[9].push_back(varpep[k]);
											(*it)[10] = "SNP_STOP";
											break;
										}
										else if (varpep[k] != refpep[k]) {
											(*it)[9].push_back(refpep[k]);
											(*it)[9].append(std::to_string(k + 1));
											(*it)[9].push_back(varpep[k]);
											(*it)[10] = "SNP";
											break;
										}
									}
								}
							}
							else {
								for (int k = 0; k < (*it)[2].size(); k++) {
									var[static_cast<size_t>(rpos - start) + k] = (*it)[2][k];
								}

								if (annot_.strand_[j] == false) {
									var = reverse(var);
								}
								std::string varpep = translate(var);

								if (varpep == refpep) {
									uint32_t mutpos = std::round((rpos - start + 1) / 3);
									if (annot_.strand_[j] == false) {
										mutpos = ((end - start + 1) / 3) - mutpos + 1;
									}
									mutpos--;

									(*it)[9].push_back(refpep[mutpos]);
									(*it)[9].append(std::to_string(mutpos + 1));
									(*it)[9].push_back(varpep[mutpos]);
									(*it)[10] = "SNP_silent";
								}
								else {
									std::string varp;
									int pos = 0;
									for (int k = 0; k < varpep.size(); k++) {
										if (varpep[k] == '*') {
											if (k + 1 == varpep.size()) {
												break;
											}

											if ((*it)[9].size() != 0) {
												(*it)[9].clear();
											}

											(*it)[9].push_back(refpep[k]);
											(*it)[9].append(std::to_string(k + 1));
											(*it)[9].push_back(varpep[k]);
											(*it)[10] = "SNP_STOP";
											break;
										}
										else if (varpep[k] != refpep[k]) {
											(*it)[9].push_back(refpep[k]);
											varp.push_back(varpep[k]);

											if ((*it)[9].size() == 1) {
												pos = k + 1;
											}
										}
									}

									if ((*it)[10].empty()) {
										(*it)[9].append(std::to_string(pos));
										(*it)[9].append(varp);
										(*it)[10] = "SNP";
									}
								}
							}
						}

						break;
					}
				}
			}
		}

		return nucmer;
	}

public:
	void operator()(Group* sample, const Gff& annot_, const std::string& refseq_) {
		annotation(merge(&(sample->del), 'D'), annot_, refseq_);
		annotation(merge(&(sample->ins), 'I'), annot_, refseq_);
		annotation(merge(&(sample->snp), 'S'), annot_, refseq_);
	};
};

class FileIO {
private:
	char buf[4096];
public:
	uint8_t ReadSnpParallel(const std::string& path, std::vector<Row>* rows_, ThreadPool* pool_) {
		std::ifstream test_(path, std::ios::binary);

		if (!test_.is_open()) {
			std::cout << "cannot open snp file" << std::endl;
			return 1;
		}

		uint8_t blocks = pool_->get_threads_total();
		std::vector<std::vector<Row>> buckets(blocks, std::vector<Row>());
		test_.seekg(0, test_.end);
		size_t length = test_.tellg();
		test_.seekg(0, test_.beg);
		test_.close();
		size_t block_size = length / blocks;

		/*std::function<void(int)> f([block_size, &buckets, path, this](int index) {
			std::ifstream snp_(path, std::ios::binary);
			std::ifstream::pos_type bytes = 0;
			std::ifstream::pos_type buf_pos = index * 512;

			if (index == 0) {
				for (int i = 0; bytes < block_size && !snp_.eof(); i++) {
					snp_.getline(&(buf[buf_pos]), 512);

					if (i > 3) {
						if (std::strlen(buf) < 10) {
							continue;
						}

						buckets[index].push_back(Row(&(buf[buf_pos])));
					}
					bytes += snp_.tellg() - bytes;
				}
			}
			else {
				snp_.seekg(index * block_size);

				for (int i = 0; bytes < block_size && !snp_.eof(); i++) {
					snp_.getline(&(buf[buf_pos]), 512);

					if (i > 0) {
						if (std::strlen(buf) < 10) {
							continue;
						}

						buckets[index].push_back(Row(&(buf[buf_pos])));
					}
					bytes += snp_.tellg() - bytes;
				}
			}
			snp_.close();
			});

		for (int i = 0; i < blocks; i++) {
			f(i);
		}*/

		for (int i = 0; i < blocks; i++) {
			pool_->push([block_size, &buckets, path, this](int index) {
				std::ifstream snp_(path, std::ios::binary);
				std::ifstream::pos_type bytes = 0;
				size_t buf_pos = index * 256;

				if (index == 0) {
					for (int i = 0; bytes < block_size && !snp_.eof(); i++) {
						snp_.getline(&(buf[buf_pos]), 256);

						if (i > 3) {
							if (std::strlen(buf) < 10) {
								continue;
							}

							buckets[index].push_back(Row(&(buf[buf_pos])));
						}
						bytes += snp_.tellg() - bytes;
					}
				}
				else {
					snp_.seekg(index * block_size);

					for (int i = 0; bytes < block_size && !snp_.eof(); i++) {
						snp_.getline(&(buf[buf_pos]), 512);

						if (i > 0) {
							if (std::strlen(buf) < 10) {
								continue;
							}

							buckets[index].push_back(Row(&(buf[buf_pos])));
						}
						bytes += snp_.tellg() - bytes - (index * block_size);
					}
				}
				snp_.close();
				}, i);
		}

		pool_->wait();

		for (int i = 0; i < blocks; i++) {
			rows_->insert(rows_->end(), buckets[i].begin(), buckets[i].end());
		}
		return 0;
	}

	/*uint8_t ReadSnp(const std::string& path, std::vector<Row>* rows_) {
		std::ifstream snp_(path, std::ios::binary);

		if (!snp_.is_open()) {
			std::cout << "cannot open snp file" << std::endl;
			return 1;
		}

		for (int i = 0; !snp_.eof(); i++) {
			snp_.getline(buf, 1024);
			if (i > 3) {
				if (std::strlen(buf) < 10) {
					continue;
				}
				rows_->push_back(Row(buf));
			}
		}

		snp_.close();
		return 0;
	}*/

	uint8_t ReadMeta(const std::string& gffpath, const std::string& refseqpath, Gff* annot, std::string* refseq) {
		std::ifstream ref_(refseqpath, std::ios::binary);

		if (!ref_.is_open()) {
			std::cout << "cannot open ref file" << std::endl;
			return 1;
		}

		for (int i = 0; !ref_.eof(); i++) {
			ref_.getline(buf, 1024);
			if (i > 0) {
				std::string temp = buf;
				refseq->append(temp);
			}
		}

		refseq->erase(std::remove(refseq->begin(), refseq->end(), '\r'), refseq->end());
		ref_.close();
		std::ifstream gff_(gffpath, std::ios::binary);

		if (!gff_.is_open()) {
			std::cout << "cannot open gff file" << std::endl;
			return 1;
		}

		memset(buf, 0, 1024);
		for (int i = 0; !gff_.eof(); i++) {
			gff_.getline(buf, 1024);
			if (i > 0) {
				uint8_t tabs = 0;
				uint16_t f = 0;
				uint16_t r = 0;

				for (uint16_t j = 0; buf[j] != '\0'; j++) {
					if (buf[j] == '\t') {
						std::string buff = buf;
						switch (tabs)
						{
						case 3:
							f = j + 1;
							break;
						case 4:
							r = j - f;
							annot->starts_.push_back(NumConvert(std::string(buff, f, r)));
							f = j + 1;
							break;
						case 5:
							r = j - f;
							annot->ends_.push_back(NumConvert(std::string(buff, f, r)));
							break;
						case 7:
							if (buf[j - 2] == '-') {
								annot->strand_.set(i - 1, false);
							}
							else {
								annot->strand_.set(i - 1, true);
							}
							break;
						case 8:
							f = j + 1;
							break;
						case 9:
							r = j - f;
							annot->varclass_.push_back(std::string(buff, f, r));
							annot->annotation_.push_back(std::string(buff, j+1));
				            annot->varclass_[i - 1].pop_back();
							annot->annotation_[i - 1].pop_back();
							annot->varclass_[i - 1].erase(0, 1);
							annot->annotation_[i - 1].erase(0, 1);
							break;
						default:
							break;
						}
						tabs++;
					}
				}
			}
		}

		return 0;
	}

	uint8_t WriteCsv(const std::unordered_map<std::string, Group>& idx) {
		std::ofstream csv_("output.csv");
		csv_ << "sample,refpos,refvar,qvar,qpos,qlength,protein,variant,varclass,annotation\n";

		for (auto it = idx.begin(); it != idx.end(); it++) {
			std::vector<Row> all = it->second.del;
			all.insert(all.end(), it->second.ins.begin(), it->second.ins.end());
			all.insert(all.end(), it->second.snp.begin(), it->second.snp.end());

			for (auto all_it = all.begin(); all_it != all.end(); all_it++) {
				std::string row = (*all_it)[7];
				row.push_back(',');
				row.append((*all_it)[0]);
				row.push_back(',');
				row.append((*all_it)[1]);
				row.push_back(',');
				row.append((*all_it)[2]);
				row.push_back(',');
				row.append((*all_it)[3]);
				row.push_back(',');
				row.append((*all_it)[5]);
				row.push_back(',');
				row.append((*all_it)[8]);
				row.push_back(',');
				row.append((*all_it)[9]);
				row.push_back(',');
				row.append((*all_it)[10]);
				row.push_back(',');
				row.append((*all_it)[11]);
				row.append("\n");
				csv_ << row;
			}
		}

		csv_.close();
		return 0;
	}
};

class MpvProcess {
private:
	ThreadPool pool_ = {};
	std::function<void(Group*, const Gff&, const std::string&)> process_ = {Process()};
	FileIO file_io_ = {};
	Gff annot_ = {};
	std::string refseq;
	std::unordered_map<std::string, Group> ids;

	void partition(std::vector<Row>& nucmer) {
		for (int i = 0; i < nucmer.size(); i++) {
			auto it = ids.find(nucmer[i][7]);
			if (it == ids.end()) {
				it = ids.insert(std::make_pair(nucmer[i][7], Group())).first;
			}

			Row* temp = &(nucmer[i]);
			if (temp->operator[](1)[0] != '.' && temp->operator[](2)[0] != '.') {
				it->second.snp.push_back(nucmer[i]);
			}
			else if (temp->operator[](1)[0] == '.') {
				it->second.ins.push_back(nucmer[i]);
			}
			else if (temp->operator[](2)[0] == '.') {
				it->second.del.push_back(nucmer[i]);
			}
		}
	}
public:
	MpvProcess(): 
	ids(std::unordered_map<std::string, Group>()){
		std::cout<<"initializing..."<<std::endl;
		char *path = nullptr;
		path = get_current_dir_name();
		std::string gff(path);
		std::string ref(path);
		gff.append("/gff_.tsv");
		ref.append("/ref_.fasta");
		file_io_.ReadMeta(gff, ref, &annot_, &refseq);
		annot_.UTR5_ = annot_.starts_[0];
		annot_.UTR3_ = annot_.ends_[annot_.ends_.size() - 1];
		std::cout<<"done"<<std::endl;
	}

	MpvProcess(const std::string& gff, const std::string& ref) :
	ids(std::unordered_map<std::string, Group>()) {
		std::cout << "initializing..." << std::endl;
		file_io_.ReadMeta(gff, ref, &annot_, &refseq);
		annot_.UTR5_ = annot_.starts_[0];
		annot_.UTR3_ = annot_.ends_[annot_.ends_.size() - 1];
		std::cout << "done" << std::endl;
	}

	void operator()(const std::string& path){
		std::vector<Row> reads_;
		std::cout<<"reading snps"<<std::endl;
		file_io_.ReadSnpParallel(path, &reads_, &pool_);
		partition(reads_);
		
		std::cout<<"start merge and annotation"<<std::endl;

		for (auto it = ids.begin(); it != ids.end();) {
			Group* nucmer = &(it->second);
			pool_.push(process_, nucmer, annot_, refseq);
			it++;
		}

		pool_.wait();
		std::cout<<"merge and annotation done, generating csv file"<<std::endl;
		file_io_.WriteCsv(ids);
		std::cout<<"done"<<std::endl;
	}

	void run(const std::string& path){
		operator()(path);
	}
};
