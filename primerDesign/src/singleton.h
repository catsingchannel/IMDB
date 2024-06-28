#pragma once

template<typename T>
class singleton {
private:
	singleton() = default;

public:
	singleton(const singleton&) = delete;
	singleton& operator=(const singleton&) = delete;
	~singleton() = default;

	static T& get_instance() {
		static T instance;
		return instance;
	}
};