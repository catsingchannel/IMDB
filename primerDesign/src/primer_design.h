// primer design.h: 标准系统包含文件的包含文件
// 或项目特定的包含文件。

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <iostream>
#include "probes.h"
#include "mopso.h"

namespace py = pybind11;

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::sort;
using std::unique;
using std::roundf;
using utils::get_coverage;
using structs::primer_pair;
using structs::assay;
using probe::design_probes;
using probe::get_probe;

#define assays singleton<std::vector<structs::assay>>::get_instance()

py::list primer_design_(const py::kwargs&);

PYBIND11_MODULE(primer_design, m) {
	m.def("primer_design", &primer_design_, py::return_value_policy::move, "A function design qPCR assays with given template");
}

// TODO: 在此处引用程序需要的其他标头。
