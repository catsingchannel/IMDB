#include <pybind11/pybind11.h>
#include "MPVProcess.h"

namespace py = pybind11;

PYBIND11_MODULE(pysnpprocess, m){
	m.doc() = "snp merge and annotation";

	py::class_<MpvProcess>(m, "MpvProcess")
	.def(py::init<const std::string&, const std::string&>())
	.def("run", &MpvProcess::run, py::call_guard<py::gil_scoped_release>())
	;
}
