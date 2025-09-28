#include "unidesign/SequenceHandle.hpp"

#ifdef UNIDESIGN_ENABLE_PYBIND11
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace {

unidesign::SequenceHandle clone_sequence(const unidesign::SequenceHandle& source) {
  return unidesign::SequenceHandle(source);
}

}  // namespace

PYBIND11_MODULE(_unidesign, m) {
  py::class_<unidesign::SequenceHandle>(m, "SequenceHandle")
      .def(py::init<>())
      .def("clone", &clone_sequence)
      .def("raw", [](unidesign::SequenceHandle& handle) { return handle.get(); },
           py::return_value_policy::reference);
}

#endif  // UNIDESIGN_ENABLE_PYBIND11

