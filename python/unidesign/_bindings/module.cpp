#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "Atom.h"
#include "AtomParamsSet.h"
#include "Chain.h"
#include "EnergyFunction.h"
#include "ProgramFunction.h"
#include "Residue.h"
#include "ResidueTopology.h"
#include "Structure.h"
#include "Utility.h"

namespace py = pybind11;

namespace {

template <typename Fn>
void check_status(Fn&& fn, const char* context) {
  int code = fn();
  if (FAILED(code)) {
    throw std::runtime_error(std::string(context) + " failed with code " + std::to_string(code));
  }
}

inline std::vector<char> make_buffer(const std::string& value) {
  std::vector<char> buffer(value.begin(), value.end());
  buffer.push_back('\0');
  return buffer;
}

struct AtomHandle {
  AtomHandle() { check_status([&]() { return AtomCreate(&value); }, "AtomCreate"); }
  AtomHandle(const AtomHandle& other) {
    check_status([&]() { return AtomCreate(&value); }, "AtomCreate");
    check_status([&]() { return AtomCopy(&value, const_cast<Atom*>(&other.value)); }, "AtomCopy");
  }
  AtomHandle& operator=(const AtomHandle& other) {
    if (this != &other) {
      check_status([&]() { return AtomCopy(&value, const_cast<Atom*>(&other.value)); }, "AtomCopy");
    }
    return *this;
  }
  ~AtomHandle() { AtomDestroy(&value); }

  Atom value{};
};

struct ResidueHandle {
  ResidueHandle() { check_status([&]() { return ResidueCreate(&value); }, "ResidueCreate"); }
  ResidueHandle(const ResidueHandle& other) {
    check_status([&]() { return ResidueCreate(&value); }, "ResidueCreate");
    check_status([&]() { return ResidueCopy(&value, const_cast<Residue*>(&other.value)); }, "ResidueCopy");
  }
  ResidueHandle& operator=(const ResidueHandle& other) {
    if (this != &other) {
      check_status([&]() { return ResidueCopy(&value, const_cast<Residue*>(&other.value)); }, "ResidueCopy");
    }
    return *this;
  }
  ~ResidueHandle() { ResidueDestroy(&value); }

  Residue value{};
};

struct ChainHandle {
  ChainHandle() { check_status([&]() { return ChainCreate(&value); }, "ChainCreate"); }
  ChainHandle(const ChainHandle& other) {
    check_status([&]() { return ChainCreate(&value); }, "ChainCreate");
    check_status([&]() { return ChainCopy(&value, const_cast<Chain*>(&other.value)); }, "ChainCopy");
  }
  ChainHandle& operator=(const ChainHandle& other) {
    if (this != &other) {
      check_status([&]() { return ChainCopy(&value, const_cast<Chain*>(&other.value)); }, "ChainCopy");
    }
    return *this;
  }
  ~ChainHandle() { ChainDestroy(&value); }

  Chain value{};
};

struct StructureHandle {
  StructureHandle() { check_status([&]() { return StructureCreate(&value); }, "StructureCreate"); }
  StructureHandle(const StructureHandle& other) {
    check_status([&]() { return StructureCreate(&value); }, "StructureCreate");
    check_status([&]() { return StructureCopy(&value, const_cast<Structure*>(&other.value)); }, "StructureCopy");
  }
  StructureHandle& operator=(const StructureHandle& other) {
    if (this != &other) {
      check_status([&]() { return StructureCopy(&value, const_cast<Structure*>(&other.value)); }, "StructureCopy");
    }
    return *this;
  }
  ~StructureHandle() { StructureDestroy(&value); }

  Structure value{};
};

std::vector<std::string> collect_atom_names(Residue* residue) {
  std::vector<std::string> result;
  int count = ResidueGetAtomCount(residue);
  result.reserve(static_cast<size_t>(count));
  for (int i = 0; i < count; ++i) {
    Atom* atom = ResidueGetAtom(residue, i);
    if (atom != nullptr) {
      result.emplace_back(AtomGetName(atom));
    }
  }
  return result;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Pybind11 bindings for the UniDesign native library (initial skeleton)";

  py::register_exception<std::runtime_error>(m, "NativeError");

  py::enum_<Type_Chain>(m, "ChainType")
      .value("PROTEIN", Type_Chain_Protein)
      .value("DNA", Type_Chain_DNA)
      .value("RNA", Type_Chain_RNA)
      .value("SMALL_MOLECULE", Type_Chain_SmallMol)
      .value("METAL_ION", Type_Chain_MetalIon)
      .value("WATER", Type_Chain_Water)
      .value("UNKNOWN", Type_Chain_Unknown);

  py::class_<AtomHandle>(m, "Atom")
      .def(py::init<>())
      .def_property(
          "name",
          [](AtomHandle& self) { return std::string(AtomGetName(&self.value)); },
          [](AtomHandle& self, const std::string& new_name) {
            auto buffer = make_buffer(new_name);
            check_status([&]() { return AtomSetName(&self.value, buffer.data()); }, "AtomSetName");
          })
      .def_property(
          "chain",
          [](AtomHandle& self) { return std::string(AtomGetChainName(&self.value)); },
          [](AtomHandle& self, const std::string& new_chain) {
            auto buffer = make_buffer(new_chain);
            check_status([&]() { return AtomSetChainName(&self.value, buffer.data()); }, "AtomSetChainName");
          })
      .def_property(
          "position",
          [](AtomHandle& self) { return AtomGetPosInChain(&self.value); },
          [](AtomHandle& self, int pos) {
            check_status([&]() { return AtomSetPosInChain(&self.value, pos); }, "AtomSetPosInChain");
          });

  py::class_<ResidueHandle>(m, "Residue")
      .def(py::init<>())
      .def_property(
          "name",
          [](ResidueHandle& self) { return std::string(ResidueGetName(&self.value)); },
          [](ResidueHandle& self, const std::string& new_name) {
            auto buffer = make_buffer(new_name);
            check_status([&]() { return ResidueSetName(&self.value, buffer.data()); }, "ResidueSetName");
          })
      .def_property(
          "chain",
          [](ResidueHandle& self) { return std::string(ResidueGetChainName(&self.value)); },
          [](ResidueHandle& self, const std::string& chain_name) {
            auto buffer = make_buffer(chain_name);
            check_status([&]() { return ResidueSetChainName(&self.value, buffer.data()); }, "ResidueSetChainName");
          })
      .def_property(
          "position",
          [](ResidueHandle& self) { return ResidueGetPosInChain(&self.value); },
          [](ResidueHandle& self, int pos) {
            check_status([&]() { return ResidueSetPosInChain(&self.value, pos); }, "ResidueSetPosInChain");
          })
      .def("atom_count", [](ResidueHandle& self) { return ResidueGetAtomCount(&self.value); })
      .def("atom_names", [](ResidueHandle& self) { return collect_atom_names(&self.value); });

  py::class_<ChainHandle>(m, "Chain")
      .def(py::init<>())
      .def_property(
          "name",
          [](ChainHandle& self) { return std::string(ChainGetName(&self.value)); },
          [](ChainHandle& self, const std::string& new_name) {
            auto buffer = make_buffer(new_name);
            check_status([&]() { return ChainSetName(&self.value, buffer.data()); }, "ChainSetName");
          })
      .def_property(
          "type",
          [](ChainHandle& self) { return ChainGetType(&self.value); },
          [](ChainHandle& self, Type_Chain type) {
            check_status([&]() { return ChainSetType(&self.value, type); }, "ChainSetType");
          })
      .def("residue_count", [](ChainHandle& self) { return ChainGetResidueCount(&self.value); })
      .def("append_residue", [](ChainHandle& self, ResidueHandle& residue) {
        check_status([&]() { return ChainAppendResidue(&self.value, &residue.value); }, "ChainAppendResidue");
      })
      .def("residue", [](ChainHandle& self, int index) {
        Residue* ptr = ChainGetResidue(&self.value, index);
        if (!ptr) {
          throw std::out_of_range("Residue index out of range");
        }
        auto copy = ResidueHandle();
        check_status([&]() { return ResidueCopy(&copy.value, ptr); }, "ResidueCopy");
        return copy;
      });

  py::class_<StructureHandle>(m, "Structure")
      .def(py::init<>())
      .def_property(
          "name",
          [](StructureHandle& self) { return std::string(StructureGetName(&self.value)); },
          [](StructureHandle& self, const std::string& new_name) {
            auto buffer = make_buffer(new_name);
            check_status([&]() { return StructureSetName(&self.value, buffer.data()); }, "StructureSetName");
          })
      .def("chain_count", [](StructureHandle& self) { return StructureGetChainCount(&self.value); })
      .def("add_chain", [](StructureHandle& self, ChainHandle& chain) {
        check_status([&]() { return StructureAddChain(&self.value, &chain.value); }, "StructureAddChain");
      })
      .def("chain", [](StructureHandle& self, int index) {
        Chain* ptr = StructureGetChain(&self.value, index);
        if (!ptr) {
          throw std::out_of_range("Chain index out of range");
        }
        auto copy = ChainHandle();
        check_status([&]() { return ChainCopy(&copy.value, ptr); }, "ChainCopy");
        return copy;
      })
      .def(
          "compute_stability",
          [](StructureHandle& self,
             const std::string& weight_file,
             const std::string& aapp_file,
             const std::string& rama_file) {
            double energy_terms[MAX_ENERGY_TERM];
            check_status([&]() { return EnergyTermInitialize(energy_terms); }, "EnergyTermInitialize");

            auto weight_buffer = make_buffer(weight_file);
            check_status([&]() { return EnergyWeightRead(weight_buffer.data()); }, "EnergyWeightRead");

            AAppTable aap_table{};
            auto aapp_buffer = make_buffer(aapp_file);
            check_status(
                [&]() { return AApropensityTableReadFromFile(&aap_table, aapp_buffer.data()); },
                "AApropensityTableReadFromFile");

            RamaTable rama_table{};
            auto rama_buffer = make_buffer(rama_file);
            check_status(
                [&]() { return RamaTableReadFromFile(&rama_table, rama_buffer.data()); },
                "RamaTableReadFromFile");

            check_status(
                [&]() { return ComputeStructureStabilitySilent(&self.value, &aap_table, &rama_table, energy_terms); },
                "ComputeStructureStabilitySilent");

            return py::cast(std::vector<double>(energy_terms, energy_terms + MAX_ENERGY_TERM));
          },
          py::arg("weight_file"),
          py::arg("aapp_file"),
          py::arg("rama_file"))
      .def(
          "read_pdb",
          [](StructureHandle& self,
             const std::string& pdb_path,
             const std::string& atom_param_path,
             const std::string& topo_path) {
            check_status([&]() { return StructureDestroy(&self.value); }, "StructureDestroy");
            check_status([&]() { return StructureCreate(&self.value); }, "StructureCreate");

            struct AtomParamGuard {
              AtomParamsSet value;
              AtomParamGuard() { check_status([&]() { return AtomParamsSetCreate(&value); }, "AtomParamsSetCreate"); }
              ~AtomParamGuard() { AtomParamsSetDestroy(&value); }
            } atom_params;

            struct ResiTopoGuard {
              ResiTopoSet value;
              ResiTopoGuard() { check_status([&]() { return ResiTopoSetCreate(&value); }, "ResiTopoSetCreate"); }
              ~ResiTopoGuard() { ResiTopoSetDestroy(&value); }
            } resi_topos;

            auto atom_param_buffer = make_buffer(atom_param_path);
            auto topo_buffer = make_buffer(topo_path);
            check_status(
                [&]() { return AtomParameterRead(&atom_params.value, atom_param_buffer.data()); }, "AtomParameterRead");
            check_status(
                [&]() { return ResiTopoSetRead(&resi_topos.value, topo_buffer.data()); }, "ResiTopoSetRead");

            auto pdb_buffer = make_buffer(pdb_path);
            check_status(
                [&]() {
                  return StructureReadPDB(&self.value, pdb_buffer.data(), &atom_params.value, &resi_topos.value);
                },
                "StructureReadPDB");
          },
          py::arg("pdb_path"),
          py::arg("atom_params"),
          py::arg("topology"));

  m.def(
      "print_version",
      []() {
        return PrintVersion();
      },
      "Proxy around the native PrintVersion helper");
}
