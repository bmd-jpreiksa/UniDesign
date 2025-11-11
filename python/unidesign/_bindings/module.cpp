#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cctype>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
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
#include "PyDesign.h"
#include "SmallMol.h"
#include "SmallMolEEF1.h"
#include "SmallMolParAndTopo.h"
#include "RotamerBuilder.h"

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

Type_ResidueDesignType parse_design_type(const py::handle& handle,
                                         Type_ResidueDesignType fallback) {
  if (handle.is_none()) {
    return fallback;
  }
  std::string raw = py::cast<std::string>(handle);
  std::string upper;
  upper.reserve(raw.size());
  for (char ch : raw) {
    upper.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
  }
  if (upper == "FIXED") {
    return Type_DesType_Fixed;
  }
  if (upper == "MUTABLE") {
    return Type_DesType_Mutable;
  }
  if (upper == "REPACKABLE") {
    return Type_DesType_Repackable;
  }
  if (upper == "SMALL_MOLECULE" || upper == "SMALLMOLECULE") {
    return Type_DesType_SmallMol;
  }
  if (upper == "CATALYTIC") {
    return Type_DesType_Catalytic;
  }
  if (upper == "NATROT") {
    return Type_DesType_NatRot;
  }
  return fallback;
}

PyDesignSiteSpec parse_site_spec(const py::handle& obj, Type_ResidueDesignType fallback) {
  PyDesignSiteSpec spec;
  py::dict mapping = py::cast<py::dict>(obj);
  if (!mapping.contains("chain") || !mapping.contains("position")) {
    throw std::invalid_argument("design/repack site requires 'chain' and 'position'");
  }
  spec.chain = py::cast<std::string>(mapping["chain"]);
  spec.position = py::cast<int>(mapping["position"]);
  spec.allowed = mapping.contains("allowed") && !mapping["allowed"].is_none()
                     ? py::cast<std::string>(mapping["allowed"])
                     : std::string();
  if (mapping.contains("design_type")) {
    spec.design_type = parse_design_type(mapping["design_type"], fallback);
  } else {
    spec.design_type = fallback;
  }
  return spec;
}

PyMonomerDesignOptions parse_monomer_options(const py::dict& kwargs) {
  auto require_string = [&](const char* key) -> std::string {
    if (!kwargs.contains(key)) {
      throw std::invalid_argument(std::string("Missing required option '") + key + "'");
    }
    return py::cast<std::string>(kwargs[key]);
  };

  auto optional_string = [&](const char* key, const std::string& fallback) -> std::string {
    if (!kwargs.contains(key) || kwargs[key].is_none()) {
      return fallback;
    }
    return py::cast<std::string>(kwargs[key]);
  };

  PyMonomerDesignOptions opts{};
  opts.program_path = optional_string("program_path", ".");
  opts.working_directory = optional_string("working_directory", opts.program_path);
  opts.atom_params_path = require_string("atom_params");
  opts.topology_path = require_string("topology");
  opts.weight_file = require_string("weight_file");
  opts.aapp_file = require_string("aapp_file");
  opts.rama_file = require_string("rama_file");
  opts.rotlib_bin = require_string("rotlib_bin");
  opts.resfile_contents = optional_string("resfile_text", "");
  opts.design_chains = optional_string("design_chains", "");

  opts.profile_weight =
      kwargs.contains("profile_weight") ? py::cast<double>(kwargs["profile_weight"]) : 1.0;
  opts.binding_weight =
      kwargs.contains("binding_weight") ? py::cast<double>(kwargs["binding_weight"]) : 1.0;
  opts.trajectories =
      kwargs.contains("trajectories") ? py::cast<int>(kwargs["trajectories"]) : 1;

  opts.interface_only =
      kwargs.contains("interface_only") ? py::cast<bool>(kwargs["interface_only"]) : false;
  opts.design_from_native =
      kwargs.contains("design_from_native") ? py::cast<bool>(kwargs["design_from_native"]) : false;
  opts.use_input_sc =
      kwargs.contains("use_input_sc") ? py::cast<bool>(kwargs["use_input_sc"]) : true;
  opts.rotate_hydroxyl =
      kwargs.contains("rotate_hydroxyl") ? py::cast<bool>(kwargs["rotate_hydroxyl"]) : true;
  opts.exclude_cys_rotamers = kwargs.contains("exclude_cys_rotamers")
                                  ? py::cast<bool>(kwargs["exclude_cys_rotamers"])
                                  : false;
  opts.wildtype_only =
      kwargs.contains("wildtype_only") ? py::cast<bool>(kwargs["wildtype_only"]) : false;
  opts.quiet_output =
      kwargs.contains("quiet") ? py::cast<bool>(kwargs["quiet"]) : false;

  opts.has_ligand = kwargs.contains("has_ligand");
  if (opts.has_ligand) {
    opts.ligand_mol2 = optional_string("ligand_mol2", "");
    opts.ligand_parameters = optional_string("ligand_params", "");
    opts.ligand_topology = optional_string("ligand_topology", "");
    opts.ligand_conformers = optional_string("ligand_conformers", "");
  }

  opts.design_sites.clear();
  if (kwargs.contains("design_sites") && !kwargs["design_sites"].is_none()) {
    for (const py::handle& entry : py::cast<py::iterable>(kwargs["design_sites"])) {
      opts.design_sites.push_back(parse_site_spec(entry, Type_DesType_Mutable));
    }
  }

  opts.repack_sites.clear();
  if (kwargs.contains("repack_sites") && !kwargs["repack_sites"].is_none()) {
    for (const py::handle& entry : py::cast<py::iterable>(kwargs["repack_sites"])) {
      opts.repack_sites.push_back(parse_site_spec(entry, Type_DesType_Repackable));
    }
  }

  return opts;
}

PyMinimizeOptions parse_minimize_options(const py::dict& kwargs) {
  auto require_string = [&](const char* key) -> std::string {
    if (!kwargs.contains(key)) {
      throw std::invalid_argument(std::string("Missing required option '") + key + "'");
    }
    return py::cast<std::string>(kwargs[key]);
  };

  auto optional_string = [&](const char* key, const std::string& fallback) -> std::string {
    if (!kwargs.contains(key) || kwargs[key].is_none()) {
      return fallback;
    }
    return py::cast<std::string>(kwargs[key]);
  };

  PyMinimizeOptions opts{};
  opts.program_path = optional_string("program_path", ".");
  opts.working_directory = optional_string("working_directory", opts.program_path);
  opts.atom_params_path = require_string("atom_params");
  opts.topology_path = require_string("topology");
  opts.rotlib_bin = require_string("rotlib_bin");
  opts.weight_file = require_string("weight_file");
  opts.use_input_sc =
      kwargs.contains("use_input_sc") ? py::cast<bool>(kwargs["use_input_sc"]) : true;
  opts.rotate_hydroxyl =
      kwargs.contains("rotate_hydroxyl") ? py::cast<bool>(kwargs["rotate_hydroxyl"]) : true;
  opts.quiet_output =
      kwargs.contains("quiet") ? py::cast<bool>(kwargs["quiet"]) : false;
  opts.respect_design_types = kwargs.contains("respect_design_types")
                                  ? py::cast<bool>(kwargs["respect_design_types"])
                                  : false;

  opts.has_ligand = kwargs.contains("has_ligand");
  if (opts.has_ligand) {
    opts.ligand_mol2 = optional_string("ligand_mol2", "");
    opts.ligand_parameters = optional_string("ligand_params", "");
    opts.ligand_topology = optional_string("ligand_topology", "");
    if (opts.ligand_mol2.empty() || opts.ligand_parameters.empty() || opts.ligand_topology.empty()) {
      throw std::invalid_argument(
          "Ligand minimization requires 'ligand_mol2', 'ligand_params', and 'ligand_topology'");
    }
  }

  opts.design_sites.clear();
  if (kwargs.contains("design_sites") && !kwargs["design_sites"].is_none()) {
    for (const py::handle& entry : py::cast<py::iterable>(kwargs["design_sites"])) {
      opts.design_sites.push_back(parse_site_spec(entry, Type_DesType_Mutable));
    }
  }

  opts.repack_sites.clear();
  if (kwargs.contains("repack_sites") && !kwargs["repack_sites"].is_none()) {
    for (const py::handle& entry : py::cast<py::iterable>(kwargs["repack_sites"])) {
      opts.repack_sites.push_back(parse_site_spec(entry, Type_DesType_Repackable));
    }
  }

  return opts;
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

void reset_residue_energy_terms(Structure* structure) {
  for (int chain_index = 0; chain_index < StructureGetChainCount(structure); ++chain_index) {
    Chain* chain = StructureGetChain(structure, chain_index);
    for (int res_index = 0; res_index < ChainGetResidueCount(chain); ++res_index) {
      Residue* residue = ChainGetResidue(chain, res_index);
      residue->aapp = 0.0;
      residue->rama = 0.0;
      ResidueSetDunbrack(residue, 0.0);
    }
  }
}

std::unordered_set<std::string> collect_chain_names(Structure* structure) {
  std::unordered_set<std::string> result;
  for (int i = 0; i < StructureGetChainCount(structure); ++i) {
    Chain* chain = StructureGetChain(structure, i);
    if (!chain) {
      continue;
    }
    const char* name = ChainGetName(chain);
    if (name != nullptr && name[0] != '\0') {
      result.emplace(name);
    }
  }
  return result;
}

void validate_chain_group(
    const std::string& group, const std::unordered_set<std::string>& available, const char* label) {
  if (group.empty()) {
    throw std::invalid_argument(std::string("Chain split '") + label + "' is empty");
  }

  std::string token;
  bool has_any = false;
  auto flush_token = [&]() {
    if (token.empty()) {
      return;
    }
    has_any = true;
    if (!available.count(token)) {
      bool ok = true;
      for (char ch : token) {
        if (std::isspace(static_cast<unsigned char>(ch)) || ch == ',') {
          continue;
        }
        std::string unit(1, ch);
        if (!available.count(unit)) {
          ok = false;
          break;
        }
      }
      if (!ok) {
        throw std::invalid_argument(
            std::string("Unknown chain identifier '") + token + "' in " + label);
      }
    }
    token.clear();
  };

  for (char ch : group) {
    if (ch == ',' || std::isspace(static_cast<unsigned char>(ch))) {
      flush_token();
    } else {
      token.push_back(ch);
    }
  }
  flush_token();

  if (!has_any) {
    throw std::invalid_argument(std::string("Chain split '") + label + "' does not reference any chains");
  }
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Pybind11 bindings for the UniDesign native library";

  py::register_exception<std::runtime_error>(m, "NativeError");

  py::enum_<Type_ResidueDesignType>(m, "ResidueDesignType")
      .value("FIXED", Type_DesType_Fixed)
      .value("MUTABLE", Type_DesType_Mutable)
      .value("REPACKABLE", Type_DesType_Repackable)
      .value("SMALL_MOLECULE", Type_DesType_SmallMol)
      .value("CATALYTIC", Type_DesType_Catalytic)
      .value("NATROT", Type_DesType_NatRot);

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
          })
      .def_property(
          "coords",
          [](AtomHandle& self) {
            return py::make_tuple(self.value.xyz.X, self.value.xyz.Y, self.value.xyz.Z);
          },
          [](AtomHandle& self, const py::tuple& coords) {
            if (coords.size() != 3) {
              throw std::invalid_argument("coords must be a 3-tuple");
            }
            self.value.xyz.X = coords[0].cast<double>();
            self.value.xyz.Y = coords[1].cast<double>();
            self.value.xyz.Z = coords[2].cast<double>();
            self.value.isXyzValid = TRUE;
          })
      .def_property(
          "bfactor",
          [](AtomHandle& self) { return self.value.bfactor; },
          [](AtomHandle& self, double value) { self.value.bfactor = value; });

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
      .def("atom_names", [](ResidueHandle& self) { return collect_atom_names(&self.value); })
      .def_property(
          "design_type",
          [](ResidueHandle& self) {
            return static_cast<Type_ResidueDesignType>(ResidueGetDesignType(&self.value));
          },
          [](ResidueHandle& self, Type_ResidueDesignType design_type) {
            check_status([&]() { return ResidueSetDesignType(&self.value, design_type); }, "ResidueSetDesignType");
          })
      .def("copy_from", [](ResidueHandle& self, ResidueHandle& other) {
        check_status([&]() { return ResidueCopy(&self.value, &other.value); }, "ResidueCopy");
      })
      .def("add_atom", [](ResidueHandle& self, AtomHandle& atom) {
        check_status([&]() { return ResidueAddAtom(&self.value, &atom.value); }, "ResidueAddAtom");
      })
      .def("atoms", [](ResidueHandle& self) {
        std::vector<AtomHandle> atoms;
        int count = ResidueGetAtomCount(&self.value);
        atoms.reserve(count);
        for (int i = 0; i < count; ++i) {
          Atom* ptr = ResidueGetAtom(&self.value, i);
          if (ptr == nullptr) {
            continue;
          }
          AtomHandle copy;
          check_status([&]() { return AtomCopy(&copy.value, ptr); }, "AtomCopy");
          atoms.push_back(copy);
        }
        return atoms;
      });

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
        ResidueHandle copy;
        check_status([&]() { return ResidueCopy(&copy.value, ptr); }, "ResidueCopy");
        return copy;
      })
      .def("copy_from", [](ChainHandle& self, ChainHandle& other) {
        check_status([&]() { return ChainCopy(&self.value, &other.value); }, "ChainCopy");
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
        ChainHandle copy;
        check_status([&]() { return ChainCopy(&copy.value, ptr); }, "ChainCopy");
        return copy;
      })
      .def("copy_from", [](StructureHandle& self, StructureHandle& other) {
        check_status([&]() { return StructureCopy(&self.value, &other.value); }, "StructureCopy");
      })
      .def(
          "compute_binding",
          [](StructureHandle& self,
             const std::string& weight_file,
             const py::object& split1,
             const py::object& split2) {
            if (StructureGetChainCount(&self.value) < 2) {
              throw std::invalid_argument("compute_binding requires at least two chains");
            }

            if (!weight_file.empty()) {
              auto weight_buffer = make_buffer(weight_file);
              check_status([&]() { return EnergyWeightRead(weight_buffer.data()); }, "EnergyWeightRead");
            }

            if (!split1.is_none() && !split2.is_none()) {
              auto part1 = split1.cast<std::string>();
              auto part2 = split2.cast<std::string>();
              auto available = collect_chain_names(&self.value);
              validate_chain_group(part1, available, "split1");
              validate_chain_group(part2, available, "split2");
              auto buffer1 = make_buffer(part1);
              auto buffer2 = make_buffer(part2);
              check_status(
                  [&]() { return ComputeBindingWithChainSplitting(&self.value, buffer1.data(), buffer2.data()); },
                  "ComputeBindingWithChainSplitting");
            } else if (split1.is_none() && split2.is_none()) {
              check_status([&]() { return ComputeBinding(&self.value); }, "ComputeBinding");
            } else {
              throw std::invalid_argument("split1 and split2 must both be provided or both omitted");
            }
          },
          py::arg("weight_file") = std::string{},
          py::arg("split1") = py::none(),
          py::arg("split2") = py::none())
      .def(
          "compute_stability",
          [](StructureHandle& self,
             const std::string& weight_file,
             const std::string& aapp_file,
             const std::string& rama_file,
             const std::string& rotlib_bin) {
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

            reset_residue_energy_terms(&self.value);

            if (!rotlib_bin.empty()) {
              auto rotlib_buffer = make_buffer(rotlib_bin);
              struct RotlibGuard {
                BBdepRotamerLib value;
                explicit RotlibGuard(char* path) {
                  check_status([&]() { return BBdepRotamerLibCreate2(&value, path); }, "BBdepRotamerLibCreate2");
                }
                ~RotlibGuard() { BBdepRotamerLibDestroy(&value); }
              } rotamer_lib(rotlib_buffer.data());

              check_status(
                  [&]() {
                    return ComputeStructureStabilityByBBdepRotLib2(
                        &self.value, &aap_table, &rama_table, rotlib_buffer.data(), energy_terms);
                  },
                  "ComputeStructureStabilityByBBdepRotLib2");
            } else {
              check_status(
                  [&]() { return ComputeStructureStabilitySilent(&self.value, &aap_table, &rama_table, energy_terms); },
                  "ComputeStructureStabilitySilent");
            }

            return py::cast(std::vector<double>(energy_terms, energy_terms + MAX_ENERGY_TERM));
          },
          py::arg("weight_file"),
          py::arg("aapp_file"),
          py::arg("rama_file"),
          py::arg("rotlib_bin") = std::string{})
      .def("calc_phi_psi", [](StructureHandle& self) {
        check_status([&]() { return StructureCalcPhiPsi(&self.value); }, "StructureCalcPhiPsi");
      })
      .def(
          "calc_propensity",
          [](StructureHandle& self, const std::string& aapp_file, const std::string& rama_file) {
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
                [&]() { return StructureCalcAminoAcidPropensityAndRamaEnergy(&self.value, &aap_table, &rama_table); },
                "StructureCalcAminoAcidPropensityAndRamaEnergy");
          },
          py::arg("aapp_file"),
          py::arg("rama_file"))
      .def(
          "calc_dunbrack",
          [](StructureHandle& self, const std::string& rotlib_bin) {
            auto rotlib_buffer = make_buffer(rotlib_bin);
            struct RotlibGuard {
              BBdepRotamerLib value;
              explicit RotlibGuard(char* path) {
                check_status([&]() { return BBdepRotamerLibCreate2(&value, path); }, "BBdepRotamerLibCreate2");
              }
              ~RotlibGuard() { BBdepRotamerLibDestroy(&value); }
            } rotamer_lib(rotlib_buffer.data());

            check_status(
                [&]() { return StructureCalcAminoAcidDunbrackEnergy(&self.value, &rotamer_lib.value); },
                "StructureCalcAminoAcidDunbrackEnergy");
          },
          py::arg("rotlib_bin"))
      .def(
          "reset_energy_terms",
          [](StructureHandle& self) {
            reset_residue_energy_terms(&self.value);
          })
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
            check_status([&]() { return StructureCalcPhiPsi(&self.value); }, "StructureCalcPhiPsi");
          },
          py::arg("pdb_path"),
          py::arg("atom_params"),
          py::arg("topology"))
      .def(
          "write_pdb",
          [](StructureHandle& self, const std::string& output_path) {
            FILE* file = fopen(output_path.c_str(), "w");
            if (file == nullptr) {
              throw std::runtime_error("Failed to open PDB path for writing: " + output_path);
            }
            std::unique_ptr<FILE, decltype(&fclose)> guard(file, fclose);
            check_status([&]() { return StructureShowInPDBFormat(&self.value, file); }, "StructureShowInPDBFormat");
          },
          py::arg("output_path"))
      .def(
          "read_mol2",
          [](StructureHandle& self,
             const std::string& mol2_path,
             const std::string& atom_param_path,
             const std::string& topo_path) {
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

            auto mol2_buffer = make_buffer(mol2_path);
            check_status(
                [&]() {
                  return StructureReadMol2(&self.value, mol2_buffer.data(), &atom_params.value, &resi_topos.value);
                },                "StructureReadMol2");
          },
          py::arg("mol2_path"),
          py::arg("atom_params"),
          py::arg("topology"))
      .def(
          "read_small_mol_rotamers",
          [](StructureHandle& self, const std::string& rotamer_path) {
            ResiTopoSet resi_topos;
            ResiTopoSetCreate(&resi_topos);
            auto rotamer_buffer = make_buffer(rotamer_path);
            check_status(
                [&]() { return StructureReadSmallMolRotamers(&self.value, &resi_topos, rotamer_buffer.data()); },
                "StructureReadSmallMolRotamers");
            ResiTopoSetDestroy(&resi_topos);
          },
          py::arg("rotamer_path"))
      .def(
          "run_monomer_design",
          [](StructureHandle& self, const py::dict& options) {
            PyMonomerDesignOptions native_options = parse_monomer_options(options);
            PyMonomerDesignResult native_result;
            int status = RunMonomerDesignWorkflow(&self.value, native_options, &native_result);
            if (FAILED(status)) {
              throw std::runtime_error("Monomer design failed with code " + std::to_string(status));
            }

            py::dict payload;
            payload["sequence_string"] = native_result.sequence_string;
            payload["trajectory_index"] = native_result.trajectory_index;
            payload["sequence_identity"] = native_result.sequence_identity;
            payload["energy_total"] = native_result.energy_total;
            payload["energy_evolution"] = native_result.energy_evolution;
            payload["energy_physical"] = native_result.energy_physical;
            payload["energy_binding"] = native_result.energy_binding;
            payload["unsatisfied_constraints"] = native_result.unsatisfied_constraints;

            if (native_result.has_best_structure) {
              StructureHandle best;
              check_status(
                  [&]() { return StructureDestroy(&best.value); }, "StructureDestroy(best_structure)");
              check_status(
                  [&]() { return StructureCreate(&best.value); }, "StructureCreate(best_structure)");
              check_status(
                  [&]() { return StructureCopy(&best.value, &native_result.best_structure); },
                  "StructureCopy(best_structure)");
              payload["best_structure"] = best;
            } else {
              payload["best_structure"] = py::none();
            }

            if (native_result.has_best_sites_structure) {
              StructureHandle best_sites;
              check_status(
                  [&]() { return StructureDestroy(&best_sites.value); },
                  "StructureDestroy(best_sites_structure)");
              check_status(
                  [&]() { return StructureCreate(&best_sites.value); },
                  "StructureCreate(best_sites_structure)");
              check_status(
                  [&]() { return StructureCopy(&best_sites.value, &native_result.best_sites_structure); },
                  "StructureCopy(best_sites_structure)");
              payload["best_sites_structure"] = best_sites;
            } else {
              payload["best_sites_structure"] = py::none();
            }

            if (native_result.has_best_mutable_sites_structure) {
              StructureHandle best_mutable_sites;
              check_status(
                  [&]() { return StructureDestroy(&best_mutable_sites.value); },
                  "StructureDestroy(best_mutable_sites_structure)");
              check_status(
                  [&]() { return StructureCreate(&best_mutable_sites.value); },
                  "StructureCreate(best_mutable_sites_structure)");
              check_status(
                  [&]() {
                    return StructureCopy(&best_mutable_sites.value,
                                         &native_result.best_mutable_sites_structure);
                  },
                  "StructureCopy(best_mutable_sites_structure)");
              payload["best_mutable_sites_structure"] = best_mutable_sites;
            } else {
              payload["best_mutable_sites_structure"] = py::none();
            }

            py::list residue_energy_list;
            for (const auto& entry : native_result.residue_self_energies) {
              py::dict row;
              row["chain"] = entry.chain_name;
              row["position"] = entry.position;
              row["self_energy"] = entry.self_energy;
              row["binding_energy"] = entry.binding_energy;
              residue_energy_list.append(row);
            }
            payload["residue_self_energies"] = std::move(residue_energy_list);

            return payload;
          },
          py::arg("options"))
      .def(
          "run_minimization",
          [](StructureHandle& self, const py::dict& options) {
            PyMinimizeOptions native_options = parse_minimize_options(options);
            PyMinimizeResult native_result;
            int status = RunMinimizeWorkflow(&self.value, native_options, &native_result);
            if (FAILED(status)) {
              throw std::runtime_error("Minimization failed with code " + std::to_string(status));
            }

            py::dict payload;
            if (native_result.has_structure) {
              StructureHandle minimized;
              check_status(
                  [&]() { return StructureCopy(&minimized.value, &native_result.minimized_structure); },
                  "StructureCopy(minimized_structure)");
              payload["structure"] = minimized;
            } else {
              payload["structure"] = py::none();
            }

            return payload;
          },
          py::arg("options"));

  m.def(
      "print_version",
      []() {
        return PrintVersion();
      },
      "Proxy around the native PrintVersion helper");

  m.def(
      "generate_small_mol_parameter_and_topology",
      [](const std::string& mol2file, const std::string& parfile, const std::string& topfile, const std::string& iniAtom1, const std::string& iniAtom2, const std::string& iniAtom3) {
        auto mol2file_buf = make_buffer(mol2file);
        auto parfile_buf = make_buffer(parfile);
        auto topfile_buf = make_buffer(topfile);
        auto iniAtom1_buf = make_buffer(iniAtom1);
        auto iniAtom2_buf = make_buffer(iniAtom2);
        auto iniAtom3_buf = make_buffer(iniAtom3);
        check_status([&]() { return GenerateSmallMolParameterAndTopologyFromMol2(mol2file_buf.data(), parfile_buf.data(), topfile_buf.data(), iniAtom1_buf.data(), iniAtom2_buf.data(), iniAtom3_buf.data()); }, "GenerateSmallMolParameterAndTopologyFromMol2");
      },
      "Generate small molecule parameter and topology files from a MOL2 file",
      py::arg("mol2file"),
      py::arg("parfile"),
      py::arg("topfile"),
      py::arg("iniAtom1"),
      py::arg("iniAtom2"),
      py::arg("iniAtom3")
  );
}
