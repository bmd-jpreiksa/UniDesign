#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ManagedTypes.hpp"
#include "ProgramFunction.h"
#include "EnergyFunction.h"
#include "Structure.h"
#include "Sequence.h"
#include "ErrorTracker.h"
#include "Utility.h"

#include <array>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace py = pybind11;

namespace
{
void ThrowIfFailed(int code, const char* functionName)
{
  if (FAILED(code))
  {
    throw std::runtime_error(std::string(functionName) + " failed with error code " + std::to_string(code));
  }
}

struct StringConfigEntry
{
  char* buffer;
  std::size_t capacity;
};

template <std::size_t N>
StringConfigEntry MakeEntry(char (&arr)[N])
{
  return {arr, N - 1};
}

struct BoolConfigEntry
{
  BOOL* value;
};

struct DoubleConfigEntry
{
  double* value;
};

struct IntConfigEntry
{
  int* value;
};

std::mutex& GetConfigMutex()
{
  static std::mutex mutex;
  return mutex;
}

// Configuration globals declared in Main.cpp
extern char PROGRAM_PATH[MAX_LEN_FILE_NAME + 1];
extern char PROGRAM_NAME[MAX_LEN_FILE_NAME + 1];
extern char PROGRAM_VERSION[MAX_LEN_FILE_NAME + 1];
extern char FILE_ATOMPARAM[MAX_LEN_FILE_NAME + 1];
extern char FILE_TOPO[MAX_LEN_FILE_NAME + 1];
extern char LIG_PARAM[MAX_LEN_FILE_NAME + 1];
extern char LIG_TOPO[MAX_LEN_FILE_NAME + 1];
extern char FILE_AAPROPENSITY[MAX_LEN_FILE_NAME + 1];
extern char FILE_RAMACHANDRAN[MAX_LEN_FILE_NAME + 1];
extern char FILE_WEIGHT_READ[MAX_LEN_FILE_NAME + 1];
extern char FILE_ROTLIB[MAX_LEN_FILE_NAME + 1];
extern char FILE_ROTLIB_BIN[MAX_LEN_FILE_NAME + 1];
extern char USER_ROTLIB_NAME[MAX_LEN_FILE_NAME + 1];
extern char TGT_PRF[MAX_LEN_FILE_NAME + 1];
extern char TGT_MSA[MAX_LEN_FILE_NAME + 1];
extern char TGT_SA[MAX_LEN_FILE_NAME + 1];
extern char TGT_SS[MAX_LEN_FILE_NAME + 1];
extern char TGT_SEQ[MAX_LEN_FILE_NAME + 1];
extern char TGT_PHIPSI[MAX_LEN_FILE_NAME + 1];
extern char FILE_CATACONS[MAX_LEN_FILE_NAME + 1];
extern char FILE_LIG_PLACEMENT[MAX_LEN_FILE_NAME + 1];
extern char FILE_SELF_ENERGY[MAX_LEN_FILE_NAME + 1];
extern char FILE_ROTLIST[MAX_LEN_FILE_NAME + 1];
extern char FILE_ROTLIST_SEC[MAX_LEN_FILE_NAME + 1];
extern char FILE_DESROT_NDX[MAX_LEN_FILE_NAME + 1];
extern char FILE_DESSEQS[MAX_LEN_FILE_NAME + 1];
extern char FILE_BESTSEQS[MAX_LEN_FILE_NAME + 1];
extern char FILE_BESTSTRUCT[MAX_LEN_FILE_NAME + 1];
extern char FILE_BEST_ALL_SITES[MAX_LEN_FILE_NAME + 1];
extern char FILE_BEST_MUT_SITES[MAX_LEN_FILE_NAME + 1];
extern char FILE_BEST_LIG_MOL2[MAX_LEN_FILE_NAME + 1];
extern char PREFIX[MAX_LEN_FILE_NAME + 1];
extern char PDB[MAX_LEN_FILE_NAME + 1];
extern char PDBPATH[MAX_LEN_FILE_NAME + 1];
extern char PDBNAME[MAX_LEN_FILE_NAME + 1];
extern char PDBID[MAX_LEN_FILE_NAME + 1];
extern char MOL2[MAX_LEN_FILE_NAME + 1];
extern char DES_CHAINS[10];
extern char INI_ATOM1[MAX_LEN_ATOM_NAME + 1];
extern char INI_ATOM2[MAX_LEN_ATOM_NAME + 1];
extern char INI_ATOM3[MAX_LEN_ATOM_NAME + 1];
extern char SPLIT_CHAINS[10];
extern char SPLIT_PART1[10];
extern char SPLIT_PART2[10];
extern char RESI[10];
extern char EXCL_RESI[100];
extern char RESI_PAIR[2 * MAX_LEN_CHAIN_NAME + 1];
extern char MUTANT_FILE[MAX_LEN_FILE_NAME + 1];
extern char FILE_LIG_POSES_IN[MAX_LEN_FILE_NAME + 1];
extern char FILE_LIG_POSES_OUT[MAX_LEN_FILE_NAME + 1];
extern char FILE_LIG_SCREEN_BY_ORITENTATION[MAX_LEN_FILE_NAME + 1];
extern char PDB2[MAX_LEN_FILE_NAME + 1];
extern char PDBLIST[MAX_LEN_FILE_NAME + 1];
extern char FILE_RESFILE[MAX_LEN_FILE_NAME + 1];
extern char REFERENCE_RESIDUES[MAX_LEN_ONE_LINE_CONTENT + 1];

extern BOOL FLAG_USER_ROTLIB;
extern BOOL FLAG_CHAIN_SPLIT;
extern BOOL FLAG_LIG_POSES;
extern BOOL FLAG_LIG_SCREEN_BY_ORIENTATION;
extern BOOL FLAG_LIG_SCREEN_BY_TOPVDW;
extern BOOL FLAG_LIG_SCREEN_BY_RMSD;
extern BOOL FLAG_PDB;
extern BOOL FLAG_MOL2;
extern BOOL FLAG_MONOMER;
extern BOOL FLAG_PPI;
extern BOOL FLAG_PROT_LIG;
extern BOOL FLAG_ENZYME;
extern BOOL FLAG_PHYSICS;
extern BOOL FLAG_EVOLUTION;
extern BOOL FLAG_EVOPHIPSI;
extern BOOL FLAG_BBDEP_ROTLIB;
extern BOOL FLAG_USE_INPUT_SC;
extern BOOL FLAG_ROTATE_HYDROXYL;
extern BOOL FLAG_WILDTYPE_ONLY;
extern BOOL FLAG_INTERFACE_ONLY;
extern BOOL FLAG_EXCL_CYS_ROTS;
extern BOOL FLAG_RESFILE;
extern BOOL FLAG_DESIGN_FROM_NATAA;
extern BOOL FLAG_READ_HYDROGEN;
extern BOOL FLAG_WRITE_HYDROGEN;

extern double CUT_EXCL_LOW_PROB_ROT;
extern double WGT_PROFILE;
extern double WGT_BIND;
extern double CUT_PPI_DIST_SHELL1;
extern double CUT_PPI_DIST_SHELL2;
extern double CUT_PLI_DIST_SHELL1;
extern double CUT_PLI_DIST_SHELL2;
extern double LIG_SCREEN_TOP_VDW_PERCENTILE;
extern double LIG_SCREEN_RMSD_CUTOFF;

extern int NTRAJ;
extern int NTRAJ_START_NDX;
extern int MAX_NUM_OF_RUNS;
extern int CUT_NUM_CB_CORE;
extern int CUT_NUM_CB_SURF;
extern int PROT_LEN_NORM;

const std::unordered_map<std::string, StringConfigEntry>& GetStringEntries()
{
  static const std::unordered_map<std::string, StringConfigEntry> entries = {
    {"PROGRAM_PATH", MakeEntry(PROGRAM_PATH)},
    {"PROGRAM_NAME", MakeEntry(PROGRAM_NAME)},
    {"PROGRAM_VERSION", MakeEntry(PROGRAM_VERSION)},
    {"FILE_ATOMPARAM", MakeEntry(FILE_ATOMPARAM)},
    {"FILE_TOPO", MakeEntry(FILE_TOPO)},
    {"LIG_PARAM", MakeEntry(LIG_PARAM)},
    {"LIG_TOPO", MakeEntry(LIG_TOPO)},
    {"FILE_AAPROPENSITY", MakeEntry(FILE_AAPROPENSITY)},
    {"FILE_RAMACHANDRAN", MakeEntry(FILE_RAMACHANDRAN)},
    {"FILE_WEIGHT_READ", MakeEntry(FILE_WEIGHT_READ)},
    {"FILE_ROTLIB", MakeEntry(FILE_ROTLIB)},
    {"FILE_ROTLIB_BIN", MakeEntry(FILE_ROTLIB_BIN)},
    {"USER_ROTLIB_NAME", MakeEntry(USER_ROTLIB_NAME)},
    {"TGT_PRF", MakeEntry(TGT_PRF)},
    {"TGT_MSA", MakeEntry(TGT_MSA)},
    {"TGT_SA", MakeEntry(TGT_SA)},
    {"TGT_SS", MakeEntry(TGT_SS)},
    {"TGT_SEQ", MakeEntry(TGT_SEQ)},
    {"TGT_PHIPSI", MakeEntry(TGT_PHIPSI)},
    {"FILE_CATACONS", MakeEntry(FILE_CATACONS)},
    {"FILE_LIG_PLACEMENT", MakeEntry(FILE_LIG_PLACEMENT)},
    {"FILE_SELF_ENERGY", MakeEntry(FILE_SELF_ENERGY)},
    {"FILE_ROTLIST", MakeEntry(FILE_ROTLIST)},
    {"FILE_ROTLIST_SEC", MakeEntry(FILE_ROTLIST_SEC)},
    {"FILE_DESROT_NDX", MakeEntry(FILE_DESROT_NDX)},
    {"FILE_DESSEQS", MakeEntry(FILE_DESSEQS)},
    {"FILE_BESTSEQS", MakeEntry(FILE_BESTSEQS)},
    {"FILE_BESTSTRUCT", MakeEntry(FILE_BESTSTRUCT)},
    {"FILE_BEST_ALL_SITES", MakeEntry(FILE_BEST_ALL_SITES)},
    {"FILE_BEST_MUT_SITES", MakeEntry(FILE_BEST_MUT_SITES)},
    {"FILE_BEST_LIG_MOL2", MakeEntry(FILE_BEST_LIG_MOL2)},
    {"PREFIX", MakeEntry(PREFIX)},
    {"PDB", MakeEntry(PDB)},
    {"PDBPATH", MakeEntry(PDBPATH)},
    {"PDBNAME", MakeEntry(PDBNAME)},
    {"PDBID", MakeEntry(PDBID)},
    {"MOL2", MakeEntry(MOL2)},
    {"DES_CHAINS", MakeEntry(DES_CHAINS)},
    {"INI_ATOM1", MakeEntry(INI_ATOM1)},
    {"INI_ATOM2", MakeEntry(INI_ATOM2)},
    {"INI_ATOM3", MakeEntry(INI_ATOM3)},
    {"SPLIT_CHAINS", MakeEntry(SPLIT_CHAINS)},
    {"SPLIT_PART1", MakeEntry(SPLIT_PART1)},
    {"SPLIT_PART2", MakeEntry(SPLIT_PART2)},
    {"RESI", MakeEntry(RESI)},
    {"EXCL_RESI", MakeEntry(EXCL_RESI)},
    {"RESI_PAIR", MakeEntry(RESI_PAIR)},
    {"MUTANT_FILE", MakeEntry(MUTANT_FILE)},
    {"FILE_LIG_POSES_IN", MakeEntry(FILE_LIG_POSES_IN)},
    {"FILE_LIG_POSES_OUT", MakeEntry(FILE_LIG_POSES_OUT)},
    {"FILE_LIG_SCREEN_BY_ORITENTATION", MakeEntry(FILE_LIG_SCREEN_BY_ORITENTATION)},
    {"PDB2", MakeEntry(PDB2)},
    {"PDBLIST", MakeEntry(PDBLIST)},
    {"FILE_RESFILE", MakeEntry(FILE_RESFILE)},
    {"REFERENCE_RESIDUES", MakeEntry(REFERENCE_RESIDUES)},
  };
  return entries;
}

const std::unordered_map<std::string, BoolConfigEntry>& GetBoolEntries()
{
  static const std::unordered_map<std::string, BoolConfigEntry> entries = {
    {"FLAG_USER_ROTLIB", {&FLAG_USER_ROTLIB}},
    {"FLAG_CHAIN_SPLIT", {&FLAG_CHAIN_SPLIT}},
    {"FLAG_LIG_POSES", {&FLAG_LIG_POSES}},
    {"FLAG_LIG_SCREEN_BY_ORIENTATION", {&FLAG_LIG_SCREEN_BY_ORIENTATION}},
    {"FLAG_LIG_SCREEN_BY_TOPVDW", {&FLAG_LIG_SCREEN_BY_TOPVDW}},
    {"FLAG_LIG_SCREEN_BY_RMSD", {&FLAG_LIG_SCREEN_BY_RMSD}},
    {"FLAG_PDB", {&FLAG_PDB}},
    {"FLAG_MOL2", {&FLAG_MOL2}},
    {"FLAG_MONOMER", {&FLAG_MONOMER}},
    {"FLAG_PPI", {&FLAG_PPI}},
    {"FLAG_PROT_LIG", {&FLAG_PROT_LIG}},
    {"FLAG_ENZYME", {&FLAG_ENZYME}},
    {"FLAG_PHYSICS", {&FLAG_PHYSICS}},
    {"FLAG_EVOLUTION", {&FLAG_EVOLUTION}},
    {"FLAG_EVOPHIPSI", {&FLAG_EVOPHIPSI}},
    {"FLAG_BBDEP_ROTLIB", {&FLAG_BBDEP_ROTLIB}},
    {"FLAG_USE_INPUT_SC", {&FLAG_USE_INPUT_SC}},
    {"FLAG_ROTATE_HYDROXYL", {&FLAG_ROTATE_HYDROXYL}},
    {"FLAG_WILDTYPE_ONLY", {&FLAG_WILDTYPE_ONLY}},
    {"FLAG_INTERFACE_ONLY", {&FLAG_INTERFACE_ONLY}},
    {"FLAG_EXCL_CYS_ROTS", {&FLAG_EXCL_CYS_ROTS}},
    {"FLAG_RESFILE", {&FLAG_RESFILE}},
    {"FLAG_DESIGN_FROM_NATAA", {&FLAG_DESIGN_FROM_NATAA}},
    {"FLAG_READ_HYDROGEN", {&FLAG_READ_HYDROGEN}},
    {"FLAG_WRITE_HYDROGEN", {&FLAG_WRITE_HYDROGEN}},
  };
  return entries;
}

const std::unordered_map<std::string, DoubleConfigEntry>& GetDoubleEntries()
{
  static const std::unordered_map<std::string, DoubleConfigEntry> entries = {
    {"CUT_EXCL_LOW_PROB_ROT", {&CUT_EXCL_LOW_PROB_ROT}},
    {"WGT_PROFILE", {&WGT_PROFILE}},
    {"WGT_BIND", {&WGT_BIND}},
    {"CUT_PPI_DIST_SHELL1", {&CUT_PPI_DIST_SHELL1}},
    {"CUT_PPI_DIST_SHELL2", {&CUT_PPI_DIST_SHELL2}},
    {"CUT_PLI_DIST_SHELL1", {&CUT_PLI_DIST_SHELL1}},
    {"CUT_PLI_DIST_SHELL2", {&CUT_PLI_DIST_SHELL2}},
    {"LIG_SCREEN_TOP_VDW_PERCENTILE", {&LIG_SCREEN_TOP_VDW_PERCENTILE}},
    {"LIG_SCREEN_RMSD_CUTOFF", {&LIG_SCREEN_RMSD_CUTOFF}},
  };
  return entries;
}

const std::unordered_map<std::string, IntConfigEntry>& GetIntEntries()
{
  static const std::unordered_map<std::string, IntConfigEntry> entries = {
    {"NTRAJ", {&NTRAJ}},
    {"NTRAJ_START_NDX", {&NTRAJ_START_NDX}},
    {"MAX_NUM_OF_RUNS", {&MAX_NUM_OF_RUNS}},
    {"CUT_NUM_CB_CORE", {&CUT_NUM_CB_CORE}},
    {"CUT_NUM_CB_SURF", {&CUT_NUM_CB_SURF}},
    {"PROT_LEN_NORM", {&PROT_LEN_NORM}},
  };
  return entries;
}

std::string GetPathConfig(const std::string& name)
{
  const auto& entries = GetStringEntries();
  auto it = entries.find(name);
  if (it == entries.end())
  {
    throw std::invalid_argument("Unknown path configuration key: " + name);
  }
  std::lock_guard<std::mutex> lock(GetConfigMutex());
  return std::string(it->second.buffer);
}

void SetPathConfig(const std::string& name, const std::string& value)
{
  const auto& entries = GetStringEntries();
  auto it = entries.find(name);
  if (it == entries.end())
  {
    throw std::invalid_argument("Unknown path configuration key: " + name);
  }
  if (value.size() > it->second.capacity)
  {
    throw std::length_error("Value for " + name + " exceeds maximum length of " + std::to_string(it->second.capacity));
  }
  std::lock_guard<std::mutex> lock(GetConfigMutex());
  std::strncpy(it->second.buffer, value.c_str(), it->second.capacity);
  it->second.buffer[it->second.capacity] = '\0';
}

bool GetFlagConfig(const std::string& name)
{
  const auto& entries = GetBoolEntries();
  auto it = entries.find(name);
  if (it == entries.end())
  {
    throw std::invalid_argument("Unknown flag configuration key: " + name);
  }
  std::lock_guard<std::mutex> lock(GetConfigMutex());
  return *(it->second.value) != FALSE;
}

void SetFlagConfig(const std::string& name, bool value)
{
  const auto& entries = GetBoolEntries();
  auto it = entries.find(name);
  if (it == entries.end())
  {
    throw std::invalid_argument("Unknown flag configuration key: " + name);
  }
  std::lock_guard<std::mutex> lock(GetConfigMutex());
  *(it->second.value) = value ? TRUE : FALSE;
}

double GetDoubleConfig(const std::string& name)
{
  const auto& entries = GetDoubleEntries();
  auto it = entries.find(name);
  if (it == entries.end())
  {
    throw std::invalid_argument("Unknown numeric configuration key: " + name);
  }
  std::lock_guard<std::mutex> lock(GetConfigMutex());
  return *(it->second.value);
}

void SetDoubleConfig(const std::string& name, double value)
{
  const auto& entries = GetDoubleEntries();
  auto it = entries.find(name);
  if (it == entries.end())
  {
    throw std::invalid_argument("Unknown numeric configuration key: " + name);
  }
  std::lock_guard<std::mutex> lock(GetConfigMutex());
  *(it->second.value) = value;
}

int GetIntConfig(const std::string& name)
{
  const auto& entries = GetIntEntries();
  auto it = entries.find(name);
  if (it == entries.end())
  {
    throw std::invalid_argument("Unknown integer configuration key: " + name);
  }
  std::lock_guard<std::mutex> lock(GetConfigMutex());
  return *(it->second.value);
}

void SetIntConfig(const std::string& name, int value)
{
  const auto& entries = GetIntEntries();
  auto it = entries.find(name);
  if (it == entries.end())
  {
    throw std::invalid_argument("Unknown integer configuration key: " + name);
  }
  std::lock_guard<std::mutex> lock(GetConfigMutex());
  *(it->second.value) = value;
}

std::vector<std::string> ListKeys(const std::unordered_map<std::string, StringConfigEntry>& entries)
{
  std::vector<std::string> keys;
  keys.reserve(entries.size());
  for (const auto& kv : entries)
  {
    keys.push_back(kv.first);
  }
  return keys;
}

std::vector<std::string> ListKeys(const std::unordered_map<std::string, BoolConfigEntry>& entries)
{
  std::vector<std::string> keys;
  keys.reserve(entries.size());
  for (const auto& kv : entries)
  {
    keys.push_back(kv.first);
  }
  return keys;
}

std::vector<std::string> ListKeys(const std::unordered_map<std::string, DoubleConfigEntry>& entries)
{
  std::vector<std::string> keys;
  keys.reserve(entries.size());
  for (const auto& kv : entries)
  {
    keys.push_back(kv.first);
  }
  return keys;
}

std::vector<std::string> ListKeys(const std::unordered_map<std::string, IntConfigEntry>& entries)
{
  std::vector<std::string> keys;
  keys.reserve(entries.size());
  for (const auto& kv : entries)
  {
    keys.push_back(kv.first);
  }
  return keys;
}

std::vector<double> ComputeStructureStabilityFromConfig(StructureHandle& structure)
{
  std::array<double, MAX_ENERGY_TERM> energyTerms{};
  AAppTable aapptable;
  RamaTable ramatable;
  ThrowIfFailed(AApropensityTableReadFromFile(&aapptable, FILE_AAPROPENSITY), "AApropensityTableReadFromFile");
  ThrowIfFailed(RamaTableReadFromFile(&ramatable, FILE_RAMACHANDRAN), "RamaTableReadFromFile");
  ThrowIfFailed(ComputeStructureStability(structure.get(), &aapptable, &ramatable, energyTerms.data()), "ComputeStructureStability");
  return std::vector<double>(energyTerms.begin(), energyTerms.end());
}

void ComputeBindingFromConfig(StructureHandle& structure)
{
  ThrowIfFailed(ComputeBinding(structure.get()), "ComputeBinding");
}

void ProteinDesignFromConfig()
{
  throw std::runtime_error("ProteinDesign workflow is not yet exposed through the Python module");
}

} // namespace

PYBIND11_MODULE(unidesign, m)
{
  m.doc() = "Python bindings for the UniDesign core library";

  py::class_<StructureHandle>(m, "StructureHandle")
    .def(py::init<>())
    .def(py::init<const StructureHandle&>())
    .def("swap", &StructureHandle::swap)
    .def("__repr__", [](const StructureHandle*) { return std::string("<unidesign.StructureHandle>"); });

  py::class_<SequenceHandle>(m, "SequenceHandle")
    .def(py::init<>())
    .def(py::init<const SequenceHandle&>())
    .def("swap", &SequenceHandle::swap)
    .def("__repr__", [](const SequenceHandle*) { return std::string("<unidesign.SequenceHandle>"); });

  m.def("compute_structure_stability", &ComputeStructureStabilityFromConfig, py::arg("structure"),
    "Compute structure stability using the global configuration files to populate tables.");

  m.def("compute_binding", &ComputeBindingFromConfig, py::arg("structure"),
    "Compute binding energy for the provided structure.");

  m.def("protein_design", &ProteinDesignFromConfig,
    "Run the protein design workflow using the current configuration.");

  m.def("get_path", &GetPathConfig, py::arg("name"), "Get a string configuration value.");
  m.def("set_path", &SetPathConfig, py::arg("name"), py::arg("value"), "Set a string configuration value.");

  m.def("get_flag", &GetFlagConfig, py::arg("name"), "Get a boolean configuration flag.");
  m.def("set_flag", &SetFlagConfig, py::arg("name"), py::arg("value"), "Set a boolean configuration flag.");

  m.def("get_cutoff", &GetDoubleConfig, py::arg("name"), "Get a floating-point configuration value.");
  m.def("set_cutoff", &SetDoubleConfig, py::arg("name"), py::arg("value"), "Set a floating-point configuration value.");

  m.def("get_integer", &GetIntConfig, py::arg("name"), "Get an integer configuration value.");
  m.def("set_integer", &SetIntConfig, py::arg("name"), py::arg("value"), "Set an integer configuration value.");

  m.def("list_paths", [] { return ListKeys(GetStringEntries()); }, "List available string configuration keys.");
  m.def("list_flags", [] { return ListKeys(GetBoolEntries()); }, "List available boolean configuration keys.");
  m.def("list_cutoffs", [] { return ListKeys(GetDoubleEntries()); }, "List available floating-point configuration keys.");
  m.def("list_integers", [] { return ListKeys(GetIntEntries()); }, "List available integer configuration keys.");
}
