#include "PyDesign.h"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "Atom.h"
#include "AtomParamsSet.h"
#include "DesignSite.h"
#include "EnergyMatrix.h"
#include "Residue.h"
#include "ResidueTopology.h"
#include "RotamerBuilder.h"
#include "Sequence.h"
#include "Structure.h"
#include "Utility.h"

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
extern BOOL FLAG_LIG_POSES;

extern int NTRAJ;
extern int NTRAJ_START_NDX;

extern double WGT_PROFILE;
extern double WGT_BIND;

extern char PROGRAM_PATH[MAX_LEN_FILE_NAME + 1];
extern char FILE_ATOMPARAM[MAX_LEN_FILE_NAME + 1];
extern char FILE_TOPO[MAX_LEN_FILE_NAME + 1];
extern char FILE_WEIGHT_READ[MAX_LEN_FILE_NAME + 1];
extern char FILE_AAPROPENSITY[MAX_LEN_FILE_NAME + 1];
extern char FILE_RAMACHANDRAN[MAX_LEN_FILE_NAME + 1];
extern char FILE_ROTLIB[MAX_LEN_FILE_NAME + 1];
extern char FILE_ROTLIB_BIN[MAX_LEN_FILE_NAME + 1];
extern char FILE_RESFILE[MAX_LEN_FILE_NAME + 1];
extern char PREFIX[MAX_LEN_FILE_NAME + 1];
extern char FILE_SELF_ENERGY[MAX_LEN_FILE_NAME + 1];
extern char FILE_ROTLIST[MAX_LEN_FILE_NAME + 1];
extern char FILE_ROTLIST_SEC[MAX_LEN_FILE_NAME + 1];
extern char FILE_BESTSEQS[MAX_LEN_FILE_NAME + 1];
extern char FILE_BESTSTRUCT[MAX_LEN_FILE_NAME + 1];
extern char FILE_BEST_ALL_SITES[MAX_LEN_FILE_NAME + 1];
extern char FILE_BEST_MUT_SITES[MAX_LEN_FILE_NAME + 1];
extern char FILE_BEST_LIG_MOL2[MAX_LEN_FILE_NAME + 1];
extern char FILE_LIG_POSES_IN[MAX_LEN_FILE_NAME + 1];
extern char FILE_LIG_POSES_OUT[MAX_LEN_FILE_NAME + 1];
extern char DES_CHAINS[10];

namespace {

class ScopedOutputSilencer {
 public:
  explicit ScopedOutputSilencer(bool enable) : enabled_(enable) {
    if (!enabled_) {
      return;
    }
    fflush(stdout);
    fflush(stderr);
    stdout_fd_ = dup(fileno(stdout));
    stderr_fd_ = dup(fileno(stderr));
#if defined(_WIN32)
    const char* null_path = "NUL";
#else
    const char* null_path = "/dev/null";
#endif
    null_out_ = fopen(null_path, "w");
    null_err_ = fopen(null_path, "w");
    if (null_out_ != nullptr) {
      dup2(fileno(null_out_), fileno(stdout));
    }
    if (null_err_ != nullptr) {
      dup2(fileno(null_err_), fileno(stderr));
    }
  }

  ~ScopedOutputSilencer() {
    if (!enabled_) {
      return;
    }
    fflush(stdout);
    fflush(stderr);
    if (stdout_fd_ >= 0) {
      dup2(stdout_fd_, fileno(stdout));
      close(stdout_fd_);
    }
    if (stderr_fd_ >= 0) {
      dup2(stderr_fd_, fileno(stderr));
      close(stderr_fd_);
    }
    if (null_out_ != nullptr) {
      fclose(null_out_);
    }
    if (null_err_ != nullptr) {
      fclose(null_err_);
    }
  }

 private:
  bool enabled_;
  int stdout_fd_{-1};
  int stderr_fd_{-1};
  FILE* null_out_{nullptr};
  FILE* null_err_{nullptr};
};

static void CaptureResidueSelfEnergies(Structure* structure,
                                       const std::string& filepath,
                                       std::vector<PyResidueSelfEnergy>* output) {
  if (structure == nullptr || output == nullptr) {
    return;
  }
  output->clear();
  std::ifstream input(filepath);
  if (!input.is_open()) {
    return;
  }
  const int site_count = StructureGetDesignSiteCount(structure);
  if (site_count <= 0) {
    return;
  }
  std::vector<double> best_self(site_count, std::numeric_limits<double>::infinity());
  std::vector<double> best_binding(site_count, 0.0);
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    std::istringstream stream(line);
    int site_index = -1;
    int rotamer_index = -1;
    double total_energy = 0.0;
    double binding_energy = 0.0;
    if (!(stream >> site_index >> rotamer_index >> total_energy)) {
      continue;
    }
    if (!(stream >> binding_energy)) {
      binding_energy = 0.0;
    }
    if (site_index < 0 || site_index >= site_count) {
      continue;
    }
    if (total_energy < best_self[site_index]) {
      best_self[site_index] = total_energy;
      best_binding[site_index] = binding_energy;
    }
  }

  output->reserve(site_count);
  for (int i = 0; i < site_count; ++i) {
    DesignSite* site = StructureGetDesignSite(structure, i);
    if (site == nullptr) {
      continue;
    }
    Chain* chain = StructureGetChain(structure, site->chnNdx);
    if (chain == nullptr) {
      continue;
    }
    Residue* residue = ChainGetResidue(chain, site->resNdx);
    if (residue == nullptr) {
      continue;
    }
    PyResidueSelfEnergy entry;
    entry.chain_name = ChainGetName(chain);
    entry.position = ResidueGetPosInChain(residue);
    double total_energy = best_self[i];
    if (!std::isfinite(total_energy)) {
      total_energy = 0.0;
    }
    entry.self_energy = total_energy;
    entry.binding_energy = best_binding[i];
    output->push_back(entry);
  }
}

struct GlobalDesignStateGuard {
  BOOL flag_monomer;
  BOOL flag_ppi;
  BOOL flag_prot_lig;
  BOOL flag_enzyme;
  BOOL flag_physics;
  BOOL flag_evolution;
  BOOL flag_evophipsi;
  BOOL flag_bbdep_rotlib;
  BOOL flag_use_input_sc;
  BOOL flag_rotate_hydroxyl;
  BOOL flag_wildtype_only;
  BOOL flag_interface_only;
  BOOL flag_excl_cys_rots;
  BOOL flag_resfile;
  BOOL flag_design_from_nataa;
  BOOL flag_read_hydrogen;
  BOOL flag_write_hydrogen;
  BOOL flag_lig_poses;

  int ntraj;
  int ntraj_start_ndx;

  double wgt_profile;
  double wgt_bind;

  char program_path_copy[MAX_LEN_FILE_NAME + 1];
  char file_atomparam_copy[MAX_LEN_FILE_NAME + 1];
  char file_topo_copy[MAX_LEN_FILE_NAME + 1];
  char file_weight_copy[MAX_LEN_FILE_NAME + 1];
  char file_aapp_copy[MAX_LEN_FILE_NAME + 1];
  char file_rama_copy[MAX_LEN_FILE_NAME + 1];
  char file_rotlib_copy[MAX_LEN_FILE_NAME + 1];
  char file_rotlib_bin_copy[MAX_LEN_FILE_NAME + 1];
  char file_resfile_copy[MAX_LEN_FILE_NAME + 1];
  char prefix_copy[MAX_LEN_FILE_NAME + 1];
  char file_self_energy_copy[MAX_LEN_FILE_NAME + 1];
  char file_rotlist_copy[MAX_LEN_FILE_NAME + 1];
  char file_rotlist_sec_copy[MAX_LEN_FILE_NAME + 1];
  char file_bestseqs_copy[MAX_LEN_FILE_NAME + 1];
  char file_beststruct_copy[MAX_LEN_FILE_NAME + 1];
  char file_bestsites_copy[MAX_LEN_FILE_NAME + 1];
  char file_bestmutsites_copy[MAX_LEN_FILE_NAME + 1];
  char file_bestlig_copy[MAX_LEN_FILE_NAME + 1];
  char file_lig_pose_in_copy[MAX_LEN_FILE_NAME + 1];
  char file_lig_pose_out_copy[MAX_LEN_FILE_NAME + 1];
  char des_chains_copy[sizeof(DES_CHAINS)];

  GlobalDesignStateGuard() {
    flag_monomer = FLAG_MONOMER;
    flag_ppi = FLAG_PPI;
    flag_prot_lig = FLAG_PROT_LIG;
    flag_enzyme = FLAG_ENZYME;
    flag_physics = FLAG_PHYSICS;
    flag_evolution = FLAG_EVOLUTION;
    flag_evophipsi = FLAG_EVOPHIPSI;
    flag_bbdep_rotlib = FLAG_BBDEP_ROTLIB;
    flag_use_input_sc = FLAG_USE_INPUT_SC;
    flag_rotate_hydroxyl = FLAG_ROTATE_HYDROXYL;
    flag_wildtype_only = FLAG_WILDTYPE_ONLY;
    flag_interface_only = FLAG_INTERFACE_ONLY;
    flag_excl_cys_rots = FLAG_EXCL_CYS_ROTS;
    flag_resfile = FLAG_RESFILE;
    flag_design_from_nataa = FLAG_DESIGN_FROM_NATAA;
    flag_read_hydrogen = FLAG_READ_HYDROGEN;
    flag_write_hydrogen = FLAG_WRITE_HYDROGEN;
    flag_lig_poses = FLAG_LIG_POSES;

    ntraj = NTRAJ;
    ntraj_start_ndx = NTRAJ_START_NDX;

    wgt_profile = WGT_PROFILE;
    wgt_bind = WGT_BIND;

    strcpy(program_path_copy, PROGRAM_PATH);
    strcpy(file_atomparam_copy, FILE_ATOMPARAM);
    strcpy(file_topo_copy, FILE_TOPO);
    strcpy(file_weight_copy, FILE_WEIGHT_READ);
    strcpy(file_aapp_copy, FILE_AAPROPENSITY);
    strcpy(file_rama_copy, FILE_RAMACHANDRAN);
    strcpy(file_rotlib_copy, FILE_ROTLIB);
    strcpy(file_rotlib_bin_copy, FILE_ROTLIB_BIN);
    strcpy(file_resfile_copy, FILE_RESFILE);
    strcpy(prefix_copy, PREFIX);
    strcpy(file_self_energy_copy, FILE_SELF_ENERGY);
    strcpy(file_rotlist_copy, FILE_ROTLIST);
    strcpy(file_rotlist_sec_copy, FILE_ROTLIST_SEC);
    strcpy(file_bestseqs_copy, FILE_BESTSEQS);
    strcpy(file_beststruct_copy, FILE_BESTSTRUCT);
    strcpy(file_bestsites_copy, FILE_BEST_ALL_SITES);
    strcpy(file_bestmutsites_copy, FILE_BEST_MUT_SITES);
    strcpy(file_bestlig_copy, FILE_BEST_LIG_MOL2);
    strcpy(file_lig_pose_in_copy, FILE_LIG_POSES_IN);
    strcpy(file_lig_pose_out_copy, FILE_LIG_POSES_OUT);
    strcpy(des_chains_copy, DES_CHAINS);
  }

  ~GlobalDesignStateGuard() {
    FLAG_MONOMER = flag_monomer;
    FLAG_PPI = flag_ppi;
    FLAG_PROT_LIG = flag_prot_lig;
    FLAG_ENZYME = flag_enzyme;
    FLAG_PHYSICS = flag_physics;
    FLAG_EVOLUTION = flag_evolution;
    FLAG_EVOPHIPSI = flag_evophipsi;
    FLAG_BBDEP_ROTLIB = flag_bbdep_rotlib;
    FLAG_USE_INPUT_SC = flag_use_input_sc;
    FLAG_ROTATE_HYDROXYL = flag_rotate_hydroxyl;
    FLAG_WILDTYPE_ONLY = flag_wildtype_only;
    FLAG_INTERFACE_ONLY = flag_interface_only;
    FLAG_EXCL_CYS_ROTS = flag_excl_cys_rots;
    FLAG_RESFILE = flag_resfile;
    FLAG_DESIGN_FROM_NATAA = flag_design_from_nataa;
    FLAG_READ_HYDROGEN = flag_read_hydrogen;
    FLAG_WRITE_HYDROGEN = flag_write_hydrogen;
    FLAG_LIG_POSES = flag_lig_poses;

    NTRAJ = ntraj;
    NTRAJ_START_NDX = ntraj_start_ndx;

    WGT_PROFILE = wgt_profile;
    WGT_BIND = wgt_bind;

    strcpy(PROGRAM_PATH, program_path_copy);
    strcpy(FILE_ATOMPARAM, file_atomparam_copy);
    strcpy(FILE_TOPO, file_topo_copy);
    strcpy(FILE_WEIGHT_READ, file_weight_copy);
    strcpy(FILE_AAPROPENSITY, file_aapp_copy);
    strcpy(FILE_RAMACHANDRAN, file_rama_copy);
    strcpy(FILE_ROTLIB, file_rotlib_copy);
    strcpy(FILE_ROTLIB_BIN, file_rotlib_bin_copy);
    strcpy(FILE_RESFILE, file_resfile_copy);
    strcpy(PREFIX, prefix_copy);
    strcpy(FILE_SELF_ENERGY, file_self_energy_copy);
    strcpy(FILE_ROTLIST, file_rotlist_copy);
    strcpy(FILE_ROTLIST_SEC, file_rotlist_sec_copy);
    strcpy(FILE_BESTSEQS, file_bestseqs_copy);
    strcpy(FILE_BESTSTRUCT, file_beststruct_copy);
    strcpy(FILE_BEST_ALL_SITES, file_bestsites_copy);
    strcpy(FILE_BEST_MUT_SITES, file_bestmutsites_copy);
    strcpy(FILE_BEST_LIG_MOL2, file_bestlig_copy);
    strcpy(FILE_LIG_POSES_IN, file_lig_pose_in_copy);
    strcpy(FILE_LIG_POSES_OUT, file_lig_pose_out_copy);
    strcpy(DES_CHAINS, des_chains_copy);
  }
};

class TempDirectory {
 public:
  TempDirectory(const std::string& base, const std::string& prefix) {
    std::string root = base.empty() ? "." : base;
    std::string templ = root;
    if (!templ.empty() && templ.back() != '/') templ.push_back('/');
    templ += prefix;
    templ += "XXXXXX";

    std::vector<char> buffer(templ.begin(), templ.end());
    buffer.push_back('\0');
    char* created = mkdtemp(buffer.data());
    if (created == nullptr) {
      throw std::runtime_error("Failed to create temporary directory");
    }
    path_ = created;
  }

  ~TempDirectory() {
    // The caller is responsible for removing files; attempt to remove directory.
    if (!path_.empty()) {
      rmdir(path_.c_str());
    }
  }

  const std::string& path() const { return path_; }

 private:
  std::string path_;
};

static void FixDesignSitePointers(Structure* structure) {
  for (int i = 0; i < StructureGetDesignSiteCount(structure); ++i) {
    DesignSite* site = StructureGetDesignSite(structure, i);
    Chain* chain = StructureGetChain(structure, site->chnNdx);
    site->pRes = ChainGetResidue(chain, site->resNdx);
  }
}

static std::string DetermineDefaultDesignChains(Structure* structure) {
  std::string result;
  for (int i = 0; i < StructureGetChainCount(structure); ++i) {
    Chain* chain = StructureGetChain(structure, i);
    if (ChainGetType(chain) == Type_Chain_Protein) {
      if (!result.empty()) {
        result += ',';
      }
      result += ChainGetName(chain);
    }
  }
  return result;
}

static int WriteFile(const std::string& path, const std::string& contents) {
  FILE* handle = fopen(path.c_str(), "w");
  if (!handle) {
    return IOError;
  }
  size_t written = fwrite(contents.data(), 1, contents.size(), handle);
  fclose(handle);
  return written == contents.size() ? Success : IOError;
}

static int CopyFileTo(const std::string& source, const std::string& destination) {
  std::ifstream input(source, std::ios::binary);
  if (!input) {
    return IOError;
  }
  std::ofstream output(destination, std::ios::binary);
  if (!output) {
    return IOError;
  }
  output << input.rdbuf();
  if (!output.good()) {
    return IOError;
  }
  return Success;
}

static void AssignPath(char* target, size_t capacity, const std::string& value) {
  if (capacity == 0) {
    return;
  }
  std::strncpy(target, value.c_str(), capacity - 1);
  target[capacity - 1] = '\0';
}

static void RemoveIfExists(const std::string& path) {
  if (path.empty()) {
    return;
  }
  remove(path.c_str());
}

static bool FileExists(const std::string& path) {
  struct stat info {};
  return stat(path.c_str(), &info) == 0 && S_ISREG(info.st_mode);
}

static int ApplyExplicitSiteSpecs(Structure* structure,
                                  const std::vector<PyDesignSiteSpec>& specs,
                                  Type_ResidueDesignType fallback_type) {
  for (const auto& spec : specs) {
    int chain_index = -1;
    for (int i = 0; i < StructureGetChainCount(structure); ++i) {
      if (strcmp(ChainGetName(StructureGetChain(structure, i)), spec.chain.c_str()) == 0) {
        chain_index = i;
        break;
      }
    }
    if (chain_index < 0) {
      return ValueError;
    }
    Chain* chain = StructureGetChain(structure, chain_index);
    int residue_index = -1;
    ChainFindResidueByPosInChain(chain, spec.position, &residue_index);
    if (residue_index < 0) {
      return ValueError;
    }
    Residue* residue = ChainGetResidue(chain, residue_index);
    Type_ResidueDesignType design_type =
        spec.design_type == Type_DesType_Fixed ? fallback_type : spec.design_type;
    ResidueSetDesignType(residue, design_type);
    if (!spec.allowed.empty()) {
      strncpy(residue->designAATypes, spec.allowed.c_str(), sizeof(residue->designAATypes) - 1);
    }
  }
  return Success;
}

static std::string ReadLastDesignLine(const std::string& path) {
  std::ifstream input(path);
  std::string line;
  std::string last;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    last = line;
  }
  return last;
}

static int ParseBestSequenceLine(const std::string& line, PyMonomerDesignResult* result) {
  if (line.empty()) {
    return ValueError;
  }
  std::istringstream stream(line);
  std::vector<std::string> parts;
  std::string token;
  while (stream >> token) {
    parts.push_back(token);
  }
  if (parts.size() < 8) {
    return ValueError;
  }

  result->sequence_string = parts[0];
  result->trajectory_index = std::atoi(parts[1].c_str());
  result->sequence_identity = std::atof(parts[2].c_str());
  result->energy_total = std::atof(parts[3].c_str());
  result->energy_evolution = std::atof(parts[4].c_str());
  result->energy_physical = std::atof(parts[5].c_str());
  result->energy_binding = std::atof(parts[6].c_str());
  result->unsatisfied_constraints = static_cast<int>(std::atoi(parts[7].c_str()));
  return Success;
}

}  // namespace

PyMonomerDesignResult::PyMonomerDesignResult()
    : sequence_string(),
      trajectory_index(0),
      sequence_identity(0.0),
      energy_total(0.0),
      energy_evolution(0.0),
      energy_physical(0.0),
      energy_binding(0.0),
      unsatisfied_constraints(0),
      has_best_structure(false),
      has_best_sites_structure(false),
      has_best_mutable_sites_structure(false) {
  StructureCreate(&best_structure);
  StructureCreate(&best_sites_structure);
  StructureCreate(&best_mutable_sites_structure);
}

PyMonomerDesignResult::~PyMonomerDesignResult() {
  StructureDestroy(&best_structure);
  StructureDestroy(&best_sites_structure);
  StructureDestroy(&best_mutable_sites_structure);
}

PyMinimizeResult::PyMinimizeResult() : has_structure(false) {
  StructureCreate(&minimized_structure);
}

PyMinimizeResult::~PyMinimizeResult() {
  StructureDestroy(&minimized_structure);
}

int RunMonomerDesignWorkflow(Structure* input_structure,
                             const PyMonomerDesignOptions& options,
                             PyMonomerDesignResult* result) {
  if (input_structure == nullptr || result == nullptr) {
    return ValueError;
  }

  ScopedOutputSilencer silencer(options.quiet_output);
  GlobalDesignStateGuard state_guard;

  Structure working_structure;
  StructureCreate(&working_structure);
  StructureCopy(&working_structure, input_structure);
  FixDesignSitePointers(&working_structure);

  AtomParamsSet atom_params;
  AtomParamsSetCreate(&atom_params);
  ResiTopoSet resi_topos;
  ResiTopoSetCreate(&resi_topos);
  AAppTable aapp_table{};
  RamaTable rama_table{};
  BBdepRotamerLib rotamer_lib;
  int code = Success;
  bool rotamer_lib_created = false;

  TempDirectory temp_dir(options.working_directory, "ud_");
  const std::string temp_path = temp_dir.path();

  std::string resfile_path;
  std::string atom_param_path;
  std::string topology_path;
  std::string weight_path;
  std::string aapp_path;
  std::string rama_path;
  std::string rotlib_path;
  std::string ligand_param_path;
  std::string ligand_topology_path;
  std::string ligand_mol2_path;
  std::string ligand_conformer_path;
  std::string ligand_pose_in_path;
  std::string ligand_pose_out_path;

  auto cleanup = [&](int status) {
    RemoveIfExists(resfile_path);
    RemoveIfExists(atom_param_path);
    RemoveIfExists(topology_path);
    RemoveIfExists(weight_path);
    RemoveIfExists(aapp_path);
    RemoveIfExists(rama_path);
    RemoveIfExists(rotlib_path);
    RemoveIfExists(ligand_param_path);
    RemoveIfExists(ligand_topology_path);
    RemoveIfExists(ligand_mol2_path);
    RemoveIfExists(ligand_conformer_path);
    RemoveIfExists(ligand_pose_in_path);
    RemoveIfExists(ligand_pose_out_path);
    RemoveIfExists(FILE_SELF_ENERGY);
    RemoveIfExists(FILE_ROTLIST);
    RemoveIfExists(FILE_ROTLIST_SEC);
    RemoveIfExists(std::string(FILE_BESTSEQS) + ".txt");
    RemoveIfExists(std::string(FILE_BESTSTRUCT) + "0001.pdb");
    RemoveIfExists(std::string(FILE_BEST_ALL_SITES) + "0001.pdb");
    RemoveIfExists(std::string(FILE_BEST_MUT_SITES) + "0001.pdb");
    if (rotamer_lib_created) {
      BBdepRotamerLibDestroy(&rotamer_lib);
    }
    AtomParamsSetDestroy(&atom_params);
    ResiTopoSetDestroy(&resi_topos);
    StructureDestroy(&working_structure);
    return status;
  };

  atom_param_path = temp_path + "/atom_params.prm";
  code = CopyFileTo(options.atom_params_path, atom_param_path);
  if (FAILED(code)) {
    return cleanup(code);
  }
  topology_path = temp_path + "/topology.inp";
  code = CopyFileTo(options.topology_path, topology_path);
  if (FAILED(code)) {
    return cleanup(code);
  }
  weight_path = temp_path + "/weights.wgt";
  code = CopyFileTo(options.weight_file, weight_path);
  if (FAILED(code)) {
    return cleanup(code);
  }
  aapp_path = temp_path + "/aapp.nrg";
  code = CopyFileTo(options.aapp_file, aapp_path);
  if (FAILED(code)) {
    return cleanup(code);
  }
  rama_path = temp_path + "/rama.nrg";
  code = CopyFileTo(options.rama_file, rama_path);
  if (FAILED(code)) {
    return cleanup(code);
  }
  rotlib_path = temp_path + "/rotlib.bin";
  code = CopyFileTo(options.rotlib_bin, rotlib_path);
  if (FAILED(code)) {
    return cleanup(code);
  }

  if (options.has_ligand) {
    ligand_mol2_path = temp_path + "/ligand.mol2";
    code = CopyFileTo(options.ligand_mol2, ligand_mol2_path);
    if (FAILED(code)) {
      return cleanup(code);
    }
    ligand_param_path = temp_path + "/ligand.prm";
    code = CopyFileTo(options.ligand_parameters, ligand_param_path);
    if (FAILED(code)) {
      return cleanup(code);
    }
    ligand_topology_path = temp_path + "/ligand.top";
    code = CopyFileTo(options.ligand_topology, ligand_topology_path);
    if (FAILED(code)) {
      return cleanup(code);
    }
    ligand_pose_in_path = temp_path + "/ligand_poses_in.pdb";
    ligand_pose_out_path = temp_path + "/ligand_poses_out.pdb";
    AssignPath(FILE_LIG_POSES_IN, sizeof(FILE_LIG_POSES_IN), ligand_pose_in_path);
    AssignPath(FILE_LIG_POSES_OUT, sizeof(FILE_LIG_POSES_OUT), ligand_pose_out_path);
    FLAG_LIG_POSES = TRUE;
    if (!options.ligand_conformers.empty()) {
      ligand_conformer_path = ligand_pose_in_path;
      code = CopyFileTo(options.ligand_conformers, ligand_conformer_path);
      if (FAILED(code)) {
        return cleanup(code);
      }
    } else {
      ligand_conformer_path.clear();
    }
  }

  AssignPath(PROGRAM_PATH, sizeof(PROGRAM_PATH), temp_path);
  AssignPath(FILE_ATOMPARAM, sizeof(FILE_ATOMPARAM), atom_param_path);
  AssignPath(FILE_TOPO, sizeof(FILE_TOPO), topology_path);
  AssignPath(FILE_WEIGHT_READ, sizeof(FILE_WEIGHT_READ), weight_path);
  AssignPath(FILE_AAPROPENSITY, sizeof(FILE_AAPROPENSITY), aapp_path);
  AssignPath(FILE_RAMACHANDRAN, sizeof(FILE_RAMACHANDRAN), rama_path);
  AssignPath(FILE_ROTLIB_BIN, sizeof(FILE_ROTLIB_BIN), rotlib_path);
  AssignPath(FILE_ROTLIB, sizeof(FILE_ROTLIB), rotlib_path);

  std::string prefix = temp_path + "/design";
  AssignPath(PREFIX, sizeof(PREFIX), prefix);
  snprintf(FILE_SELF_ENERGY, sizeof(FILE_SELF_ENERGY), "%s_selfenergy.txt", prefix.c_str());
  snprintf(FILE_ROTLIST, sizeof(FILE_ROTLIST), "%s_rotlist.txt", prefix.c_str());
  snprintf(FILE_ROTLIST_SEC, sizeof(FILE_ROTLIST_SEC), "%s_rotlistSEC.txt", prefix.c_str());
  snprintf(FILE_BESTSEQS, sizeof(FILE_BESTSEQS), "%s_bestseqs", prefix.c_str());
  snprintf(FILE_BESTSTRUCT, sizeof(FILE_BESTSTRUCT), "%s_beststruct", prefix.c_str());
  snprintf(FILE_BEST_ALL_SITES, sizeof(FILE_BEST_ALL_SITES), "%s_bestsites", prefix.c_str());
  snprintf(FILE_BEST_MUT_SITES, sizeof(FILE_BEST_MUT_SITES), "%s_bestmutsites", prefix.c_str());
  snprintf(FILE_BEST_LIG_MOL2, sizeof(FILE_BEST_LIG_MOL2), "%s_bestlig", prefix.c_str());

  FLAG_PROT_LIG = options.has_ligand ? TRUE : FALSE;
  FLAG_ENZYME = FALSE;
  FLAG_PPI = FALSE;
  FLAG_MONOMER = options.has_ligand ? FALSE : TRUE;
  FLAG_PHYSICS = TRUE;
  FLAG_EVOLUTION = FALSE;
  FLAG_EVOPHIPSI = FALSE;
  FLAG_BBDEP_ROTLIB = TRUE;
  FLAG_USE_INPUT_SC = options.use_input_sc ? TRUE : FALSE;
  FLAG_ROTATE_HYDROXYL = options.rotate_hydroxyl ? TRUE : FALSE;
  FLAG_WILDTYPE_ONLY = options.wildtype_only ? TRUE : FALSE;
  FLAG_INTERFACE_ONLY = options.interface_only ? TRUE : FALSE;
  FLAG_EXCL_CYS_ROTS = options.exclude_cys_rotamers ? TRUE : FALSE;
  FLAG_RESFILE = FALSE;
  FLAG_DESIGN_FROM_NATAA = options.design_from_native ? TRUE : FALSE;
  FLAG_READ_HYDROGEN = TRUE;
  FLAG_WRITE_HYDROGEN = TRUE;

  WGT_PROFILE = options.profile_weight;
  WGT_BIND = options.binding_weight;

  NTRAJ = options.trajectories > 0 ? options.trajectories : 1;
  NTRAJ_START_NDX = 1;

  std::string design_chains = options.design_chains.empty()
                                  ? DetermineDefaultDesignChains(&working_structure)
                                  : options.design_chains;
  if (design_chains.empty()) {
    return cleanup(ValueError);
  }
  strncpy(DES_CHAINS, design_chains.c_str(), sizeof(DES_CHAINS) - 1);

  code = AtomParameterRead(&atom_params, FILE_ATOMPARAM);
  if (options.has_ligand) {
    code = AtomParameterRead(&atom_params, const_cast<char*>(ligand_param_path.c_str()));
    if (FAILED(code)) {
      return cleanup(code);
    }
  }
  if (FAILED(code)) {
    return cleanup(code);
  }
  code = ResiTopoSetRead(&resi_topos, FILE_TOPO);
  if (options.has_ligand) {
    code = ResiTopoSetRead(&resi_topos, const_cast<char*>(ligand_topology_path.c_str()));
    if (FAILED(code)) {
      return cleanup(code);
    }
  }
  if (FAILED(code)) {
    return cleanup(code);
  }

  code = EnergyWeightRead(FILE_WEIGHT_READ);
  if (FAILED(code)) {
    // Weight file missing; continue with defaults.
    code = Success;
  }

  code = AApropensityTableReadFromFile(&aapp_table, FILE_AAPROPENSITY);
  if (FAILED(code)) {
    return cleanup(code);
  }
  code = RamaTableReadFromFile(&rama_table, FILE_RAMACHANDRAN);
  if (FAILED(code)) {
    return cleanup(code);
  }
  code = BBdepRotamerLibCreate2(&rotamer_lib, FILE_ROTLIB_BIN);
  if (FAILED(code)) {
    return cleanup(code);
  }
  rotamer_lib_created = true;

  code = StructureCalcAminoAcidPropensityAndRamaEnergy(&working_structure, &aapp_table, &rama_table);
  if (FAILED(code)) {
    return cleanup(code);
  }
  code = StructureCalcAminoAcidDunbrackEnergy(&working_structure, &rotamer_lib);
  if (FAILED(code)) {
    return cleanup(code);
  }

  if (!options.resfile_contents.empty()) {
    resfile_path = temp_path + "/design.res";
    code = WriteFile(resfile_path, options.resfile_contents);
    if (FAILED(code)) {
      return cleanup(code);
    }
    strncpy(FILE_RESFILE, resfile_path.c_str(), sizeof(FILE_RESFILE) - 1);
    FLAG_RESFILE = TRUE;
  } else {
    code = ApplyExplicitSiteSpecs(&working_structure, options.design_sites, Type_DesType_Mutable);
    if (FAILED(code)) {
      return cleanup(code);
    }
    code = ApplyExplicitSiteSpecs(&working_structure, options.repack_sites, Type_DesType_Repackable);
    if (FAILED(code)) {
      return cleanup(code);
    }
  }

  if (FLAG_INTERFACE_ONLY && FLAG_PPI) {
    code = StructureBuildPPIRotamersByBBdepRotLib(&working_structure, &rotamer_lib, &atom_params,
                                                  &resi_topos);
  } else if (FLAG_PROT_LIG) {
    code = StructureBuildPLIShell1RotamersByBBdepRotLib(&working_structure, &rotamer_lib,
                                                        &atom_params, &resi_topos, FILE_RESFILE);
    if (!FAILED(code)) {
      code = StructureBuildPLIShell2RotamersByBBdepRotLib(&working_structure, &rotamer_lib,
                                                          &atom_params, &resi_topos, FILE_RESFILE);
    }
  } else if (FLAG_ENZYME) {
    code = StructureBuildCatalyticRotamersByBBdepRotLib(&working_structure, &rotamer_lib,
                                                        &atom_params, &resi_topos, FILE_RESFILE);
    if (!FAILED(code)) {
      code = StructureBuildPLIShell1RotamersByBBdepRotLib(&working_structure, &rotamer_lib,
                                                          &atom_params, &resi_topos, FILE_RESFILE);
    }
    if (!FAILED(code)) {
      code = StructureBuildPLIShell2RotamersByBBdepRotLib(&working_structure, &rotamer_lib,
                                                          &atom_params, &resi_topos, FILE_RESFILE);
    }
  } else if (FLAG_RESFILE) {
    code = StructureBuildResfileRotamersByBBdepRotLib(&working_structure, &rotamer_lib, &atom_params,
                                                      &resi_topos, FILE_RESFILE);
  } else {
    code = StructureBuildAllRotamersByBBdepRotLib(&working_structure, &rotamer_lib, &atom_params,
                                                  &resi_topos);
  }
  if (FAILED(code)) {
    return cleanup(code);
  }

  if (options.has_ligand && (FLAG_PROT_LIG == TRUE || FLAG_ENZYME == TRUE)) {
    bool ligand_rotamers_loaded = false;
    if (FLAG_LIG_POSES && !ligand_pose_in_path.empty()) {
      code = StructureReadSmallMolRotamers(&working_structure, &resi_topos,
                                           const_cast<char*>(ligand_pose_in_path.c_str()));
      if (!FAILED(code)) {
        printf("read ligand poses from %s\n", ligand_pose_in_path.c_str());
        ligand_rotamers_loaded = true;
      }
    }
    if (!ligand_rotamers_loaded) {
      printf("use the ligand pose in mol2 file for design\n");
      FILE* ligand_pose_out = fopen(ligand_pose_out_path.c_str(), "w");
      if (ligand_pose_out == nullptr) {
        return cleanup(IOError);
      }
      Model(1, ligand_pose_out);
      Residue* small_molecule = nullptr;
      int sm_status = StructureFindSmallMol(&working_structure, &small_molecule);
      if (FAILED(sm_status) || small_molecule == nullptr) {
        fclose(ligand_pose_out);
        return cleanup(sm_status);
      }
      char atom_header[] = "ATOM";
      AtomArrayShowInPDBFormat(ResidueGetAllAtoms(small_molecule), atom_header,
                               ResidueGetName(small_molecule), ResidueGetChainName(small_molecule), 1,
                               ResidueGetPosInChain(small_molecule), ligand_pose_out);
      EndModel(ligand_pose_out);
      fclose(ligand_pose_out);
      code = StructureReadSmallMolRotamers(&working_structure, &resi_topos,
                                           const_cast<char*>(ligand_pose_out_path.c_str()));
      if (FAILED(code)) {
        return cleanup(code);
      }
    }
  }

  StructureShowDesignSites(&working_structure);

  std::string self_energy_file = FILE_SELF_ENERGY;
  code = SelfEnergyGenerate2(&working_structure, &aapp_table, &rama_table,
                             const_cast<char*>(self_energy_file.c_str()));
  if (FAILED(code)) {
    return cleanup(code);
  }

  RotamerList rotamer_list;
  RotamerListCreateFromStructure(&rotamer_list, &working_structure);

  std::string rotlist_file = FILE_ROTLIST;
  std::string rotlist_sec_file = FILE_ROTLIST_SEC;
  RotamerListWrite(&rotamer_list, const_cast<char*>(rotlist_file.c_str()));
  code = SelfEnergyReadAndCheck(&working_structure, &rotamer_list,
                                const_cast<char*>(self_energy_file.c_str()));
  if (FAILED(code)) {
    return cleanup(code);
  }
  CaptureResidueSelfEnergies(&working_structure, self_energy_file, &result->residue_self_energies);
  RotamerListWrite(&rotamer_list, const_cast<char*>(rotlist_sec_file.c_str()));
  RotamerListRead(&rotamer_list, const_cast<char*>(rotlist_sec_file.c_str()));
  StructureShowDesignSitesAfterRotamerDelete(&working_structure, &rotamer_list);

  code = SimulatedAnnealing(&working_structure, &rotamer_list);
  RotamerListDestroy(&rotamer_list);
  if (FAILED(code)) {
    return cleanup(code);
  }

  std::string best_seq_file = std::string(FILE_BESTSEQS) + ".txt";
  std::string line = ReadLastDesignLine(best_seq_file);
  code = ParseBestSequenceLine(line, result);
  if (FAILED(code)) {
    return cleanup(code);
  }

  std::string best_struct_file = std::string(FILE_BESTSTRUCT) + "0001.pdb";
  if (FileExists(best_struct_file)) {
    code = StructureReadPDB(&result->best_structure, const_cast<char*>(best_struct_file.c_str()),
                            &atom_params, &resi_topos);
    if (!FAILED(code)) {
      result->has_best_structure = true;
    }
  }

  std::string best_sites_file = std::string(FILE_BEST_ALL_SITES) + "0001.pdb";
  if (FileExists(best_sites_file)) {
    code = StructureReadPDB(&result->best_sites_structure,
                            const_cast<char*>(best_sites_file.c_str()), &atom_params, &resi_topos);
    if (!FAILED(code)) {
      result->has_best_sites_structure = true;
    }
  }

  std::string best_mut_file = std::string(FILE_BEST_MUT_SITES) + "0001.pdb";
  if (FileExists(best_mut_file)) {
    code = StructureReadPDB(&result->best_mutable_sites_structure,
                            const_cast<char*>(best_mut_file.c_str()), &atom_params, &resi_topos);
    if (!FAILED(code)) {
      result->has_best_mutable_sites_structure = true;
    }
  }

  return cleanup(Success);
}

int RunMinimizeWorkflow(Structure* input_structure,
                        const PyMinimizeOptions& options,
                        PyMinimizeResult* result) {
  if (input_structure == nullptr || result == nullptr) {
    return ValueError;
  }

  ScopedOutputSilencer silencer(options.quiet_output);
  GlobalDesignStateGuard state_guard;

  Structure working_structure;
  StructureCreate(&working_structure);
  StructureCopy(&working_structure, input_structure);
  FixDesignSitePointers(&working_structure);

  AtomParamsSet atom_params;
  AtomParamsSetCreate(&atom_params);
  ResiTopoSet resi_topos;
  ResiTopoSetCreate(&resi_topos);
  BBdepRotamerLib rotamer_lib;
  bool rotamer_lib_created = false;

  TempDirectory temp_dir(options.working_directory, "ud_");
  const std::string temp_path = temp_dir.path();

  std::string atom_param_path;
  std::string topology_path;
  std::string rotlib_path;
  std::string weight_path;
  std::string ligand_param_path;
  std::string ligand_topology_path;

  auto cleanup = [&](int status) {
    RemoveIfExists(atom_param_path);
    RemoveIfExists(topology_path);
    RemoveIfExists(rotlib_path);
    RemoveIfExists(weight_path);
    RemoveIfExists(ligand_param_path);
    RemoveIfExists(ligand_topology_path);
    if (rotamer_lib_created) {
      BBdepRotamerLibDestroy(&rotamer_lib);
    }
    AtomParamsSetDestroy(&atom_params);
    ResiTopoSetDestroy(&resi_topos);
    StructureDestroy(&working_structure);
    return status;
  };

  atom_param_path = temp_path + "/atom_params.prm";
  int code = CopyFileTo(options.atom_params_path, atom_param_path);
  if (FAILED(code)) {
    return cleanup(code);
  }
  topology_path = temp_path + "/topology.inp";
  code = CopyFileTo(options.topology_path, topology_path);
  if (FAILED(code)) {
    return cleanup(code);
  }
  rotlib_path = temp_path + "/rotlib.bin";
  code = CopyFileTo(options.rotlib_bin, rotlib_path);
  if (FAILED(code)) {
    return cleanup(code);
  }
  weight_path = temp_path + "/weights.wgt";
  code = CopyFileTo(options.weight_file, weight_path);
  if (FAILED(code)) {
    return cleanup(code);
  }

  if (options.has_ligand) {
    ligand_param_path = temp_path + "/ligand.prm";
    code = CopyFileTo(options.ligand_parameters, ligand_param_path);
    if (FAILED(code)) {
      return cleanup(code);
    }
    ligand_topology_path = temp_path + "/ligand.top";
    code = CopyFileTo(options.ligand_topology, ligand_topology_path);
    if (FAILED(code)) {
      return cleanup(code);
    }
  }

  AssignPath(PROGRAM_PATH, sizeof(PROGRAM_PATH), temp_path);
  AssignPath(FILE_ATOMPARAM, sizeof(FILE_ATOMPARAM), atom_param_path);
  AssignPath(FILE_TOPO, sizeof(FILE_TOPO), topology_path);
  AssignPath(FILE_WEIGHT_READ, sizeof(FILE_WEIGHT_READ), weight_path);
  AssignPath(FILE_ROTLIB_BIN, sizeof(FILE_ROTLIB_BIN), rotlib_path);
  AssignPath(FILE_ROTLIB, sizeof(FILE_ROTLIB), rotlib_path);

  FLAG_PROT_LIG = options.has_ligand ? TRUE : FALSE;
  FLAG_ENZYME = FALSE;
  FLAG_PPI = FALSE;
  FLAG_MONOMER = TRUE;
  FLAG_BBDEP_ROTLIB = TRUE;
  FLAG_USE_INPUT_SC = options.use_input_sc ? TRUE : FALSE;
  FLAG_ROTATE_HYDROXYL = options.rotate_hydroxyl ? TRUE : FALSE;
  FLAG_READ_HYDROGEN = TRUE;
  FLAG_WRITE_HYDROGEN = TRUE;
  FLAG_LIG_POSES = FALSE;

  code = EnergyWeightRead(FILE_WEIGHT_READ);
  if (FAILED(code)) {
    code = Success;
  }

  code = AtomParameterRead(&atom_params, FILE_ATOMPARAM);
  if (FAILED(code)) {
    return cleanup(code);
  }
  if (options.has_ligand) {
    code = AtomParameterRead(&atom_params, const_cast<char*>(ligand_param_path.c_str()));
    if (FAILED(code)) {
      return cleanup(code);
    }
  }

  code = ResiTopoSetRead(&resi_topos, FILE_TOPO);
  if (FAILED(code)) {
    return cleanup(code);
  }
  if (options.has_ligand) {
    code = ResiTopoSetRead(&resi_topos, const_cast<char*>(ligand_topology_path.c_str()));
    if (FAILED(code)) {
      return cleanup(code);
    }
  }

  code = StructureCalcPhiPsi(&working_structure);
  if (FAILED(code)) {
    return cleanup(code);
  }

  code = BBdepRotamerLibCreate2(&rotamer_lib, FILE_ROTLIB_BIN);
  if (FAILED(code)) {
    return cleanup(code);
  }
  rotamer_lib_created = true;

  code = StructureCalcAminoAcidDunbrackEnergy(&working_structure, &rotamer_lib);
  if (FAILED(code)) {
    return cleanup(code);
  }

  if (!options.design_sites.empty() || !options.repack_sites.empty()) {
    code = ApplyExplicitSiteSpecs(&working_structure, options.design_sites, Type_DesType_Mutable);
    if (FAILED(code)) {
      return cleanup(code);
    }
    code = ApplyExplicitSiteSpecs(&working_structure, options.repack_sites, Type_DesType_Repackable);
    if (FAILED(code)) {
      return cleanup(code);
    }
  }

  code = EnergyMinimizationByBBdepRotLib(&working_structure, &rotamer_lib, &atom_params, &resi_topos,
                                         NULL, options.respect_design_types);
  if (FAILED(code)) {
    return cleanup(code);
  }

  result->has_structure = true;
  StructureCopy(&result->minimized_structure, &working_structure);

  return cleanup(Success);
}
