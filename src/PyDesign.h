/*******************************************************************************************************************************
Copyright (c) Xiaoqiang Huang

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

#ifndef PY_DESIGN_H
#define PY_DESIGN_H

#include <string>
#include <vector>

#include "ProgramFunction.h"
#include "ProteinDesign.h"
#include "Structure.h"

struct PyDesignSiteSpec {
  std::string chain;
  int position;
  std::string allowed;
  Type_ResidueDesignType design_type;
};

struct PyMonomerDesignOptions {
  std::string program_path;
  std::string working_directory;
  std::string atom_params_path;
  std::string topology_path;
  std::string weight_file;
  std::string aapp_file;
  std::string rama_file;
  std::string rotlib_bin;
  std::string resfile_contents;
  std::string design_chains;

  double profile_weight;
  double binding_weight;
  bool debug_only;
  double rotamer_probability_cutoff;

  int trajectories;

  bool interface_only;
  bool design_from_native;
  bool use_input_sc;
  bool rotate_hydroxyl;
  bool exclude_cys_rotamers;
  bool wildtype_only;
  bool quiet_output;

  bool has_ligand = false;
  std::string ligand_mol2;
  std::string ligand_parameters;
  std::string ligand_topology;
  std::string ligand_conformers;

  std::vector<PyDesignSiteSpec> design_sites;
  std::vector<PyDesignSiteSpec> repack_sites;
};

struct PyResidueSelfEnergy {
  std::string chain_name;
  int position;
  double self_energy;
  double binding_energy;
};

struct PyMonomerDesignResult {
  std::string sequence_string;
  int trajectory_index;
  double sequence_identity;
  double energy_total;
  double energy_evolution;
  double energy_physical;
  double energy_binding;
  int unsatisfied_constraints;

  bool has_best_structure;
  bool has_best_sites_structure;
  bool has_best_mutable_sites_structure;
  bool has_energy_terms_initial;
  bool has_energy_terms_final;
  double energy_terms_initial[MAX_ENERGY_TERM];
  double energy_terms_final[MAX_ENERGY_TERM];

  Structure best_structure;
  Structure best_sites_structure;
  Structure best_mutable_sites_structure;
  std::vector<PyResidueSelfEnergy> residue_self_energies;

  PyMonomerDesignResult();
  ~PyMonomerDesignResult();
};

int RunMonomerDesignWorkflow(Structure* input_structure,
                             const PyMonomerDesignOptions& options,
                             PyMonomerDesignResult* result);

struct PyMinimizeOptions {
  std::string program_path;
  std::string working_directory;
  std::string atom_params_path;
  std::string topology_path;
  std::string rotlib_bin;
  std::string weight_file;

  bool use_input_sc;
  bool rotate_hydroxyl;
  bool quiet_output;
  bool respect_design_types;

  bool has_ligand = false;
  std::string ligand_mol2;
  std::string ligand_parameters;
  std::string ligand_topology;

  std::vector<PyDesignSiteSpec> design_sites;
  std::vector<PyDesignSiteSpec> repack_sites;
};

struct PyMinimizeResult {
  bool has_structure;
  Structure minimized_structure;

  PyMinimizeResult();
  ~PyMinimizeResult();
};

int RunMinimizeWorkflow(Structure* input_structure,
                        const PyMinimizeOptions& options,
                        PyMinimizeResult* result);

#endif /* PY_DESIGN_H */
