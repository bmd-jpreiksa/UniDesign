# UniDesign

## Introduction
UniDesign is a computational framework for protein design, targeting a diversity of protein design and engineering tasks, and it can also be used for protein structure modeling and scoring. UniDesign is extended mainly from the EvoEF2 protein design program with the evolutionary feature taken from EvoDesign, but it has many new features. More details can be found in <a href=https://github.com/tommyhuangthu/UniDesign/blob/master/manual.docx>manual.docx</a> within this package.
<br>
<p align="center">
  <img src="img/EnzymeDesign.png" width="80%" height="80%">
  <br>
  Fig 1. UniDesign's capacities and its workflow for enzyme design
</p>

## Main Applications

#### Protein Design

• Design monomer protein <br>
• Design protein-protein interaction <br>
• Design protein-ligand interaction <br>
• Design protein-nucleic acid interaction <br>
• Design enzymes <br>
<br>
All protein design tasks are conducted with UniDesign command <i><b>ProteinDesign</b></i> with specific options. See <a href=https://github.com/tommyhuangthu/UniDesign/blob/master/manual.docx>manual.docx</a> for details and examples.

#### Protein Structure Modeling

• Protein side-chain packing  (command <i><b>ProteinDesign</b></i> with option <i><b>--wildtype</b></i>) <br>
• Repair incomplete protein sidechains  (command <i><b>RepairStructure</b></i>) <br>
• Protein minimization to remove sidechain clashes  (command <i><b>Minimize</b></i>) <br>
• Build mutant structural models  (command <i><b>BuildMutant</b></i>) <br>
• Add polar hydrogen atoms  (command <i><b>AddPolarHydrogen</b></i>) <br>
• Optimize hydrogen atom's position  (command <i><b>OptimizeHydrogen</b></i>) <br>

#### Protein Scoring

• Compute protein stability (command <i><b>ComputeStability</b></i>) <br>
• Compute protein interchain binding interaction (command <i><b>ComputeBinding</b></i>) <br>

## Installation, Usage and Tutorial
Please refer to the <a href=https://github.com/tommyhuangthu/UniDesign/blob/master/manual.docx>manual.docx</a> for details.

## Managed C++ helpers for Structure and Sequence
Repeated manual calls to `StructureCreate`/`StructureDestroy` and `SequenceCreate`/`SequenceDestroy` (for example, the loops in `src/ProteinDesign.cpp` that rebuild temporary sequences) can now be replaced with RAII helpers declared in `src/ManagedTypes.hpp`. The `StructureHandle` and `SequenceHandle` classes automatically invoke the appropriate `Create` routine on construction, perform `Copy` on copy/move, and guarantee `Destroy` is executed exactly once, preventing double frees.

```cpp
#include "ManagedTypes.hpp"

void OptimizeDesign()
{
  StructureHandle structure;           // StructureCreate + StructureDestroy are automatic
  SequenceHandle workingSequence;

  StructureHandle reference = structure; // uses StructureCopy under the hood
  workingSequence = SequenceHandle{};    // resets via copy-and-swap semantics

  // Access the underlying C structs when calling existing APIs.
  PerformDesign(structure.get(), workingSequence.get());
}
```

For code exposed through pybind11, the handles can be bound directly while still returning raw pointers to existing C APIs:

```cpp
py::class_<StructureHandle>(m, "StructureHandle")
  .def(py::init<>())
  .def("raw", [](StructureHandle& self) { return self.get(); },
       py::return_value_policy::reference);
```

Developers that prefer pointer semantics can also use `StructurePtr`/`SequencePtr`, which are `std::unique_ptr` aliases with custom deleters that call the corresponding `Destroy` routine before deleting the raw pointer.

## Contact
For suggestions, please contact xiaoqiah@umich.edu or xiaoqiah@outlook.com.

## References
Huang, X., Zhou, J., Yang, D., Zhang, J., Xia, X., Chen, Y. E., and Xu, J. Decoding CRISPR–Cas PAM recognition with UniDesign. Briefings in Bioinformatics. 2023, 24(3):bbad133. doi: 10.1093/bib/bbad133 <br><br>
Huang, X., Sun, Y., Osawa, Y., Chen, Y. E., and Zhang, H. Computational Redesign of Cytochrome P450 CYP102A1 for Highly Stereoselective Omeprazole Hydroxylation by UniDesign. Journal of Biological Chemistry. 2023, 299(8):105050.

