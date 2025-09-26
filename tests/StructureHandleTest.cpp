#include <cassert>
#include <cstring>
#include <utility>

#include "Structure.h"
#include "unidesign/StructureHandle.hpp"

int main() {
  using unidesign::StructureHandle;

  StructureHandle handle;
  char name[] = "TEST";
  assert(StructureSetName(handle.get(), name) == Success);
  assert(std::strcmp(StructureGetName(handle.get()), "TEST") == 0);

  StructureHandle copy(handle);
  assert(std::strcmp(StructureGetName(copy.get()), "TEST") == 0);

  StructureHandle assigned;
  assigned = copy;
  assert(std::strcmp(StructureGetName(assigned.get()), "TEST") == 0);

  StructureHandle moved(std::move(copy));
  assert(std::strcmp(StructureGetName(moved.get()), "TEST") == 0);

  StructureHandle movedAssign;
  movedAssign = std::move(moved);
  assert(std::strcmp(StructureGetName(movedAssign.get()), "TEST") == 0);

  Structure* raw = movedAssign.release();
  assert(raw != nullptr);
  assert(StructureDestroy(raw) == Success);

  return 0;
}
