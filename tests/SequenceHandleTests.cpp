#include "unidesign/SequenceHandle.hpp"

#include <cassert>
#include <utility>

double CUT_EXCL_LOW_PROB_ROT = 0.03;
double CUT_TORSION_DEVIATION = 20.0;
BOOL FLAG_READ_HYDROGEN = TRUE;
BOOL FLAG_WRITE_HYDROGEN = TRUE;
char MOL2[MAX_LEN_FILE_NAME + 1] = "";
char DES_CHAINS[MAX_LEN_ONE_LINE_CONTENT + 1] = "";

int main() {
  unidesign::SequenceHandle handle;
  handle->desSiteCount = 3;
  handle->etot = 4.2;

  unidesign::SequenceHandle copy(handle);
  assert(copy->desSiteCount == 3);
  assert(copy->etot == 4.2);

  unidesign::SequenceHandle assigned;
  assigned = handle;
  assert(assigned->desSiteCount == 3);
  assert(assigned->etot == 4.2);

  unidesign::SequenceHandle moved(std::move(handle));
  assert(moved->desSiteCount == 3);
  assert(moved->etot == 4.2);
  assert(handle.get() == nullptr);

  unidesign::SequenceHandle moveAssigned;
  moveAssigned = std::move(copy);
  assert(moveAssigned->desSiteCount == 3);
  assert(moveAssigned->etot == 4.2);
  assert(copy.get() == nullptr);

  return 0;
}

