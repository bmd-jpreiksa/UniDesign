#pragma once

#include <cstring>
#include <stdexcept>
#include <utility>

#include "ErrorTracker.h"
#include "Structure.h"

namespace unidesign {

class StructureHandle {
public:
  StructureHandle() {
    if (StructureCreate(&structure_) != Success) {
      throw std::runtime_error("StructureCreate failed");
    }
    owns_ = true;
  }

  StructureHandle(const StructureHandle& other) : StructureHandle() {
    if (StructureCopy(&structure_, const_cast<Structure*>(&other.structure_)) != Success) {
      reset();
      throw std::runtime_error("StructureCopy failed");
    }
  }

  StructureHandle& operator=(const StructureHandle& other) {
    if (this != &other) {
      StructureHandle temp(other);
      swap(temp);
    }
    return *this;
  }

  StructureHandle(StructureHandle&& other) noexcept
      : structure_(other.structure_), owns_(other.owns_) {
    other.owns_ = false;
    std::memset(&other.structure_, 0, sizeof(Structure));
  }

  StructureHandle& operator=(StructureHandle&& other) noexcept {
    if (this != &other) {
      reset();
      structure_ = other.structure_;
      owns_ = other.owns_;
      other.owns_ = false;
      std::memset(&other.structure_, 0, sizeof(Structure));
    }
    return *this;
  }

  ~StructureHandle() {
    reset();
  }

  Structure* get() noexcept { return &structure_; }
  const Structure* get() const noexcept { return &structure_; }

  Structure& operator*() noexcept { return structure_; }
  const Structure& operator*() const noexcept { return structure_; }

  Structure* operator->() noexcept { return &structure_; }
  const Structure* operator->() const noexcept { return &structure_; }

  Structure* release() noexcept {
    owns_ = false;
    return &structure_;
  }

  void swap(StructureHandle& other) noexcept {
    using std::swap;
    swap(structure_, other.structure_);
    swap(owns_, other.owns_);
  }

private:
  void reset() noexcept {
    if (owns_) {
      if (StructureDestroy(&structure_) == Success) {
        owns_ = false;
        std::memset(&structure_, 0, sizeof(Structure));
      }
    }
  }

  Structure structure_{};
  bool owns_{false};
};

inline void swap(StructureHandle& lhs, StructureHandle& rhs) noexcept {
  lhs.swap(rhs);
}

}  // namespace unidesign

