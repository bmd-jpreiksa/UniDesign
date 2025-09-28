#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include "Sequence.h"

namespace unidesign {

class SequenceError : public std::runtime_error {
public:
  SequenceError(const std::string& context, int code)
      : std::runtime_error(context + " failed with error code " + std::to_string(code)), code_(code) {}

  int code() const noexcept { return code_; }

private:
  int code_;
};

class SequenceHandle {
public:
  SequenceHandle()
      : seq_(allocateSequence()) {}

  SequenceHandle(const SequenceHandle& other)
      : seq_(nullptr) {
    copyFrom(other);
  }

  SequenceHandle& operator=(const SequenceHandle& other) {
    if (this != &other) {
      copyFrom(other);
    }
    return *this;
  }

  SequenceHandle(SequenceHandle&& other) noexcept = default;
  SequenceHandle& operator=(SequenceHandle&& other) noexcept = default;

  ~SequenceHandle() { reset(); }

  Sequence* get() noexcept { return seq_.get(); }
  const Sequence* get() const noexcept { return seq_.get(); }

  Sequence& operator*() {
    if (!seq_) {
      throw std::runtime_error("SequenceHandle does not own a Sequence");
    }
    return *seq_;
  }

  const Sequence& operator*() const {
    if (!seq_) {
      throw std::runtime_error("SequenceHandle does not own a Sequence");
    }
    return *seq_;
  }

  Sequence* operator->() noexcept { return seq_.get(); }
  const Sequence* operator->() const noexcept { return seq_.get(); }

private:
  struct SequenceStorageDeleter {
    void operator()(Sequence* seq) const noexcept { delete seq; }
  };

  using Storage = std::unique_ptr<Sequence, SequenceStorageDeleter>;

  static Storage allocateSequence() {
    Storage seq(new Sequence{});
    const int result = SequenceCreate(seq.get());
    if (result != Success) {
      SequenceDestroy(seq.get());
      throw SequenceError("SequenceCreate", result);
    }
    return seq;
  }

  void ensureAllocated() {
    if (!seq_) {
      seq_ = allocateSequence();
    }
  }

  void reset() noexcept {
    if (seq_) {
      SequenceDestroy(seq_.get());
      seq_.reset();
    }
  }

  void copyFrom(const SequenceHandle& other) {
    if (!other.seq_) {
      reset();
      ensureAllocated();
      return;
    }

    ensureAllocated();
    const int result = SequenceCopy(seq_.get(), const_cast<Sequence*>(other.seq_.get()));
    if (result != Success) {
      SequenceDestroy(seq_.get());
      const int recreate = SequenceCreate(seq_.get());
      if (recreate != Success) {
        throw SequenceError("SequenceCreate", recreate);
      }
      throw SequenceError("SequenceCopy", result);
    }
  }

  Storage seq_;
};

}  // namespace unidesign

