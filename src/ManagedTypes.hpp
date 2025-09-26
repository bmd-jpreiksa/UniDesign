#ifndef MANAGED_TYPES_HPP
#define MANAGED_TYPES_HPP

#include "Structure.h"
#include "Sequence.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace managed_types_detail
{
inline void ThrowIfFailed(int code, const char* functionName)
{
  if (FAILED(code))
  {
    throw std::runtime_error(std::string(functionName) + " failed with error code " + std::to_string(code));
  }
}
} // namespace managed_types_detail

class StructureHandle
{
public:
  StructureHandle()
  {
    managed_types_detail::ThrowIfFailed(StructureCreate(&mValue), "StructureCreate");
  }

  ~StructureHandle()
  {
    StructureDestroy(&mValue);
  }

  StructureHandle(const StructureHandle& other)
  {
    managed_types_detail::ThrowIfFailed(StructureCreate(&mValue), "StructureCreate");
    managed_types_detail::ThrowIfFailed(StructureCopy(&mValue, const_cast<Structure*>(&other.mValue)), "StructureCopy");
  }

  StructureHandle(StructureHandle&& other) noexcept
  {
    managed_types_detail::ThrowIfFailed(StructureCreate(&mValue), "StructureCreate");
    swap(other);
  }

  StructureHandle& operator=(StructureHandle other) noexcept
  {
    swap(other);
    return *this;
  }

  void swap(StructureHandle& other) noexcept
  {
    using std::swap;
    swap(mValue, other.mValue);
  }

  Structure* get() noexcept { return &mValue; }
  const Structure* get() const noexcept { return &mValue; }

  Structure& operator*() noexcept { return mValue; }
  const Structure& operator*() const noexcept { return mValue; }

  Structure* operator->() noexcept { return &mValue; }
  const Structure* operator->() const noexcept { return &mValue; }

private:
  Structure mValue{};
};

inline void swap(StructureHandle& lhs, StructureHandle& rhs) noexcept
{
  lhs.swap(rhs);
}

class SequenceHandle
{
public:
  SequenceHandle()
  {
    managed_types_detail::ThrowIfFailed(SequenceCreate(&mValue), "SequenceCreate");
  }

  ~SequenceHandle()
  {
    SequenceDestroy(&mValue);
  }

  SequenceHandle(const SequenceHandle& other)
  {
    managed_types_detail::ThrowIfFailed(SequenceCreate(&mValue), "SequenceCreate");
    managed_types_detail::ThrowIfFailed(SequenceCopy(&mValue, const_cast<Sequence*>(&other.mValue)), "SequenceCopy");
  }

  SequenceHandle(SequenceHandle&& other) noexcept
  {
    managed_types_detail::ThrowIfFailed(SequenceCreate(&mValue), "SequenceCreate");
    swap(other);
  }

  SequenceHandle& operator=(SequenceHandle other) noexcept
  {
    swap(other);
    return *this;
  }

  void swap(SequenceHandle& other) noexcept
  {
    using std::swap;
    swap(mValue, other.mValue);
  }

  Sequence* get() noexcept { return &mValue; }
  const Sequence* get() const noexcept { return &mValue; }

  Sequence& operator*() noexcept { return mValue; }
  const Sequence& operator*() const noexcept { return mValue; }

  Sequence* operator->() noexcept { return &mValue; }
  const Sequence* operator->() const noexcept { return &mValue; }

private:
  Sequence mValue{};
};

inline void swap(SequenceHandle& lhs, SequenceHandle& rhs) noexcept
{
  lhs.swap(rhs);
}

struct StructureDeleter
{
  void operator()(Structure* ptr) const noexcept
  {
    if (ptr != nullptr)
    {
      StructureDestroy(ptr);
      delete ptr;
    }
  }
};

using StructurePtr = std::unique_ptr<Structure, StructureDeleter>;

inline StructureHandle MakeStructure()
{
  return StructureHandle{};
}

inline StructurePtr MakeStructurePtr()
{
  StructurePtr ptr(new Structure{});
  managed_types_detail::ThrowIfFailed(StructureCreate(ptr.get()), "StructureCreate");
  return ptr;
}

struct SequenceDeleter
{
  void operator()(Sequence* ptr) const noexcept
  {
    if (ptr != nullptr)
    {
      SequenceDestroy(ptr);
      delete ptr;
    }
  }
};

using SequencePtr = std::unique_ptr<Sequence, SequenceDeleter>;

inline SequenceHandle MakeSequence()
{
  return SequenceHandle{};
}

inline SequencePtr MakeSequencePtr()
{
  SequencePtr ptr(new Sequence{});
  managed_types_detail::ThrowIfFailed(SequenceCreate(ptr.get()), "SequenceCreate");
  return ptr;
}

#endif // MANAGED_TYPES_HPP
