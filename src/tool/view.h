
#ifndef TOOL_VIEW_H
#define TOOL_VIEW_H

#include "utils.h"
#include <array>

enum class ContiguousLayout {
  Fortran, C
};

/// @brief Represents a non-owning view around a 3D array
template<typename T>
class View{
  T* ptr_;
  std::array<int, 3> shape_;
  std::array<int, 3> strides_;

  void check_invariants_() const {
    if (ptr_ == nullptr) {
      for (int i = 0; i < 3; i++){
        GRCLI_REQUIRE(shape_[i] == 0 && strides_[i] == 0, 
                     "shape and strides must be 0 for a nullptr");
      }
    } else {
      for (int i = 0; i < 3; i++){
        GRCLI_REQUIRE(shape_[i] > 0 && strides_[i] >= 0,
                      "shape must be positive & strides must be non-negative");
      }
    }
  }

public:

  View()
    : ptr_(nullptr), shape_{0,0,0}, strides_{0,0,0}
  {}

  View(T* ptr, std::array<int, 3> shape, std::array<int,3> strides)
    : ptr_(ptr), shape_{shape}, strides_{strides}
  { }

  /// Construct a new View from a contiguous region of memory where
  /// the shape is structured such that axis 0 is fast (like in fortran)
  View(T* ptr, int nx, int ny, int nz, ContiguousLayout layout)
    : ptr_(ptr), shape_{nx, ny, nz}
  {
    if (layout == ContiguousLayout::Fortran) {
      strides_[0] = 1;
      strides_[1] = strides_[0] * shape_[0];
      strides_[2] = strides_[1] * shape_[1];
    } else if (layout == ContiguousLayout::C) {
      strides_[2] = 1;
      strides_[1] = strides_[2] * shape_[2];
      strides_[0] = strides_[1] * shape_[1];
    } else {
      GRCLI_ERROR("Unhandled case!");
    }
  }

  View(T* ptr, std::array<int,3> shape, ContiguousLayout layout)
    : View(ptr, shape[0], shape[1], shape[2], layout)
  {}

  T* data() const noexcept {
    return ptr_;
  }

  int stride(int i) const {
    GRCLI_REQUIRE(i >= 0 && i <= 2, "i must be 0, 1, or 2");
    return strides_[i];
  }

  int extent(int i) const {
    GRCLI_REQUIRE(i >= 0 && i <= 2, "i must be 0, 1, or 2");
    return shape_[i];
  }

  inline T& operator()(int i, int j, int k) const noexcept {
    return ptr_[i*strides_[0] + j*strides_[1] + k*strides_[2]];
  }

};

#endif /* TOOL_VIEW_H */
