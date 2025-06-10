// See LICENSE file for license and copyright information

/// @file printer.hpp
/// @brief declare/implement the Printer visitor

#ifndef VISITOR_PRINTER_HPP
#define VISITOR_PRINTER_HPP

#include "common.hpp"
#include "status_reporting.h"  // GR_INTERNAL_REQUIRE
#include <cstdio>              // printf

namespace grackle::impl::visitor {

// there's probably an opportunity to reduce some code duplication with
// grackle::impl::print_contiguous_row_

/// A visitor that prints the contents of a type (intended for debugging
/// purposes)
///
/// @note
/// This is implemented very inefficiently!
///
/// @note
/// If we want to make it possible to customize the way we write stuff or the
/// destination, or the format, that's all fine as long as we retain a version
/// using vanilla printf, since vanilla printf is important for debugging on
/// GPUs (with that said, we would probably need to dramatically cut down on
/// the number of calls to printf since many GPU threads execute at one time)
class Printer {
  VisitorCtx ctx;
  int indent_width;
  /// modified by previsit_struct_member and begin_visit
  const char* name_of_struct_;

  template <typename Visitor>
  friend void previsit_struct_member(const char*, Visitor&);

  template <typename Visitor>
  friend void begin_visit(const char*, Visitor&);

  template <typename Visitor>
  friend void end_visit(Visitor&);

  void print_indent_() const {
    if (indent_width > 0) {
      std::printf("%*c", indent_width - 1, ' ');
    }
  }

  void print_prefix_(const char* name) const {
    print_indent_();
    if (name != nullptr) {
      std::printf(".%s = ", name);
    }
  }

public:
  Printer() = delete;
  Printer(VisitorCtx ctx, int indent_width)
      : ctx(ctx), indent_width(indent_width), name_of_struct_(nullptr) {}

  void operator()(const char* name, int*& val, const BufLenSpec& spec) const {
    print_prefix_(name);
    if (val == nullptr) {
      std::printf("nullptr");
    } else {
      std::printf("[");
      std::printf("%d", val[0]);
      std::size_t n_elem = get_buf_len(spec, ctx);
      for (std::size_t i = 1; i < n_elem; i++) {
        std::printf(", %d", val[i]);
      }
      std::printf("]");
    }
    std::printf(",\n");
  }

  void operator()(const char* name, long long*& val,
                  const BufLenSpec& spec) const {
    print_prefix_(name);
    if (val == nullptr) {
      std::printf("nullptr");
    } else {
      std::printf("[");
      std::printf("%lld", val[0]);
      std::size_t n_elem = get_buf_len(spec, ctx);
      for (std::size_t i = 1; i < n_elem; i++) {
        std::printf(", %lld", val[i]);
      }
      std::printf("]");
    }
    std::printf(",\n");
  }

  void operator()(const char* name, double*& val,
                  const BufLenSpec& spec) const {
    print_prefix_(name);
    if (val == nullptr) {
      std::printf("nullptr");
    } else {
      std::printf("[");
      std::printf("%e", val[0]);
      std::size_t n_elem = get_buf_len(spec, ctx);
      for (std::size_t i = 1; i < n_elem; i++) {
        std::printf(", %e", val[i]);
      }
      std::printf("]");
    }
    std::printf(",\n");
  }
};

template <>
inline void previsit_struct_member(const char* name, Printer& visitor) {
  GR_INTERNAL_REQUIRE(visitor.name_of_struct_ == nullptr,
                      "sanity check failed!");
  visitor.name_of_struct_ = name;
}

template <>
inline void begin_visit(const char* type_name, Printer& visitor) {
  visitor.print_prefix_(visitor.name_of_struct_);
  if (type_name) {
    std::printf("%s {\n", type_name);
  } else {
    std::printf("{\n");
  }
  visitor.indent_width += 2;
  visitor.name_of_struct_ = nullptr;
}

template <>
inline void end_visit(Printer& visitor) {
  visitor.indent_width -= 2;
  visitor.print_indent_();
  std::printf("}\n");
}

}  // namespace grackle::impl::visitor

#endif  // VISITOR_PRINTER_HPP
