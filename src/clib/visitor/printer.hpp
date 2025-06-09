// See LICENSE file for license and copyright information

/// @file printer.hpp
/// @brief declare/implement the Printer visitor

#ifndef VISITOR_PRINTER_HPP
#define VISITOR_PRINTER_HPP

#include "common.hpp"
#include <cstdio>  // printf

namespace grackle::impl::visitor {

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

  template <typename Visitor>
  friend void begin_visit(const char* type_name, Visitor& visitor);

  template <typename Visitor>
  friend void end_visit(Visitor& visitor);

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
      : ctx(ctx), indent_width(indent_width) {}

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
inline void begin_visit(const char* type_name, Printer& visitor) {
  visitor.print_indent_();
  std::printf("{\n");
  visitor.indent_width += 2;
}

template <>
inline void end_visit(Printer& visitor) {
  visitor.indent_width -= 2;
  visitor.print_indent_();
  std::printf("}\n");
}

}  // namespace grackle::impl::visitor

#endif  // VISITOR_PRINTER_HPP
