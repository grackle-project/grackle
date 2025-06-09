// See LICENSE file for license and copyright information

#include "gtest/gtest.h"
#include "visitor/common.hpp"
#include "visitor/copy.hpp"
#include "visitor/memory.hpp"

// the idea here is to set up a very simple struct and then test that the
// visitors do what they are supposed to do

struct MyDummyStruct {
  double* attr1;
  long long* attr2;
};

namespace grackle {
/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template <class BinaryVisitor>
void visit_member_pair(MyDummyStruct& obj0, MyDummyStruct& obj1,
                       BinaryVisitor f) {
  impl::visitor::begin_visit("MyDummyStruct", f);
  f(VIS_MEMBER_NAME("attr1"), obj0.attr1, obj1.attr1,
    grackle::impl::visitor::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("attr2"), obj0.attr2, obj1.attr2,
    grackle::impl::visitor::idx_range_len_multiple(1));
  impl::visitor::end_visit(f);
}

template <class UnaryFn>
void visit_member(MyDummyStruct* obj, UnaryFn fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, MyDummyStruct, obj, fn)
}

}  // namespace grackle

TEST(VisitorTest, SimpleTest) {
  grackle::impl::visitor::VisitorCtx ctx{2u};

  MyDummyStruct tmp{nullptr, nullptr};

  grackle::visit_member(&tmp, grackle::impl::visitor::AllocateMembers{ctx});
  ASSERT_NE(tmp.attr1, nullptr);
  ASSERT_NE(tmp.attr2, nullptr);

  tmp.attr1[0] = 1.0;
  tmp.attr1[1] = 2.0;
  tmp.attr2[0] = -2;
  tmp.attr2[1] = -1;

  MyDummyStruct tmp2{nullptr, nullptr};
  grackle::visit_member(&tmp2, grackle::impl::visitor::AllocateMembers{ctx});
  grackle::visit_member_pair(tmp, tmp2,
                             grackle::impl::visitor::CopyMembers{ctx});

  ASSERT_EQ(tmp2.attr1[0], 1.0);
  ASSERT_EQ(tmp2.attr1[1], 2.0);
  ASSERT_EQ(tmp2.attr2[0], -2);
  ASSERT_EQ(tmp2.attr2[1], -1);

  grackle::visit_member(&tmp, grackle::impl::visitor::FreeMembers{});
  grackle::visit_member(&tmp2, grackle::impl::visitor::FreeMembers{});
}