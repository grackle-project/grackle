// See LICENSE file for license and copyright information

#include "gtest/gtest.h"
#include "visitor/common.hpp"
#include "visitor/copy.hpp"
#include "visitor/memory.hpp"

// uncomment if you want to test printing
// #include "visitor/printer.hpp"

// the idea here is to set up a very simple struct and then test that the
// visitors do what they are supposed to do

struct MyDummyStruct {
  double* attr1;
  long long* attr2;
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template <class BinaryVisitor>
void visit_member_pair(MyDummyStruct& obj0, MyDummyStruct& obj1,
                       BinaryVisitor f) {
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("MyDummyStruct", f);
  f(VIS_MEMBER_NAME("attr1"), obj0.attr1, obj1.attr1,
    vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("attr2"), obj0.attr2, obj1.attr2,
    vis::idx_range_len_multiple(1));
  vis::end_visit(f);
}

template <class UnaryFn>
void visit_member(MyDummyStruct* obj, UnaryFn fn){
    GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, MyDummyStruct, obj, fn)}

TEST(VisitorTest, Simple) {
  grackle::impl::visitor::VisitorCtx ctx{2u};

  MyDummyStruct tmp{nullptr, nullptr};

  visit_member(&tmp, grackle::impl::visitor::AllocateMembers{ctx});
  ASSERT_NE(tmp.attr1, nullptr);
  ASSERT_NE(tmp.attr2, nullptr);

  tmp.attr1[0] = 1.0;
  tmp.attr1[1] = 2.0;
  tmp.attr2[0] = -2;
  tmp.attr2[1] = -1;

  MyDummyStruct tmp2{nullptr, nullptr};
  visit_member(&tmp2, grackle::impl::visitor::AllocateMembers{ctx});
  visit_member_pair(tmp, tmp2, grackle::impl::visitor::CopyMembers{ctx});

  ASSERT_EQ(tmp2.attr1[0], 1.0);
  ASSERT_EQ(tmp2.attr1[1], 2.0);
  ASSERT_EQ(tmp2.attr2[0], -2);
  ASSERT_EQ(tmp2.attr2[1], -1);

  visit_member(&tmp, grackle::impl::visitor::FreeMembers{});
  visit_member(&tmp2, grackle::impl::visitor::FreeMembers{});
}

struct MyNestedStruct {
  int* attr0;
  MyDummyStruct dummy;
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template <class BinaryVisitor>
void visit_member_pair(MyNestedStruct& obj0, MyNestedStruct& obj1,
                       BinaryVisitor f) {
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("MyNestedStruct", f);
  f(VIS_MEMBER_NAME("attr0"), obj0.attr0, obj1.attr0,
    vis::idx_range_len_multiple(1));
  vis::previsit_struct_member(VIS_MEMBER_NAME("dummy"), f);
  visit_member_pair(obj0.dummy, obj1.dummy, f);
  vis::end_visit(f);
}

template <class UnaryFn>
void visit_member(MyNestedStruct* obj, UnaryFn fn){
    GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, MyNestedStruct, obj, fn)}

TEST(VisitorTest, Nested) {
  grackle::impl::visitor::VisitorCtx ctx{2u};

  MyNestedStruct tmp{nullptr, {nullptr, nullptr}};

  visit_member(&tmp, grackle::impl::visitor::AllocateMembers{ctx});
  ASSERT_NE(tmp.attr0, nullptr);
  ASSERT_NE(tmp.dummy.attr1, nullptr);
  ASSERT_NE(tmp.dummy.attr2, nullptr);

  tmp.attr0[0] = 100;
  tmp.attr0[1] = -100;
  tmp.dummy.attr1[0] = 1.0;
  tmp.dummy.attr1[1] = 2.0;
  tmp.dummy.attr2[0] = -2;
  tmp.dummy.attr2[1] = -1;

  // to print this out, you just need to include "visitor/printer.hpp" and
  // uncomment the following line
  // grackle::visit_member(&tmp, grackle::impl::visitor::Printer(ctx, 2));

  MyNestedStruct tmp2{nullptr, {nullptr, nullptr}};
  visit_member(&tmp2, grackle::impl::visitor::AllocateMembers{ctx});
  visit_member_pair(tmp, tmp2, grackle::impl::visitor::CopyMembers{ctx});

  ASSERT_EQ(tmp2.attr0[0], 100);
  ASSERT_EQ(tmp2.attr0[1], -100);
  ASSERT_EQ(tmp2.dummy.attr1[0], 1.0);
  ASSERT_EQ(tmp2.dummy.attr1[1], 2.0);
  ASSERT_EQ(tmp2.dummy.attr2[0], -2);
  ASSERT_EQ(tmp2.dummy.attr2[1], -1);

  visit_member(&tmp, grackle::impl::visitor::FreeMembers{});
  visit_member(&tmp2, grackle::impl::visitor::FreeMembers{});
}
