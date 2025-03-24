#include <memory>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "grackle.h"
#include "grackle_experimental_api.h"

// these headers are included just so we can use grFields_from_field_names_
// (since the public API's initializer hasn't been created yet)
#include "fields/gr_fields_private.h" // grFields_from_field_names_
#include "FrozenKeyIdxBiMap.h"

/// helper function to initialize a map from a vector
FrozenKeyIdxBiMap* FrozenKeyIdxBiMap_from_vec(
  const std::vector<std::string>& vec_
)
{
  std::size_t key_count = vec_.size();

  // create a vector of pointers
  std::vector<const char*> key_ptr_l(key_count, nullptr);
  for (std::size_t i = 0; i < key_count; i++) {
    key_ptr_l[i] = vec_[i].c_str();
  }

  FrozenKeyIdxBiMap* out = nullptr;
  new_FrozenKeyIdxBiMap(&out, key_ptr_l.data(), key_count,
                        BIMAP_COPIES_KEYDATA);
  return out;
}


namespace wrapper {

// this internally deals with the issues of creating FrozenKeyIdxBiMap and
// deallocating if we don't succeed. This is the best way to keep the test
// cases clean without introducing nice C++ wrapper classes
int grFields_from_field_names_(
  gr_fields** out, const std::vector<std::string>& field_names_vec,
  uint64_t flags, uint64_t nelements_per_field
)
{
  FrozenKeyIdxBiMap* field_names = FrozenKeyIdxBiMap_from_vec(field_names_vec);
  int ret = grFields_from_field_names_(out, field_names, flags,
                                       nelements_per_field);
  if (ret != GR_SUCCESS) { drop_FrozenKeyIdxBiMap(field_names); }
  return ret;
}

} // namespace wrapper


// Here we define some tests for initialization
// ============================================
// -> technically, we are testing grFields_from_field_names_, which is an
//    internal (private function).
// -> this is because, the public API function (that wraps this function) is
//    not yet complete yet

TEST(GrFieldInitTest, SimpleCase) {
  gr_fields* fields = nullptr;
  ASSERT_EQ(
    wrapper::grFields_from_field_names_(
      &fields, {"density", "internal_energy"}, 0, 0
    ),
    GR_SUCCESS
  );
  ASSERT_NE(fields, nullptr);

  grunstableFields_free(fields);
}

TEST(GrFieldInitTest, SimpleManagedData) {
  gr_fields* fields = nullptr;
  ASSERT_EQ(
    wrapper::grFields_from_field_names_(
      &fields, {"density", "internal_energy"}, GR_FCREATE_MANAGED_FIELDDATA, 64
    ),
    GR_SUCCESS
  );
  ASSERT_NE(fields, nullptr);

  grunstableFields_free(fields);
}

TEST(GrFieldInitTest, BadLenArgUnamanged) {
  gr_fields* fields = nullptr;
  EXPECT_NE(
    wrapper::grFields_from_field_names_(
      &fields, {"density", "internal_energy"}, 0, 64
    ),
    GR_SUCCESS
  ) << "specifying a non-zero value for nelements_per_field should cause a "
    << "a failure when flags = 0";
  EXPECT_EQ(fields, nullptr);

  if (fields != nullptr) { grunstableFields_free(fields); }
}

TEST(GrFieldInitTest, BadLenArgManaged) {
  gr_fields* fields = nullptr;
  EXPECT_NE(
    wrapper::grFields_from_field_names_(
      &fields, {"density", "internal_energy"}, GR_FCREATE_MANAGED_FIELDDATA, 0
    ), GR_SUCCESS
  ) << "when the GR_FCREATE_MANAGED_FIELDDATA flag is provided, then this "
    << "function should require nelements_per_field > 0";
  EXPECT_EQ(fields, nullptr);

  if (fields != nullptr) { grunstableFields_free(fields); }
}

// when we finish implementing grunstableFields_init_from_local and start
// testing it (rather than grFields_from_field_names_), it will not be able to
// fail in the following manner. The following test will still be important for
// the more distant future when we add support for dynamically specified fields
TEST(GrFieldInitTest, InvalidFieldName) {
  gr_fields* fields = nullptr;
  EXPECT_NE(
    wrapper::grFields_from_field_names_(
      &fields, {"density", "not-a-field"}, 0, 0
    ),
    GR_SUCCESS
  ) << "specifying an unknown field-name should make initialization fail";
  EXPECT_EQ(fields, nullptr);

  if (fields != nullptr) { grunstableFields_free(fields); }
}

// =====================================================

struct gr_fields_deleter_ {
  void operator() (gr_fields* p) { if (p) { grunstableFields_free(p); } }
};
using unique_gr_field_ptr = std::unique_ptr<gr_fields, gr_fields_deleter_>;

// this is used within a fixture
class BasicGrFieldManager {
  unique_gr_field_ptr ptr_ = nullptr;

public:
  static std::vector<std::string> get_field_names()
  { return {"density", "metal_density", "internal_energy"}; }

  gr_fields* get_fields() { return ptr_.get(); }

  BasicGrFieldManager() = default;
  BasicGrFieldManager(uint64_t nelements_per_field)
    : ptr_{nullptr}
  {
    uint64_t flag = 0;
    if (nelements_per_field > 0) { flag = GR_FCREATE_MANAGED_FIELDDATA; }

    gr_fields* ptr = nullptr;
    wrapper::grFields_from_field_names_(
      &ptr, this->get_field_names(), flag, nelements_per_field
    );
    ptr_ = unique_gr_field_ptr(ptr, gr_fields_deleter_{});
  }
};

/// This is a parameterized test fixture
/// - when the parameter is 0, then the field data is unmanaged
/// - when the parameter is anything else, that specifies the number of
///   entries per field
class GrFieldBasicTest : public testing::TestWithParam<uint64_t> {
protected:

  BasicGrFieldManager manager_;
  gr_fields* fields;

  void SetUp() override {
    manager_ = BasicGrFieldManager(this->GetParam());
    fields = manager_.get_fields();
  }

  // no need to define TearDown(), thanks to BasicGrFieldManager's destructor

  std::vector<std::string> get_field_name_vec()
  { return manager_.get_field_names(); }
};

TEST_P(GrFieldBasicTest, NumFields)
{
  int64_t num_fields=-1;
  ASSERT_EQ(grunstableFields_num_fields(fields, &num_fields), GR_SUCCESS);
  ASSERT_EQ(num_fields, (int64_t)get_field_name_vec().size());
}

TEST_P(GrFieldBasicTest, NameByIdx)
{
  std::vector<std::string> field_name_vec = get_field_name_vec();
  int64_t num_fields = (int64_t)field_name_vec.size();
  for (int64_t i = 0; i < num_fields; i++) {
    const char* name = nullptr;
    ASSERT_EQ(grunstableFields_name_by_idx(fields, i, &name), GR_SUCCESS);
    EXPECT_EQ(std::string(name), field_name_vec[i]);
  }
}

TEST_P(GrFieldBasicTest, NameByIdxNegative)
{
  const char* name = nullptr;
  ASSERT_NE(grunstableFields_name_by_idx(fields, -1, &name), GR_SUCCESS);
}

TEST_P(GrFieldBasicTest, NameByIdxTooLarge)
{
  std::vector<std::string> field_name_vec = get_field_name_vec();
  int64_t bad_idx = (int64_t)field_name_vec.size();
  
  const char* name = nullptr;
  ASSERT_NE(grunstableFields_name_by_idx(fields, bad_idx, &name), GR_SUCCESS);
}

TEST_P(GrFieldBasicTest, IdxByName)
{
  std::vector<std::string> field_name_vec = get_field_name_vec();
  int64_t num_fields = (int64_t)field_name_vec.size();
  for (int64_t i = 0; i < num_fields; i++) {
    const char* name = field_name_vec[i].c_str();
    int64_t idx = i + 1;
    EXPECT_EQ(grunstableFields_idx_by_name(fields, name, &idx), GR_SUCCESS);
    EXPECT_EQ(idx, i);
  }
}

TEST_P(GrFieldBasicTest, IdxByNameNotAField)
{
  int64_t idx = 0;
  EXPECT_EQ(
    grunstableFields_idx_by_name(fields, "not-a-field", &idx), GR_SUCCESS
  );
  EXPECT_EQ(idx, -1);
}

TEST_P(GrFieldBasicTest, SetGridWidth)
{
  EXPECT_EQ(grunstableFields_set_grid_dx(fields, 1.0), GR_SUCCESS);
}

TEST_P(GrFieldBasicTest, SetGridLayout)
{
  int64_t extent[3] = {1, 1, 1};
  int64_t start[3]  = {0, 0, 0};
  int64_t stop[3]   = {1, 1, 1};

  for (int grid_rank = 1; grid_rank <= 3; grid_rank++) {
    EXPECT_EQ(
      grunstableFields_set_grid_layout(
        fields, grid_rank, extent, start, stop
      ), GR_SUCCESS
    );
  }
}

TEST_P(GrFieldBasicTest, SetGridLayoutWithNULL)
{
  int64_t extent[3] = {1, 1, 1};

  for (int grid_rank = 1; grid_rank <= 3; grid_rank++) {
    EXPECT_EQ(
      grunstableFields_set_grid_layout(
        fields, grid_rank, extent, nullptr, nullptr
      ), GR_SUCCESS
    );
  }
}

TEST_P(GrFieldBasicTest, SetGridLayoutWithInvalidNULL)
{
  int64_t extent[3] = {1, 1, 1};
  int64_t start[3]  = {0, 0, 0};
  int64_t stop[3]   = {1, 1, 1};

  for (int grid_rank = 1; grid_rank <= 3; grid_rank++) {
    EXPECT_NE(
      grunstableFields_set_grid_layout(
        fields, grid_rank, extent, start, nullptr
      ), GR_SUCCESS
    );

    EXPECT_NE(
      grunstableFields_set_grid_layout(
        fields, grid_rank, extent, nullptr, stop
      ), GR_SUCCESS
    );

    EXPECT_NE(
      grunstableFields_set_grid_layout(
        fields, grid_rank, nullptr, start, stop
      ), GR_SUCCESS
    );

    EXPECT_NE(
      grunstableFields_set_grid_layout(
        fields, grid_rank, nullptr, nullptr, nullptr
      ), GR_SUCCESS
    );

  }
}


TEST_P(GrFieldBasicTest, SetGridLayoutBadRank)
{
  int64_t extent[4] = {1, 1, 1, 1};

  EXPECT_NE(
    grunstableFields_set_grid_layout(fields, 0, extent, nullptr, nullptr),
    GR_SUCCESS
  );

  EXPECT_NE(
    grunstableFields_set_grid_layout(fields, -1, extent, nullptr, nullptr),
    GR_SUCCESS
  );

  EXPECT_NE(
    grunstableFields_set_grid_layout(fields, 4, extent, nullptr, nullptr),
    GR_SUCCESS
  );
}

TEST_P(GrFieldBasicTest, GetSlotByName)
{
  std::vector<std::string> field_name_vec = get_field_name_vec();
  int64_t num_fields = (int64_t)field_name_vec.size();
  std::vector<gr_float**> pointer_vec{};
  for (int64_t i = 0; i < num_fields; i++) {
    const char* name = field_name_vec[i].c_str();
    gr_float** curslot = grunstableFields_get_slot_by_name(fields, name);
    ASSERT_NE(curslot, nullptr);

    // confirm that we haven't seen this before
    for (int64_t j = 0; j < i; j++) {
      ASSERT_NE(curslot, pointer_vec[j])
        << "the slots for the `" << name << "` and `" << field_name_vec[j]
        << "should specify different pointers";
    }
    pointer_vec.push_back(curslot);

    if (GetParam() == 0) {
      ASSERT_EQ(*curslot, nullptr)
        << "pointers initially held by slots (for unmanaged field data) "
        << "should all be nullptr";
    } else {
      ASSERT_NE(*curslot, nullptr)
        << "pointers held by slots (for managed field data) should not be "
        << "nullptr";
    }
  }
}

TEST_P(GrFieldBasicTest, GetSlotByNameNonField)
{
  ASSERT_EQ(
    grunstableFields_get_slot_by_name(fields, "not-a-field"), nullptr
  );
}

INSTANTIATE_TEST_SUITE_P(
  , /* <- leaving Instantiation name empty */
  GrFieldBasicTest,
  testing::Values(0, 64),
  [](const testing::TestParamInfo<GrFieldBasicTest::ParamType>& info) {
    if (info.param == 0) {
      return std::string("UnmanagedFieldData");
    } else {
      return std::to_string(info.param) + "ElemPerString";
    }
  }
);

// ====================================================

TEST(GrFieldAssorted, StorePtr)
{
  // initialize an instance of gr_fields, where it doesn't manage field data
  BasicGrFieldManager manager(0);
  std::vector<std::string> field_name_vec = manager.get_field_names();
  gr_fields* fields = manager.get_fields();

  // set up our external field data
  int nelements_per_field = 3;
  int64_t num_fields = (int64_t)field_name_vec.size();
  std::vector<gr_float> external_field_data(nelements_per_field * num_fields);

  // this is a lambda function to check that the slot has the expected value
  auto check_slot_val = [&](const char* name, gr_float* expected_ptr) -> bool
  {
    gr_float** slot = grunstableFields_get_slot_by_name(fields, name);
    return (*slot) == expected_ptr;
  };

  // now we will go through and store pointers to our fields
  for (int64_t i = 0; i < num_fields; i++) {
    const char* name = field_name_vec[i].c_str();

    // sanity check:
    ASSERT_TRUE(check_slot_val(name, nullptr));

    gr_float* ptr = external_field_data.data() + i * nelements_per_field;

    ASSERT_EQ(
      grunstableFields_store_by_name(fields, name, ptr), GR_SUCCESS
    );

    ASSERT_TRUE(check_slot_val(name, ptr)) << "the update didn't work";

  }
}

TEST(GrFieldAssorted, StorePtrNonField)
{
  // initialize an instance of gr_fields, where it doesn't manage field data
  BasicGrFieldManager manager(0);
  gr_fields* fields = manager.get_fields();

  std::vector<gr_float> data(3);
  ASSERT_NE(
    grunstableFields_store_by_name(fields, "not-a-field", data.data()),
    GR_SUCCESS
  ) << "storing a pointer for an unknown field should be forbidden";
}

TEST(GrFieldAssorted, StorePtrManaged)
{
  // initialize an instance of gr_fields, where it manages field data
  BasicGrFieldManager manager(64);
  std::vector<std::string> field_name_vec = manager.get_field_names();
  gr_fields* fields = manager.get_fields();

  std::vector<gr_float> data(3);
  const char* field_name = field_name_vec[0].c_str();
  ASSERT_NE(
    grunstableFields_store_by_name(fields, field_name, data.data()),
    GR_SUCCESS
  ) << "an instance of gr_fields, with managed field_data, should not be "
    << "allowed to store pointers in any slots";
}

TEST(GrFieldAssorted, GetSlotByNameManaged)
{
  // initialize an instance of gr_fields, where it manages field data
  BasicGrFieldManager manager(64);
  std::vector<std::string> field_name_vec = manager.get_field_names();
  gr_fields* fields = manager.get_fields();

  // now we will go through and store pointers to our fields
  for (int64_t i = 0; i < field_name_vec.size(); i++) {
    const char* name = field_name_vec[i].c_str();
    gr_float** slot = grunstableFields_get_slot_by_name(fields, name);
    ASSERT_NE(slot, nullptr);
    ASSERT_NE(*slot, nullptr);
  }
}

// this is a really long drawn out test...
TEST(GrFieldAssorted, FieldListLocalPrimordialChem2)
{
  // Setup 1: set up the units system.
  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24;
  my_units.length_units = 1.0;
  my_units.time_units = 1.0e12;
  my_units.a_units = 1.0;
  my_units.a_value = 1.0;

  // Setup 2: create a chemistry object to hold configuration parameters.
  chemistry_data *my_chemistry = new chemistry_data;
  {
    int tmp = set_default_chemistry_parameters(my_chemistry);
    if (tmp != GR_SUCCESS) { delete my_chemistry; }
    ASSERT_EQ(tmp, GR_SUCCESS);
  }

  // Setup 3: set the parameters
  my_chemistry->use_grackle = 1; // chemistry on
  my_chemistry->with_radiative_cooling = 1; // cooling on
  my_chemistry->primordial_chemistry = 2; // molecular network with H & He
  // for simplicity, we disable all features using the datafile. After we merge
  // in changes from #254, we can start using datafiles in unit tests.
  my_chemistry->metal_cooling = 0;          // metal cooling off
  my_chemistry->UVbackground = 0;           // UV background off
  my_chemistry->grackle_data_file = nullptr;

  // Setup 4: Create chemistry data storage object to store rates.
  chemistry_data_storage* my_rates = new chemistry_data_storage;

  // Setup 5: initialize the chemistry object.
  {
    int tmp = local_initialize_chemistry_data(my_chemistry, my_rates,
                                              &my_units);
    if (tmp != GR_SUCCESS) {
      delete my_rates;
      delete my_chemistry;
    }
    ASSERT_EQ(tmp, GR_SUCCESS);
  }

  // Invoke the function we are actually testing
  // ===========================================
  gr_fields* fields = nullptr;
  int ret = grunstableFields_init_from_local(&fields, my_chemistry, my_rates,
                                             0, 0);
  if (ret != GR_SUCCESS) {
    delete my_rates;
    delete my_chemistry;
  }
  ASSERT_EQ(ret, GR_SUCCESS);
  ASSERT_NE(fields, nullptr);

  // get the number of fields
  // =======================================================
  int64_t num_fields = -1;
  EXPECT_EQ(grunstableFields_num_fields(fields, &num_fields), GR_SUCCESS);
  if (num_fields < 1) {
    grunstableFields_free(fields);
    local_free_chemistry_data(my_chemistry, my_rates);
    delete my_rates;
    delete my_chemistry;
  }
  ASSERT_GE(num_fields, 1)
    << "it shouldn't be possible to have a non-positive field count";

  // Now, confirm that some fields are in the field list
  // ===================================================
  std::vector<std::string> expected_field_names = {
    "density", "internal_energy", "HI_density", "HII_density", "HeI_density",
    "HeII_density", "HeIII_density", "e_density", "H2I_density",
    "H2II_density", "HM_density"
  };
  for (const std::string& name : expected_field_names) { // range-based loop
    int64_t idx = 0;
    EXPECT_EQ(
      grunstableFields_idx_by_name(fields, name.c_str(), &idx), GR_SUCCESS
    ) << "encountered a problem while querying the idx associated with the `"
      << name << "` field";
    EXPECT_GT(idx, -1) << '`' << name << "` field seems to be missing";
    EXPECT_LT(idx, num_fields)
      << "index associated with `" << name << "` field doesn't make sense";
  }

  // Now, check that some fields are missing from the field list
  // ===========================================================
  // (this is not exhaustive)
  std::vector<std::string> unexpected_field_names = {
    "DII_density", "not-a-field", "metal_density", "dust_density",
    "volumetric_heating_rate", "temperature_floor", "RT_HeI_ionization_rate",
    "H2_self_shielding_length", "isrf_habing"
  };
  for (const std::string& name : unexpected_field_names) { // range-based loop
    int64_t idx = 0;
    EXPECT_EQ(
      grunstableFields_idx_by_name(fields, name.c_str(), &idx), GR_SUCCESS
    ) << "encountered a problem while querying the idx associated with the `"
      << name << "` field";
    EXPECT_EQ(idx, -1) << '`' << name << "` field should be missing";
  }

  // CLEANUP
  // =======
  grunstableFields_free(fields);
  local_free_chemistry_data(my_chemistry, my_rates);
  delete my_rates;
  delete my_chemistry;
}
