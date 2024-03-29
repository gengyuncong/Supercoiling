// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/types/TrajectoryBarrier.proto

#include "lm/types/TrajectoryBarrier.pb.h"

#include <algorithm>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/wire_format_lite.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
extern PROTOBUF_INTERNAL_EXPORT_lm_2ftypes_2fTrajectoryBarrier_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_TrajectoryBarrier_lm_2ftypes_2fTrajectoryBarrier_2eproto;
namespace lm {
namespace types {
class TrajectoryBarrierDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<TrajectoryBarrier> _instance;
} _TrajectoryBarrier_default_instance_;
class TrajectoryBarriersDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<TrajectoryBarriers> _instance;
} _TrajectoryBarriers_default_instance_;
}  // namespace types
}  // namespace lm
static void InitDefaultsscc_info_TrajectoryBarrier_lm_2ftypes_2fTrajectoryBarrier_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::types::_TrajectoryBarrier_default_instance_;
    new (ptr) ::lm::types::TrajectoryBarrier();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::types::TrajectoryBarrier::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_TrajectoryBarrier_lm_2ftypes_2fTrajectoryBarrier_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_TrajectoryBarrier_lm_2ftypes_2fTrajectoryBarrier_2eproto}, {}};

static void InitDefaultsscc_info_TrajectoryBarriers_lm_2ftypes_2fTrajectoryBarrier_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::types::_TrajectoryBarriers_default_instance_;
    new (ptr) ::lm::types::TrajectoryBarriers();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::types::TrajectoryBarriers::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_TrajectoryBarriers_lm_2ftypes_2fTrajectoryBarrier_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_TrajectoryBarriers_lm_2ftypes_2fTrajectoryBarrier_2eproto}, {
      &scc_info_TrajectoryBarrier_lm_2ftypes_2fTrajectoryBarrier_2eproto.base,}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2ftypes_2fTrajectoryBarrier_2eproto[2];
static const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* file_level_enum_descriptors_lm_2ftypes_2fTrajectoryBarrier_2eproto[2];
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2ftypes_2fTrajectoryBarrier_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2ftypes_2fTrajectoryBarrier_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::types::TrajectoryBarrier, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::types::TrajectoryBarrier, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::types::TrajectoryBarrier, type_),
  PROTOBUF_FIELD_OFFSET(::lm::types::TrajectoryBarrier, type_arg_),
  PROTOBUF_FIELD_OFFSET(::lm::types::TrajectoryBarrier, behavior_),
  PROTOBUF_FIELD_OFFSET(::lm::types::TrajectoryBarrier, behavior_value_double_),
  PROTOBUF_FIELD_OFFSET(::lm::types::TrajectoryBarrier, behavior_value_int_),
  0,
  1,
  2,
  3,
  4,
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::lm::types::TrajectoryBarriers, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::types::TrajectoryBarriers, barrier_),
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 10, sizeof(::lm::types::TrajectoryBarrier)},
  { 15, -1, sizeof(::lm::types::TrajectoryBarriers)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::types::_TrajectoryBarrier_default_instance_),
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::types::_TrajectoryBarriers_default_instance_),
};

const char descriptor_table_protodef_lm_2ftypes_2fTrajectoryBarrier_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n lm/types/TrajectoryBarrier.proto\022\010lm.t"
  "ypes\"\201\004\n\021TrajectoryBarrier\0225\n\004type\030\001 \002(\016"
  "2\'.lm.types.TrajectoryBarrier.BarrierTyp"
  "e\022\020\n\010type_arg\030\002 \001(\r\022=\n\010behavior\030\013 \002(\0162+."
  "lm.types.TrajectoryBarrier.BarrierBehavi"
  "or\022 \n\025behavior_value_double\030\014 \001(\001:\0010\022\035\n\022"
  "behavior_value_int\030\r \001(\003:\0010\"/\n\013BarrierTy"
  "pe\022\013\n\007SPECIES\020\000\022\023\n\017ORDER_PARAMETER\020\001\"\361\001\n"
  "\017BarrierBehavior\022\016\n\nREFLECTING\020\000\022\031\n\025REFL"
  "ECTING_DECREASING\020\001\022\031\n\025REFLECTING_INCREA"
  "SING\020\002\022\014\n\010TRACKING\020\003\022!\n\035TRACKING_DECREAS"
  "ING_INCLUSIVE\020\004\022!\n\035TRACKING_INCREASING_I"
  "NCLUSIVE\020\005\022!\n\035TRACKING_DECREASING_EXCLUS"
  "IVE\020\006\022!\n\035TRACKING_INCREASING_EXCLUSIVE\020\007"
  "\"B\n\022TrajectoryBarriers\022,\n\007barrier\030\001 \003(\0132"
  "\033.lm.types.TrajectoryBarrier"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto_sccs[2] = {
  &scc_info_TrajectoryBarrier_lm_2ftypes_2fTrajectoryBarrier_2eproto.base,
  &scc_info_TrajectoryBarriers_lm_2ftypes_2fTrajectoryBarrier_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto = {
  false, false, descriptor_table_protodef_lm_2ftypes_2fTrajectoryBarrier_2eproto, "lm/types/TrajectoryBarrier.proto", 628,
  &descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto_once, descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto_sccs, descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto_deps, 2, 0,
  schemas, file_default_instances, TableStruct_lm_2ftypes_2fTrajectoryBarrier_2eproto::offsets,
  file_level_metadata_lm_2ftypes_2fTrajectoryBarrier_2eproto, 2, file_level_enum_descriptors_lm_2ftypes_2fTrajectoryBarrier_2eproto, file_level_service_descriptors_lm_2ftypes_2fTrajectoryBarrier_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2ftypes_2fTrajectoryBarrier_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto)), true);
namespace lm {
namespace types {
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* TrajectoryBarrier_BarrierType_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto);
  return file_level_enum_descriptors_lm_2ftypes_2fTrajectoryBarrier_2eproto[0];
}
bool TrajectoryBarrier_BarrierType_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)
constexpr TrajectoryBarrier_BarrierType TrajectoryBarrier::SPECIES;
constexpr TrajectoryBarrier_BarrierType TrajectoryBarrier::ORDER_PARAMETER;
constexpr TrajectoryBarrier_BarrierType TrajectoryBarrier::BarrierType_MIN;
constexpr TrajectoryBarrier_BarrierType TrajectoryBarrier::BarrierType_MAX;
constexpr int TrajectoryBarrier::BarrierType_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* TrajectoryBarrier_BarrierBehavior_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_lm_2ftypes_2fTrajectoryBarrier_2eproto);
  return file_level_enum_descriptors_lm_2ftypes_2fTrajectoryBarrier_2eproto[1];
}
bool TrajectoryBarrier_BarrierBehavior_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::REFLECTING;
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::REFLECTING_DECREASING;
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::REFLECTING_INCREASING;
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::TRACKING;
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::TRACKING_DECREASING_INCLUSIVE;
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::TRACKING_INCREASING_INCLUSIVE;
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::TRACKING_DECREASING_EXCLUSIVE;
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::TRACKING_INCREASING_EXCLUSIVE;
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::BarrierBehavior_MIN;
constexpr TrajectoryBarrier_BarrierBehavior TrajectoryBarrier::BarrierBehavior_MAX;
constexpr int TrajectoryBarrier::BarrierBehavior_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)

// ===================================================================

void TrajectoryBarrier::InitAsDefaultInstance() {
}
class TrajectoryBarrier::_Internal {
 public:
  using HasBits = decltype(std::declval<TrajectoryBarrier>()._has_bits_);
  static void set_has_type(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_type_arg(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static void set_has_behavior(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_behavior_value_double(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
  static void set_has_behavior_value_int(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000005) ^ 0x00000005) != 0;
  }
};

TrajectoryBarrier::TrajectoryBarrier(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.types.TrajectoryBarrier)
}
TrajectoryBarrier::TrajectoryBarrier(const TrajectoryBarrier& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&type_, &from.type_,
    static_cast<size_t>(reinterpret_cast<char*>(&behavior_value_int_) -
    reinterpret_cast<char*>(&type_)) + sizeof(behavior_value_int_));
  // @@protoc_insertion_point(copy_constructor:lm.types.TrajectoryBarrier)
}

void TrajectoryBarrier::SharedCtor() {
  ::memset(&type_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&behavior_value_int_) -
      reinterpret_cast<char*>(&type_)) + sizeof(behavior_value_int_));
}

TrajectoryBarrier::~TrajectoryBarrier() {
  // @@protoc_insertion_point(destructor:lm.types.TrajectoryBarrier)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void TrajectoryBarrier::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void TrajectoryBarrier::ArenaDtor(void* object) {
  TrajectoryBarrier* _this = reinterpret_cast< TrajectoryBarrier* >(object);
  (void)_this;
}
void TrajectoryBarrier::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void TrajectoryBarrier::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const TrajectoryBarrier& TrajectoryBarrier::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_TrajectoryBarrier_lm_2ftypes_2fTrajectoryBarrier_2eproto.base);
  return *internal_default_instance();
}


void TrajectoryBarrier::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.types.TrajectoryBarrier)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    ::memset(&type_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&behavior_value_int_) -
        reinterpret_cast<char*>(&type_)) + sizeof(behavior_value_int_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* TrajectoryBarrier::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // required .lm.types.TrajectoryBarrier.BarrierType type = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          ::PROTOBUF_NAMESPACE_ID::uint64 val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::lm::types::TrajectoryBarrier_BarrierType_IsValid(val))) {
            _internal_set_type(static_cast<::lm::types::TrajectoryBarrier_BarrierType>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(1, val, mutable_unknown_fields());
          }
        } else goto handle_unusual;
        continue;
      // optional uint32 type_arg = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 16)) {
          _Internal::set_has_type_arg(&has_bits);
          type_arg_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // required .lm.types.TrajectoryBarrier.BarrierBehavior behavior = 11;
      case 11:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 88)) {
          ::PROTOBUF_NAMESPACE_ID::uint64 val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::lm::types::TrajectoryBarrier_BarrierBehavior_IsValid(val))) {
            _internal_set_behavior(static_cast<::lm::types::TrajectoryBarrier_BarrierBehavior>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(11, val, mutable_unknown_fields());
          }
        } else goto handle_unusual;
        continue;
      // optional double behavior_value_double = 12 [default = 0];
      case 12:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 97)) {
          _Internal::set_has_behavior_value_double(&has_bits);
          behavior_value_double_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else goto handle_unusual;
        continue;
      // optional int64 behavior_value_int = 13 [default = 0];
      case 13:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 104)) {
          _Internal::set_has_behavior_value_int(&has_bits);
          behavior_value_int_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      default: {
      handle_unusual:
        if ((tag & 7) == 4 || tag == 0) {
          ctx->SetLastTag(tag);
          goto success;
        }
        ptr = UnknownFieldParse(tag,
            _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
            ptr, ctx);
        CHK_(ptr != nullptr);
        continue;
      }
    }  // switch
  }  // while
success:
  _has_bits_.Or(has_bits);
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* TrajectoryBarrier::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.types.TrajectoryBarrier)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required .lm.types.TrajectoryBarrier.BarrierType type = 1;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteEnumToArray(
      1, this->_internal_type(), target);
  }

  // optional uint32 type_arg = 2;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt32ToArray(2, this->_internal_type_arg(), target);
  }

  // required .lm.types.TrajectoryBarrier.BarrierBehavior behavior = 11;
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteEnumToArray(
      11, this->_internal_behavior(), target);
  }

  // optional double behavior_value_double = 12 [default = 0];
  if (cached_has_bits & 0x00000008u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(12, this->_internal_behavior_value_double(), target);
  }

  // optional int64 behavior_value_int = 13 [default = 0];
  if (cached_has_bits & 0x00000010u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt64ToArray(13, this->_internal_behavior_value_int(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.types.TrajectoryBarrier)
  return target;
}

size_t TrajectoryBarrier::RequiredFieldsByteSizeFallback() const {
// @@protoc_insertion_point(required_fields_byte_size_fallback_start:lm.types.TrajectoryBarrier)
  size_t total_size = 0;

  if (_internal_has_type()) {
    // required .lm.types.TrajectoryBarrier.BarrierType type = 1;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::EnumSize(this->_internal_type());
  }

  if (_internal_has_behavior()) {
    // required .lm.types.TrajectoryBarrier.BarrierBehavior behavior = 11;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::EnumSize(this->_internal_behavior());
  }

  return total_size;
}
size_t TrajectoryBarrier::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.types.TrajectoryBarrier)
  size_t total_size = 0;

  if (((_has_bits_[0] & 0x00000005) ^ 0x00000005) == 0) {  // All required fields are present.
    // required .lm.types.TrajectoryBarrier.BarrierType type = 1;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::EnumSize(this->_internal_type());

    // required .lm.types.TrajectoryBarrier.BarrierBehavior behavior = 11;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::EnumSize(this->_internal_behavior());

  } else {
    total_size += RequiredFieldsByteSizeFallback();
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // optional uint32 type_arg = 2;
  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000002u) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt32Size(
        this->_internal_type_arg());
  }

  if (cached_has_bits & 0x00000018u) {
    // optional double behavior_value_double = 12 [default = 0];
    if (cached_has_bits & 0x00000008u) {
      total_size += 1 + 8;
    }

    // optional int64 behavior_value_int = 13 [default = 0];
    if (cached_has_bits & 0x00000010u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int64Size(
          this->_internal_behavior_value_int());
    }

  }
  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void TrajectoryBarrier::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.types.TrajectoryBarrier)
  GOOGLE_DCHECK_NE(&from, this);
  const TrajectoryBarrier* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<TrajectoryBarrier>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.types.TrajectoryBarrier)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.types.TrajectoryBarrier)
    MergeFrom(*source);
  }
}

void TrajectoryBarrier::MergeFrom(const TrajectoryBarrier& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.types.TrajectoryBarrier)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    if (cached_has_bits & 0x00000001u) {
      type_ = from.type_;
    }
    if (cached_has_bits & 0x00000002u) {
      type_arg_ = from.type_arg_;
    }
    if (cached_has_bits & 0x00000004u) {
      behavior_ = from.behavior_;
    }
    if (cached_has_bits & 0x00000008u) {
      behavior_value_double_ = from.behavior_value_double_;
    }
    if (cached_has_bits & 0x00000010u) {
      behavior_value_int_ = from.behavior_value_int_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void TrajectoryBarrier::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.types.TrajectoryBarrier)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void TrajectoryBarrier::CopyFrom(const TrajectoryBarrier& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.types.TrajectoryBarrier)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool TrajectoryBarrier::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  return true;
}

void TrajectoryBarrier::InternalSwap(TrajectoryBarrier* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(TrajectoryBarrier, behavior_value_int_)
      + sizeof(TrajectoryBarrier::behavior_value_int_)
      - PROTOBUF_FIELD_OFFSET(TrajectoryBarrier, type_)>(
          reinterpret_cast<char*>(&type_),
          reinterpret_cast<char*>(&other->type_));
}

::PROTOBUF_NAMESPACE_ID::Metadata TrajectoryBarrier::GetMetadata() const {
  return GetMetadataStatic();
}


// ===================================================================

void TrajectoryBarriers::InitAsDefaultInstance() {
}
class TrajectoryBarriers::_Internal {
 public:
};

TrajectoryBarriers::TrajectoryBarriers(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  barrier_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.types.TrajectoryBarriers)
}
TrajectoryBarriers::TrajectoryBarriers(const TrajectoryBarriers& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      barrier_(from.barrier_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:lm.types.TrajectoryBarriers)
}

void TrajectoryBarriers::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_TrajectoryBarriers_lm_2ftypes_2fTrajectoryBarrier_2eproto.base);
}

TrajectoryBarriers::~TrajectoryBarriers() {
  // @@protoc_insertion_point(destructor:lm.types.TrajectoryBarriers)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void TrajectoryBarriers::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void TrajectoryBarriers::ArenaDtor(void* object) {
  TrajectoryBarriers* _this = reinterpret_cast< TrajectoryBarriers* >(object);
  (void)_this;
}
void TrajectoryBarriers::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void TrajectoryBarriers::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const TrajectoryBarriers& TrajectoryBarriers::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_TrajectoryBarriers_lm_2ftypes_2fTrajectoryBarrier_2eproto.base);
  return *internal_default_instance();
}


void TrajectoryBarriers::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.types.TrajectoryBarriers)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  barrier_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* TrajectoryBarriers::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated .lm.types.TrajectoryBarrier barrier = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_barrier(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<10>(ptr));
        } else goto handle_unusual;
        continue;
      default: {
      handle_unusual:
        if ((tag & 7) == 4 || tag == 0) {
          ctx->SetLastTag(tag);
          goto success;
        }
        ptr = UnknownFieldParse(tag,
            _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
            ptr, ctx);
        CHK_(ptr != nullptr);
        continue;
      }
    }  // switch
  }  // while
success:
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* TrajectoryBarriers::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.types.TrajectoryBarriers)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated .lm.types.TrajectoryBarrier barrier = 1;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_barrier_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(1, this->_internal_barrier(i), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.types.TrajectoryBarriers)
  return target;
}

size_t TrajectoryBarriers::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.types.TrajectoryBarriers)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .lm.types.TrajectoryBarrier barrier = 1;
  total_size += 1UL * this->_internal_barrier_size();
  for (const auto& msg : this->barrier_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void TrajectoryBarriers::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.types.TrajectoryBarriers)
  GOOGLE_DCHECK_NE(&from, this);
  const TrajectoryBarriers* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<TrajectoryBarriers>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.types.TrajectoryBarriers)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.types.TrajectoryBarriers)
    MergeFrom(*source);
  }
}

void TrajectoryBarriers::MergeFrom(const TrajectoryBarriers& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.types.TrajectoryBarriers)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  barrier_.MergeFrom(from.barrier_);
}

void TrajectoryBarriers::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.types.TrajectoryBarriers)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void TrajectoryBarriers::CopyFrom(const TrajectoryBarriers& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.types.TrajectoryBarriers)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool TrajectoryBarriers::IsInitialized() const {
  if (!::PROTOBUF_NAMESPACE_ID::internal::AllAreInitialized(barrier_)) return false;
  return true;
}

void TrajectoryBarriers::InternalSwap(TrajectoryBarriers* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  barrier_.InternalSwap(&other->barrier_);
}

::PROTOBUF_NAMESPACE_ID::Metadata TrajectoryBarriers::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace types
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::types::TrajectoryBarrier* Arena::CreateMaybeMessage< ::lm::types::TrajectoryBarrier >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::types::TrajectoryBarrier >(arena);
}
template<> PROTOBUF_NOINLINE ::lm::types::TrajectoryBarriers* Arena::CreateMaybeMessage< ::lm::types::TrajectoryBarriers >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::types::TrajectoryBarriers >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
