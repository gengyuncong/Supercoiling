// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/input/ffpilot/FFPilotPhaseLimit.proto

#include "lm/input/ffpilot/FFPilotPhaseLimit.pb.h"

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
namespace lm {
namespace input {
namespace ffpilot {
class FFPilotPhaseLimitDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<FFPilotPhaseLimit> _instance;
} _FFPilotPhaseLimit_default_instance_;
}  // namespace ffpilot
}  // namespace input
}  // namespace lm
static void InitDefaultsscc_info_FFPilotPhaseLimit_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::input::ffpilot::_FFPilotPhaseLimit_default_instance_;
    new (ptr) ::lm::input::ffpilot::FFPilotPhaseLimit();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::input::ffpilot::FFPilotPhaseLimit::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_FFPilotPhaseLimit_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_FFPilotPhaseLimit_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto[1];
static const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* file_level_enum_descriptors_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::input::ffpilot::FFPilotPhaseLimit, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::input::ffpilot::FFPilotPhaseLimit, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::input::ffpilot::FFPilotPhaseLimit, stop_condition_),
  PROTOBUF_FIELD_OFFSET(::lm::input::ffpilot::FFPilotPhaseLimit, events_per_trajectory_),
  PROTOBUF_FIELD_OFFSET(::lm::input::ffpilot::FFPilotPhaseLimit, trajectories_per_phase_),
  PROTOBUF_FIELD_OFFSET(::lm::input::ffpilot::FFPilotPhaseLimit, dvalue_),
  PROTOBUF_FIELD_OFFSET(::lm::input::ffpilot::FFPilotPhaseLimit, uvalue_),
  0,
  1,
  2,
  3,
  4,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 10, sizeof(::lm::input::ffpilot::FFPilotPhaseLimit)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::input::ffpilot::_FFPilotPhaseLimit_default_instance_),
};

const char descriptor_table_protodef_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n(lm/input/ffpilot/FFPilotPhaseLimit.pro"
  "to\022\020lm.input.ffpilot\"\222\002\n\021FFPilotPhaseLim"
  "it\022Y\n\016stop_condition\030\013 \001(\01621.lm.input.ff"
  "pilot.FFPilotPhaseLimit.StopCondition:\016F"
  "ORWARD_FLUXES\022\035\n\025events_per_trajectory\030\025"
  " \001(\004\022\036\n\026trajectories_per_phase\030\026 \001(\004\022\016\n\006"
  "dvalue\030\037 \001(\001\022\016\n\006uvalue\030  \001(\004\"C\n\rStopCond"
  "ition\022\022\n\016FORWARD_FLUXES\020\000\022\024\n\020TRAJECTORY_"
  "COUNT\020\001\022\010\n\004TIME\020\002"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto_sccs[1] = {
  &scc_info_FFPilotPhaseLimit_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto = {
  false, false, descriptor_table_protodef_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto, "lm/input/ffpilot/FFPilotPhaseLimit.proto", 337,
  &descriptor_table_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto_once, descriptor_table_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto_sccs, descriptor_table_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto_deps, 1, 0,
  schemas, file_default_instances, TableStruct_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto::offsets,
  file_level_metadata_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto, 1, file_level_enum_descriptors_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto, file_level_service_descriptors_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto)), true);
namespace lm {
namespace input {
namespace ffpilot {
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* FFPilotPhaseLimit_StopCondition_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto);
  return file_level_enum_descriptors_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto[0];
}
bool FFPilotPhaseLimit_StopCondition_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
    case 2:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)
constexpr FFPilotPhaseLimit_StopCondition FFPilotPhaseLimit::FORWARD_FLUXES;
constexpr FFPilotPhaseLimit_StopCondition FFPilotPhaseLimit::TRAJECTORY_COUNT;
constexpr FFPilotPhaseLimit_StopCondition FFPilotPhaseLimit::TIME;
constexpr FFPilotPhaseLimit_StopCondition FFPilotPhaseLimit::StopCondition_MIN;
constexpr FFPilotPhaseLimit_StopCondition FFPilotPhaseLimit::StopCondition_MAX;
constexpr int FFPilotPhaseLimit::StopCondition_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)

// ===================================================================

void FFPilotPhaseLimit::InitAsDefaultInstance() {
}
class FFPilotPhaseLimit::_Internal {
 public:
  using HasBits = decltype(std::declval<FFPilotPhaseLimit>()._has_bits_);
  static void set_has_stop_condition(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_events_per_trajectory(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static void set_has_trajectories_per_phase(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_dvalue(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
  static void set_has_uvalue(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
};

FFPilotPhaseLimit::FFPilotPhaseLimit(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.input.ffpilot.FFPilotPhaseLimit)
}
FFPilotPhaseLimit::FFPilotPhaseLimit(const FFPilotPhaseLimit& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&stop_condition_, &from.stop_condition_,
    static_cast<size_t>(reinterpret_cast<char*>(&uvalue_) -
    reinterpret_cast<char*>(&stop_condition_)) + sizeof(uvalue_));
  // @@protoc_insertion_point(copy_constructor:lm.input.ffpilot.FFPilotPhaseLimit)
}

void FFPilotPhaseLimit::SharedCtor() {
  ::memset(&stop_condition_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&uvalue_) -
      reinterpret_cast<char*>(&stop_condition_)) + sizeof(uvalue_));
}

FFPilotPhaseLimit::~FFPilotPhaseLimit() {
  // @@protoc_insertion_point(destructor:lm.input.ffpilot.FFPilotPhaseLimit)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void FFPilotPhaseLimit::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void FFPilotPhaseLimit::ArenaDtor(void* object) {
  FFPilotPhaseLimit* _this = reinterpret_cast< FFPilotPhaseLimit* >(object);
  (void)_this;
}
void FFPilotPhaseLimit::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void FFPilotPhaseLimit::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const FFPilotPhaseLimit& FFPilotPhaseLimit::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_FFPilotPhaseLimit_lm_2finput_2fffpilot_2fFFPilotPhaseLimit_2eproto.base);
  return *internal_default_instance();
}


void FFPilotPhaseLimit::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.input.ffpilot.FFPilotPhaseLimit)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    ::memset(&stop_condition_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&uvalue_) -
        reinterpret_cast<char*>(&stop_condition_)) + sizeof(uvalue_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* FFPilotPhaseLimit::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // optional .lm.input.ffpilot.FFPilotPhaseLimit.StopCondition stop_condition = 11 [default = FORWARD_FLUXES];
      case 11:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 88)) {
          ::PROTOBUF_NAMESPACE_ID::uint64 val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::lm::input::ffpilot::FFPilotPhaseLimit_StopCondition_IsValid(val))) {
            _internal_set_stop_condition(static_cast<::lm::input::ffpilot::FFPilotPhaseLimit_StopCondition>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(11, val, mutable_unknown_fields());
          }
        } else goto handle_unusual;
        continue;
      // optional uint64 events_per_trajectory = 21;
      case 21:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 168)) {
          _Internal::set_has_events_per_trajectory(&has_bits);
          events_per_trajectory_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional uint64 trajectories_per_phase = 22;
      case 22:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 176)) {
          _Internal::set_has_trajectories_per_phase(&has_bits);
          trajectories_per_phase_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional double dvalue = 31;
      case 31:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 249)) {
          _Internal::set_has_dvalue(&has_bits);
          dvalue_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else goto handle_unusual;
        continue;
      // optional uint64 uvalue = 32;
      case 32:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 0)) {
          _Internal::set_has_uvalue(&has_bits);
          uvalue_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
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

::PROTOBUF_NAMESPACE_ID::uint8* FFPilotPhaseLimit::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.input.ffpilot.FFPilotPhaseLimit)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // optional .lm.input.ffpilot.FFPilotPhaseLimit.StopCondition stop_condition = 11 [default = FORWARD_FLUXES];
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteEnumToArray(
      11, this->_internal_stop_condition(), target);
  }

  // optional uint64 events_per_trajectory = 21;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt64ToArray(21, this->_internal_events_per_trajectory(), target);
  }

  // optional uint64 trajectories_per_phase = 22;
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt64ToArray(22, this->_internal_trajectories_per_phase(), target);
  }

  // optional double dvalue = 31;
  if (cached_has_bits & 0x00000008u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(31, this->_internal_dvalue(), target);
  }

  // optional uint64 uvalue = 32;
  if (cached_has_bits & 0x00000010u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt64ToArray(32, this->_internal_uvalue(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.input.ffpilot.FFPilotPhaseLimit)
  return target;
}

size_t FFPilotPhaseLimit::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.input.ffpilot.FFPilotPhaseLimit)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    // optional .lm.input.ffpilot.FFPilotPhaseLimit.StopCondition stop_condition = 11 [default = FORWARD_FLUXES];
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::EnumSize(this->_internal_stop_condition());
    }

    // optional uint64 events_per_trajectory = 21;
    if (cached_has_bits & 0x00000002u) {
      total_size += 2 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt64Size(
          this->_internal_events_per_trajectory());
    }

    // optional uint64 trajectories_per_phase = 22;
    if (cached_has_bits & 0x00000004u) {
      total_size += 2 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt64Size(
          this->_internal_trajectories_per_phase());
    }

    // optional double dvalue = 31;
    if (cached_has_bits & 0x00000008u) {
      total_size += 2 + 8;
    }

    // optional uint64 uvalue = 32;
    if (cached_has_bits & 0x00000010u) {
      total_size += 2 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt64Size(
          this->_internal_uvalue());
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

void FFPilotPhaseLimit::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.input.ffpilot.FFPilotPhaseLimit)
  GOOGLE_DCHECK_NE(&from, this);
  const FFPilotPhaseLimit* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<FFPilotPhaseLimit>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.input.ffpilot.FFPilotPhaseLimit)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.input.ffpilot.FFPilotPhaseLimit)
    MergeFrom(*source);
  }
}

void FFPilotPhaseLimit::MergeFrom(const FFPilotPhaseLimit& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.input.ffpilot.FFPilotPhaseLimit)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x0000001fu) {
    if (cached_has_bits & 0x00000001u) {
      stop_condition_ = from.stop_condition_;
    }
    if (cached_has_bits & 0x00000002u) {
      events_per_trajectory_ = from.events_per_trajectory_;
    }
    if (cached_has_bits & 0x00000004u) {
      trajectories_per_phase_ = from.trajectories_per_phase_;
    }
    if (cached_has_bits & 0x00000008u) {
      dvalue_ = from.dvalue_;
    }
    if (cached_has_bits & 0x00000010u) {
      uvalue_ = from.uvalue_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void FFPilotPhaseLimit::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.input.ffpilot.FFPilotPhaseLimit)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void FFPilotPhaseLimit::CopyFrom(const FFPilotPhaseLimit& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.input.ffpilot.FFPilotPhaseLimit)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool FFPilotPhaseLimit::IsInitialized() const {
  return true;
}

void FFPilotPhaseLimit::InternalSwap(FFPilotPhaseLimit* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(FFPilotPhaseLimit, uvalue_)
      + sizeof(FFPilotPhaseLimit::uvalue_)
      - PROTOBUF_FIELD_OFFSET(FFPilotPhaseLimit, stop_condition_)>(
          reinterpret_cast<char*>(&stop_condition_),
          reinterpret_cast<char*>(&other->stop_condition_));
}

::PROTOBUF_NAMESPACE_ID::Metadata FFPilotPhaseLimit::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace ffpilot
}  // namespace input
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::input::ffpilot::FFPilotPhaseLimit* Arena::CreateMaybeMessage< ::lm::input::ffpilot::FFPilotPhaseLimit >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::input::ffpilot::FFPilotPhaseLimit >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>