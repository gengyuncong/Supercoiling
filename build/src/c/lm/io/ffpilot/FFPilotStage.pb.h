// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/io/ffpilot/FFPilotStage.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fffpilot_2fFFPilotStage_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fffpilot_2fFFPilotStage_2eproto

#include <limits>
#include <string>

#include <google/protobuf/port_def.inc>
#if PROTOBUF_VERSION < 3013000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers. Please update
#error your headers.
#endif
#if 3013000 < PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers. Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/port_undef.inc>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/arena.h>
#include <google/protobuf/arenastring.h>
#include <google/protobuf/generated_message_table_driven.h>
#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/inlined_string_field.h>
#include <google/protobuf/metadata_lite.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>  // IWYU pragma: export
#include <google/protobuf/extension_set.h>  // IWYU pragma: export
#include <google/protobuf/unknown_field_set.h>
#include "lm/io/ffpilot/FFPilotPhase.pb.h"
#include "lm/io/ffpilot/FFPilotStageOutput.pb.h"
#include "lm/types/OrderParameters.pb.h"
#include "lm/types/Tilings.pb.h"
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_lm_2fio_2fffpilot_2fFFPilotStage_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2fio_2fffpilot_2fFFPilotStage_2eproto {
  static const ::PROTOBUF_NAMESPACE_ID::internal::ParseTableField entries[]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::AuxiliaryParseTableField aux[]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::ParseTable schema[1]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::FieldMetadata field_metadata[];
  static const ::PROTOBUF_NAMESPACE_ID::internal::SerializationTable serialization_table[];
  static const ::PROTOBUF_NAMESPACE_ID::uint32 offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto;
namespace lm {
namespace io {
namespace ffpilot {
class FFPilotStage;
class FFPilotStageDefaultTypeInternal;
extern FFPilotStageDefaultTypeInternal _FFPilotStage_default_instance_;
}  // namespace ffpilot
}  // namespace io
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::io::ffpilot::FFPilotStage* Arena::CreateMaybeMessage<::lm::io::ffpilot::FFPilotStage>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace io {
namespace ffpilot {

// ===================================================================

class FFPilotStage PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.io.ffpilot.FFPilotStage) */ {
 public:
  inline FFPilotStage() : FFPilotStage(nullptr) {}
  virtual ~FFPilotStage();

  FFPilotStage(const FFPilotStage& from);
  FFPilotStage(FFPilotStage&& from) noexcept
    : FFPilotStage() {
    *this = ::std::move(from);
  }

  inline FFPilotStage& operator=(const FFPilotStage& from) {
    CopyFrom(from);
    return *this;
  }
  inline FFPilotStage& operator=(FFPilotStage&& from) noexcept {
    if (GetArena() == from.GetArena()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }

  inline const ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet& unknown_fields() const {
    return _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance);
  }
  inline ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet* mutable_unknown_fields() {
    return _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
  }

  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* descriptor() {
    return GetDescriptor();
  }
  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* GetDescriptor() {
    return GetMetadataStatic().descriptor;
  }
  static const ::PROTOBUF_NAMESPACE_ID::Reflection* GetReflection() {
    return GetMetadataStatic().reflection;
  }
  static const FFPilotStage& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const FFPilotStage* internal_default_instance() {
    return reinterpret_cast<const FFPilotStage*>(
               &_FFPilotStage_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(FFPilotStage& a, FFPilotStage& b) {
    a.Swap(&b);
  }
  inline void Swap(FFPilotStage* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(FFPilotStage* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline FFPilotStage* New() const final {
    return CreateMaybeMessage<FFPilotStage>(nullptr);
  }

  FFPilotStage* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<FFPilotStage>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const FFPilotStage& from);
  void MergeFrom(const FFPilotStage& from);
  PROTOBUF_ATTRIBUTE_REINITIALIZES void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  const char* _InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) final;
  ::PROTOBUF_NAMESPACE_ID::uint8* _InternalSerialize(
      ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  inline void SharedCtor();
  inline void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(FFPilotStage* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.io.ffpilot.FFPilotStage";
  }
  protected:
  explicit FFPilotStage(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto);
    return ::descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kFfpilotPhasesFieldNumber = 201,
    kTilingFieldNumber = 3,
    kOrderParameterFieldNumber = 5,
    kPilotStageOutputFieldNumber = 301,
    kBasinIdFieldNumber = 1,
    kTilingIdFieldNumber = 2,
    kReplicateIdFieldNumber = 4,
    kIsPilotStageFieldNumber = 21,
    kNeedsPilotStageFieldNumber = 22,
  };
  // repeated .lm.io.ffpilot.FFPilotPhase ffpilot_phases = 201;
  int ffpilot_phases_size() const;
  private:
  int _internal_ffpilot_phases_size() const;
  public:
  void clear_ffpilot_phases();
  ::lm::io::ffpilot::FFPilotPhase* mutable_ffpilot_phases(int index);
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotPhase >*
      mutable_ffpilot_phases();
  private:
  const ::lm::io::ffpilot::FFPilotPhase& _internal_ffpilot_phases(int index) const;
  ::lm::io::ffpilot::FFPilotPhase* _internal_add_ffpilot_phases();
  public:
  const ::lm::io::ffpilot::FFPilotPhase& ffpilot_phases(int index) const;
  ::lm::io::ffpilot::FFPilotPhase* add_ffpilot_phases();
  const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotPhase >&
      ffpilot_phases() const;

  // optional .lm.types.Tiling tiling = 3;
  bool has_tiling() const;
  private:
  bool _internal_has_tiling() const;
  public:
  void clear_tiling();
  const ::lm::types::Tiling& tiling() const;
  ::lm::types::Tiling* release_tiling();
  ::lm::types::Tiling* mutable_tiling();
  void set_allocated_tiling(::lm::types::Tiling* tiling);
  private:
  const ::lm::types::Tiling& _internal_tiling() const;
  ::lm::types::Tiling* _internal_mutable_tiling();
  public:
  void unsafe_arena_set_allocated_tiling(
      ::lm::types::Tiling* tiling);
  ::lm::types::Tiling* unsafe_arena_release_tiling();

  // optional .lm.types.OrderParameters order_parameter = 5;
  bool has_order_parameter() const;
  private:
  bool _internal_has_order_parameter() const;
  public:
  void clear_order_parameter();
  const ::lm::types::OrderParameters& order_parameter() const;
  ::lm::types::OrderParameters* release_order_parameter();
  ::lm::types::OrderParameters* mutable_order_parameter();
  void set_allocated_order_parameter(::lm::types::OrderParameters* order_parameter);
  private:
  const ::lm::types::OrderParameters& _internal_order_parameter() const;
  ::lm::types::OrderParameters* _internal_mutable_order_parameter();
  public:
  void unsafe_arena_set_allocated_order_parameter(
      ::lm::types::OrderParameters* order_parameter);
  ::lm::types::OrderParameters* unsafe_arena_release_order_parameter();

  // optional .lm.io.ffpilot.FFPilotStageOutput pilot_stage_output = 301;
  bool has_pilot_stage_output() const;
  private:
  bool _internal_has_pilot_stage_output() const;
  public:
  void clear_pilot_stage_output();
  const ::lm::io::ffpilot::FFPilotStageOutput& pilot_stage_output() const;
  ::lm::io::ffpilot::FFPilotStageOutput* release_pilot_stage_output();
  ::lm::io::ffpilot::FFPilotStageOutput* mutable_pilot_stage_output();
  void set_allocated_pilot_stage_output(::lm::io::ffpilot::FFPilotStageOutput* pilot_stage_output);
  private:
  const ::lm::io::ffpilot::FFPilotStageOutput& _internal_pilot_stage_output() const;
  ::lm::io::ffpilot::FFPilotStageOutput* _internal_mutable_pilot_stage_output();
  public:
  void unsafe_arena_set_allocated_pilot_stage_output(
      ::lm::io::ffpilot::FFPilotStageOutput* pilot_stage_output);
  ::lm::io::ffpilot::FFPilotStageOutput* unsafe_arena_release_pilot_stage_output();

  // optional int64 basin_id = 1;
  bool has_basin_id() const;
  private:
  bool _internal_has_basin_id() const;
  public:
  void clear_basin_id();
  ::PROTOBUF_NAMESPACE_ID::int64 basin_id() const;
  void set_basin_id(::PROTOBUF_NAMESPACE_ID::int64 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::int64 _internal_basin_id() const;
  void _internal_set_basin_id(::PROTOBUF_NAMESPACE_ID::int64 value);
  public:

  // optional uint64 tiling_id = 2;
  bool has_tiling_id() const;
  private:
  bool _internal_has_tiling_id() const;
  public:
  void clear_tiling_id();
  ::PROTOBUF_NAMESPACE_ID::uint64 tiling_id() const;
  void set_tiling_id(::PROTOBUF_NAMESPACE_ID::uint64 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::uint64 _internal_tiling_id() const;
  void _internal_set_tiling_id(::PROTOBUF_NAMESPACE_ID::uint64 value);
  public:

  // optional uint64 replicate_id = 4 [default = 0];
  bool has_replicate_id() const;
  private:
  bool _internal_has_replicate_id() const;
  public:
  void clear_replicate_id();
  ::PROTOBUF_NAMESPACE_ID::uint64 replicate_id() const;
  void set_replicate_id(::PROTOBUF_NAMESPACE_ID::uint64 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::uint64 _internal_replicate_id() const;
  void _internal_set_replicate_id(::PROTOBUF_NAMESPACE_ID::uint64 value);
  public:

  // optional bool is_pilot_stage = 21 [default = false];
  bool has_is_pilot_stage() const;
  private:
  bool _internal_has_is_pilot_stage() const;
  public:
  void clear_is_pilot_stage();
  bool is_pilot_stage() const;
  void set_is_pilot_stage(bool value);
  private:
  bool _internal_is_pilot_stage() const;
  void _internal_set_is_pilot_stage(bool value);
  public:

  // optional bool needs_pilot_stage = 22 [default = false];
  bool has_needs_pilot_stage() const;
  private:
  bool _internal_has_needs_pilot_stage() const;
  public:
  void clear_needs_pilot_stage();
  bool needs_pilot_stage() const;
  void set_needs_pilot_stage(bool value);
  private:
  bool _internal_needs_pilot_stage() const;
  void _internal_set_needs_pilot_stage(bool value);
  public:

  // @@protoc_insertion_point(class_scope:lm.io.ffpilot.FFPilotStage)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotPhase > ffpilot_phases_;
  ::lm::types::Tiling* tiling_;
  ::lm::types::OrderParameters* order_parameter_;
  ::lm::io::ffpilot::FFPilotStageOutput* pilot_stage_output_;
  ::PROTOBUF_NAMESPACE_ID::int64 basin_id_;
  ::PROTOBUF_NAMESPACE_ID::uint64 tiling_id_;
  ::PROTOBUF_NAMESPACE_ID::uint64 replicate_id_;
  bool is_pilot_stage_;
  bool needs_pilot_stage_;
  friend struct ::TableStruct_lm_2fio_2fffpilot_2fFFPilotStage_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// FFPilotStage

// optional int64 basin_id = 1;
inline bool FFPilotStage::_internal_has_basin_id() const {
  bool value = (_has_bits_[0] & 0x00000008u) != 0;
  return value;
}
inline bool FFPilotStage::has_basin_id() const {
  return _internal_has_basin_id();
}
inline void FFPilotStage::clear_basin_id() {
  basin_id_ = PROTOBUF_LONGLONG(0);
  _has_bits_[0] &= ~0x00000008u;
}
inline ::PROTOBUF_NAMESPACE_ID::int64 FFPilotStage::_internal_basin_id() const {
  return basin_id_;
}
inline ::PROTOBUF_NAMESPACE_ID::int64 FFPilotStage::basin_id() const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotStage.basin_id)
  return _internal_basin_id();
}
inline void FFPilotStage::_internal_set_basin_id(::PROTOBUF_NAMESPACE_ID::int64 value) {
  _has_bits_[0] |= 0x00000008u;
  basin_id_ = value;
}
inline void FFPilotStage::set_basin_id(::PROTOBUF_NAMESPACE_ID::int64 value) {
  _internal_set_basin_id(value);
  // @@protoc_insertion_point(field_set:lm.io.ffpilot.FFPilotStage.basin_id)
}

// optional uint64 tiling_id = 2;
inline bool FFPilotStage::_internal_has_tiling_id() const {
  bool value = (_has_bits_[0] & 0x00000010u) != 0;
  return value;
}
inline bool FFPilotStage::has_tiling_id() const {
  return _internal_has_tiling_id();
}
inline void FFPilotStage::clear_tiling_id() {
  tiling_id_ = PROTOBUF_ULONGLONG(0);
  _has_bits_[0] &= ~0x00000010u;
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 FFPilotStage::_internal_tiling_id() const {
  return tiling_id_;
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 FFPilotStage::tiling_id() const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotStage.tiling_id)
  return _internal_tiling_id();
}
inline void FFPilotStage::_internal_set_tiling_id(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  _has_bits_[0] |= 0x00000010u;
  tiling_id_ = value;
}
inline void FFPilotStage::set_tiling_id(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  _internal_set_tiling_id(value);
  // @@protoc_insertion_point(field_set:lm.io.ffpilot.FFPilotStage.tiling_id)
}

// optional uint64 replicate_id = 4 [default = 0];
inline bool FFPilotStage::_internal_has_replicate_id() const {
  bool value = (_has_bits_[0] & 0x00000020u) != 0;
  return value;
}
inline bool FFPilotStage::has_replicate_id() const {
  return _internal_has_replicate_id();
}
inline void FFPilotStage::clear_replicate_id() {
  replicate_id_ = PROTOBUF_ULONGLONG(0);
  _has_bits_[0] &= ~0x00000020u;
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 FFPilotStage::_internal_replicate_id() const {
  return replicate_id_;
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 FFPilotStage::replicate_id() const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotStage.replicate_id)
  return _internal_replicate_id();
}
inline void FFPilotStage::_internal_set_replicate_id(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  _has_bits_[0] |= 0x00000020u;
  replicate_id_ = value;
}
inline void FFPilotStage::set_replicate_id(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  _internal_set_replicate_id(value);
  // @@protoc_insertion_point(field_set:lm.io.ffpilot.FFPilotStage.replicate_id)
}

// optional .lm.types.Tiling tiling = 3;
inline bool FFPilotStage::_internal_has_tiling() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  PROTOBUF_ASSUME(!value || tiling_ != nullptr);
  return value;
}
inline bool FFPilotStage::has_tiling() const {
  return _internal_has_tiling();
}
inline const ::lm::types::Tiling& FFPilotStage::_internal_tiling() const {
  const ::lm::types::Tiling* p = tiling_;
  return p != nullptr ? *p : *reinterpret_cast<const ::lm::types::Tiling*>(
      &::lm::types::_Tiling_default_instance_);
}
inline const ::lm::types::Tiling& FFPilotStage::tiling() const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotStage.tiling)
  return _internal_tiling();
}
inline void FFPilotStage::unsafe_arena_set_allocated_tiling(
    ::lm::types::Tiling* tiling) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(tiling_);
  }
  tiling_ = tiling;
  if (tiling) {
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.io.ffpilot.FFPilotStage.tiling)
}
inline ::lm::types::Tiling* FFPilotStage::release_tiling() {
  _has_bits_[0] &= ~0x00000001u;
  ::lm::types::Tiling* temp = tiling_;
  tiling_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::lm::types::Tiling* FFPilotStage::unsafe_arena_release_tiling() {
  // @@protoc_insertion_point(field_release:lm.io.ffpilot.FFPilotStage.tiling)
  _has_bits_[0] &= ~0x00000001u;
  ::lm::types::Tiling* temp = tiling_;
  tiling_ = nullptr;
  return temp;
}
inline ::lm::types::Tiling* FFPilotStage::_internal_mutable_tiling() {
  _has_bits_[0] |= 0x00000001u;
  if (tiling_ == nullptr) {
    auto* p = CreateMaybeMessage<::lm::types::Tiling>(GetArena());
    tiling_ = p;
  }
  return tiling_;
}
inline ::lm::types::Tiling* FFPilotStage::mutable_tiling() {
  // @@protoc_insertion_point(field_mutable:lm.io.ffpilot.FFPilotStage.tiling)
  return _internal_mutable_tiling();
}
inline void FFPilotStage::set_allocated_tiling(::lm::types::Tiling* tiling) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(tiling_);
  }
  if (tiling) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(tiling)->GetArena();
    if (message_arena != submessage_arena) {
      tiling = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, tiling, submessage_arena);
    }
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  tiling_ = tiling;
  // @@protoc_insertion_point(field_set_allocated:lm.io.ffpilot.FFPilotStage.tiling)
}

// optional .lm.types.OrderParameters order_parameter = 5;
inline bool FFPilotStage::_internal_has_order_parameter() const {
  bool value = (_has_bits_[0] & 0x00000002u) != 0;
  PROTOBUF_ASSUME(!value || order_parameter_ != nullptr);
  return value;
}
inline bool FFPilotStage::has_order_parameter() const {
  return _internal_has_order_parameter();
}
inline const ::lm::types::OrderParameters& FFPilotStage::_internal_order_parameter() const {
  const ::lm::types::OrderParameters* p = order_parameter_;
  return p != nullptr ? *p : *reinterpret_cast<const ::lm::types::OrderParameters*>(
      &::lm::types::_OrderParameters_default_instance_);
}
inline const ::lm::types::OrderParameters& FFPilotStage::order_parameter() const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotStage.order_parameter)
  return _internal_order_parameter();
}
inline void FFPilotStage::unsafe_arena_set_allocated_order_parameter(
    ::lm::types::OrderParameters* order_parameter) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(order_parameter_);
  }
  order_parameter_ = order_parameter;
  if (order_parameter) {
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.io.ffpilot.FFPilotStage.order_parameter)
}
inline ::lm::types::OrderParameters* FFPilotStage::release_order_parameter() {
  _has_bits_[0] &= ~0x00000002u;
  ::lm::types::OrderParameters* temp = order_parameter_;
  order_parameter_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::lm::types::OrderParameters* FFPilotStage::unsafe_arena_release_order_parameter() {
  // @@protoc_insertion_point(field_release:lm.io.ffpilot.FFPilotStage.order_parameter)
  _has_bits_[0] &= ~0x00000002u;
  ::lm::types::OrderParameters* temp = order_parameter_;
  order_parameter_ = nullptr;
  return temp;
}
inline ::lm::types::OrderParameters* FFPilotStage::_internal_mutable_order_parameter() {
  _has_bits_[0] |= 0x00000002u;
  if (order_parameter_ == nullptr) {
    auto* p = CreateMaybeMessage<::lm::types::OrderParameters>(GetArena());
    order_parameter_ = p;
  }
  return order_parameter_;
}
inline ::lm::types::OrderParameters* FFPilotStage::mutable_order_parameter() {
  // @@protoc_insertion_point(field_mutable:lm.io.ffpilot.FFPilotStage.order_parameter)
  return _internal_mutable_order_parameter();
}
inline void FFPilotStage::set_allocated_order_parameter(::lm::types::OrderParameters* order_parameter) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(order_parameter_);
  }
  if (order_parameter) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(order_parameter)->GetArena();
    if (message_arena != submessage_arena) {
      order_parameter = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, order_parameter, submessage_arena);
    }
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  order_parameter_ = order_parameter;
  // @@protoc_insertion_point(field_set_allocated:lm.io.ffpilot.FFPilotStage.order_parameter)
}

// optional bool is_pilot_stage = 21 [default = false];
inline bool FFPilotStage::_internal_has_is_pilot_stage() const {
  bool value = (_has_bits_[0] & 0x00000040u) != 0;
  return value;
}
inline bool FFPilotStage::has_is_pilot_stage() const {
  return _internal_has_is_pilot_stage();
}
inline void FFPilotStage::clear_is_pilot_stage() {
  is_pilot_stage_ = false;
  _has_bits_[0] &= ~0x00000040u;
}
inline bool FFPilotStage::_internal_is_pilot_stage() const {
  return is_pilot_stage_;
}
inline bool FFPilotStage::is_pilot_stage() const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotStage.is_pilot_stage)
  return _internal_is_pilot_stage();
}
inline void FFPilotStage::_internal_set_is_pilot_stage(bool value) {
  _has_bits_[0] |= 0x00000040u;
  is_pilot_stage_ = value;
}
inline void FFPilotStage::set_is_pilot_stage(bool value) {
  _internal_set_is_pilot_stage(value);
  // @@protoc_insertion_point(field_set:lm.io.ffpilot.FFPilotStage.is_pilot_stage)
}

// optional bool needs_pilot_stage = 22 [default = false];
inline bool FFPilotStage::_internal_has_needs_pilot_stage() const {
  bool value = (_has_bits_[0] & 0x00000080u) != 0;
  return value;
}
inline bool FFPilotStage::has_needs_pilot_stage() const {
  return _internal_has_needs_pilot_stage();
}
inline void FFPilotStage::clear_needs_pilot_stage() {
  needs_pilot_stage_ = false;
  _has_bits_[0] &= ~0x00000080u;
}
inline bool FFPilotStage::_internal_needs_pilot_stage() const {
  return needs_pilot_stage_;
}
inline bool FFPilotStage::needs_pilot_stage() const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotStage.needs_pilot_stage)
  return _internal_needs_pilot_stage();
}
inline void FFPilotStage::_internal_set_needs_pilot_stage(bool value) {
  _has_bits_[0] |= 0x00000080u;
  needs_pilot_stage_ = value;
}
inline void FFPilotStage::set_needs_pilot_stage(bool value) {
  _internal_set_needs_pilot_stage(value);
  // @@protoc_insertion_point(field_set:lm.io.ffpilot.FFPilotStage.needs_pilot_stage)
}

// repeated .lm.io.ffpilot.FFPilotPhase ffpilot_phases = 201;
inline int FFPilotStage::_internal_ffpilot_phases_size() const {
  return ffpilot_phases_.size();
}
inline int FFPilotStage::ffpilot_phases_size() const {
  return _internal_ffpilot_phases_size();
}
inline ::lm::io::ffpilot::FFPilotPhase* FFPilotStage::mutable_ffpilot_phases(int index) {
  // @@protoc_insertion_point(field_mutable:lm.io.ffpilot.FFPilotStage.ffpilot_phases)
  return ffpilot_phases_.Mutable(index);
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotPhase >*
FFPilotStage::mutable_ffpilot_phases() {
  // @@protoc_insertion_point(field_mutable_list:lm.io.ffpilot.FFPilotStage.ffpilot_phases)
  return &ffpilot_phases_;
}
inline const ::lm::io::ffpilot::FFPilotPhase& FFPilotStage::_internal_ffpilot_phases(int index) const {
  return ffpilot_phases_.Get(index);
}
inline const ::lm::io::ffpilot::FFPilotPhase& FFPilotStage::ffpilot_phases(int index) const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotStage.ffpilot_phases)
  return _internal_ffpilot_phases(index);
}
inline ::lm::io::ffpilot::FFPilotPhase* FFPilotStage::_internal_add_ffpilot_phases() {
  return ffpilot_phases_.Add();
}
inline ::lm::io::ffpilot::FFPilotPhase* FFPilotStage::add_ffpilot_phases() {
  // @@protoc_insertion_point(field_add:lm.io.ffpilot.FFPilotStage.ffpilot_phases)
  return _internal_add_ffpilot_phases();
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotPhase >&
FFPilotStage::ffpilot_phases() const {
  // @@protoc_insertion_point(field_list:lm.io.ffpilot.FFPilotStage.ffpilot_phases)
  return ffpilot_phases_;
}

// optional .lm.io.ffpilot.FFPilotStageOutput pilot_stage_output = 301;
inline bool FFPilotStage::_internal_has_pilot_stage_output() const {
  bool value = (_has_bits_[0] & 0x00000004u) != 0;
  PROTOBUF_ASSUME(!value || pilot_stage_output_ != nullptr);
  return value;
}
inline bool FFPilotStage::has_pilot_stage_output() const {
  return _internal_has_pilot_stage_output();
}
inline const ::lm::io::ffpilot::FFPilotStageOutput& FFPilotStage::_internal_pilot_stage_output() const {
  const ::lm::io::ffpilot::FFPilotStageOutput* p = pilot_stage_output_;
  return p != nullptr ? *p : *reinterpret_cast<const ::lm::io::ffpilot::FFPilotStageOutput*>(
      &::lm::io::ffpilot::_FFPilotStageOutput_default_instance_);
}
inline const ::lm::io::ffpilot::FFPilotStageOutput& FFPilotStage::pilot_stage_output() const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotStage.pilot_stage_output)
  return _internal_pilot_stage_output();
}
inline void FFPilotStage::unsafe_arena_set_allocated_pilot_stage_output(
    ::lm::io::ffpilot::FFPilotStageOutput* pilot_stage_output) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(pilot_stage_output_);
  }
  pilot_stage_output_ = pilot_stage_output;
  if (pilot_stage_output) {
    _has_bits_[0] |= 0x00000004u;
  } else {
    _has_bits_[0] &= ~0x00000004u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.io.ffpilot.FFPilotStage.pilot_stage_output)
}
inline ::lm::io::ffpilot::FFPilotStageOutput* FFPilotStage::release_pilot_stage_output() {
  _has_bits_[0] &= ~0x00000004u;
  ::lm::io::ffpilot::FFPilotStageOutput* temp = pilot_stage_output_;
  pilot_stage_output_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::lm::io::ffpilot::FFPilotStageOutput* FFPilotStage::unsafe_arena_release_pilot_stage_output() {
  // @@protoc_insertion_point(field_release:lm.io.ffpilot.FFPilotStage.pilot_stage_output)
  _has_bits_[0] &= ~0x00000004u;
  ::lm::io::ffpilot::FFPilotStageOutput* temp = pilot_stage_output_;
  pilot_stage_output_ = nullptr;
  return temp;
}
inline ::lm::io::ffpilot::FFPilotStageOutput* FFPilotStage::_internal_mutable_pilot_stage_output() {
  _has_bits_[0] |= 0x00000004u;
  if (pilot_stage_output_ == nullptr) {
    auto* p = CreateMaybeMessage<::lm::io::ffpilot::FFPilotStageOutput>(GetArena());
    pilot_stage_output_ = p;
  }
  return pilot_stage_output_;
}
inline ::lm::io::ffpilot::FFPilotStageOutput* FFPilotStage::mutable_pilot_stage_output() {
  // @@protoc_insertion_point(field_mutable:lm.io.ffpilot.FFPilotStage.pilot_stage_output)
  return _internal_mutable_pilot_stage_output();
}
inline void FFPilotStage::set_allocated_pilot_stage_output(::lm::io::ffpilot::FFPilotStageOutput* pilot_stage_output) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(pilot_stage_output_);
  }
  if (pilot_stage_output) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(pilot_stage_output)->GetArena();
    if (message_arena != submessage_arena) {
      pilot_stage_output = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, pilot_stage_output, submessage_arena);
    }
    _has_bits_[0] |= 0x00000004u;
  } else {
    _has_bits_[0] &= ~0x00000004u;
  }
  pilot_stage_output_ = pilot_stage_output;
  // @@protoc_insertion_point(field_set_allocated:lm.io.ffpilot.FFPilotStage.pilot_stage_output)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace ffpilot
}  // namespace io
}  // namespace lm

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fffpilot_2fFFPilotStage_2eproto
