// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/message/WorkUnitStatus.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2fmessage_2fWorkUnitStatus_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2fmessage_2fWorkUnitStatus_2eproto

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
#include <google/protobuf/generated_enum_reflection.h>
#include <google/protobuf/unknown_field_set.h>
#include "lm/io/TrajectoryState.pb.h"
#include "lm/types/TrajectoryLimits.pb.h"
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_lm_2fmessage_2fWorkUnitStatus_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2fmessage_2fWorkUnitStatus_2eproto {
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
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fmessage_2fWorkUnitStatus_2eproto;
namespace lm {
namespace message {
class WorkUnitStatus;
class WorkUnitStatusDefaultTypeInternal;
extern WorkUnitStatusDefaultTypeInternal _WorkUnitStatus_default_instance_;
}  // namespace message
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::message::WorkUnitStatus* Arena::CreateMaybeMessage<::lm::message::WorkUnitStatus>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace message {

enum WorkUnitStatus_Status : int {
  WorkUnitStatus_Status_NONE = 0,
  WorkUnitStatus_Status_STEPS_FINISHED = 1,
  WorkUnitStatus_Status_LIMIT_REACHED = 2,
  WorkUnitStatus_Status_ERROR = 3
};
bool WorkUnitStatus_Status_IsValid(int value);
constexpr WorkUnitStatus_Status WorkUnitStatus_Status_Status_MIN = WorkUnitStatus_Status_NONE;
constexpr WorkUnitStatus_Status WorkUnitStatus_Status_Status_MAX = WorkUnitStatus_Status_ERROR;
constexpr int WorkUnitStatus_Status_Status_ARRAYSIZE = WorkUnitStatus_Status_Status_MAX + 1;

const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* WorkUnitStatus_Status_descriptor();
template<typename T>
inline const std::string& WorkUnitStatus_Status_Name(T enum_t_value) {
  static_assert(::std::is_same<T, WorkUnitStatus_Status>::value ||
    ::std::is_integral<T>::value,
    "Incorrect type passed to function WorkUnitStatus_Status_Name.");
  return ::PROTOBUF_NAMESPACE_ID::internal::NameOfEnum(
    WorkUnitStatus_Status_descriptor(), enum_t_value);
}
inline bool WorkUnitStatus_Status_Parse(
    ::PROTOBUF_NAMESPACE_ID::ConstStringParam name, WorkUnitStatus_Status* value) {
  return ::PROTOBUF_NAMESPACE_ID::internal::ParseNamedEnum<WorkUnitStatus_Status>(
    WorkUnitStatus_Status_descriptor(), name, value);
}
// ===================================================================

class WorkUnitStatus PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.message.WorkUnitStatus) */ {
 public:
  inline WorkUnitStatus() : WorkUnitStatus(nullptr) {}
  virtual ~WorkUnitStatus();

  WorkUnitStatus(const WorkUnitStatus& from);
  WorkUnitStatus(WorkUnitStatus&& from) noexcept
    : WorkUnitStatus() {
    *this = ::std::move(from);
  }

  inline WorkUnitStatus& operator=(const WorkUnitStatus& from) {
    CopyFrom(from);
    return *this;
  }
  inline WorkUnitStatus& operator=(WorkUnitStatus&& from) noexcept {
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
  static const WorkUnitStatus& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const WorkUnitStatus* internal_default_instance() {
    return reinterpret_cast<const WorkUnitStatus*>(
               &_WorkUnitStatus_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(WorkUnitStatus& a, WorkUnitStatus& b) {
    a.Swap(&b);
  }
  inline void Swap(WorkUnitStatus* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(WorkUnitStatus* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline WorkUnitStatus* New() const final {
    return CreateMaybeMessage<WorkUnitStatus>(nullptr);
  }

  WorkUnitStatus* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<WorkUnitStatus>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const WorkUnitStatus& from);
  void MergeFrom(const WorkUnitStatus& from);
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
  void InternalSwap(WorkUnitStatus* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.message.WorkUnitStatus";
  }
  protected:
  explicit WorkUnitStatus(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2fmessage_2fWorkUnitStatus_2eproto);
    return ::descriptor_table_lm_2fmessage_2fWorkUnitStatus_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  typedef WorkUnitStatus_Status Status;
  static constexpr Status NONE =
    WorkUnitStatus_Status_NONE;
  static constexpr Status STEPS_FINISHED =
    WorkUnitStatus_Status_STEPS_FINISHED;
  static constexpr Status LIMIT_REACHED =
    WorkUnitStatus_Status_LIMIT_REACHED;
  static constexpr Status ERROR =
    WorkUnitStatus_Status_ERROR;
  static inline bool Status_IsValid(int value) {
    return WorkUnitStatus_Status_IsValid(value);
  }
  static constexpr Status Status_MIN =
    WorkUnitStatus_Status_Status_MIN;
  static constexpr Status Status_MAX =
    WorkUnitStatus_Status_Status_MAX;
  static constexpr int Status_ARRAYSIZE =
    WorkUnitStatus_Status_Status_ARRAYSIZE;
  static inline const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor*
  Status_descriptor() {
    return WorkUnitStatus_Status_descriptor();
  }
  template<typename T>
  static inline const std::string& Status_Name(T enum_t_value) {
    static_assert(::std::is_same<T, Status>::value ||
      ::std::is_integral<T>::value,
      "Incorrect type passed to function Status_Name.");
    return WorkUnitStatus_Status_Name(enum_t_value);
  }
  static inline bool Status_Parse(::PROTOBUF_NAMESPACE_ID::ConstStringParam name,
      Status* value) {
    return WorkUnitStatus_Status_Parse(name, value);
  }

  // accessors -------------------------------------------------------

  enum : int {
    kErrorMessageFieldNumber = 2,
    kFinalStateFieldNumber = 3,
    kLimitReachedFieldNumber = 4,
    kStatusFieldNumber = 1,
  };
  // optional string error_message = 2;
  bool has_error_message() const;
  private:
  bool _internal_has_error_message() const;
  public:
  void clear_error_message();
  const std::string& error_message() const;
  void set_error_message(const std::string& value);
  void set_error_message(std::string&& value);
  void set_error_message(const char* value);
  void set_error_message(const char* value, size_t size);
  std::string* mutable_error_message();
  std::string* release_error_message();
  void set_allocated_error_message(std::string* error_message);
  private:
  const std::string& _internal_error_message() const;
  void _internal_set_error_message(const std::string& value);
  std::string* _internal_mutable_error_message();
  public:

  // required .lm.io.TrajectoryState final_state = 3;
  bool has_final_state() const;
  private:
  bool _internal_has_final_state() const;
  public:
  void clear_final_state();
  const ::lm::io::TrajectoryState& final_state() const;
  ::lm::io::TrajectoryState* release_final_state();
  ::lm::io::TrajectoryState* mutable_final_state();
  void set_allocated_final_state(::lm::io::TrajectoryState* final_state);
  private:
  const ::lm::io::TrajectoryState& _internal_final_state() const;
  ::lm::io::TrajectoryState* _internal_mutable_final_state();
  public:
  void unsafe_arena_set_allocated_final_state(
      ::lm::io::TrajectoryState* final_state);
  ::lm::io::TrajectoryState* unsafe_arena_release_final_state();

  // optional .lm.types.TrajectoryLimit limit_reached = 4;
  bool has_limit_reached() const;
  private:
  bool _internal_has_limit_reached() const;
  public:
  void clear_limit_reached();
  const ::lm::types::TrajectoryLimit& limit_reached() const;
  ::lm::types::TrajectoryLimit* release_limit_reached();
  ::lm::types::TrajectoryLimit* mutable_limit_reached();
  void set_allocated_limit_reached(::lm::types::TrajectoryLimit* limit_reached);
  private:
  const ::lm::types::TrajectoryLimit& _internal_limit_reached() const;
  ::lm::types::TrajectoryLimit* _internal_mutable_limit_reached();
  public:
  void unsafe_arena_set_allocated_limit_reached(
      ::lm::types::TrajectoryLimit* limit_reached);
  ::lm::types::TrajectoryLimit* unsafe_arena_release_limit_reached();

  // required .lm.message.WorkUnitStatus.Status status = 1;
  bool has_status() const;
  private:
  bool _internal_has_status() const;
  public:
  void clear_status();
  ::lm::message::WorkUnitStatus_Status status() const;
  void set_status(::lm::message::WorkUnitStatus_Status value);
  private:
  ::lm::message::WorkUnitStatus_Status _internal_status() const;
  void _internal_set_status(::lm::message::WorkUnitStatus_Status value);
  public:

  // @@protoc_insertion_point(class_scope:lm.message.WorkUnitStatus)
 private:
  class _Internal;

  // helper for ByteSizeLong()
  size_t RequiredFieldsByteSizeFallback() const;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr error_message_;
  ::lm::io::TrajectoryState* final_state_;
  ::lm::types::TrajectoryLimit* limit_reached_;
  int status_;
  friend struct ::TableStruct_lm_2fmessage_2fWorkUnitStatus_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// WorkUnitStatus

// required .lm.message.WorkUnitStatus.Status status = 1;
inline bool WorkUnitStatus::_internal_has_status() const {
  bool value = (_has_bits_[0] & 0x00000008u) != 0;
  return value;
}
inline bool WorkUnitStatus::has_status() const {
  return _internal_has_status();
}
inline void WorkUnitStatus::clear_status() {
  status_ = 0;
  _has_bits_[0] &= ~0x00000008u;
}
inline ::lm::message::WorkUnitStatus_Status WorkUnitStatus::_internal_status() const {
  return static_cast< ::lm::message::WorkUnitStatus_Status >(status_);
}
inline ::lm::message::WorkUnitStatus_Status WorkUnitStatus::status() const {
  // @@protoc_insertion_point(field_get:lm.message.WorkUnitStatus.status)
  return _internal_status();
}
inline void WorkUnitStatus::_internal_set_status(::lm::message::WorkUnitStatus_Status value) {
  assert(::lm::message::WorkUnitStatus_Status_IsValid(value));
  _has_bits_[0] |= 0x00000008u;
  status_ = value;
}
inline void WorkUnitStatus::set_status(::lm::message::WorkUnitStatus_Status value) {
  _internal_set_status(value);
  // @@protoc_insertion_point(field_set:lm.message.WorkUnitStatus.status)
}

// optional string error_message = 2;
inline bool WorkUnitStatus::_internal_has_error_message() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool WorkUnitStatus::has_error_message() const {
  return _internal_has_error_message();
}
inline void WorkUnitStatus::clear_error_message() {
  error_message_.ClearToEmpty(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
  _has_bits_[0] &= ~0x00000001u;
}
inline const std::string& WorkUnitStatus::error_message() const {
  // @@protoc_insertion_point(field_get:lm.message.WorkUnitStatus.error_message)
  return _internal_error_message();
}
inline void WorkUnitStatus::set_error_message(const std::string& value) {
  _internal_set_error_message(value);
  // @@protoc_insertion_point(field_set:lm.message.WorkUnitStatus.error_message)
}
inline std::string* WorkUnitStatus::mutable_error_message() {
  // @@protoc_insertion_point(field_mutable:lm.message.WorkUnitStatus.error_message)
  return _internal_mutable_error_message();
}
inline const std::string& WorkUnitStatus::_internal_error_message() const {
  return error_message_.Get();
}
inline void WorkUnitStatus::_internal_set_error_message(const std::string& value) {
  _has_bits_[0] |= 0x00000001u;
  error_message_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), value, GetArena());
}
inline void WorkUnitStatus::set_error_message(std::string&& value) {
  _has_bits_[0] |= 0x00000001u;
  error_message_.Set(
    &::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), ::std::move(value), GetArena());
  // @@protoc_insertion_point(field_set_rvalue:lm.message.WorkUnitStatus.error_message)
}
inline void WorkUnitStatus::set_error_message(const char* value) {
  GOOGLE_DCHECK(value != nullptr);
  _has_bits_[0] |= 0x00000001u;
  error_message_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), ::std::string(value),
              GetArena());
  // @@protoc_insertion_point(field_set_char:lm.message.WorkUnitStatus.error_message)
}
inline void WorkUnitStatus::set_error_message(const char* value,
    size_t size) {
  _has_bits_[0] |= 0x00000001u;
  error_message_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), ::std::string(
      reinterpret_cast<const char*>(value), size), GetArena());
  // @@protoc_insertion_point(field_set_pointer:lm.message.WorkUnitStatus.error_message)
}
inline std::string* WorkUnitStatus::_internal_mutable_error_message() {
  _has_bits_[0] |= 0x00000001u;
  return error_message_.Mutable(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
}
inline std::string* WorkUnitStatus::release_error_message() {
  // @@protoc_insertion_point(field_release:lm.message.WorkUnitStatus.error_message)
  if (!_internal_has_error_message()) {
    return nullptr;
  }
  _has_bits_[0] &= ~0x00000001u;
  return error_message_.ReleaseNonDefault(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
}
inline void WorkUnitStatus::set_allocated_error_message(std::string* error_message) {
  if (error_message != nullptr) {
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  error_message_.SetAllocated(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), error_message,
      GetArena());
  // @@protoc_insertion_point(field_set_allocated:lm.message.WorkUnitStatus.error_message)
}

// required .lm.io.TrajectoryState final_state = 3;
inline bool WorkUnitStatus::_internal_has_final_state() const {
  bool value = (_has_bits_[0] & 0x00000002u) != 0;
  PROTOBUF_ASSUME(!value || final_state_ != nullptr);
  return value;
}
inline bool WorkUnitStatus::has_final_state() const {
  return _internal_has_final_state();
}
inline const ::lm::io::TrajectoryState& WorkUnitStatus::_internal_final_state() const {
  const ::lm::io::TrajectoryState* p = final_state_;
  return p != nullptr ? *p : *reinterpret_cast<const ::lm::io::TrajectoryState*>(
      &::lm::io::_TrajectoryState_default_instance_);
}
inline const ::lm::io::TrajectoryState& WorkUnitStatus::final_state() const {
  // @@protoc_insertion_point(field_get:lm.message.WorkUnitStatus.final_state)
  return _internal_final_state();
}
inline void WorkUnitStatus::unsafe_arena_set_allocated_final_state(
    ::lm::io::TrajectoryState* final_state) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(final_state_);
  }
  final_state_ = final_state;
  if (final_state) {
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.message.WorkUnitStatus.final_state)
}
inline ::lm::io::TrajectoryState* WorkUnitStatus::release_final_state() {
  _has_bits_[0] &= ~0x00000002u;
  ::lm::io::TrajectoryState* temp = final_state_;
  final_state_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::lm::io::TrajectoryState* WorkUnitStatus::unsafe_arena_release_final_state() {
  // @@protoc_insertion_point(field_release:lm.message.WorkUnitStatus.final_state)
  _has_bits_[0] &= ~0x00000002u;
  ::lm::io::TrajectoryState* temp = final_state_;
  final_state_ = nullptr;
  return temp;
}
inline ::lm::io::TrajectoryState* WorkUnitStatus::_internal_mutable_final_state() {
  _has_bits_[0] |= 0x00000002u;
  if (final_state_ == nullptr) {
    auto* p = CreateMaybeMessage<::lm::io::TrajectoryState>(GetArena());
    final_state_ = p;
  }
  return final_state_;
}
inline ::lm::io::TrajectoryState* WorkUnitStatus::mutable_final_state() {
  // @@protoc_insertion_point(field_mutable:lm.message.WorkUnitStatus.final_state)
  return _internal_mutable_final_state();
}
inline void WorkUnitStatus::set_allocated_final_state(::lm::io::TrajectoryState* final_state) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(final_state_);
  }
  if (final_state) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(final_state)->GetArena();
    if (message_arena != submessage_arena) {
      final_state = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, final_state, submessage_arena);
    }
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  final_state_ = final_state;
  // @@protoc_insertion_point(field_set_allocated:lm.message.WorkUnitStatus.final_state)
}

// optional .lm.types.TrajectoryLimit limit_reached = 4;
inline bool WorkUnitStatus::_internal_has_limit_reached() const {
  bool value = (_has_bits_[0] & 0x00000004u) != 0;
  PROTOBUF_ASSUME(!value || limit_reached_ != nullptr);
  return value;
}
inline bool WorkUnitStatus::has_limit_reached() const {
  return _internal_has_limit_reached();
}
inline const ::lm::types::TrajectoryLimit& WorkUnitStatus::_internal_limit_reached() const {
  const ::lm::types::TrajectoryLimit* p = limit_reached_;
  return p != nullptr ? *p : *reinterpret_cast<const ::lm::types::TrajectoryLimit*>(
      &::lm::types::_TrajectoryLimit_default_instance_);
}
inline const ::lm::types::TrajectoryLimit& WorkUnitStatus::limit_reached() const {
  // @@protoc_insertion_point(field_get:lm.message.WorkUnitStatus.limit_reached)
  return _internal_limit_reached();
}
inline void WorkUnitStatus::unsafe_arena_set_allocated_limit_reached(
    ::lm::types::TrajectoryLimit* limit_reached) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(limit_reached_);
  }
  limit_reached_ = limit_reached;
  if (limit_reached) {
    _has_bits_[0] |= 0x00000004u;
  } else {
    _has_bits_[0] &= ~0x00000004u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.message.WorkUnitStatus.limit_reached)
}
inline ::lm::types::TrajectoryLimit* WorkUnitStatus::release_limit_reached() {
  _has_bits_[0] &= ~0x00000004u;
  ::lm::types::TrajectoryLimit* temp = limit_reached_;
  limit_reached_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::lm::types::TrajectoryLimit* WorkUnitStatus::unsafe_arena_release_limit_reached() {
  // @@protoc_insertion_point(field_release:lm.message.WorkUnitStatus.limit_reached)
  _has_bits_[0] &= ~0x00000004u;
  ::lm::types::TrajectoryLimit* temp = limit_reached_;
  limit_reached_ = nullptr;
  return temp;
}
inline ::lm::types::TrajectoryLimit* WorkUnitStatus::_internal_mutable_limit_reached() {
  _has_bits_[0] |= 0x00000004u;
  if (limit_reached_ == nullptr) {
    auto* p = CreateMaybeMessage<::lm::types::TrajectoryLimit>(GetArena());
    limit_reached_ = p;
  }
  return limit_reached_;
}
inline ::lm::types::TrajectoryLimit* WorkUnitStatus::mutable_limit_reached() {
  // @@protoc_insertion_point(field_mutable:lm.message.WorkUnitStatus.limit_reached)
  return _internal_mutable_limit_reached();
}
inline void WorkUnitStatus::set_allocated_limit_reached(::lm::types::TrajectoryLimit* limit_reached) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(limit_reached_);
  }
  if (limit_reached) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(limit_reached)->GetArena();
    if (message_arena != submessage_arena) {
      limit_reached = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, limit_reached, submessage_arena);
    }
    _has_bits_[0] |= 0x00000004u;
  } else {
    _has_bits_[0] &= ~0x00000004u;
  }
  limit_reached_ = limit_reached;
  // @@protoc_insertion_point(field_set_allocated:lm.message.WorkUnitStatus.limit_reached)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace message
}  // namespace lm

PROTOBUF_NAMESPACE_OPEN

template <> struct is_proto_enum< ::lm::message::WorkUnitStatus_Status> : ::std::true_type {};
template <>
inline const EnumDescriptor* GetEnumDescriptor< ::lm::message::WorkUnitStatus_Status>() {
  return ::lm::message::WorkUnitStatus_Status_descriptor();
}

PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2fmessage_2fWorkUnitStatus_2eproto
