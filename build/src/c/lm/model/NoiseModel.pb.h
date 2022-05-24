// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/model/NoiseModel.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2fmodel_2fNoiseModel_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2fmodel_2fNoiseModel_2eproto

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
#include "robertslab/pbuf/NDArray.pb.h"
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_lm_2fmodel_2fNoiseModel_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2fmodel_2fNoiseModel_2eproto {
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
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fmodel_2fNoiseModel_2eproto;
namespace lm {
namespace model {
class NoiseModel;
class NoiseModelDefaultTypeInternal;
extern NoiseModelDefaultTypeInternal _NoiseModel_default_instance_;
}  // namespace model
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::model::NoiseModel* Arena::CreateMaybeMessage<::lm::model::NoiseModel>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace model {

// ===================================================================

class NoiseModel PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.model.NoiseModel) */ {
 public:
  inline NoiseModel() : NoiseModel(nullptr) {}
  virtual ~NoiseModel();

  NoiseModel(const NoiseModel& from);
  NoiseModel(NoiseModel&& from) noexcept
    : NoiseModel() {
    *this = ::std::move(from);
  }

  inline NoiseModel& operator=(const NoiseModel& from) {
    CopyFrom(from);
    return *this;
  }
  inline NoiseModel& operator=(NoiseModel&& from) noexcept {
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
  static const NoiseModel& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const NoiseModel* internal_default_instance() {
    return reinterpret_cast<const NoiseModel*>(
               &_NoiseModel_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(NoiseModel& a, NoiseModel& b) {
    a.Swap(&b);
  }
  inline void Swap(NoiseModel* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(NoiseModel* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline NoiseModel* New() const final {
    return CreateMaybeMessage<NoiseModel>(nullptr);
  }

  NoiseModel* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<NoiseModel>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const NoiseModel& from);
  void MergeFrom(const NoiseModel& from);
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
  void InternalSwap(NoiseModel* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.model.NoiseModel";
  }
  protected:
  explicit NoiseModel(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2fmodel_2fNoiseModel_2eproto);
    return ::descriptor_table_lm_2fmodel_2fNoiseModel_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kProcessTypesFieldNumber = 3,
    kProcessParametersFieldNumber = 4,
    kReactionDependenciesFieldNumber = 5,
    kProcessUpdateIntervalFieldNumber = 2,
    kNumberProcessesFieldNumber = 1,
  };
  // required .robertslab.pbuf.NDArray process_types = 3;
  bool has_process_types() const;
  private:
  bool _internal_has_process_types() const;
  public:
  void clear_process_types();
  const ::robertslab::pbuf::NDArray& process_types() const;
  ::robertslab::pbuf::NDArray* release_process_types();
  ::robertslab::pbuf::NDArray* mutable_process_types();
  void set_allocated_process_types(::robertslab::pbuf::NDArray* process_types);
  private:
  const ::robertslab::pbuf::NDArray& _internal_process_types() const;
  ::robertslab::pbuf::NDArray* _internal_mutable_process_types();
  public:
  void unsafe_arena_set_allocated_process_types(
      ::robertslab::pbuf::NDArray* process_types);
  ::robertslab::pbuf::NDArray* unsafe_arena_release_process_types();

  // required .robertslab.pbuf.NDArray process_parameters = 4;
  bool has_process_parameters() const;
  private:
  bool _internal_has_process_parameters() const;
  public:
  void clear_process_parameters();
  const ::robertslab::pbuf::NDArray& process_parameters() const;
  ::robertslab::pbuf::NDArray* release_process_parameters();
  ::robertslab::pbuf::NDArray* mutable_process_parameters();
  void set_allocated_process_parameters(::robertslab::pbuf::NDArray* process_parameters);
  private:
  const ::robertslab::pbuf::NDArray& _internal_process_parameters() const;
  ::robertslab::pbuf::NDArray* _internal_mutable_process_parameters();
  public:
  void unsafe_arena_set_allocated_process_parameters(
      ::robertslab::pbuf::NDArray* process_parameters);
  ::robertslab::pbuf::NDArray* unsafe_arena_release_process_parameters();

  // required .robertslab.pbuf.NDArray reaction_dependencies = 5;
  bool has_reaction_dependencies() const;
  private:
  bool _internal_has_reaction_dependencies() const;
  public:
  void clear_reaction_dependencies();
  const ::robertslab::pbuf::NDArray& reaction_dependencies() const;
  ::robertslab::pbuf::NDArray* release_reaction_dependencies();
  ::robertslab::pbuf::NDArray* mutable_reaction_dependencies();
  void set_allocated_reaction_dependencies(::robertslab::pbuf::NDArray* reaction_dependencies);
  private:
  const ::robertslab::pbuf::NDArray& _internal_reaction_dependencies() const;
  ::robertslab::pbuf::NDArray* _internal_mutable_reaction_dependencies();
  public:
  void unsafe_arena_set_allocated_reaction_dependencies(
      ::robertslab::pbuf::NDArray* reaction_dependencies);
  ::robertslab::pbuf::NDArray* unsafe_arena_release_reaction_dependencies();

  // required double process_update_interval = 2;
  bool has_process_update_interval() const;
  private:
  bool _internal_has_process_update_interval() const;
  public:
  void clear_process_update_interval();
  double process_update_interval() const;
  void set_process_update_interval(double value);
  private:
  double _internal_process_update_interval() const;
  void _internal_set_process_update_interval(double value);
  public:

  // required uint32 number_processes = 1;
  bool has_number_processes() const;
  private:
  bool _internal_has_number_processes() const;
  public:
  void clear_number_processes();
  ::PROTOBUF_NAMESPACE_ID::uint32 number_processes() const;
  void set_number_processes(::PROTOBUF_NAMESPACE_ID::uint32 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::uint32 _internal_number_processes() const;
  void _internal_set_number_processes(::PROTOBUF_NAMESPACE_ID::uint32 value);
  public:

  // @@protoc_insertion_point(class_scope:lm.model.NoiseModel)
 private:
  class _Internal;

  // helper for ByteSizeLong()
  size_t RequiredFieldsByteSizeFallback() const;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  ::robertslab::pbuf::NDArray* process_types_;
  ::robertslab::pbuf::NDArray* process_parameters_;
  ::robertslab::pbuf::NDArray* reaction_dependencies_;
  double process_update_interval_;
  ::PROTOBUF_NAMESPACE_ID::uint32 number_processes_;
  friend struct ::TableStruct_lm_2fmodel_2fNoiseModel_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// NoiseModel

// required uint32 number_processes = 1;
inline bool NoiseModel::_internal_has_number_processes() const {
  bool value = (_has_bits_[0] & 0x00000010u) != 0;
  return value;
}
inline bool NoiseModel::has_number_processes() const {
  return _internal_has_number_processes();
}
inline void NoiseModel::clear_number_processes() {
  number_processes_ = 0u;
  _has_bits_[0] &= ~0x00000010u;
}
inline ::PROTOBUF_NAMESPACE_ID::uint32 NoiseModel::_internal_number_processes() const {
  return number_processes_;
}
inline ::PROTOBUF_NAMESPACE_ID::uint32 NoiseModel::number_processes() const {
  // @@protoc_insertion_point(field_get:lm.model.NoiseModel.number_processes)
  return _internal_number_processes();
}
inline void NoiseModel::_internal_set_number_processes(::PROTOBUF_NAMESPACE_ID::uint32 value) {
  _has_bits_[0] |= 0x00000010u;
  number_processes_ = value;
}
inline void NoiseModel::set_number_processes(::PROTOBUF_NAMESPACE_ID::uint32 value) {
  _internal_set_number_processes(value);
  // @@protoc_insertion_point(field_set:lm.model.NoiseModel.number_processes)
}

// required double process_update_interval = 2;
inline bool NoiseModel::_internal_has_process_update_interval() const {
  bool value = (_has_bits_[0] & 0x00000008u) != 0;
  return value;
}
inline bool NoiseModel::has_process_update_interval() const {
  return _internal_has_process_update_interval();
}
inline void NoiseModel::clear_process_update_interval() {
  process_update_interval_ = 0;
  _has_bits_[0] &= ~0x00000008u;
}
inline double NoiseModel::_internal_process_update_interval() const {
  return process_update_interval_;
}
inline double NoiseModel::process_update_interval() const {
  // @@protoc_insertion_point(field_get:lm.model.NoiseModel.process_update_interval)
  return _internal_process_update_interval();
}
inline void NoiseModel::_internal_set_process_update_interval(double value) {
  _has_bits_[0] |= 0x00000008u;
  process_update_interval_ = value;
}
inline void NoiseModel::set_process_update_interval(double value) {
  _internal_set_process_update_interval(value);
  // @@protoc_insertion_point(field_set:lm.model.NoiseModel.process_update_interval)
}

// required .robertslab.pbuf.NDArray process_types = 3;
inline bool NoiseModel::_internal_has_process_types() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  PROTOBUF_ASSUME(!value || process_types_ != nullptr);
  return value;
}
inline bool NoiseModel::has_process_types() const {
  return _internal_has_process_types();
}
inline const ::robertslab::pbuf::NDArray& NoiseModel::_internal_process_types() const {
  const ::robertslab::pbuf::NDArray* p = process_types_;
  return p != nullptr ? *p : *reinterpret_cast<const ::robertslab::pbuf::NDArray*>(
      &::robertslab::pbuf::_NDArray_default_instance_);
}
inline const ::robertslab::pbuf::NDArray& NoiseModel::process_types() const {
  // @@protoc_insertion_point(field_get:lm.model.NoiseModel.process_types)
  return _internal_process_types();
}
inline void NoiseModel::unsafe_arena_set_allocated_process_types(
    ::robertslab::pbuf::NDArray* process_types) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(process_types_);
  }
  process_types_ = process_types;
  if (process_types) {
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.model.NoiseModel.process_types)
}
inline ::robertslab::pbuf::NDArray* NoiseModel::release_process_types() {
  _has_bits_[0] &= ~0x00000001u;
  ::robertslab::pbuf::NDArray* temp = process_types_;
  process_types_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::robertslab::pbuf::NDArray* NoiseModel::unsafe_arena_release_process_types() {
  // @@protoc_insertion_point(field_release:lm.model.NoiseModel.process_types)
  _has_bits_[0] &= ~0x00000001u;
  ::robertslab::pbuf::NDArray* temp = process_types_;
  process_types_ = nullptr;
  return temp;
}
inline ::robertslab::pbuf::NDArray* NoiseModel::_internal_mutable_process_types() {
  _has_bits_[0] |= 0x00000001u;
  if (process_types_ == nullptr) {
    auto* p = CreateMaybeMessage<::robertslab::pbuf::NDArray>(GetArena());
    process_types_ = p;
  }
  return process_types_;
}
inline ::robertslab::pbuf::NDArray* NoiseModel::mutable_process_types() {
  // @@protoc_insertion_point(field_mutable:lm.model.NoiseModel.process_types)
  return _internal_mutable_process_types();
}
inline void NoiseModel::set_allocated_process_types(::robertslab::pbuf::NDArray* process_types) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(process_types_);
  }
  if (process_types) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(process_types)->GetArena();
    if (message_arena != submessage_arena) {
      process_types = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, process_types, submessage_arena);
    }
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  process_types_ = process_types;
  // @@protoc_insertion_point(field_set_allocated:lm.model.NoiseModel.process_types)
}

// required .robertslab.pbuf.NDArray process_parameters = 4;
inline bool NoiseModel::_internal_has_process_parameters() const {
  bool value = (_has_bits_[0] & 0x00000002u) != 0;
  PROTOBUF_ASSUME(!value || process_parameters_ != nullptr);
  return value;
}
inline bool NoiseModel::has_process_parameters() const {
  return _internal_has_process_parameters();
}
inline const ::robertslab::pbuf::NDArray& NoiseModel::_internal_process_parameters() const {
  const ::robertslab::pbuf::NDArray* p = process_parameters_;
  return p != nullptr ? *p : *reinterpret_cast<const ::robertslab::pbuf::NDArray*>(
      &::robertslab::pbuf::_NDArray_default_instance_);
}
inline const ::robertslab::pbuf::NDArray& NoiseModel::process_parameters() const {
  // @@protoc_insertion_point(field_get:lm.model.NoiseModel.process_parameters)
  return _internal_process_parameters();
}
inline void NoiseModel::unsafe_arena_set_allocated_process_parameters(
    ::robertslab::pbuf::NDArray* process_parameters) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(process_parameters_);
  }
  process_parameters_ = process_parameters;
  if (process_parameters) {
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.model.NoiseModel.process_parameters)
}
inline ::robertslab::pbuf::NDArray* NoiseModel::release_process_parameters() {
  _has_bits_[0] &= ~0x00000002u;
  ::robertslab::pbuf::NDArray* temp = process_parameters_;
  process_parameters_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::robertslab::pbuf::NDArray* NoiseModel::unsafe_arena_release_process_parameters() {
  // @@protoc_insertion_point(field_release:lm.model.NoiseModel.process_parameters)
  _has_bits_[0] &= ~0x00000002u;
  ::robertslab::pbuf::NDArray* temp = process_parameters_;
  process_parameters_ = nullptr;
  return temp;
}
inline ::robertslab::pbuf::NDArray* NoiseModel::_internal_mutable_process_parameters() {
  _has_bits_[0] |= 0x00000002u;
  if (process_parameters_ == nullptr) {
    auto* p = CreateMaybeMessage<::robertslab::pbuf::NDArray>(GetArena());
    process_parameters_ = p;
  }
  return process_parameters_;
}
inline ::robertslab::pbuf::NDArray* NoiseModel::mutable_process_parameters() {
  // @@protoc_insertion_point(field_mutable:lm.model.NoiseModel.process_parameters)
  return _internal_mutable_process_parameters();
}
inline void NoiseModel::set_allocated_process_parameters(::robertslab::pbuf::NDArray* process_parameters) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(process_parameters_);
  }
  if (process_parameters) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(process_parameters)->GetArena();
    if (message_arena != submessage_arena) {
      process_parameters = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, process_parameters, submessage_arena);
    }
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  process_parameters_ = process_parameters;
  // @@protoc_insertion_point(field_set_allocated:lm.model.NoiseModel.process_parameters)
}

// required .robertslab.pbuf.NDArray reaction_dependencies = 5;
inline bool NoiseModel::_internal_has_reaction_dependencies() const {
  bool value = (_has_bits_[0] & 0x00000004u) != 0;
  PROTOBUF_ASSUME(!value || reaction_dependencies_ != nullptr);
  return value;
}
inline bool NoiseModel::has_reaction_dependencies() const {
  return _internal_has_reaction_dependencies();
}
inline const ::robertslab::pbuf::NDArray& NoiseModel::_internal_reaction_dependencies() const {
  const ::robertslab::pbuf::NDArray* p = reaction_dependencies_;
  return p != nullptr ? *p : *reinterpret_cast<const ::robertslab::pbuf::NDArray*>(
      &::robertslab::pbuf::_NDArray_default_instance_);
}
inline const ::robertslab::pbuf::NDArray& NoiseModel::reaction_dependencies() const {
  // @@protoc_insertion_point(field_get:lm.model.NoiseModel.reaction_dependencies)
  return _internal_reaction_dependencies();
}
inline void NoiseModel::unsafe_arena_set_allocated_reaction_dependencies(
    ::robertslab::pbuf::NDArray* reaction_dependencies) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(reaction_dependencies_);
  }
  reaction_dependencies_ = reaction_dependencies;
  if (reaction_dependencies) {
    _has_bits_[0] |= 0x00000004u;
  } else {
    _has_bits_[0] &= ~0x00000004u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.model.NoiseModel.reaction_dependencies)
}
inline ::robertslab::pbuf::NDArray* NoiseModel::release_reaction_dependencies() {
  _has_bits_[0] &= ~0x00000004u;
  ::robertslab::pbuf::NDArray* temp = reaction_dependencies_;
  reaction_dependencies_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::robertslab::pbuf::NDArray* NoiseModel::unsafe_arena_release_reaction_dependencies() {
  // @@protoc_insertion_point(field_release:lm.model.NoiseModel.reaction_dependencies)
  _has_bits_[0] &= ~0x00000004u;
  ::robertslab::pbuf::NDArray* temp = reaction_dependencies_;
  reaction_dependencies_ = nullptr;
  return temp;
}
inline ::robertslab::pbuf::NDArray* NoiseModel::_internal_mutable_reaction_dependencies() {
  _has_bits_[0] |= 0x00000004u;
  if (reaction_dependencies_ == nullptr) {
    auto* p = CreateMaybeMessage<::robertslab::pbuf::NDArray>(GetArena());
    reaction_dependencies_ = p;
  }
  return reaction_dependencies_;
}
inline ::robertslab::pbuf::NDArray* NoiseModel::mutable_reaction_dependencies() {
  // @@protoc_insertion_point(field_mutable:lm.model.NoiseModel.reaction_dependencies)
  return _internal_mutable_reaction_dependencies();
}
inline void NoiseModel::set_allocated_reaction_dependencies(::robertslab::pbuf::NDArray* reaction_dependencies) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(reaction_dependencies_);
  }
  if (reaction_dependencies) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(reaction_dependencies)->GetArena();
    if (message_arena != submessage_arena) {
      reaction_dependencies = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, reaction_dependencies, submessage_arena);
    }
    _has_bits_[0] |= 0x00000004u;
  } else {
    _has_bits_[0] &= ~0x00000004u;
  }
  reaction_dependencies_ = reaction_dependencies;
  // @@protoc_insertion_point(field_set_allocated:lm.model.NoiseModel.reaction_dependencies)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace model
}  // namespace lm

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2fmodel_2fNoiseModel_2eproto