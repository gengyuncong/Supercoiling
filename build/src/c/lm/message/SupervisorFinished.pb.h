// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/message/SupervisorFinished.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2fmessage_2fSupervisorFinished_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2fmessage_2fSupervisorFinished_2eproto

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
#include "lm/message/Endpoint.pb.h"
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_lm_2fmessage_2fSupervisorFinished_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2fmessage_2fSupervisorFinished_2eproto {
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
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fmessage_2fSupervisorFinished_2eproto;
namespace lm {
namespace message {
class SupervisorFinished;
class SupervisorFinishedDefaultTypeInternal;
extern SupervisorFinishedDefaultTypeInternal _SupervisorFinished_default_instance_;
}  // namespace message
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::message::SupervisorFinished* Arena::CreateMaybeMessage<::lm::message::SupervisorFinished>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace message {

// ===================================================================

class SupervisorFinished PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.message.SupervisorFinished) */ {
 public:
  inline SupervisorFinished() : SupervisorFinished(nullptr) {}
  virtual ~SupervisorFinished();

  SupervisorFinished(const SupervisorFinished& from);
  SupervisorFinished(SupervisorFinished&& from) noexcept
    : SupervisorFinished() {
    *this = ::std::move(from);
  }

  inline SupervisorFinished& operator=(const SupervisorFinished& from) {
    CopyFrom(from);
    return *this;
  }
  inline SupervisorFinished& operator=(SupervisorFinished&& from) noexcept {
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
  static const SupervisorFinished& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const SupervisorFinished* internal_default_instance() {
    return reinterpret_cast<const SupervisorFinished*>(
               &_SupervisorFinished_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(SupervisorFinished& a, SupervisorFinished& b) {
    a.Swap(&b);
  }
  inline void Swap(SupervisorFinished* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(SupervisorFinished* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline SupervisorFinished* New() const final {
    return CreateMaybeMessage<SupervisorFinished>(nullptr);
  }

  SupervisorFinished* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<SupervisorFinished>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const SupervisorFinished& from);
  void MergeFrom(const SupervisorFinished& from);
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
  void InternalSwap(SupervisorFinished* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.message.SupervisorFinished";
  }
  protected:
  explicit SupervisorFinished(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2fmessage_2fSupervisorFinished_2eproto);
    return ::descriptor_table_lm_2fmessage_2fSupervisorFinished_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kAddressFieldNumber = 1,
    kIdFieldNumber = 2,
  };
  // required .lm.message.Endpoint address = 1;
  bool has_address() const;
  private:
  bool _internal_has_address() const;
  public:
  void clear_address();
  const ::lm::message::Endpoint& address() const;
  ::lm::message::Endpoint* release_address();
  ::lm::message::Endpoint* mutable_address();
  void set_allocated_address(::lm::message::Endpoint* address);
  private:
  const ::lm::message::Endpoint& _internal_address() const;
  ::lm::message::Endpoint* _internal_mutable_address();
  public:
  void unsafe_arena_set_allocated_address(
      ::lm::message::Endpoint* address);
  ::lm::message::Endpoint* unsafe_arena_release_address();

  // required int64 id = 2;
  bool has_id() const;
  private:
  bool _internal_has_id() const;
  public:
  void clear_id();
  ::PROTOBUF_NAMESPACE_ID::int64 id() const;
  void set_id(::PROTOBUF_NAMESPACE_ID::int64 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::int64 _internal_id() const;
  void _internal_set_id(::PROTOBUF_NAMESPACE_ID::int64 value);
  public:

  // @@protoc_insertion_point(class_scope:lm.message.SupervisorFinished)
 private:
  class _Internal;

  // helper for ByteSizeLong()
  size_t RequiredFieldsByteSizeFallback() const;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  ::lm::message::Endpoint* address_;
  ::PROTOBUF_NAMESPACE_ID::int64 id_;
  friend struct ::TableStruct_lm_2fmessage_2fSupervisorFinished_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// SupervisorFinished

// required .lm.message.Endpoint address = 1;
inline bool SupervisorFinished::_internal_has_address() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  PROTOBUF_ASSUME(!value || address_ != nullptr);
  return value;
}
inline bool SupervisorFinished::has_address() const {
  return _internal_has_address();
}
inline const ::lm::message::Endpoint& SupervisorFinished::_internal_address() const {
  const ::lm::message::Endpoint* p = address_;
  return p != nullptr ? *p : *reinterpret_cast<const ::lm::message::Endpoint*>(
      &::lm::message::_Endpoint_default_instance_);
}
inline const ::lm::message::Endpoint& SupervisorFinished::address() const {
  // @@protoc_insertion_point(field_get:lm.message.SupervisorFinished.address)
  return _internal_address();
}
inline void SupervisorFinished::unsafe_arena_set_allocated_address(
    ::lm::message::Endpoint* address) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(address_);
  }
  address_ = address;
  if (address) {
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.message.SupervisorFinished.address)
}
inline ::lm::message::Endpoint* SupervisorFinished::release_address() {
  _has_bits_[0] &= ~0x00000001u;
  ::lm::message::Endpoint* temp = address_;
  address_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::lm::message::Endpoint* SupervisorFinished::unsafe_arena_release_address() {
  // @@protoc_insertion_point(field_release:lm.message.SupervisorFinished.address)
  _has_bits_[0] &= ~0x00000001u;
  ::lm::message::Endpoint* temp = address_;
  address_ = nullptr;
  return temp;
}
inline ::lm::message::Endpoint* SupervisorFinished::_internal_mutable_address() {
  _has_bits_[0] |= 0x00000001u;
  if (address_ == nullptr) {
    auto* p = CreateMaybeMessage<::lm::message::Endpoint>(GetArena());
    address_ = p;
  }
  return address_;
}
inline ::lm::message::Endpoint* SupervisorFinished::mutable_address() {
  // @@protoc_insertion_point(field_mutable:lm.message.SupervisorFinished.address)
  return _internal_mutable_address();
}
inline void SupervisorFinished::set_allocated_address(::lm::message::Endpoint* address) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(address_);
  }
  if (address) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(address)->GetArena();
    if (message_arena != submessage_arena) {
      address = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, address, submessage_arena);
    }
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  address_ = address;
  // @@protoc_insertion_point(field_set_allocated:lm.message.SupervisorFinished.address)
}

// required int64 id = 2;
inline bool SupervisorFinished::_internal_has_id() const {
  bool value = (_has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool SupervisorFinished::has_id() const {
  return _internal_has_id();
}
inline void SupervisorFinished::clear_id() {
  id_ = PROTOBUF_LONGLONG(0);
  _has_bits_[0] &= ~0x00000002u;
}
inline ::PROTOBUF_NAMESPACE_ID::int64 SupervisorFinished::_internal_id() const {
  return id_;
}
inline ::PROTOBUF_NAMESPACE_ID::int64 SupervisorFinished::id() const {
  // @@protoc_insertion_point(field_get:lm.message.SupervisorFinished.id)
  return _internal_id();
}
inline void SupervisorFinished::_internal_set_id(::PROTOBUF_NAMESPACE_ID::int64 value) {
  _has_bits_[0] |= 0x00000002u;
  id_ = value;
}
inline void SupervisorFinished::set_id(::PROTOBUF_NAMESPACE_ID::int64 value) {
  _internal_set_id(value);
  // @@protoc_insertion_point(field_set:lm.message.SupervisorFinished.id)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace message
}  // namespace lm

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2fmessage_2fSupervisorFinished_2eproto
