// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/types/Tilings.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2ftypes_2fTilings_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2ftypes_2fTilings_2eproto

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
#define PROTOBUF_INTERNAL_EXPORT_lm_2ftypes_2fTilings_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2ftypes_2fTilings_2eproto {
  static const ::PROTOBUF_NAMESPACE_ID::internal::ParseTableField entries[]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::AuxiliaryParseTableField aux[]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::ParseTable schema[2]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::FieldMetadata field_metadata[];
  static const ::PROTOBUF_NAMESPACE_ID::internal::SerializationTable serialization_table[];
  static const ::PROTOBUF_NAMESPACE_ID::uint32 offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2ftypes_2fTilings_2eproto;
namespace lm {
namespace types {
class Tiling;
class TilingDefaultTypeInternal;
extern TilingDefaultTypeInternal _Tiling_default_instance_;
class Tilings;
class TilingsDefaultTypeInternal;
extern TilingsDefaultTypeInternal _Tilings_default_instance_;
}  // namespace types
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::types::Tiling* Arena::CreateMaybeMessage<::lm::types::Tiling>(Arena*);
template<> ::lm::types::Tilings* Arena::CreateMaybeMessage<::lm::types::Tilings>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace types {

// ===================================================================

class Tiling PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.types.Tiling) */ {
 public:
  inline Tiling() : Tiling(nullptr) {}
  virtual ~Tiling();

  Tiling(const Tiling& from);
  Tiling(Tiling&& from) noexcept
    : Tiling() {
    *this = ::std::move(from);
  }

  inline Tiling& operator=(const Tiling& from) {
    CopyFrom(from);
    return *this;
  }
  inline Tiling& operator=(Tiling&& from) noexcept {
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
  static const Tiling& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const Tiling* internal_default_instance() {
    return reinterpret_cast<const Tiling*>(
               &_Tiling_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(Tiling& a, Tiling& b) {
    a.Swap(&b);
  }
  inline void Swap(Tiling* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(Tiling* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline Tiling* New() const final {
    return CreateMaybeMessage<Tiling>(nullptr);
  }

  Tiling* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<Tiling>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const Tiling& from);
  void MergeFrom(const Tiling& from);
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
  void InternalSwap(Tiling* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.types.Tiling";
  }
  protected:
  explicit Tiling(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2ftypes_2fTilings_2eproto);
    return ::descriptor_table_lm_2ftypes_2fTilings_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kOrderParameterIndicesFieldNumber = 1,
    kEdgesFieldNumber = 2,
  };
  // required .robertslab.pbuf.NDArray order_parameter_indices = 1;
  bool has_order_parameter_indices() const;
  private:
  bool _internal_has_order_parameter_indices() const;
  public:
  void clear_order_parameter_indices();
  const ::robertslab::pbuf::NDArray& order_parameter_indices() const;
  ::robertslab::pbuf::NDArray* release_order_parameter_indices();
  ::robertslab::pbuf::NDArray* mutable_order_parameter_indices();
  void set_allocated_order_parameter_indices(::robertslab::pbuf::NDArray* order_parameter_indices);
  private:
  const ::robertslab::pbuf::NDArray& _internal_order_parameter_indices() const;
  ::robertslab::pbuf::NDArray* _internal_mutable_order_parameter_indices();
  public:
  void unsafe_arena_set_allocated_order_parameter_indices(
      ::robertslab::pbuf::NDArray* order_parameter_indices);
  ::robertslab::pbuf::NDArray* unsafe_arena_release_order_parameter_indices();

  // required .robertslab.pbuf.NDArray edges = 2;
  bool has_edges() const;
  private:
  bool _internal_has_edges() const;
  public:
  void clear_edges();
  const ::robertslab::pbuf::NDArray& edges() const;
  ::robertslab::pbuf::NDArray* release_edges();
  ::robertslab::pbuf::NDArray* mutable_edges();
  void set_allocated_edges(::robertslab::pbuf::NDArray* edges);
  private:
  const ::robertslab::pbuf::NDArray& _internal_edges() const;
  ::robertslab::pbuf::NDArray* _internal_mutable_edges();
  public:
  void unsafe_arena_set_allocated_edges(
      ::robertslab::pbuf::NDArray* edges);
  ::robertslab::pbuf::NDArray* unsafe_arena_release_edges();

  // @@protoc_insertion_point(class_scope:lm.types.Tiling)
 private:
  class _Internal;

  // helper for ByteSizeLong()
  size_t RequiredFieldsByteSizeFallback() const;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  ::robertslab::pbuf::NDArray* order_parameter_indices_;
  ::robertslab::pbuf::NDArray* edges_;
  friend struct ::TableStruct_lm_2ftypes_2fTilings_2eproto;
};
// -------------------------------------------------------------------

class Tilings PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.types.Tilings) */ {
 public:
  inline Tilings() : Tilings(nullptr) {}
  virtual ~Tilings();

  Tilings(const Tilings& from);
  Tilings(Tilings&& from) noexcept
    : Tilings() {
    *this = ::std::move(from);
  }

  inline Tilings& operator=(const Tilings& from) {
    CopyFrom(from);
    return *this;
  }
  inline Tilings& operator=(Tilings&& from) noexcept {
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
  static const Tilings& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const Tilings* internal_default_instance() {
    return reinterpret_cast<const Tilings*>(
               &_Tilings_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    1;

  friend void swap(Tilings& a, Tilings& b) {
    a.Swap(&b);
  }
  inline void Swap(Tilings* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(Tilings* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline Tilings* New() const final {
    return CreateMaybeMessage<Tilings>(nullptr);
  }

  Tilings* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<Tilings>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const Tilings& from);
  void MergeFrom(const Tilings& from);
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
  void InternalSwap(Tilings* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.types.Tilings";
  }
  protected:
  explicit Tilings(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2ftypes_2fTilings_2eproto);
    return ::descriptor_table_lm_2ftypes_2fTilings_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kTilingFieldNumber = 1,
  };
  // repeated .lm.types.Tiling tiling = 1;
  int tiling_size() const;
  private:
  int _internal_tiling_size() const;
  public:
  void clear_tiling();
  ::lm::types::Tiling* mutable_tiling(int index);
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::types::Tiling >*
      mutable_tiling();
  private:
  const ::lm::types::Tiling& _internal_tiling(int index) const;
  ::lm::types::Tiling* _internal_add_tiling();
  public:
  const ::lm::types::Tiling& tiling(int index) const;
  ::lm::types::Tiling* add_tiling();
  const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::types::Tiling >&
      tiling() const;

  // @@protoc_insertion_point(class_scope:lm.types.Tilings)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::types::Tiling > tiling_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  friend struct ::TableStruct_lm_2ftypes_2fTilings_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// Tiling

// required .robertslab.pbuf.NDArray order_parameter_indices = 1;
inline bool Tiling::_internal_has_order_parameter_indices() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  PROTOBUF_ASSUME(!value || order_parameter_indices_ != nullptr);
  return value;
}
inline bool Tiling::has_order_parameter_indices() const {
  return _internal_has_order_parameter_indices();
}
inline const ::robertslab::pbuf::NDArray& Tiling::_internal_order_parameter_indices() const {
  const ::robertslab::pbuf::NDArray* p = order_parameter_indices_;
  return p != nullptr ? *p : *reinterpret_cast<const ::robertslab::pbuf::NDArray*>(
      &::robertslab::pbuf::_NDArray_default_instance_);
}
inline const ::robertslab::pbuf::NDArray& Tiling::order_parameter_indices() const {
  // @@protoc_insertion_point(field_get:lm.types.Tiling.order_parameter_indices)
  return _internal_order_parameter_indices();
}
inline void Tiling::unsafe_arena_set_allocated_order_parameter_indices(
    ::robertslab::pbuf::NDArray* order_parameter_indices) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(order_parameter_indices_);
  }
  order_parameter_indices_ = order_parameter_indices;
  if (order_parameter_indices) {
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.types.Tiling.order_parameter_indices)
}
inline ::robertslab::pbuf::NDArray* Tiling::release_order_parameter_indices() {
  _has_bits_[0] &= ~0x00000001u;
  ::robertslab::pbuf::NDArray* temp = order_parameter_indices_;
  order_parameter_indices_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::robertslab::pbuf::NDArray* Tiling::unsafe_arena_release_order_parameter_indices() {
  // @@protoc_insertion_point(field_release:lm.types.Tiling.order_parameter_indices)
  _has_bits_[0] &= ~0x00000001u;
  ::robertslab::pbuf::NDArray* temp = order_parameter_indices_;
  order_parameter_indices_ = nullptr;
  return temp;
}
inline ::robertslab::pbuf::NDArray* Tiling::_internal_mutable_order_parameter_indices() {
  _has_bits_[0] |= 0x00000001u;
  if (order_parameter_indices_ == nullptr) {
    auto* p = CreateMaybeMessage<::robertslab::pbuf::NDArray>(GetArena());
    order_parameter_indices_ = p;
  }
  return order_parameter_indices_;
}
inline ::robertslab::pbuf::NDArray* Tiling::mutable_order_parameter_indices() {
  // @@protoc_insertion_point(field_mutable:lm.types.Tiling.order_parameter_indices)
  return _internal_mutable_order_parameter_indices();
}
inline void Tiling::set_allocated_order_parameter_indices(::robertslab::pbuf::NDArray* order_parameter_indices) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(order_parameter_indices_);
  }
  if (order_parameter_indices) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(order_parameter_indices)->GetArena();
    if (message_arena != submessage_arena) {
      order_parameter_indices = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, order_parameter_indices, submessage_arena);
    }
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  order_parameter_indices_ = order_parameter_indices;
  // @@protoc_insertion_point(field_set_allocated:lm.types.Tiling.order_parameter_indices)
}

// required .robertslab.pbuf.NDArray edges = 2;
inline bool Tiling::_internal_has_edges() const {
  bool value = (_has_bits_[0] & 0x00000002u) != 0;
  PROTOBUF_ASSUME(!value || edges_ != nullptr);
  return value;
}
inline bool Tiling::has_edges() const {
  return _internal_has_edges();
}
inline const ::robertslab::pbuf::NDArray& Tiling::_internal_edges() const {
  const ::robertslab::pbuf::NDArray* p = edges_;
  return p != nullptr ? *p : *reinterpret_cast<const ::robertslab::pbuf::NDArray*>(
      &::robertslab::pbuf::_NDArray_default_instance_);
}
inline const ::robertslab::pbuf::NDArray& Tiling::edges() const {
  // @@protoc_insertion_point(field_get:lm.types.Tiling.edges)
  return _internal_edges();
}
inline void Tiling::unsafe_arena_set_allocated_edges(
    ::robertslab::pbuf::NDArray* edges) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(edges_);
  }
  edges_ = edges;
  if (edges) {
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.types.Tiling.edges)
}
inline ::robertslab::pbuf::NDArray* Tiling::release_edges() {
  _has_bits_[0] &= ~0x00000002u;
  ::robertslab::pbuf::NDArray* temp = edges_;
  edges_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::robertslab::pbuf::NDArray* Tiling::unsafe_arena_release_edges() {
  // @@protoc_insertion_point(field_release:lm.types.Tiling.edges)
  _has_bits_[0] &= ~0x00000002u;
  ::robertslab::pbuf::NDArray* temp = edges_;
  edges_ = nullptr;
  return temp;
}
inline ::robertslab::pbuf::NDArray* Tiling::_internal_mutable_edges() {
  _has_bits_[0] |= 0x00000002u;
  if (edges_ == nullptr) {
    auto* p = CreateMaybeMessage<::robertslab::pbuf::NDArray>(GetArena());
    edges_ = p;
  }
  return edges_;
}
inline ::robertslab::pbuf::NDArray* Tiling::mutable_edges() {
  // @@protoc_insertion_point(field_mutable:lm.types.Tiling.edges)
  return _internal_mutable_edges();
}
inline void Tiling::set_allocated_edges(::robertslab::pbuf::NDArray* edges) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(edges_);
  }
  if (edges) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(edges)->GetArena();
    if (message_arena != submessage_arena) {
      edges = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, edges, submessage_arena);
    }
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  edges_ = edges;
  // @@protoc_insertion_point(field_set_allocated:lm.types.Tiling.edges)
}

// -------------------------------------------------------------------

// Tilings

// repeated .lm.types.Tiling tiling = 1;
inline int Tilings::_internal_tiling_size() const {
  return tiling_.size();
}
inline int Tilings::tiling_size() const {
  return _internal_tiling_size();
}
inline void Tilings::clear_tiling() {
  tiling_.Clear();
}
inline ::lm::types::Tiling* Tilings::mutable_tiling(int index) {
  // @@protoc_insertion_point(field_mutable:lm.types.Tilings.tiling)
  return tiling_.Mutable(index);
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::types::Tiling >*
Tilings::mutable_tiling() {
  // @@protoc_insertion_point(field_mutable_list:lm.types.Tilings.tiling)
  return &tiling_;
}
inline const ::lm::types::Tiling& Tilings::_internal_tiling(int index) const {
  return tiling_.Get(index);
}
inline const ::lm::types::Tiling& Tilings::tiling(int index) const {
  // @@protoc_insertion_point(field_get:lm.types.Tilings.tiling)
  return _internal_tiling(index);
}
inline ::lm::types::Tiling* Tilings::_internal_add_tiling() {
  return tiling_.Add();
}
inline ::lm::types::Tiling* Tilings::add_tiling() {
  // @@protoc_insertion_point(field_add:lm.types.Tilings.tiling)
  return _internal_add_tiling();
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::types::Tiling >&
Tilings::tiling() const {
  // @@protoc_insertion_point(field_list:lm.types.Tilings.tiling)
  return tiling_;
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__
// -------------------------------------------------------------------


// @@protoc_insertion_point(namespace_scope)

}  // namespace types
}  // namespace lm

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2ftypes_2fTilings_2eproto
