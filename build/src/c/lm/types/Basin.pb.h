// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/types/Basin.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2ftypes_2fBasin_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2ftypes_2fBasin_2eproto

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
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_lm_2ftypes_2fBasin_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2ftypes_2fBasin_2eproto {
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
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2ftypes_2fBasin_2eproto;
namespace lm {
namespace types {
class Basin;
class BasinDefaultTypeInternal;
extern BasinDefaultTypeInternal _Basin_default_instance_;
}  // namespace types
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::types::Basin* Arena::CreateMaybeMessage<::lm::types::Basin>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace types {

// ===================================================================

class Basin PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.types.Basin) */ {
 public:
  inline Basin() : Basin(nullptr) {}
  virtual ~Basin();

  Basin(const Basin& from);
  Basin(Basin&& from) noexcept
    : Basin() {
    *this = ::std::move(from);
  }

  inline Basin& operator=(const Basin& from) {
    CopyFrom(from);
    return *this;
  }
  inline Basin& operator=(Basin&& from) noexcept {
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
  static const Basin& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const Basin* internal_default_instance() {
    return reinterpret_cast<const Basin*>(
               &_Basin_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(Basin& a, Basin& b) {
    a.Swap(&b);
  }
  inline void Swap(Basin* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(Basin* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline Basin* New() const final {
    return CreateMaybeMessage<Basin>(nullptr);
  }

  Basin* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<Basin>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const Basin& from);
  void MergeFrom(const Basin& from);
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
  void InternalSwap(Basin* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.types.Basin";
  }
  protected:
  explicit Basin(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2ftypes_2fBasin_2eproto);
    return ::descriptor_table_lm_2ftypes_2fBasin_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kSpeciesCountsFieldNumber = 1,
  };
  // repeated int32 species_counts = 1 [packed = true];
  int species_counts_size() const;
  private:
  int _internal_species_counts_size() const;
  public:
  void clear_species_counts();
  private:
  ::PROTOBUF_NAMESPACE_ID::int32 _internal_species_counts(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
      _internal_species_counts() const;
  void _internal_add_species_counts(::PROTOBUF_NAMESPACE_ID::int32 value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
      _internal_mutable_species_counts();
  public:
  ::PROTOBUF_NAMESPACE_ID::int32 species_counts(int index) const;
  void set_species_counts(int index, ::PROTOBUF_NAMESPACE_ID::int32 value);
  void add_species_counts(::PROTOBUF_NAMESPACE_ID::int32 value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
      species_counts() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
      mutable_species_counts();

  // @@protoc_insertion_point(class_scope:lm.types.Basin)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 > species_counts_;
  mutable std::atomic<int> _species_counts_cached_byte_size_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  friend struct ::TableStruct_lm_2ftypes_2fBasin_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// Basin

// repeated int32 species_counts = 1 [packed = true];
inline int Basin::_internal_species_counts_size() const {
  return species_counts_.size();
}
inline int Basin::species_counts_size() const {
  return _internal_species_counts_size();
}
inline void Basin::clear_species_counts() {
  species_counts_.Clear();
}
inline ::PROTOBUF_NAMESPACE_ID::int32 Basin::_internal_species_counts(int index) const {
  return species_counts_.Get(index);
}
inline ::PROTOBUF_NAMESPACE_ID::int32 Basin::species_counts(int index) const {
  // @@protoc_insertion_point(field_get:lm.types.Basin.species_counts)
  return _internal_species_counts(index);
}
inline void Basin::set_species_counts(int index, ::PROTOBUF_NAMESPACE_ID::int32 value) {
  species_counts_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.types.Basin.species_counts)
}
inline void Basin::_internal_add_species_counts(::PROTOBUF_NAMESPACE_ID::int32 value) {
  species_counts_.Add(value);
}
inline void Basin::add_species_counts(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_add_species_counts(value);
  // @@protoc_insertion_point(field_add:lm.types.Basin.species_counts)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
Basin::_internal_species_counts() const {
  return species_counts_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
Basin::species_counts() const {
  // @@protoc_insertion_point(field_list:lm.types.Basin.species_counts)
  return _internal_species_counts();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
Basin::_internal_mutable_species_counts() {
  return &species_counts_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
Basin::mutable_species_counts() {
  // @@protoc_insertion_point(field_mutable_list:lm.types.Basin.species_counts)
  return _internal_mutable_species_counts();
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace types
}  // namespace lm

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2ftypes_2fBasin_2eproto
