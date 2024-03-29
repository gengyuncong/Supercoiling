// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/io/SpeciesCounts.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fSpeciesCounts_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fSpeciesCounts_2eproto

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
#define PROTOBUF_INTERNAL_EXPORT_lm_2fio_2fSpeciesCounts_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2fio_2fSpeciesCounts_2eproto {
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
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fio_2fSpeciesCounts_2eproto;
namespace lm {
namespace io {
class SpeciesCounts;
class SpeciesCountsDefaultTypeInternal;
extern SpeciesCountsDefaultTypeInternal _SpeciesCounts_default_instance_;
}  // namespace io
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::io::SpeciesCounts* Arena::CreateMaybeMessage<::lm::io::SpeciesCounts>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace io {

// ===================================================================

class SpeciesCounts PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.io.SpeciesCounts) */ {
 public:
  inline SpeciesCounts() : SpeciesCounts(nullptr) {}
  virtual ~SpeciesCounts();

  SpeciesCounts(const SpeciesCounts& from);
  SpeciesCounts(SpeciesCounts&& from) noexcept
    : SpeciesCounts() {
    *this = ::std::move(from);
  }

  inline SpeciesCounts& operator=(const SpeciesCounts& from) {
    CopyFrom(from);
    return *this;
  }
  inline SpeciesCounts& operator=(SpeciesCounts&& from) noexcept {
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
  static const SpeciesCounts& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const SpeciesCounts* internal_default_instance() {
    return reinterpret_cast<const SpeciesCounts*>(
               &_SpeciesCounts_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(SpeciesCounts& a, SpeciesCounts& b) {
    a.Swap(&b);
  }
  inline void Swap(SpeciesCounts* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(SpeciesCounts* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline SpeciesCounts* New() const final {
    return CreateMaybeMessage<SpeciesCounts>(nullptr);
  }

  SpeciesCounts* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<SpeciesCounts>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const SpeciesCounts& from);
  void MergeFrom(const SpeciesCounts& from);
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
  void InternalSwap(SpeciesCounts* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.io.SpeciesCounts";
  }
  protected:
  explicit SpeciesCounts(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2fio_2fSpeciesCounts_2eproto);
    return ::descriptor_table_lm_2fio_2fSpeciesCounts_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kSpeciesCountFieldNumber = 4,
    kTimeFieldNumber = 5,
    kSpeciesCountPreviousFieldNumber = 6,
    kTrajectoryIdFieldNumber = 1,
    kNumberSpeciesFieldNumber = 2,
    kNumberEntriesFieldNumber = 3,
  };
  // repeated int32 species_count = 4 [packed = true];
  int species_count_size() const;
  private:
  int _internal_species_count_size() const;
  public:
  void clear_species_count();
  private:
  ::PROTOBUF_NAMESPACE_ID::int32 _internal_species_count(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
      _internal_species_count() const;
  void _internal_add_species_count(::PROTOBUF_NAMESPACE_ID::int32 value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
      _internal_mutable_species_count();
  public:
  ::PROTOBUF_NAMESPACE_ID::int32 species_count(int index) const;
  void set_species_count(int index, ::PROTOBUF_NAMESPACE_ID::int32 value);
  void add_species_count(::PROTOBUF_NAMESPACE_ID::int32 value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
      species_count() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
      mutable_species_count();

  // repeated double time = 5 [packed = true];
  int time_size() const;
  private:
  int _internal_time_size() const;
  public:
  void clear_time();
  private:
  double _internal_time(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      _internal_time() const;
  void _internal_add_time(double value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      _internal_mutable_time();
  public:
  double time(int index) const;
  void set_time(int index, double value);
  void add_time(double value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      time() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      mutable_time();

  // repeated int32 species_count_previous = 6 [packed = true];
  int species_count_previous_size() const;
  private:
  int _internal_species_count_previous_size() const;
  public:
  void clear_species_count_previous();
  private:
  ::PROTOBUF_NAMESPACE_ID::int32 _internal_species_count_previous(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
      _internal_species_count_previous() const;
  void _internal_add_species_count_previous(::PROTOBUF_NAMESPACE_ID::int32 value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
      _internal_mutable_species_count_previous();
  public:
  ::PROTOBUF_NAMESPACE_ID::int32 species_count_previous(int index) const;
  void set_species_count_previous(int index, ::PROTOBUF_NAMESPACE_ID::int32 value);
  void add_species_count_previous(::PROTOBUF_NAMESPACE_ID::int32 value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
      species_count_previous() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
      mutable_species_count_previous();

  // required uint64 trajectory_id = 1;
  bool has_trajectory_id() const;
  private:
  bool _internal_has_trajectory_id() const;
  public:
  void clear_trajectory_id();
  ::PROTOBUF_NAMESPACE_ID::uint64 trajectory_id() const;
  void set_trajectory_id(::PROTOBUF_NAMESPACE_ID::uint64 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::uint64 _internal_trajectory_id() const;
  void _internal_set_trajectory_id(::PROTOBUF_NAMESPACE_ID::uint64 value);
  public:

  // required int32 number_species = 2;
  bool has_number_species() const;
  private:
  bool _internal_has_number_species() const;
  public:
  void clear_number_species();
  ::PROTOBUF_NAMESPACE_ID::int32 number_species() const;
  void set_number_species(::PROTOBUF_NAMESPACE_ID::int32 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::int32 _internal_number_species() const;
  void _internal_set_number_species(::PROTOBUF_NAMESPACE_ID::int32 value);
  public:

  // required int32 number_entries = 3;
  bool has_number_entries() const;
  private:
  bool _internal_has_number_entries() const;
  public:
  void clear_number_entries();
  ::PROTOBUF_NAMESPACE_ID::int32 number_entries() const;
  void set_number_entries(::PROTOBUF_NAMESPACE_ID::int32 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::int32 _internal_number_entries() const;
  void _internal_set_number_entries(::PROTOBUF_NAMESPACE_ID::int32 value);
  public:

  // @@protoc_insertion_point(class_scope:lm.io.SpeciesCounts)
 private:
  class _Internal;

  // helper for ByteSizeLong()
  size_t RequiredFieldsByteSizeFallback() const;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 > species_count_;
  mutable std::atomic<int> _species_count_cached_byte_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double > time_;
  mutable std::atomic<int> _time_cached_byte_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 > species_count_previous_;
  mutable std::atomic<int> _species_count_previous_cached_byte_size_;
  ::PROTOBUF_NAMESPACE_ID::uint64 trajectory_id_;
  ::PROTOBUF_NAMESPACE_ID::int32 number_species_;
  ::PROTOBUF_NAMESPACE_ID::int32 number_entries_;
  friend struct ::TableStruct_lm_2fio_2fSpeciesCounts_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// SpeciesCounts

// required uint64 trajectory_id = 1;
inline bool SpeciesCounts::_internal_has_trajectory_id() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool SpeciesCounts::has_trajectory_id() const {
  return _internal_has_trajectory_id();
}
inline void SpeciesCounts::clear_trajectory_id() {
  trajectory_id_ = PROTOBUF_ULONGLONG(0);
  _has_bits_[0] &= ~0x00000001u;
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 SpeciesCounts::_internal_trajectory_id() const {
  return trajectory_id_;
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 SpeciesCounts::trajectory_id() const {
  // @@protoc_insertion_point(field_get:lm.io.SpeciesCounts.trajectory_id)
  return _internal_trajectory_id();
}
inline void SpeciesCounts::_internal_set_trajectory_id(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  _has_bits_[0] |= 0x00000001u;
  trajectory_id_ = value;
}
inline void SpeciesCounts::set_trajectory_id(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  _internal_set_trajectory_id(value);
  // @@protoc_insertion_point(field_set:lm.io.SpeciesCounts.trajectory_id)
}

// required int32 number_species = 2;
inline bool SpeciesCounts::_internal_has_number_species() const {
  bool value = (_has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool SpeciesCounts::has_number_species() const {
  return _internal_has_number_species();
}
inline void SpeciesCounts::clear_number_species() {
  number_species_ = 0;
  _has_bits_[0] &= ~0x00000002u;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 SpeciesCounts::_internal_number_species() const {
  return number_species_;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 SpeciesCounts::number_species() const {
  // @@protoc_insertion_point(field_get:lm.io.SpeciesCounts.number_species)
  return _internal_number_species();
}
inline void SpeciesCounts::_internal_set_number_species(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _has_bits_[0] |= 0x00000002u;
  number_species_ = value;
}
inline void SpeciesCounts::set_number_species(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_set_number_species(value);
  // @@protoc_insertion_point(field_set:lm.io.SpeciesCounts.number_species)
}

// required int32 number_entries = 3;
inline bool SpeciesCounts::_internal_has_number_entries() const {
  bool value = (_has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool SpeciesCounts::has_number_entries() const {
  return _internal_has_number_entries();
}
inline void SpeciesCounts::clear_number_entries() {
  number_entries_ = 0;
  _has_bits_[0] &= ~0x00000004u;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 SpeciesCounts::_internal_number_entries() const {
  return number_entries_;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 SpeciesCounts::number_entries() const {
  // @@protoc_insertion_point(field_get:lm.io.SpeciesCounts.number_entries)
  return _internal_number_entries();
}
inline void SpeciesCounts::_internal_set_number_entries(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _has_bits_[0] |= 0x00000004u;
  number_entries_ = value;
}
inline void SpeciesCounts::set_number_entries(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_set_number_entries(value);
  // @@protoc_insertion_point(field_set:lm.io.SpeciesCounts.number_entries)
}

// repeated int32 species_count = 4 [packed = true];
inline int SpeciesCounts::_internal_species_count_size() const {
  return species_count_.size();
}
inline int SpeciesCounts::species_count_size() const {
  return _internal_species_count_size();
}
inline void SpeciesCounts::clear_species_count() {
  species_count_.Clear();
}
inline ::PROTOBUF_NAMESPACE_ID::int32 SpeciesCounts::_internal_species_count(int index) const {
  return species_count_.Get(index);
}
inline ::PROTOBUF_NAMESPACE_ID::int32 SpeciesCounts::species_count(int index) const {
  // @@protoc_insertion_point(field_get:lm.io.SpeciesCounts.species_count)
  return _internal_species_count(index);
}
inline void SpeciesCounts::set_species_count(int index, ::PROTOBUF_NAMESPACE_ID::int32 value) {
  species_count_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.io.SpeciesCounts.species_count)
}
inline void SpeciesCounts::_internal_add_species_count(::PROTOBUF_NAMESPACE_ID::int32 value) {
  species_count_.Add(value);
}
inline void SpeciesCounts::add_species_count(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_add_species_count(value);
  // @@protoc_insertion_point(field_add:lm.io.SpeciesCounts.species_count)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
SpeciesCounts::_internal_species_count() const {
  return species_count_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
SpeciesCounts::species_count() const {
  // @@protoc_insertion_point(field_list:lm.io.SpeciesCounts.species_count)
  return _internal_species_count();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
SpeciesCounts::_internal_mutable_species_count() {
  return &species_count_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
SpeciesCounts::mutable_species_count() {
  // @@protoc_insertion_point(field_mutable_list:lm.io.SpeciesCounts.species_count)
  return _internal_mutable_species_count();
}

// repeated int32 species_count_previous = 6 [packed = true];
inline int SpeciesCounts::_internal_species_count_previous_size() const {
  return species_count_previous_.size();
}
inline int SpeciesCounts::species_count_previous_size() const {
  return _internal_species_count_previous_size();
}
inline void SpeciesCounts::clear_species_count_previous() {
  species_count_previous_.Clear();
}
inline ::PROTOBUF_NAMESPACE_ID::int32 SpeciesCounts::_internal_species_count_previous(int index) const {
  return species_count_previous_.Get(index);
}
inline ::PROTOBUF_NAMESPACE_ID::int32 SpeciesCounts::species_count_previous(int index) const {
  // @@protoc_insertion_point(field_get:lm.io.SpeciesCounts.species_count_previous)
  return _internal_species_count_previous(index);
}
inline void SpeciesCounts::set_species_count_previous(int index, ::PROTOBUF_NAMESPACE_ID::int32 value) {
  species_count_previous_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.io.SpeciesCounts.species_count_previous)
}
inline void SpeciesCounts::_internal_add_species_count_previous(::PROTOBUF_NAMESPACE_ID::int32 value) {
  species_count_previous_.Add(value);
}
inline void SpeciesCounts::add_species_count_previous(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_add_species_count_previous(value);
  // @@protoc_insertion_point(field_add:lm.io.SpeciesCounts.species_count_previous)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
SpeciesCounts::_internal_species_count_previous() const {
  return species_count_previous_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >&
SpeciesCounts::species_count_previous() const {
  // @@protoc_insertion_point(field_list:lm.io.SpeciesCounts.species_count_previous)
  return _internal_species_count_previous();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
SpeciesCounts::_internal_mutable_species_count_previous() {
  return &species_count_previous_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::int32 >*
SpeciesCounts::mutable_species_count_previous() {
  // @@protoc_insertion_point(field_mutable_list:lm.io.SpeciesCounts.species_count_previous)
  return _internal_mutable_species_count_previous();
}

// repeated double time = 5 [packed = true];
inline int SpeciesCounts::_internal_time_size() const {
  return time_.size();
}
inline int SpeciesCounts::time_size() const {
  return _internal_time_size();
}
inline void SpeciesCounts::clear_time() {
  time_.Clear();
}
inline double SpeciesCounts::_internal_time(int index) const {
  return time_.Get(index);
}
inline double SpeciesCounts::time(int index) const {
  // @@protoc_insertion_point(field_get:lm.io.SpeciesCounts.time)
  return _internal_time(index);
}
inline void SpeciesCounts::set_time(int index, double value) {
  time_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.io.SpeciesCounts.time)
}
inline void SpeciesCounts::_internal_add_time(double value) {
  time_.Add(value);
}
inline void SpeciesCounts::add_time(double value) {
  _internal_add_time(value);
  // @@protoc_insertion_point(field_add:lm.io.SpeciesCounts.time)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
SpeciesCounts::_internal_time() const {
  return time_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
SpeciesCounts::time() const {
  // @@protoc_insertion_point(field_list:lm.io.SpeciesCounts.time)
  return _internal_time();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
SpeciesCounts::_internal_mutable_time() {
  return &time_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
SpeciesCounts::mutable_time() {
  // @@protoc_insertion_point(field_mutable_list:lm.io.SpeciesCounts.time)
  return _internal_mutable_time();
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace io
}  // namespace lm

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fSpeciesCounts_2eproto
