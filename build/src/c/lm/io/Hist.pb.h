// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/io/Hist.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fHist_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fHist_2eproto

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
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_lm_2fio_2fHist_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2fio_2fHist_2eproto {
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
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fio_2fHist_2eproto;
namespace lm {
namespace io {
class Hist;
class HistDefaultTypeInternal;
extern HistDefaultTypeInternal _Hist_default_instance_;
}  // namespace io
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::io::Hist* Arena::CreateMaybeMessage<::lm::io::Hist>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace io {

enum Hist_EdgeType : int {
  Hist_EdgeType_CLOSED = 0,
  Hist_EdgeType_OPEN = 1,
  Hist_EdgeType_IGNORED = 2
};
bool Hist_EdgeType_IsValid(int value);
constexpr Hist_EdgeType Hist_EdgeType_EdgeType_MIN = Hist_EdgeType_CLOSED;
constexpr Hist_EdgeType Hist_EdgeType_EdgeType_MAX = Hist_EdgeType_IGNORED;
constexpr int Hist_EdgeType_EdgeType_ARRAYSIZE = Hist_EdgeType_EdgeType_MAX + 1;

const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* Hist_EdgeType_descriptor();
template<typename T>
inline const std::string& Hist_EdgeType_Name(T enum_t_value) {
  static_assert(::std::is_same<T, Hist_EdgeType>::value ||
    ::std::is_integral<T>::value,
    "Incorrect type passed to function Hist_EdgeType_Name.");
  return ::PROTOBUF_NAMESPACE_ID::internal::NameOfEnum(
    Hist_EdgeType_descriptor(), enum_t_value);
}
inline bool Hist_EdgeType_Parse(
    ::PROTOBUF_NAMESPACE_ID::ConstStringParam name, Hist_EdgeType* value) {
  return ::PROTOBUF_NAMESPACE_ID::internal::ParseNamedEnum<Hist_EdgeType>(
    Hist_EdgeType_descriptor(), name, value);
}
// ===================================================================

class Hist PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.io.Hist) */ {
 public:
  inline Hist() : Hist(nullptr) {}
  virtual ~Hist();

  Hist(const Hist& from);
  Hist(Hist&& from) noexcept
    : Hist() {
    *this = ::std::move(from);
  }

  inline Hist& operator=(const Hist& from) {
    CopyFrom(from);
    return *this;
  }
  inline Hist& operator=(Hist&& from) noexcept {
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
  static const Hist& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const Hist* internal_default_instance() {
    return reinterpret_cast<const Hist*>(
               &_Hist_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(Hist& a, Hist& b) {
    a.Swap(&b);
  }
  inline void Swap(Hist* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(Hist* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline Hist* New() const final {
    return CreateMaybeMessage<Hist>(nullptr);
  }

  Hist* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<Hist>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const Hist& from);
  void MergeFrom(const Hist& from);
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
  void InternalSwap(Hist* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.io.Hist";
  }
  protected:
  explicit Hist(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2fio_2fHist_2eproto);
    return ::descriptor_table_lm_2fio_2fHist_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  typedef Hist_EdgeType EdgeType;
  static constexpr EdgeType CLOSED =
    Hist_EdgeType_CLOSED;
  static constexpr EdgeType OPEN =
    Hist_EdgeType_OPEN;
  static constexpr EdgeType IGNORED =
    Hist_EdgeType_IGNORED;
  static inline bool EdgeType_IsValid(int value) {
    return Hist_EdgeType_IsValid(value);
  }
  static constexpr EdgeType EdgeType_MIN =
    Hist_EdgeType_EdgeType_MIN;
  static constexpr EdgeType EdgeType_MAX =
    Hist_EdgeType_EdgeType_MAX;
  static constexpr int EdgeType_ARRAYSIZE =
    Hist_EdgeType_EdgeType_ARRAYSIZE;
  static inline const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor*
  EdgeType_descriptor() {
    return Hist_EdgeType_descriptor();
  }
  template<typename T>
  static inline const std::string& EdgeType_Name(T enum_t_value) {
    static_assert(::std::is_same<T, EdgeType>::value ||
      ::std::is_integral<T>::value,
      "Incorrect type passed to function EdgeType_Name.");
    return Hist_EdgeType_Name(enum_t_value);
  }
  static inline bool EdgeType_Parse(::PROTOBUF_NAMESPACE_ID::ConstStringParam name,
      EdgeType* value) {
    return Hist_EdgeType_Parse(name, value);
  }

  // accessors -------------------------------------------------------

  enum : int {
    kDimsFieldNumber = 2,
    kEdgesFieldNumber = 100,
    kValsFieldNumber = 101,
    kEdgeTypesFieldNumber = 200,
    kRankFieldNumber = 1,
  };
  // repeated uint64 dims = 2 [packed = true];
  int dims_size() const;
  private:
  int _internal_dims_size() const;
  public:
  void clear_dims();
  private:
  ::PROTOBUF_NAMESPACE_ID::uint64 _internal_dims(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint64 >&
      _internal_dims() const;
  void _internal_add_dims(::PROTOBUF_NAMESPACE_ID::uint64 value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint64 >*
      _internal_mutable_dims();
  public:
  ::PROTOBUF_NAMESPACE_ID::uint64 dims(int index) const;
  void set_dims(int index, ::PROTOBUF_NAMESPACE_ID::uint64 value);
  void add_dims(::PROTOBUF_NAMESPACE_ID::uint64 value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint64 >&
      dims() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint64 >*
      mutable_dims();

  // repeated double edges = 100 [packed = true];
  int edges_size() const;
  private:
  int _internal_edges_size() const;
  public:
  void clear_edges();
  private:
  double _internal_edges(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      _internal_edges() const;
  void _internal_add_edges(double value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      _internal_mutable_edges();
  public:
  double edges(int index) const;
  void set_edges(int index, double value);
  void add_edges(double value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      edges() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      mutable_edges();

  // repeated double vals = 101 [packed = true];
  int vals_size() const;
  private:
  int _internal_vals_size() const;
  public:
  void clear_vals();
  private:
  double _internal_vals(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      _internal_vals() const;
  void _internal_add_vals(double value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      _internal_mutable_vals();
  public:
  double vals(int index) const;
  void set_vals(int index, double value);
  void add_vals(double value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      vals() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      mutable_vals();

  // repeated .lm.io.Hist.EdgeType edge_types = 200 [packed = true];
  int edge_types_size() const;
  private:
  int _internal_edge_types_size() const;
  public:
  void clear_edge_types();
  private:
  ::lm::io::Hist_EdgeType _internal_edge_types(int index) const;
  void _internal_add_edge_types(::lm::io::Hist_EdgeType value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField<int>* _internal_mutable_edge_types();
  public:
  ::lm::io::Hist_EdgeType edge_types(int index) const;
  void set_edge_types(int index, ::lm::io::Hist_EdgeType value);
  void add_edge_types(::lm::io::Hist_EdgeType value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField<int>& edge_types() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField<int>* mutable_edge_types();

  // optional uint64 rank = 1;
  bool has_rank() const;
  private:
  bool _internal_has_rank() const;
  public:
  void clear_rank();
  ::PROTOBUF_NAMESPACE_ID::uint64 rank() const;
  void set_rank(::PROTOBUF_NAMESPACE_ID::uint64 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::uint64 _internal_rank() const;
  void _internal_set_rank(::PROTOBUF_NAMESPACE_ID::uint64 value);
  public:

  // @@protoc_insertion_point(class_scope:lm.io.Hist)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint64 > dims_;
  mutable std::atomic<int> _dims_cached_byte_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double > edges_;
  mutable std::atomic<int> _edges_cached_byte_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double > vals_;
  mutable std::atomic<int> _vals_cached_byte_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField<int> edge_types_;
  mutable std::atomic<int> _edge_types_cached_byte_size_;
  ::PROTOBUF_NAMESPACE_ID::uint64 rank_;
  friend struct ::TableStruct_lm_2fio_2fHist_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// Hist

// optional uint64 rank = 1;
inline bool Hist::_internal_has_rank() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool Hist::has_rank() const {
  return _internal_has_rank();
}
inline void Hist::clear_rank() {
  rank_ = PROTOBUF_ULONGLONG(0);
  _has_bits_[0] &= ~0x00000001u;
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 Hist::_internal_rank() const {
  return rank_;
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 Hist::rank() const {
  // @@protoc_insertion_point(field_get:lm.io.Hist.rank)
  return _internal_rank();
}
inline void Hist::_internal_set_rank(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  _has_bits_[0] |= 0x00000001u;
  rank_ = value;
}
inline void Hist::set_rank(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  _internal_set_rank(value);
  // @@protoc_insertion_point(field_set:lm.io.Hist.rank)
}

// repeated uint64 dims = 2 [packed = true];
inline int Hist::_internal_dims_size() const {
  return dims_.size();
}
inline int Hist::dims_size() const {
  return _internal_dims_size();
}
inline void Hist::clear_dims() {
  dims_.Clear();
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 Hist::_internal_dims(int index) const {
  return dims_.Get(index);
}
inline ::PROTOBUF_NAMESPACE_ID::uint64 Hist::dims(int index) const {
  // @@protoc_insertion_point(field_get:lm.io.Hist.dims)
  return _internal_dims(index);
}
inline void Hist::set_dims(int index, ::PROTOBUF_NAMESPACE_ID::uint64 value) {
  dims_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.io.Hist.dims)
}
inline void Hist::_internal_add_dims(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  dims_.Add(value);
}
inline void Hist::add_dims(::PROTOBUF_NAMESPACE_ID::uint64 value) {
  _internal_add_dims(value);
  // @@protoc_insertion_point(field_add:lm.io.Hist.dims)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint64 >&
Hist::_internal_dims() const {
  return dims_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint64 >&
Hist::dims() const {
  // @@protoc_insertion_point(field_list:lm.io.Hist.dims)
  return _internal_dims();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint64 >*
Hist::_internal_mutable_dims() {
  return &dims_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint64 >*
Hist::mutable_dims() {
  // @@protoc_insertion_point(field_mutable_list:lm.io.Hist.dims)
  return _internal_mutable_dims();
}

// repeated double edges = 100 [packed = true];
inline int Hist::_internal_edges_size() const {
  return edges_.size();
}
inline int Hist::edges_size() const {
  return _internal_edges_size();
}
inline void Hist::clear_edges() {
  edges_.Clear();
}
inline double Hist::_internal_edges(int index) const {
  return edges_.Get(index);
}
inline double Hist::edges(int index) const {
  // @@protoc_insertion_point(field_get:lm.io.Hist.edges)
  return _internal_edges(index);
}
inline void Hist::set_edges(int index, double value) {
  edges_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.io.Hist.edges)
}
inline void Hist::_internal_add_edges(double value) {
  edges_.Add(value);
}
inline void Hist::add_edges(double value) {
  _internal_add_edges(value);
  // @@protoc_insertion_point(field_add:lm.io.Hist.edges)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
Hist::_internal_edges() const {
  return edges_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
Hist::edges() const {
  // @@protoc_insertion_point(field_list:lm.io.Hist.edges)
  return _internal_edges();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
Hist::_internal_mutable_edges() {
  return &edges_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
Hist::mutable_edges() {
  // @@protoc_insertion_point(field_mutable_list:lm.io.Hist.edges)
  return _internal_mutable_edges();
}

// repeated double vals = 101 [packed = true];
inline int Hist::_internal_vals_size() const {
  return vals_.size();
}
inline int Hist::vals_size() const {
  return _internal_vals_size();
}
inline void Hist::clear_vals() {
  vals_.Clear();
}
inline double Hist::_internal_vals(int index) const {
  return vals_.Get(index);
}
inline double Hist::vals(int index) const {
  // @@protoc_insertion_point(field_get:lm.io.Hist.vals)
  return _internal_vals(index);
}
inline void Hist::set_vals(int index, double value) {
  vals_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.io.Hist.vals)
}
inline void Hist::_internal_add_vals(double value) {
  vals_.Add(value);
}
inline void Hist::add_vals(double value) {
  _internal_add_vals(value);
  // @@protoc_insertion_point(field_add:lm.io.Hist.vals)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
Hist::_internal_vals() const {
  return vals_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
Hist::vals() const {
  // @@protoc_insertion_point(field_list:lm.io.Hist.vals)
  return _internal_vals();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
Hist::_internal_mutable_vals() {
  return &vals_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
Hist::mutable_vals() {
  // @@protoc_insertion_point(field_mutable_list:lm.io.Hist.vals)
  return _internal_mutable_vals();
}

// repeated .lm.io.Hist.EdgeType edge_types = 200 [packed = true];
inline int Hist::_internal_edge_types_size() const {
  return edge_types_.size();
}
inline int Hist::edge_types_size() const {
  return _internal_edge_types_size();
}
inline void Hist::clear_edge_types() {
  edge_types_.Clear();
}
inline ::lm::io::Hist_EdgeType Hist::_internal_edge_types(int index) const {
  return static_cast< ::lm::io::Hist_EdgeType >(edge_types_.Get(index));
}
inline ::lm::io::Hist_EdgeType Hist::edge_types(int index) const {
  // @@protoc_insertion_point(field_get:lm.io.Hist.edge_types)
  return _internal_edge_types(index);
}
inline void Hist::set_edge_types(int index, ::lm::io::Hist_EdgeType value) {
  assert(::lm::io::Hist_EdgeType_IsValid(value));
  edge_types_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.io.Hist.edge_types)
}
inline void Hist::_internal_add_edge_types(::lm::io::Hist_EdgeType value) {
  assert(::lm::io::Hist_EdgeType_IsValid(value));
  edge_types_.Add(value);
}
inline void Hist::add_edge_types(::lm::io::Hist_EdgeType value) {
  // @@protoc_insertion_point(field_add:lm.io.Hist.edge_types)
  _internal_add_edge_types(value);
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField<int>&
Hist::edge_types() const {
  // @@protoc_insertion_point(field_list:lm.io.Hist.edge_types)
  return edge_types_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField<int>*
Hist::_internal_mutable_edge_types() {
  return &edge_types_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField<int>*
Hist::mutable_edge_types() {
  // @@protoc_insertion_point(field_mutable_list:lm.io.Hist.edge_types)
  return _internal_mutable_edge_types();
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace io
}  // namespace lm

PROTOBUF_NAMESPACE_OPEN

template <> struct is_proto_enum< ::lm::io::Hist_EdgeType> : ::std::true_type {};
template <>
inline const EnumDescriptor* GetEnumDescriptor< ::lm::io::Hist_EdgeType>() {
  return ::lm::io::Hist_EdgeType_descriptor();
}

PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fHist_2eproto