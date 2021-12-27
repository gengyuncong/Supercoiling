// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/input/DiffusionModel.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2finput_2fDiffusionModel_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2finput_2fDiffusionModel_2eproto

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
#include "lm/types/BoundaryConditions.pb.h"
#include "lm/types/Lattice.pb.h"
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_lm_2finput_2fDiffusionModel_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2finput_2fDiffusionModel_2eproto {
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
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2finput_2fDiffusionModel_2eproto;
namespace lm {
namespace input {
class DiffusionModel;
class DiffusionModelDefaultTypeInternal;
extern DiffusionModelDefaultTypeInternal _DiffusionModel_default_instance_;
}  // namespace input
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::input::DiffusionModel* Arena::CreateMaybeMessage<::lm::input::DiffusionModel>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace input {

// ===================================================================

class DiffusionModel PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.input.DiffusionModel) */ {
 public:
  inline DiffusionModel() : DiffusionModel(nullptr) {}
  virtual ~DiffusionModel();

  DiffusionModel(const DiffusionModel& from);
  DiffusionModel(DiffusionModel&& from) noexcept
    : DiffusionModel() {
    *this = ::std::move(from);
  }

  inline DiffusionModel& operator=(const DiffusionModel& from) {
    CopyFrom(from);
    return *this;
  }
  inline DiffusionModel& operator=(DiffusionModel&& from) noexcept {
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
  static const DiffusionModel& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const DiffusionModel* internal_default_instance() {
    return reinterpret_cast<const DiffusionModel*>(
               &_DiffusionModel_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(DiffusionModel& a, DiffusionModel& b) {
    a.Swap(&b);
  }
  inline void Swap(DiffusionModel* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(DiffusionModel* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline DiffusionModel* New() const final {
    return CreateMaybeMessage<DiffusionModel>(nullptr);
  }

  DiffusionModel* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<DiffusionModel>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const DiffusionModel& from);
  void MergeFrom(const DiffusionModel& from);
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
  void InternalSwap(DiffusionModel* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.input.DiffusionModel";
  }
  protected:
  explicit DiffusionModel(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2finput_2fDiffusionModel_2eproto);
    return ::descriptor_table_lm_2finput_2fDiffusionModel_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kDiffusionMatrixFieldNumber = 4,
    kReactionLocationMatrixFieldNumber = 5,
    kInitialLatticeFieldNumber = 7,
    kBoundaryConditionsFieldNumber = 8,
    kNumberSpeciesFieldNumber = 1,
    kNumberReactionsFieldNumber = 2,
    kLatticeSpacingFieldNumber = 6,
    kNumberSiteTypesFieldNumber = 3,
  };
  // repeated double diffusion_matrix = 4 [packed = true];
  int diffusion_matrix_size() const;
  private:
  int _internal_diffusion_matrix_size() const;
  public:
  void clear_diffusion_matrix();
  private:
  double _internal_diffusion_matrix(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      _internal_diffusion_matrix() const;
  void _internal_add_diffusion_matrix(double value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      _internal_mutable_diffusion_matrix();
  public:
  double diffusion_matrix(int index) const;
  void set_diffusion_matrix(int index, double value);
  void add_diffusion_matrix(double value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      diffusion_matrix() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      mutable_diffusion_matrix();

  // repeated uint32 reaction_location_matrix = 5 [packed = true];
  int reaction_location_matrix_size() const;
  private:
  int _internal_reaction_location_matrix_size() const;
  public:
  void clear_reaction_location_matrix();
  private:
  ::PROTOBUF_NAMESPACE_ID::uint32 _internal_reaction_location_matrix(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint32 >&
      _internal_reaction_location_matrix() const;
  void _internal_add_reaction_location_matrix(::PROTOBUF_NAMESPACE_ID::uint32 value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint32 >*
      _internal_mutable_reaction_location_matrix();
  public:
  ::PROTOBUF_NAMESPACE_ID::uint32 reaction_location_matrix(int index) const;
  void set_reaction_location_matrix(int index, ::PROTOBUF_NAMESPACE_ID::uint32 value);
  void add_reaction_location_matrix(::PROTOBUF_NAMESPACE_ID::uint32 value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint32 >&
      reaction_location_matrix() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint32 >*
      mutable_reaction_location_matrix();

  // required .lm.types.Lattice initial_lattice = 7;
  bool has_initial_lattice() const;
  private:
  bool _internal_has_initial_lattice() const;
  public:
  void clear_initial_lattice();
  const ::lm::types::Lattice& initial_lattice() const;
  ::lm::types::Lattice* release_initial_lattice();
  ::lm::types::Lattice* mutable_initial_lattice();
  void set_allocated_initial_lattice(::lm::types::Lattice* initial_lattice);
  private:
  const ::lm::types::Lattice& _internal_initial_lattice() const;
  ::lm::types::Lattice* _internal_mutable_initial_lattice();
  public:
  void unsafe_arena_set_allocated_initial_lattice(
      ::lm::types::Lattice* initial_lattice);
  ::lm::types::Lattice* unsafe_arena_release_initial_lattice();

  // optional .lm.types.BoundaryConditions boundary_conditions = 8;
  bool has_boundary_conditions() const;
  private:
  bool _internal_has_boundary_conditions() const;
  public:
  void clear_boundary_conditions();
  const ::lm::types::BoundaryConditions& boundary_conditions() const;
  ::lm::types::BoundaryConditions* release_boundary_conditions();
  ::lm::types::BoundaryConditions* mutable_boundary_conditions();
  void set_allocated_boundary_conditions(::lm::types::BoundaryConditions* boundary_conditions);
  private:
  const ::lm::types::BoundaryConditions& _internal_boundary_conditions() const;
  ::lm::types::BoundaryConditions* _internal_mutable_boundary_conditions();
  public:
  void unsafe_arena_set_allocated_boundary_conditions(
      ::lm::types::BoundaryConditions* boundary_conditions);
  ::lm::types::BoundaryConditions* unsafe_arena_release_boundary_conditions();

  // required int32 number_species = 1;
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

  // required int32 number_reactions = 2;
  bool has_number_reactions() const;
  private:
  bool _internal_has_number_reactions() const;
  public:
  void clear_number_reactions();
  ::PROTOBUF_NAMESPACE_ID::int32 number_reactions() const;
  void set_number_reactions(::PROTOBUF_NAMESPACE_ID::int32 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::int32 _internal_number_reactions() const;
  void _internal_set_number_reactions(::PROTOBUF_NAMESPACE_ID::int32 value);
  public:

  // required double lattice_spacing = 6;
  bool has_lattice_spacing() const;
  private:
  bool _internal_has_lattice_spacing() const;
  public:
  void clear_lattice_spacing();
  double lattice_spacing() const;
  void set_lattice_spacing(double value);
  private:
  double _internal_lattice_spacing() const;
  void _internal_set_lattice_spacing(double value);
  public:

  // required int32 number_site_types = 3;
  bool has_number_site_types() const;
  private:
  bool _internal_has_number_site_types() const;
  public:
  void clear_number_site_types();
  ::PROTOBUF_NAMESPACE_ID::int32 number_site_types() const;
  void set_number_site_types(::PROTOBUF_NAMESPACE_ID::int32 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::int32 _internal_number_site_types() const;
  void _internal_set_number_site_types(::PROTOBUF_NAMESPACE_ID::int32 value);
  public:

  // @@protoc_insertion_point(class_scope:lm.input.DiffusionModel)
 private:
  class _Internal;

  // helper for ByteSizeLong()
  size_t RequiredFieldsByteSizeFallback() const;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double > diffusion_matrix_;
  mutable std::atomic<int> _diffusion_matrix_cached_byte_size_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint32 > reaction_location_matrix_;
  mutable std::atomic<int> _reaction_location_matrix_cached_byte_size_;
  ::lm::types::Lattice* initial_lattice_;
  ::lm::types::BoundaryConditions* boundary_conditions_;
  ::PROTOBUF_NAMESPACE_ID::int32 number_species_;
  ::PROTOBUF_NAMESPACE_ID::int32 number_reactions_;
  double lattice_spacing_;
  ::PROTOBUF_NAMESPACE_ID::int32 number_site_types_;
  friend struct ::TableStruct_lm_2finput_2fDiffusionModel_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// DiffusionModel

// required int32 number_species = 1;
inline bool DiffusionModel::_internal_has_number_species() const {
  bool value = (_has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool DiffusionModel::has_number_species() const {
  return _internal_has_number_species();
}
inline void DiffusionModel::clear_number_species() {
  number_species_ = 0;
  _has_bits_[0] &= ~0x00000004u;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 DiffusionModel::_internal_number_species() const {
  return number_species_;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 DiffusionModel::number_species() const {
  // @@protoc_insertion_point(field_get:lm.input.DiffusionModel.number_species)
  return _internal_number_species();
}
inline void DiffusionModel::_internal_set_number_species(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _has_bits_[0] |= 0x00000004u;
  number_species_ = value;
}
inline void DiffusionModel::set_number_species(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_set_number_species(value);
  // @@protoc_insertion_point(field_set:lm.input.DiffusionModel.number_species)
}

// required int32 number_reactions = 2;
inline bool DiffusionModel::_internal_has_number_reactions() const {
  bool value = (_has_bits_[0] & 0x00000008u) != 0;
  return value;
}
inline bool DiffusionModel::has_number_reactions() const {
  return _internal_has_number_reactions();
}
inline void DiffusionModel::clear_number_reactions() {
  number_reactions_ = 0;
  _has_bits_[0] &= ~0x00000008u;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 DiffusionModel::_internal_number_reactions() const {
  return number_reactions_;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 DiffusionModel::number_reactions() const {
  // @@protoc_insertion_point(field_get:lm.input.DiffusionModel.number_reactions)
  return _internal_number_reactions();
}
inline void DiffusionModel::_internal_set_number_reactions(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _has_bits_[0] |= 0x00000008u;
  number_reactions_ = value;
}
inline void DiffusionModel::set_number_reactions(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_set_number_reactions(value);
  // @@protoc_insertion_point(field_set:lm.input.DiffusionModel.number_reactions)
}

// required int32 number_site_types = 3;
inline bool DiffusionModel::_internal_has_number_site_types() const {
  bool value = (_has_bits_[0] & 0x00000020u) != 0;
  return value;
}
inline bool DiffusionModel::has_number_site_types() const {
  return _internal_has_number_site_types();
}
inline void DiffusionModel::clear_number_site_types() {
  number_site_types_ = 0;
  _has_bits_[0] &= ~0x00000020u;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 DiffusionModel::_internal_number_site_types() const {
  return number_site_types_;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 DiffusionModel::number_site_types() const {
  // @@protoc_insertion_point(field_get:lm.input.DiffusionModel.number_site_types)
  return _internal_number_site_types();
}
inline void DiffusionModel::_internal_set_number_site_types(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _has_bits_[0] |= 0x00000020u;
  number_site_types_ = value;
}
inline void DiffusionModel::set_number_site_types(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_set_number_site_types(value);
  // @@protoc_insertion_point(field_set:lm.input.DiffusionModel.number_site_types)
}

// repeated double diffusion_matrix = 4 [packed = true];
inline int DiffusionModel::_internal_diffusion_matrix_size() const {
  return diffusion_matrix_.size();
}
inline int DiffusionModel::diffusion_matrix_size() const {
  return _internal_diffusion_matrix_size();
}
inline void DiffusionModel::clear_diffusion_matrix() {
  diffusion_matrix_.Clear();
}
inline double DiffusionModel::_internal_diffusion_matrix(int index) const {
  return diffusion_matrix_.Get(index);
}
inline double DiffusionModel::diffusion_matrix(int index) const {
  // @@protoc_insertion_point(field_get:lm.input.DiffusionModel.diffusion_matrix)
  return _internal_diffusion_matrix(index);
}
inline void DiffusionModel::set_diffusion_matrix(int index, double value) {
  diffusion_matrix_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.input.DiffusionModel.diffusion_matrix)
}
inline void DiffusionModel::_internal_add_diffusion_matrix(double value) {
  diffusion_matrix_.Add(value);
}
inline void DiffusionModel::add_diffusion_matrix(double value) {
  _internal_add_diffusion_matrix(value);
  // @@protoc_insertion_point(field_add:lm.input.DiffusionModel.diffusion_matrix)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
DiffusionModel::_internal_diffusion_matrix() const {
  return diffusion_matrix_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
DiffusionModel::diffusion_matrix() const {
  // @@protoc_insertion_point(field_list:lm.input.DiffusionModel.diffusion_matrix)
  return _internal_diffusion_matrix();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
DiffusionModel::_internal_mutable_diffusion_matrix() {
  return &diffusion_matrix_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
DiffusionModel::mutable_diffusion_matrix() {
  // @@protoc_insertion_point(field_mutable_list:lm.input.DiffusionModel.diffusion_matrix)
  return _internal_mutable_diffusion_matrix();
}

// repeated uint32 reaction_location_matrix = 5 [packed = true];
inline int DiffusionModel::_internal_reaction_location_matrix_size() const {
  return reaction_location_matrix_.size();
}
inline int DiffusionModel::reaction_location_matrix_size() const {
  return _internal_reaction_location_matrix_size();
}
inline void DiffusionModel::clear_reaction_location_matrix() {
  reaction_location_matrix_.Clear();
}
inline ::PROTOBUF_NAMESPACE_ID::uint32 DiffusionModel::_internal_reaction_location_matrix(int index) const {
  return reaction_location_matrix_.Get(index);
}
inline ::PROTOBUF_NAMESPACE_ID::uint32 DiffusionModel::reaction_location_matrix(int index) const {
  // @@protoc_insertion_point(field_get:lm.input.DiffusionModel.reaction_location_matrix)
  return _internal_reaction_location_matrix(index);
}
inline void DiffusionModel::set_reaction_location_matrix(int index, ::PROTOBUF_NAMESPACE_ID::uint32 value) {
  reaction_location_matrix_.Set(index, value);
  // @@protoc_insertion_point(field_set:lm.input.DiffusionModel.reaction_location_matrix)
}
inline void DiffusionModel::_internal_add_reaction_location_matrix(::PROTOBUF_NAMESPACE_ID::uint32 value) {
  reaction_location_matrix_.Add(value);
}
inline void DiffusionModel::add_reaction_location_matrix(::PROTOBUF_NAMESPACE_ID::uint32 value) {
  _internal_add_reaction_location_matrix(value);
  // @@protoc_insertion_point(field_add:lm.input.DiffusionModel.reaction_location_matrix)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint32 >&
DiffusionModel::_internal_reaction_location_matrix() const {
  return reaction_location_matrix_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint32 >&
DiffusionModel::reaction_location_matrix() const {
  // @@protoc_insertion_point(field_list:lm.input.DiffusionModel.reaction_location_matrix)
  return _internal_reaction_location_matrix();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint32 >*
DiffusionModel::_internal_mutable_reaction_location_matrix() {
  return &reaction_location_matrix_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< ::PROTOBUF_NAMESPACE_ID::uint32 >*
DiffusionModel::mutable_reaction_location_matrix() {
  // @@protoc_insertion_point(field_mutable_list:lm.input.DiffusionModel.reaction_location_matrix)
  return _internal_mutable_reaction_location_matrix();
}

// required double lattice_spacing = 6;
inline bool DiffusionModel::_internal_has_lattice_spacing() const {
  bool value = (_has_bits_[0] & 0x00000010u) != 0;
  return value;
}
inline bool DiffusionModel::has_lattice_spacing() const {
  return _internal_has_lattice_spacing();
}
inline void DiffusionModel::clear_lattice_spacing() {
  lattice_spacing_ = 0;
  _has_bits_[0] &= ~0x00000010u;
}
inline double DiffusionModel::_internal_lattice_spacing() const {
  return lattice_spacing_;
}
inline double DiffusionModel::lattice_spacing() const {
  // @@protoc_insertion_point(field_get:lm.input.DiffusionModel.lattice_spacing)
  return _internal_lattice_spacing();
}
inline void DiffusionModel::_internal_set_lattice_spacing(double value) {
  _has_bits_[0] |= 0x00000010u;
  lattice_spacing_ = value;
}
inline void DiffusionModel::set_lattice_spacing(double value) {
  _internal_set_lattice_spacing(value);
  // @@protoc_insertion_point(field_set:lm.input.DiffusionModel.lattice_spacing)
}

// required .lm.types.Lattice initial_lattice = 7;
inline bool DiffusionModel::_internal_has_initial_lattice() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  PROTOBUF_ASSUME(!value || initial_lattice_ != nullptr);
  return value;
}
inline bool DiffusionModel::has_initial_lattice() const {
  return _internal_has_initial_lattice();
}
inline const ::lm::types::Lattice& DiffusionModel::_internal_initial_lattice() const {
  const ::lm::types::Lattice* p = initial_lattice_;
  return p != nullptr ? *p : *reinterpret_cast<const ::lm::types::Lattice*>(
      &::lm::types::_Lattice_default_instance_);
}
inline const ::lm::types::Lattice& DiffusionModel::initial_lattice() const {
  // @@protoc_insertion_point(field_get:lm.input.DiffusionModel.initial_lattice)
  return _internal_initial_lattice();
}
inline void DiffusionModel::unsafe_arena_set_allocated_initial_lattice(
    ::lm::types::Lattice* initial_lattice) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(initial_lattice_);
  }
  initial_lattice_ = initial_lattice;
  if (initial_lattice) {
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.input.DiffusionModel.initial_lattice)
}
inline ::lm::types::Lattice* DiffusionModel::release_initial_lattice() {
  _has_bits_[0] &= ~0x00000001u;
  ::lm::types::Lattice* temp = initial_lattice_;
  initial_lattice_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::lm::types::Lattice* DiffusionModel::unsafe_arena_release_initial_lattice() {
  // @@protoc_insertion_point(field_release:lm.input.DiffusionModel.initial_lattice)
  _has_bits_[0] &= ~0x00000001u;
  ::lm::types::Lattice* temp = initial_lattice_;
  initial_lattice_ = nullptr;
  return temp;
}
inline ::lm::types::Lattice* DiffusionModel::_internal_mutable_initial_lattice() {
  _has_bits_[0] |= 0x00000001u;
  if (initial_lattice_ == nullptr) {
    auto* p = CreateMaybeMessage<::lm::types::Lattice>(GetArena());
    initial_lattice_ = p;
  }
  return initial_lattice_;
}
inline ::lm::types::Lattice* DiffusionModel::mutable_initial_lattice() {
  // @@protoc_insertion_point(field_mutable:lm.input.DiffusionModel.initial_lattice)
  return _internal_mutable_initial_lattice();
}
inline void DiffusionModel::set_allocated_initial_lattice(::lm::types::Lattice* initial_lattice) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(initial_lattice_);
  }
  if (initial_lattice) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(initial_lattice)->GetArena();
    if (message_arena != submessage_arena) {
      initial_lattice = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, initial_lattice, submessage_arena);
    }
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  initial_lattice_ = initial_lattice;
  // @@protoc_insertion_point(field_set_allocated:lm.input.DiffusionModel.initial_lattice)
}

// optional .lm.types.BoundaryConditions boundary_conditions = 8;
inline bool DiffusionModel::_internal_has_boundary_conditions() const {
  bool value = (_has_bits_[0] & 0x00000002u) != 0;
  PROTOBUF_ASSUME(!value || boundary_conditions_ != nullptr);
  return value;
}
inline bool DiffusionModel::has_boundary_conditions() const {
  return _internal_has_boundary_conditions();
}
inline const ::lm::types::BoundaryConditions& DiffusionModel::_internal_boundary_conditions() const {
  const ::lm::types::BoundaryConditions* p = boundary_conditions_;
  return p != nullptr ? *p : *reinterpret_cast<const ::lm::types::BoundaryConditions*>(
      &::lm::types::_BoundaryConditions_default_instance_);
}
inline const ::lm::types::BoundaryConditions& DiffusionModel::boundary_conditions() const {
  // @@protoc_insertion_point(field_get:lm.input.DiffusionModel.boundary_conditions)
  return _internal_boundary_conditions();
}
inline void DiffusionModel::unsafe_arena_set_allocated_boundary_conditions(
    ::lm::types::BoundaryConditions* boundary_conditions) {
  if (GetArena() == nullptr) {
    delete reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(boundary_conditions_);
  }
  boundary_conditions_ = boundary_conditions;
  if (boundary_conditions) {
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  // @@protoc_insertion_point(field_unsafe_arena_set_allocated:lm.input.DiffusionModel.boundary_conditions)
}
inline ::lm::types::BoundaryConditions* DiffusionModel::release_boundary_conditions() {
  _has_bits_[0] &= ~0x00000002u;
  ::lm::types::BoundaryConditions* temp = boundary_conditions_;
  boundary_conditions_ = nullptr;
  if (GetArena() != nullptr) {
    temp = ::PROTOBUF_NAMESPACE_ID::internal::DuplicateIfNonNull(temp);
  }
  return temp;
}
inline ::lm::types::BoundaryConditions* DiffusionModel::unsafe_arena_release_boundary_conditions() {
  // @@protoc_insertion_point(field_release:lm.input.DiffusionModel.boundary_conditions)
  _has_bits_[0] &= ~0x00000002u;
  ::lm::types::BoundaryConditions* temp = boundary_conditions_;
  boundary_conditions_ = nullptr;
  return temp;
}
inline ::lm::types::BoundaryConditions* DiffusionModel::_internal_mutable_boundary_conditions() {
  _has_bits_[0] |= 0x00000002u;
  if (boundary_conditions_ == nullptr) {
    auto* p = CreateMaybeMessage<::lm::types::BoundaryConditions>(GetArena());
    boundary_conditions_ = p;
  }
  return boundary_conditions_;
}
inline ::lm::types::BoundaryConditions* DiffusionModel::mutable_boundary_conditions() {
  // @@protoc_insertion_point(field_mutable:lm.input.DiffusionModel.boundary_conditions)
  return _internal_mutable_boundary_conditions();
}
inline void DiffusionModel::set_allocated_boundary_conditions(::lm::types::BoundaryConditions* boundary_conditions) {
  ::PROTOBUF_NAMESPACE_ID::Arena* message_arena = GetArena();
  if (message_arena == nullptr) {
    delete reinterpret_cast< ::PROTOBUF_NAMESPACE_ID::MessageLite*>(boundary_conditions_);
  }
  if (boundary_conditions) {
    ::PROTOBUF_NAMESPACE_ID::Arena* submessage_arena =
      reinterpret_cast<::PROTOBUF_NAMESPACE_ID::MessageLite*>(boundary_conditions)->GetArena();
    if (message_arena != submessage_arena) {
      boundary_conditions = ::PROTOBUF_NAMESPACE_ID::internal::GetOwnedMessage(
          message_arena, boundary_conditions, submessage_arena);
    }
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  boundary_conditions_ = boundary_conditions;
  // @@protoc_insertion_point(field_set_allocated:lm.input.DiffusionModel.boundary_conditions)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace input
}  // namespace lm

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2finput_2fDiffusionModel_2eproto
