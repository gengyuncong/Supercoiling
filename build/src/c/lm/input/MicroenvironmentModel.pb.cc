// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/input/MicroenvironmentModel.proto

#include "lm/input/MicroenvironmentModel.pb.h"

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
extern PROTOBUF_INTERNAL_EXPORT_lm_2ftypes_2fBoundaryConditions_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_BoundaryConditions_lm_2ftypes_2fBoundaryConditions_2eproto;
extern PROTOBUF_INTERNAL_EXPORT_robertslab_2fpbuf_2fNDArray_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto;
namespace lm {
namespace input {
class MicroenvironmentModelDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<MicroenvironmentModel> _instance;
} _MicroenvironmentModel_default_instance_;
}  // namespace input
}  // namespace lm
static void InitDefaultsscc_info_MicroenvironmentModel_lm_2finput_2fMicroenvironmentModel_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::input::_MicroenvironmentModel_default_instance_;
    new (ptr) ::lm::input::MicroenvironmentModel();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::input::MicroenvironmentModel::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<2> scc_info_MicroenvironmentModel_lm_2finput_2fMicroenvironmentModel_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 2, 0, InitDefaultsscc_info_MicroenvironmentModel_lm_2finput_2fMicroenvironmentModel_2eproto}, {
      &scc_info_BoundaryConditions_lm_2ftypes_2fBoundaryConditions_2eproto.base,
      &scc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto.base,}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2finput_2fMicroenvironmentModel_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2finput_2fMicroenvironmentModel_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2finput_2fMicroenvironmentModel_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2finput_2fMicroenvironmentModel_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, initial_time_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, grid_shape_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, grid_spacing_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, boundaries_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, species_index_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, diffusion_coefficients_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, initial_concentrations_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, synchronization_timestep_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, number_cells_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, cell_initial_species_counts_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, cell_coordinates_),
  PROTOBUF_FIELD_OFFSET(::lm::input::MicroenvironmentModel, cell_volume_),
  6,
  ~0u,
  4,
  0,
  ~0u,
  ~0u,
  ~0u,
  5,
  7,
  1,
  2,
  3,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 17, sizeof(::lm::input::MicroenvironmentModel)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::input::_MicroenvironmentModel_default_instance_),
};

const char descriptor_table_protodef_lm_2finput_2fMicroenvironmentModel_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n$lm/input/MicroenvironmentModel.proto\022\010"
  "lm.input\032!lm/types/BoundaryConditions.pr"
  "oto\032\035robertslab/pbuf/NDArray.proto\"\332\003\n\025M"
  "icroenvironmentModel\022\027\n\014initial_time\030\014 \001"
  "(\001:\0010\022\022\n\ngrid_shape\030\001 \003(\005\022\024\n\014grid_spacin"
  "g\030\002 \002(\001\0220\n\nboundaries\030\003 \002(\0132\034.lm.types.B"
  "oundaryConditions\022\025\n\rspecies_index\030\004 \003(\005"
  "\022\036\n\026diffusion_coefficients\030\005 \003(\001\0228\n\026init"
  "ial_concentrations\030\006 \003(\0132\030.robertslab.pb"
  "uf.NDArray\022 \n\030synchronization_timestep\030\007"
  " \002(\001\022\027\n\014number_cells\030\010 \001(\r:\0010\022=\n\033cell_in"
  "itial_species_counts\030\t \001(\0132\030.robertslab."
  "pbuf.NDArray\0222\n\020cell_coordinates\030\n \001(\0132\030"
  ".robertslab.pbuf.NDArray\022-\n\013cell_volume\030"
  "\013 \001(\0132\030.robertslab.pbuf.NDArray"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2finput_2fMicroenvironmentModel_2eproto_deps[2] = {
  &::descriptor_table_lm_2ftypes_2fBoundaryConditions_2eproto,
  &::descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto,
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2finput_2fMicroenvironmentModel_2eproto_sccs[1] = {
  &scc_info_MicroenvironmentModel_lm_2finput_2fMicroenvironmentModel_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2finput_2fMicroenvironmentModel_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2finput_2fMicroenvironmentModel_2eproto = {
  false, false, descriptor_table_protodef_lm_2finput_2fMicroenvironmentModel_2eproto, "lm/input/MicroenvironmentModel.proto", 591,
  &descriptor_table_lm_2finput_2fMicroenvironmentModel_2eproto_once, descriptor_table_lm_2finput_2fMicroenvironmentModel_2eproto_sccs, descriptor_table_lm_2finput_2fMicroenvironmentModel_2eproto_deps, 1, 2,
  schemas, file_default_instances, TableStruct_lm_2finput_2fMicroenvironmentModel_2eproto::offsets,
  file_level_metadata_lm_2finput_2fMicroenvironmentModel_2eproto, 1, file_level_enum_descriptors_lm_2finput_2fMicroenvironmentModel_2eproto, file_level_service_descriptors_lm_2finput_2fMicroenvironmentModel_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2finput_2fMicroenvironmentModel_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2finput_2fMicroenvironmentModel_2eproto)), true);
namespace lm {
namespace input {

// ===================================================================

void MicroenvironmentModel::InitAsDefaultInstance() {
  ::lm::input::_MicroenvironmentModel_default_instance_._instance.get_mutable()->boundaries_ = const_cast< ::lm::types::BoundaryConditions*>(
      ::lm::types::BoundaryConditions::internal_default_instance());
  ::lm::input::_MicroenvironmentModel_default_instance_._instance.get_mutable()->cell_initial_species_counts_ = const_cast< ::robertslab::pbuf::NDArray*>(
      ::robertslab::pbuf::NDArray::internal_default_instance());
  ::lm::input::_MicroenvironmentModel_default_instance_._instance.get_mutable()->cell_coordinates_ = const_cast< ::robertslab::pbuf::NDArray*>(
      ::robertslab::pbuf::NDArray::internal_default_instance());
  ::lm::input::_MicroenvironmentModel_default_instance_._instance.get_mutable()->cell_volume_ = const_cast< ::robertslab::pbuf::NDArray*>(
      ::robertslab::pbuf::NDArray::internal_default_instance());
}
class MicroenvironmentModel::_Internal {
 public:
  using HasBits = decltype(std::declval<MicroenvironmentModel>()._has_bits_);
  static void set_has_initial_time(HasBits* has_bits) {
    (*has_bits)[0] |= 64u;
  }
  static void set_has_grid_spacing(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
  static const ::lm::types::BoundaryConditions& boundaries(const MicroenvironmentModel* msg);
  static void set_has_boundaries(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_synchronization_timestep(HasBits* has_bits) {
    (*has_bits)[0] |= 32u;
  }
  static void set_has_number_cells(HasBits* has_bits) {
    (*has_bits)[0] |= 128u;
  }
  static const ::robertslab::pbuf::NDArray& cell_initial_species_counts(const MicroenvironmentModel* msg);
  static void set_has_cell_initial_species_counts(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static const ::robertslab::pbuf::NDArray& cell_coordinates(const MicroenvironmentModel* msg);
  static void set_has_cell_coordinates(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static const ::robertslab::pbuf::NDArray& cell_volume(const MicroenvironmentModel* msg);
  static void set_has_cell_volume(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000031) ^ 0x00000031) != 0;
  }
};

const ::lm::types::BoundaryConditions&
MicroenvironmentModel::_Internal::boundaries(const MicroenvironmentModel* msg) {
  return *msg->boundaries_;
}
const ::robertslab::pbuf::NDArray&
MicroenvironmentModel::_Internal::cell_initial_species_counts(const MicroenvironmentModel* msg) {
  return *msg->cell_initial_species_counts_;
}
const ::robertslab::pbuf::NDArray&
MicroenvironmentModel::_Internal::cell_coordinates(const MicroenvironmentModel* msg) {
  return *msg->cell_coordinates_;
}
const ::robertslab::pbuf::NDArray&
MicroenvironmentModel::_Internal::cell_volume(const MicroenvironmentModel* msg) {
  return *msg->cell_volume_;
}
void MicroenvironmentModel::clear_boundaries() {
  if (boundaries_ != nullptr) boundaries_->Clear();
  _has_bits_[0] &= ~0x00000001u;
}
void MicroenvironmentModel::clear_initial_concentrations() {
  initial_concentrations_.Clear();
}
void MicroenvironmentModel::clear_cell_initial_species_counts() {
  if (cell_initial_species_counts_ != nullptr) cell_initial_species_counts_->Clear();
  _has_bits_[0] &= ~0x00000002u;
}
void MicroenvironmentModel::clear_cell_coordinates() {
  if (cell_coordinates_ != nullptr) cell_coordinates_->Clear();
  _has_bits_[0] &= ~0x00000004u;
}
void MicroenvironmentModel::clear_cell_volume() {
  if (cell_volume_ != nullptr) cell_volume_->Clear();
  _has_bits_[0] &= ~0x00000008u;
}
MicroenvironmentModel::MicroenvironmentModel(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  grid_shape_(arena),
  species_index_(arena),
  diffusion_coefficients_(arena),
  initial_concentrations_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.input.MicroenvironmentModel)
}
MicroenvironmentModel::MicroenvironmentModel(const MicroenvironmentModel& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_),
      grid_shape_(from.grid_shape_),
      species_index_(from.species_index_),
      diffusion_coefficients_(from.diffusion_coefficients_),
      initial_concentrations_(from.initial_concentrations_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  if (from._internal_has_boundaries()) {
    boundaries_ = new ::lm::types::BoundaryConditions(*from.boundaries_);
  } else {
    boundaries_ = nullptr;
  }
  if (from._internal_has_cell_initial_species_counts()) {
    cell_initial_species_counts_ = new ::robertslab::pbuf::NDArray(*from.cell_initial_species_counts_);
  } else {
    cell_initial_species_counts_ = nullptr;
  }
  if (from._internal_has_cell_coordinates()) {
    cell_coordinates_ = new ::robertslab::pbuf::NDArray(*from.cell_coordinates_);
  } else {
    cell_coordinates_ = nullptr;
  }
  if (from._internal_has_cell_volume()) {
    cell_volume_ = new ::robertslab::pbuf::NDArray(*from.cell_volume_);
  } else {
    cell_volume_ = nullptr;
  }
  ::memcpy(&grid_spacing_, &from.grid_spacing_,
    static_cast<size_t>(reinterpret_cast<char*>(&number_cells_) -
    reinterpret_cast<char*>(&grid_spacing_)) + sizeof(number_cells_));
  // @@protoc_insertion_point(copy_constructor:lm.input.MicroenvironmentModel)
}

void MicroenvironmentModel::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_MicroenvironmentModel_lm_2finput_2fMicroenvironmentModel_2eproto.base);
  ::memset(&boundaries_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&number_cells_) -
      reinterpret_cast<char*>(&boundaries_)) + sizeof(number_cells_));
}

MicroenvironmentModel::~MicroenvironmentModel() {
  // @@protoc_insertion_point(destructor:lm.input.MicroenvironmentModel)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void MicroenvironmentModel::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
  if (this != internal_default_instance()) delete boundaries_;
  if (this != internal_default_instance()) delete cell_initial_species_counts_;
  if (this != internal_default_instance()) delete cell_coordinates_;
  if (this != internal_default_instance()) delete cell_volume_;
}

void MicroenvironmentModel::ArenaDtor(void* object) {
  MicroenvironmentModel* _this = reinterpret_cast< MicroenvironmentModel* >(object);
  (void)_this;
}
void MicroenvironmentModel::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void MicroenvironmentModel::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const MicroenvironmentModel& MicroenvironmentModel::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_MicroenvironmentModel_lm_2finput_2fMicroenvironmentModel_2eproto.base);
  return *internal_default_instance();
}


void MicroenvironmentModel::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.input.MicroenvironmentModel)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  grid_shape_.Clear();
  species_index_.Clear();
  diffusion_coefficients_.Clear();
  initial_concentrations_.Clear();
  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x0000000fu) {
    if (cached_has_bits & 0x00000001u) {
      GOOGLE_DCHECK(boundaries_ != nullptr);
      boundaries_->Clear();
    }
    if (cached_has_bits & 0x00000002u) {
      GOOGLE_DCHECK(cell_initial_species_counts_ != nullptr);
      cell_initial_species_counts_->Clear();
    }
    if (cached_has_bits & 0x00000004u) {
      GOOGLE_DCHECK(cell_coordinates_ != nullptr);
      cell_coordinates_->Clear();
    }
    if (cached_has_bits & 0x00000008u) {
      GOOGLE_DCHECK(cell_volume_ != nullptr);
      cell_volume_->Clear();
    }
  }
  if (cached_has_bits & 0x000000f0u) {
    ::memset(&grid_spacing_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&number_cells_) -
        reinterpret_cast<char*>(&grid_spacing_)) + sizeof(number_cells_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* MicroenvironmentModel::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated int32 grid_shape = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_grid_shape(::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr));
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<8>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedInt32Parser(_internal_mutable_grid_shape(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // required double grid_spacing = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 17)) {
          _Internal::set_has_grid_spacing(&has_bits);
          grid_spacing_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else goto handle_unusual;
        continue;
      // required .lm.types.BoundaryConditions boundaries = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 26)) {
          ptr = ctx->ParseMessage(_internal_mutable_boundaries(), ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated int32 species_index = 4;
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 32)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_species_index(::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr));
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<32>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 34) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedInt32Parser(_internal_mutable_species_index(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated double diffusion_coefficients = 5;
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 41)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_diffusion_coefficients(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
            ptr += sizeof(double);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<41>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 42) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_diffusion_coefficients(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated .robertslab.pbuf.NDArray initial_concentrations = 6;
      case 6:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 50)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_initial_concentrations(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<50>(ptr));
        } else goto handle_unusual;
        continue;
      // required double synchronization_timestep = 7;
      case 7:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 57)) {
          _Internal::set_has_synchronization_timestep(&has_bits);
          synchronization_timestep_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else goto handle_unusual;
        continue;
      // optional uint32 number_cells = 8 [default = 0];
      case 8:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 64)) {
          _Internal::set_has_number_cells(&has_bits);
          number_cells_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional .robertslab.pbuf.NDArray cell_initial_species_counts = 9;
      case 9:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 74)) {
          ptr = ctx->ParseMessage(_internal_mutable_cell_initial_species_counts(), ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional .robertslab.pbuf.NDArray cell_coordinates = 10;
      case 10:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 82)) {
          ptr = ctx->ParseMessage(_internal_mutable_cell_coordinates(), ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional .robertslab.pbuf.NDArray cell_volume = 11;
      case 11:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 90)) {
          ptr = ctx->ParseMessage(_internal_mutable_cell_volume(), ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional double initial_time = 12 [default = 0];
      case 12:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 97)) {
          _Internal::set_has_initial_time(&has_bits);
          initial_time_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
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

::PROTOBUF_NAMESPACE_ID::uint8* MicroenvironmentModel::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.input.MicroenvironmentModel)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated int32 grid_shape = 1;
  for (int i = 0, n = this->_internal_grid_shape_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt32ToArray(1, this->_internal_grid_shape(i), target);
  }

  cached_has_bits = _has_bits_[0];
  // required double grid_spacing = 2;
  if (cached_has_bits & 0x00000010u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(2, this->_internal_grid_spacing(), target);
  }

  // required .lm.types.BoundaryConditions boundaries = 3;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        3, _Internal::boundaries(this), target, stream);
  }

  // repeated int32 species_index = 4;
  for (int i = 0, n = this->_internal_species_index_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt32ToArray(4, this->_internal_species_index(i), target);
  }

  // repeated double diffusion_coefficients = 5;
  for (int i = 0, n = this->_internal_diffusion_coefficients_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(5, this->_internal_diffusion_coefficients(i), target);
  }

  // repeated .robertslab.pbuf.NDArray initial_concentrations = 6;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_initial_concentrations_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(6, this->_internal_initial_concentrations(i), target, stream);
  }

  // required double synchronization_timestep = 7;
  if (cached_has_bits & 0x00000020u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(7, this->_internal_synchronization_timestep(), target);
  }

  // optional uint32 number_cells = 8 [default = 0];
  if (cached_has_bits & 0x00000080u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt32ToArray(8, this->_internal_number_cells(), target);
  }

  // optional .robertslab.pbuf.NDArray cell_initial_species_counts = 9;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        9, _Internal::cell_initial_species_counts(this), target, stream);
  }

  // optional .robertslab.pbuf.NDArray cell_coordinates = 10;
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        10, _Internal::cell_coordinates(this), target, stream);
  }

  // optional .robertslab.pbuf.NDArray cell_volume = 11;
  if (cached_has_bits & 0x00000008u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        11, _Internal::cell_volume(this), target, stream);
  }

  // optional double initial_time = 12 [default = 0];
  if (cached_has_bits & 0x00000040u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(12, this->_internal_initial_time(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.input.MicroenvironmentModel)
  return target;
}

size_t MicroenvironmentModel::RequiredFieldsByteSizeFallback() const {
// @@protoc_insertion_point(required_fields_byte_size_fallback_start:lm.input.MicroenvironmentModel)
  size_t total_size = 0;

  if (_internal_has_boundaries()) {
    // required .lm.types.BoundaryConditions boundaries = 3;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
        *boundaries_);
  }

  if (_internal_has_grid_spacing()) {
    // required double grid_spacing = 2;
    total_size += 1 + 8;
  }

  if (_internal_has_synchronization_timestep()) {
    // required double synchronization_timestep = 7;
    total_size += 1 + 8;
  }

  return total_size;
}
size_t MicroenvironmentModel::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.input.MicroenvironmentModel)
  size_t total_size = 0;

  if (((_has_bits_[0] & 0x00000031) ^ 0x00000031) == 0) {  // All required fields are present.
    // required .lm.types.BoundaryConditions boundaries = 3;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
        *boundaries_);

    // required double grid_spacing = 2;
    total_size += 1 + 8;

    // required double synchronization_timestep = 7;
    total_size += 1 + 8;

  } else {
    total_size += RequiredFieldsByteSizeFallback();
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated int32 grid_shape = 1;
  {
    size_t data_size = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      Int32Size(this->grid_shape_);
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_grid_shape_size());
    total_size += data_size;
  }

  // repeated int32 species_index = 4;
  {
    size_t data_size = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      Int32Size(this->species_index_);
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_species_index_size());
    total_size += data_size;
  }

  // repeated double diffusion_coefficients = 5;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_diffusion_coefficients_size());
    size_t data_size = 8UL * count;
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_diffusion_coefficients_size());
    total_size += data_size;
  }

  // repeated .robertslab.pbuf.NDArray initial_concentrations = 6;
  total_size += 1UL * this->_internal_initial_concentrations_size();
  for (const auto& msg : this->initial_concentrations_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x0000000eu) {
    // optional .robertslab.pbuf.NDArray cell_initial_species_counts = 9;
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *cell_initial_species_counts_);
    }

    // optional .robertslab.pbuf.NDArray cell_coordinates = 10;
    if (cached_has_bits & 0x00000004u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *cell_coordinates_);
    }

    // optional .robertslab.pbuf.NDArray cell_volume = 11;
    if (cached_has_bits & 0x00000008u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *cell_volume_);
    }

  }
  if (cached_has_bits & 0x000000c0u) {
    // optional double initial_time = 12 [default = 0];
    if (cached_has_bits & 0x00000040u) {
      total_size += 1 + 8;
    }

    // optional uint32 number_cells = 8 [default = 0];
    if (cached_has_bits & 0x00000080u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt32Size(
          this->_internal_number_cells());
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

void MicroenvironmentModel::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.input.MicroenvironmentModel)
  GOOGLE_DCHECK_NE(&from, this);
  const MicroenvironmentModel* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<MicroenvironmentModel>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.input.MicroenvironmentModel)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.input.MicroenvironmentModel)
    MergeFrom(*source);
  }
}

void MicroenvironmentModel::MergeFrom(const MicroenvironmentModel& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.input.MicroenvironmentModel)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  grid_shape_.MergeFrom(from.grid_shape_);
  species_index_.MergeFrom(from.species_index_);
  diffusion_coefficients_.MergeFrom(from.diffusion_coefficients_);
  initial_concentrations_.MergeFrom(from.initial_concentrations_);
  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x000000ffu) {
    if (cached_has_bits & 0x00000001u) {
      _internal_mutable_boundaries()->::lm::types::BoundaryConditions::MergeFrom(from._internal_boundaries());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_mutable_cell_initial_species_counts()->::robertslab::pbuf::NDArray::MergeFrom(from._internal_cell_initial_species_counts());
    }
    if (cached_has_bits & 0x00000004u) {
      _internal_mutable_cell_coordinates()->::robertslab::pbuf::NDArray::MergeFrom(from._internal_cell_coordinates());
    }
    if (cached_has_bits & 0x00000008u) {
      _internal_mutable_cell_volume()->::robertslab::pbuf::NDArray::MergeFrom(from._internal_cell_volume());
    }
    if (cached_has_bits & 0x00000010u) {
      grid_spacing_ = from.grid_spacing_;
    }
    if (cached_has_bits & 0x00000020u) {
      synchronization_timestep_ = from.synchronization_timestep_;
    }
    if (cached_has_bits & 0x00000040u) {
      initial_time_ = from.initial_time_;
    }
    if (cached_has_bits & 0x00000080u) {
      number_cells_ = from.number_cells_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void MicroenvironmentModel::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.input.MicroenvironmentModel)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void MicroenvironmentModel::CopyFrom(const MicroenvironmentModel& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.input.MicroenvironmentModel)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool MicroenvironmentModel::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  if (!::PROTOBUF_NAMESPACE_ID::internal::AllAreInitialized(initial_concentrations_)) return false;
  if (_internal_has_boundaries()) {
    if (!boundaries_->IsInitialized()) return false;
  }
  if (_internal_has_cell_initial_species_counts()) {
    if (!cell_initial_species_counts_->IsInitialized()) return false;
  }
  if (_internal_has_cell_coordinates()) {
    if (!cell_coordinates_->IsInitialized()) return false;
  }
  if (_internal_has_cell_volume()) {
    if (!cell_volume_->IsInitialized()) return false;
  }
  return true;
}

void MicroenvironmentModel::InternalSwap(MicroenvironmentModel* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  grid_shape_.InternalSwap(&other->grid_shape_);
  species_index_.InternalSwap(&other->species_index_);
  diffusion_coefficients_.InternalSwap(&other->diffusion_coefficients_);
  initial_concentrations_.InternalSwap(&other->initial_concentrations_);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(MicroenvironmentModel, number_cells_)
      + sizeof(MicroenvironmentModel::number_cells_)
      - PROTOBUF_FIELD_OFFSET(MicroenvironmentModel, boundaries_)>(
          reinterpret_cast<char*>(&boundaries_),
          reinterpret_cast<char*>(&other->boundaries_));
}

::PROTOBUF_NAMESPACE_ID::Metadata MicroenvironmentModel::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace input
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::input::MicroenvironmentModel* Arena::CreateMaybeMessage< ::lm::input::MicroenvironmentModel >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::input::MicroenvironmentModel >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>