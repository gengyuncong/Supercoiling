// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/types/Lattice.proto

#include "lm/types/Lattice.pb.h"

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
extern PROTOBUF_INTERNAL_EXPORT_robertslab_2fpbuf_2fNDArray_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto;
namespace lm {
namespace types {
class LatticeDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<Lattice> _instance;
} _Lattice_default_instance_;
}  // namespace types
}  // namespace lm
static void InitDefaultsscc_info_Lattice_lm_2ftypes_2fLattice_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::types::_Lattice_default_instance_;
    new (ptr) ::lm::types::Lattice();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::types::Lattice::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_Lattice_lm_2ftypes_2fLattice_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_Lattice_lm_2ftypes_2fLattice_2eproto}, {
      &scc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto.base,}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2ftypes_2fLattice_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2ftypes_2fLattice_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2ftypes_2fLattice_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2ftypes_2fLattice_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::types::Lattice, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::types::Lattice, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::types::Lattice, particles_),
  PROTOBUF_FIELD_OFFSET(::lm::types::Lattice, sites_),
  0,
  1,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 7, sizeof(::lm::types::Lattice)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::types::_Lattice_default_instance_),
};

const char descriptor_table_protodef_lm_2ftypes_2fLattice_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\026lm/types/Lattice.proto\022\010lm.types\032\035robe"
  "rtslab/pbuf/NDArray.proto\"_\n\007Lattice\022+\n\t"
  "particles\030\013 \002(\0132\030.robertslab.pbuf.NDArra"
  "y\022\'\n\005sites\030\014 \001(\0132\030.robertslab.pbuf.NDArr"
  "ay"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2ftypes_2fLattice_2eproto_deps[1] = {
  &::descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto,
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2ftypes_2fLattice_2eproto_sccs[1] = {
  &scc_info_Lattice_lm_2ftypes_2fLattice_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2ftypes_2fLattice_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2ftypes_2fLattice_2eproto = {
  false, false, descriptor_table_protodef_lm_2ftypes_2fLattice_2eproto, "lm/types/Lattice.proto", 162,
  &descriptor_table_lm_2ftypes_2fLattice_2eproto_once, descriptor_table_lm_2ftypes_2fLattice_2eproto_sccs, descriptor_table_lm_2ftypes_2fLattice_2eproto_deps, 1, 1,
  schemas, file_default_instances, TableStruct_lm_2ftypes_2fLattice_2eproto::offsets,
  file_level_metadata_lm_2ftypes_2fLattice_2eproto, 1, file_level_enum_descriptors_lm_2ftypes_2fLattice_2eproto, file_level_service_descriptors_lm_2ftypes_2fLattice_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2ftypes_2fLattice_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2ftypes_2fLattice_2eproto)), true);
namespace lm {
namespace types {

// ===================================================================

void Lattice::InitAsDefaultInstance() {
  ::lm::types::_Lattice_default_instance_._instance.get_mutable()->particles_ = const_cast< ::robertslab::pbuf::NDArray*>(
      ::robertslab::pbuf::NDArray::internal_default_instance());
  ::lm::types::_Lattice_default_instance_._instance.get_mutable()->sites_ = const_cast< ::robertslab::pbuf::NDArray*>(
      ::robertslab::pbuf::NDArray::internal_default_instance());
}
class Lattice::_Internal {
 public:
  using HasBits = decltype(std::declval<Lattice>()._has_bits_);
  static const ::robertslab::pbuf::NDArray& particles(const Lattice* msg);
  static void set_has_particles(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static const ::robertslab::pbuf::NDArray& sites(const Lattice* msg);
  static void set_has_sites(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000001) ^ 0x00000001) != 0;
  }
};

const ::robertslab::pbuf::NDArray&
Lattice::_Internal::particles(const Lattice* msg) {
  return *msg->particles_;
}
const ::robertslab::pbuf::NDArray&
Lattice::_Internal::sites(const Lattice* msg) {
  return *msg->sites_;
}
void Lattice::clear_particles() {
  if (particles_ != nullptr) particles_->Clear();
  _has_bits_[0] &= ~0x00000001u;
}
void Lattice::clear_sites() {
  if (sites_ != nullptr) sites_->Clear();
  _has_bits_[0] &= ~0x00000002u;
}
Lattice::Lattice(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.types.Lattice)
}
Lattice::Lattice(const Lattice& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  if (from._internal_has_particles()) {
    particles_ = new ::robertslab::pbuf::NDArray(*from.particles_);
  } else {
    particles_ = nullptr;
  }
  if (from._internal_has_sites()) {
    sites_ = new ::robertslab::pbuf::NDArray(*from.sites_);
  } else {
    sites_ = nullptr;
  }
  // @@protoc_insertion_point(copy_constructor:lm.types.Lattice)
}

void Lattice::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_Lattice_lm_2ftypes_2fLattice_2eproto.base);
  ::memset(&particles_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&sites_) -
      reinterpret_cast<char*>(&particles_)) + sizeof(sites_));
}

Lattice::~Lattice() {
  // @@protoc_insertion_point(destructor:lm.types.Lattice)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void Lattice::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
  if (this != internal_default_instance()) delete particles_;
  if (this != internal_default_instance()) delete sites_;
}

void Lattice::ArenaDtor(void* object) {
  Lattice* _this = reinterpret_cast< Lattice* >(object);
  (void)_this;
}
void Lattice::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void Lattice::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const Lattice& Lattice::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_Lattice_lm_2ftypes_2fLattice_2eproto.base);
  return *internal_default_instance();
}


void Lattice::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.types.Lattice)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      GOOGLE_DCHECK(particles_ != nullptr);
      particles_->Clear();
    }
    if (cached_has_bits & 0x00000002u) {
      GOOGLE_DCHECK(sites_ != nullptr);
      sites_->Clear();
    }
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* Lattice::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // required .robertslab.pbuf.NDArray particles = 11;
      case 11:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 90)) {
          ptr = ctx->ParseMessage(_internal_mutable_particles(), ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional .robertslab.pbuf.NDArray sites = 12;
      case 12:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 98)) {
          ptr = ctx->ParseMessage(_internal_mutable_sites(), ptr);
          CHK_(ptr);
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

::PROTOBUF_NAMESPACE_ID::uint8* Lattice::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.types.Lattice)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required .robertslab.pbuf.NDArray particles = 11;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        11, _Internal::particles(this), target, stream);
  }

  // optional .robertslab.pbuf.NDArray sites = 12;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        12, _Internal::sites(this), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.types.Lattice)
  return target;
}

size_t Lattice::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.types.Lattice)
  size_t total_size = 0;

  // required .robertslab.pbuf.NDArray particles = 11;
  if (_internal_has_particles()) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
        *particles_);
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // optional .robertslab.pbuf.NDArray sites = 12;
  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000002u) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
        *sites_);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void Lattice::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.types.Lattice)
  GOOGLE_DCHECK_NE(&from, this);
  const Lattice* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<Lattice>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.types.Lattice)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.types.Lattice)
    MergeFrom(*source);
  }
}

void Lattice::MergeFrom(const Lattice& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.types.Lattice)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      _internal_mutable_particles()->::robertslab::pbuf::NDArray::MergeFrom(from._internal_particles());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_mutable_sites()->::robertslab::pbuf::NDArray::MergeFrom(from._internal_sites());
    }
  }
}

void Lattice::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.types.Lattice)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void Lattice::CopyFrom(const Lattice& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.types.Lattice)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool Lattice::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  if (_internal_has_particles()) {
    if (!particles_->IsInitialized()) return false;
  }
  if (_internal_has_sites()) {
    if (!sites_->IsInitialized()) return false;
  }
  return true;
}

void Lattice::InternalSwap(Lattice* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(Lattice, sites_)
      + sizeof(Lattice::sites_)
      - PROTOBUF_FIELD_OFFSET(Lattice, particles_)>(
          reinterpret_cast<char*>(&particles_),
          reinterpret_cast<char*>(&other->particles_));
}

::PROTOBUF_NAMESPACE_ID::Metadata Lattice::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace types
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::types::Lattice* Arena::CreateMaybeMessage< ::lm::types::Lattice >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::types::Lattice >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>