// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/io/DiffusionPDEState.proto

#include "lm/io/DiffusionPDEState.pb.h"

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
namespace io {
class DiffusionPDEStateDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<DiffusionPDEState> _instance;
} _DiffusionPDEState_default_instance_;
}  // namespace io
}  // namespace lm
static void InitDefaultsscc_info_DiffusionPDEState_lm_2fio_2fDiffusionPDEState_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::io::_DiffusionPDEState_default_instance_;
    new (ptr) ::lm::io::DiffusionPDEState();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::io::DiffusionPDEState::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_DiffusionPDEState_lm_2fio_2fDiffusionPDEState_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_DiffusionPDEState_lm_2fio_2fDiffusionPDEState_2eproto}, {
      &scc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto.base,}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2fio_2fDiffusionPDEState_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2fio_2fDiffusionPDEState_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2fio_2fDiffusionPDEState_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2fio_2fDiffusionPDEState_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::io::DiffusionPDEState, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::io::DiffusionPDEState, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::io::DiffusionPDEState, time_),
  PROTOBUF_FIELD_OFFSET(::lm::io::DiffusionPDEState, concentrations_),
  0,
  ~0u,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 7, sizeof(::lm::io::DiffusionPDEState)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::io::_DiffusionPDEState_default_instance_),
};

const char descriptor_table_protodef_lm_2fio_2fDiffusionPDEState_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\035lm/io/DiffusionPDEState.proto\022\005lm.io\032\035"
  "robertslab/pbuf/NDArray.proto\"S\n\021Diffusi"
  "onPDEState\022\014\n\004time\030\001 \002(\001\0220\n\016concentratio"
  "ns\030\002 \003(\0132\030.robertslab.pbuf.NDArray"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2fio_2fDiffusionPDEState_2eproto_deps[1] = {
  &::descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto,
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2fio_2fDiffusionPDEState_2eproto_sccs[1] = {
  &scc_info_DiffusionPDEState_lm_2fio_2fDiffusionPDEState_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2fio_2fDiffusionPDEState_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fio_2fDiffusionPDEState_2eproto = {
  false, false, descriptor_table_protodef_lm_2fio_2fDiffusionPDEState_2eproto, "lm/io/DiffusionPDEState.proto", 154,
  &descriptor_table_lm_2fio_2fDiffusionPDEState_2eproto_once, descriptor_table_lm_2fio_2fDiffusionPDEState_2eproto_sccs, descriptor_table_lm_2fio_2fDiffusionPDEState_2eproto_deps, 1, 1,
  schemas, file_default_instances, TableStruct_lm_2fio_2fDiffusionPDEState_2eproto::offsets,
  file_level_metadata_lm_2fio_2fDiffusionPDEState_2eproto, 1, file_level_enum_descriptors_lm_2fio_2fDiffusionPDEState_2eproto, file_level_service_descriptors_lm_2fio_2fDiffusionPDEState_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2fio_2fDiffusionPDEState_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2fio_2fDiffusionPDEState_2eproto)), true);
namespace lm {
namespace io {

// ===================================================================

void DiffusionPDEState::InitAsDefaultInstance() {
}
class DiffusionPDEState::_Internal {
 public:
  using HasBits = decltype(std::declval<DiffusionPDEState>()._has_bits_);
  static void set_has_time(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000001) ^ 0x00000001) != 0;
  }
};

void DiffusionPDEState::clear_concentrations() {
  concentrations_.Clear();
}
DiffusionPDEState::DiffusionPDEState(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  concentrations_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.io.DiffusionPDEState)
}
DiffusionPDEState::DiffusionPDEState(const DiffusionPDEState& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_),
      concentrations_(from.concentrations_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  time_ = from.time_;
  // @@protoc_insertion_point(copy_constructor:lm.io.DiffusionPDEState)
}

void DiffusionPDEState::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_DiffusionPDEState_lm_2fio_2fDiffusionPDEState_2eproto.base);
  time_ = 0;
}

DiffusionPDEState::~DiffusionPDEState() {
  // @@protoc_insertion_point(destructor:lm.io.DiffusionPDEState)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void DiffusionPDEState::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void DiffusionPDEState::ArenaDtor(void* object) {
  DiffusionPDEState* _this = reinterpret_cast< DiffusionPDEState* >(object);
  (void)_this;
}
void DiffusionPDEState::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void DiffusionPDEState::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const DiffusionPDEState& DiffusionPDEState::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_DiffusionPDEState_lm_2fio_2fDiffusionPDEState_2eproto.base);
  return *internal_default_instance();
}


void DiffusionPDEState::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.io.DiffusionPDEState)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  concentrations_.Clear();
  time_ = 0;
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* DiffusionPDEState::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // required double time = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 9)) {
          _Internal::set_has_time(&has_bits);
          time_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else goto handle_unusual;
        continue;
      // repeated .robertslab.pbuf.NDArray concentrations = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 18)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_concentrations(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<18>(ptr));
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

::PROTOBUF_NAMESPACE_ID::uint8* DiffusionPDEState::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.io.DiffusionPDEState)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required double time = 1;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(1, this->_internal_time(), target);
  }

  // repeated .robertslab.pbuf.NDArray concentrations = 2;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_concentrations_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(2, this->_internal_concentrations(i), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.io.DiffusionPDEState)
  return target;
}

size_t DiffusionPDEState::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.io.DiffusionPDEState)
  size_t total_size = 0;

  // required double time = 1;
  if (_internal_has_time()) {
    total_size += 1 + 8;
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .robertslab.pbuf.NDArray concentrations = 2;
  total_size += 1UL * this->_internal_concentrations_size();
  for (const auto& msg : this->concentrations_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void DiffusionPDEState::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.io.DiffusionPDEState)
  GOOGLE_DCHECK_NE(&from, this);
  const DiffusionPDEState* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<DiffusionPDEState>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.io.DiffusionPDEState)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.io.DiffusionPDEState)
    MergeFrom(*source);
  }
}

void DiffusionPDEState::MergeFrom(const DiffusionPDEState& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.io.DiffusionPDEState)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  concentrations_.MergeFrom(from.concentrations_);
  if (from._internal_has_time()) {
    _internal_set_time(from._internal_time());
  }
}

void DiffusionPDEState::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.io.DiffusionPDEState)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void DiffusionPDEState::CopyFrom(const DiffusionPDEState& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.io.DiffusionPDEState)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool DiffusionPDEState::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  if (!::PROTOBUF_NAMESPACE_ID::internal::AllAreInitialized(concentrations_)) return false;
  return true;
}

void DiffusionPDEState::InternalSwap(DiffusionPDEState* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  concentrations_.InternalSwap(&other->concentrations_);
  swap(time_, other->time_);
}

::PROTOBUF_NAMESPACE_ID::Metadata DiffusionPDEState::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace io
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::io::DiffusionPDEState* Arena::CreateMaybeMessage< ::lm::io::DiffusionPDEState >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::io::DiffusionPDEState >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
