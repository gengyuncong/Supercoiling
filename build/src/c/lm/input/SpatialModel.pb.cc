// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/input/SpatialModel.proto

#include "lm/input/SpatialModel.pb.h"

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
extern PROTOBUF_INTERNAL_EXPORT_lm_2finput_2fSpatialModel_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_SpatialModel_Shape_lm_2finput_2fSpatialModel_2eproto;
namespace lm {
namespace input {
class SpatialModel_ShapeDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<SpatialModel_Shape> _instance;
} _SpatialModel_Shape_default_instance_;
class SpatialModelDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<SpatialModel> _instance;
} _SpatialModel_default_instance_;
}  // namespace input
}  // namespace lm
static void InitDefaultsscc_info_SpatialModel_lm_2finput_2fSpatialModel_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::input::_SpatialModel_default_instance_;
    new (ptr) ::lm::input::SpatialModel();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::input::SpatialModel::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_SpatialModel_lm_2finput_2fSpatialModel_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_SpatialModel_lm_2finput_2fSpatialModel_2eproto}, {
      &scc_info_SpatialModel_Shape_lm_2finput_2fSpatialModel_2eproto.base,}};

static void InitDefaultsscc_info_SpatialModel_Shape_lm_2finput_2fSpatialModel_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::input::_SpatialModel_Shape_default_instance_;
    new (ptr) ::lm::input::SpatialModel_Shape();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::input::SpatialModel_Shape::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_SpatialModel_Shape_lm_2finput_2fSpatialModel_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_SpatialModel_Shape_lm_2finput_2fSpatialModel_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2finput_2fSpatialModel_2eproto[2];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2finput_2fSpatialModel_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2finput_2fSpatialModel_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2finput_2fSpatialModel_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::input::SpatialModel_Shape, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::input::SpatialModel_Shape, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::input::SpatialModel_Shape, shape_),
  PROTOBUF_FIELD_OFFSET(::lm::input::SpatialModel_Shape, site_type_),
  PROTOBUF_FIELD_OFFSET(::lm::input::SpatialModel_Shape, shape_parameter_),
  0,
  1,
  ~0u,
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::lm::input::SpatialModel, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::input::SpatialModel, region_),
  PROTOBUF_FIELD_OFFSET(::lm::input::SpatialModel, obstacle_),
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 8, sizeof(::lm::input::SpatialModel_Shape)},
  { 11, -1, sizeof(::lm::input::SpatialModel)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::input::_SpatialModel_Shape_default_instance_),
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::input::_SpatialModel_default_instance_),
};

const char descriptor_table_protodef_lm_2finput_2fSpatialModel_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\033lm/input/SpatialModel.proto\022\010lm.input\""
  "\264\001\n\014SpatialModel\022,\n\006region\030\001 \003(\0132\034.lm.in"
  "put.SpatialModel.Shape\022.\n\010obstacle\030\002 \003(\013"
  "2\034.lm.input.SpatialModel.Shape\032F\n\005Shape\022"
  "\r\n\005shape\030\001 \002(\r\022\021\n\tsite_type\030\002 \002(\r\022\033\n\017sha"
  "pe_parameter\030\003 \003(\001B\002\020\001"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2finput_2fSpatialModel_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2finput_2fSpatialModel_2eproto_sccs[2] = {
  &scc_info_SpatialModel_lm_2finput_2fSpatialModel_2eproto.base,
  &scc_info_SpatialModel_Shape_lm_2finput_2fSpatialModel_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2finput_2fSpatialModel_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2finput_2fSpatialModel_2eproto = {
  false, false, descriptor_table_protodef_lm_2finput_2fSpatialModel_2eproto, "lm/input/SpatialModel.proto", 222,
  &descriptor_table_lm_2finput_2fSpatialModel_2eproto_once, descriptor_table_lm_2finput_2fSpatialModel_2eproto_sccs, descriptor_table_lm_2finput_2fSpatialModel_2eproto_deps, 2, 0,
  schemas, file_default_instances, TableStruct_lm_2finput_2fSpatialModel_2eproto::offsets,
  file_level_metadata_lm_2finput_2fSpatialModel_2eproto, 2, file_level_enum_descriptors_lm_2finput_2fSpatialModel_2eproto, file_level_service_descriptors_lm_2finput_2fSpatialModel_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2finput_2fSpatialModel_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2finput_2fSpatialModel_2eproto)), true);
namespace lm {
namespace input {

// ===================================================================

void SpatialModel_Shape::InitAsDefaultInstance() {
}
class SpatialModel_Shape::_Internal {
 public:
  using HasBits = decltype(std::declval<SpatialModel_Shape>()._has_bits_);
  static void set_has_shape(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_site_type(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000003) ^ 0x00000003) != 0;
  }
};

SpatialModel_Shape::SpatialModel_Shape(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  shape_parameter_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.input.SpatialModel.Shape)
}
SpatialModel_Shape::SpatialModel_Shape(const SpatialModel_Shape& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_),
      shape_parameter_(from.shape_parameter_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&shape_, &from.shape_,
    static_cast<size_t>(reinterpret_cast<char*>(&site_type_) -
    reinterpret_cast<char*>(&shape_)) + sizeof(site_type_));
  // @@protoc_insertion_point(copy_constructor:lm.input.SpatialModel.Shape)
}

void SpatialModel_Shape::SharedCtor() {
  ::memset(&shape_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&site_type_) -
      reinterpret_cast<char*>(&shape_)) + sizeof(site_type_));
}

SpatialModel_Shape::~SpatialModel_Shape() {
  // @@protoc_insertion_point(destructor:lm.input.SpatialModel.Shape)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void SpatialModel_Shape::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void SpatialModel_Shape::ArenaDtor(void* object) {
  SpatialModel_Shape* _this = reinterpret_cast< SpatialModel_Shape* >(object);
  (void)_this;
}
void SpatialModel_Shape::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void SpatialModel_Shape::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const SpatialModel_Shape& SpatialModel_Shape::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_SpatialModel_Shape_lm_2finput_2fSpatialModel_2eproto.base);
  return *internal_default_instance();
}


void SpatialModel_Shape::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.input.SpatialModel.Shape)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  shape_parameter_.Clear();
  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    ::memset(&shape_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&site_type_) -
        reinterpret_cast<char*>(&shape_)) + sizeof(site_type_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* SpatialModel_Shape::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // required uint32 shape = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          _Internal::set_has_shape(&has_bits);
          shape_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // required uint32 site_type = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 16)) {
          _Internal::set_has_site_type(&has_bits);
          site_type_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated double shape_parameter = 3 [packed = true];
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 26)) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_shape_parameter(), ptr, ctx);
          CHK_(ptr);
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 25) {
          _internal_add_shape_parameter(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
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

::PROTOBUF_NAMESPACE_ID::uint8* SpatialModel_Shape::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.input.SpatialModel.Shape)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required uint32 shape = 1;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt32ToArray(1, this->_internal_shape(), target);
  }

  // required uint32 site_type = 2;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt32ToArray(2, this->_internal_site_type(), target);
  }

  // repeated double shape_parameter = 3 [packed = true];
  if (this->_internal_shape_parameter_size() > 0) {
    target = stream->WriteFixedPacked(3, _internal_shape_parameter(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.input.SpatialModel.Shape)
  return target;
}

size_t SpatialModel_Shape::RequiredFieldsByteSizeFallback() const {
// @@protoc_insertion_point(required_fields_byte_size_fallback_start:lm.input.SpatialModel.Shape)
  size_t total_size = 0;

  if (_internal_has_shape()) {
    // required uint32 shape = 1;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt32Size(
        this->_internal_shape());
  }

  if (_internal_has_site_type()) {
    // required uint32 site_type = 2;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt32Size(
        this->_internal_site_type());
  }

  return total_size;
}
size_t SpatialModel_Shape::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.input.SpatialModel.Shape)
  size_t total_size = 0;

  if (((_has_bits_[0] & 0x00000003) ^ 0x00000003) == 0) {  // All required fields are present.
    // required uint32 shape = 1;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt32Size(
        this->_internal_shape());

    // required uint32 site_type = 2;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt32Size(
        this->_internal_site_type());

  } else {
    total_size += RequiredFieldsByteSizeFallback();
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated double shape_parameter = 3 [packed = true];
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_shape_parameter_size());
    size_t data_size = 8UL * count;
    if (data_size > 0) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
            static_cast<::PROTOBUF_NAMESPACE_ID::int32>(data_size));
    }
    int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(data_size);
    _shape_parameter_cached_byte_size_.store(cached_size,
                                    std::memory_order_relaxed);
    total_size += data_size;
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void SpatialModel_Shape::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.input.SpatialModel.Shape)
  GOOGLE_DCHECK_NE(&from, this);
  const SpatialModel_Shape* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<SpatialModel_Shape>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.input.SpatialModel.Shape)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.input.SpatialModel.Shape)
    MergeFrom(*source);
  }
}

void SpatialModel_Shape::MergeFrom(const SpatialModel_Shape& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.input.SpatialModel.Shape)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  shape_parameter_.MergeFrom(from.shape_parameter_);
  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      shape_ = from.shape_;
    }
    if (cached_has_bits & 0x00000002u) {
      site_type_ = from.site_type_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void SpatialModel_Shape::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.input.SpatialModel.Shape)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void SpatialModel_Shape::CopyFrom(const SpatialModel_Shape& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.input.SpatialModel.Shape)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool SpatialModel_Shape::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  return true;
}

void SpatialModel_Shape::InternalSwap(SpatialModel_Shape* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  shape_parameter_.InternalSwap(&other->shape_parameter_);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(SpatialModel_Shape, site_type_)
      + sizeof(SpatialModel_Shape::site_type_)
      - PROTOBUF_FIELD_OFFSET(SpatialModel_Shape, shape_)>(
          reinterpret_cast<char*>(&shape_),
          reinterpret_cast<char*>(&other->shape_));
}

::PROTOBUF_NAMESPACE_ID::Metadata SpatialModel_Shape::GetMetadata() const {
  return GetMetadataStatic();
}


// ===================================================================

void SpatialModel::InitAsDefaultInstance() {
}
class SpatialModel::_Internal {
 public:
};

SpatialModel::SpatialModel(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  region_(arena),
  obstacle_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.input.SpatialModel)
}
SpatialModel::SpatialModel(const SpatialModel& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      region_(from.region_),
      obstacle_(from.obstacle_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:lm.input.SpatialModel)
}

void SpatialModel::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_SpatialModel_lm_2finput_2fSpatialModel_2eproto.base);
}

SpatialModel::~SpatialModel() {
  // @@protoc_insertion_point(destructor:lm.input.SpatialModel)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void SpatialModel::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void SpatialModel::ArenaDtor(void* object) {
  SpatialModel* _this = reinterpret_cast< SpatialModel* >(object);
  (void)_this;
}
void SpatialModel::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void SpatialModel::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const SpatialModel& SpatialModel::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_SpatialModel_lm_2finput_2fSpatialModel_2eproto.base);
  return *internal_default_instance();
}


void SpatialModel::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.input.SpatialModel)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  region_.Clear();
  obstacle_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* SpatialModel::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated .lm.input.SpatialModel.Shape region = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_region(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<10>(ptr));
        } else goto handle_unusual;
        continue;
      // repeated .lm.input.SpatialModel.Shape obstacle = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 18)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_obstacle(), ptr);
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
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* SpatialModel::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.input.SpatialModel)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated .lm.input.SpatialModel.Shape region = 1;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_region_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(1, this->_internal_region(i), target, stream);
  }

  // repeated .lm.input.SpatialModel.Shape obstacle = 2;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_obstacle_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(2, this->_internal_obstacle(i), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.input.SpatialModel)
  return target;
}

size_t SpatialModel::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.input.SpatialModel)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .lm.input.SpatialModel.Shape region = 1;
  total_size += 1UL * this->_internal_region_size();
  for (const auto& msg : this->region_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  // repeated .lm.input.SpatialModel.Shape obstacle = 2;
  total_size += 1UL * this->_internal_obstacle_size();
  for (const auto& msg : this->obstacle_) {
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

void SpatialModel::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.input.SpatialModel)
  GOOGLE_DCHECK_NE(&from, this);
  const SpatialModel* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<SpatialModel>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.input.SpatialModel)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.input.SpatialModel)
    MergeFrom(*source);
  }
}

void SpatialModel::MergeFrom(const SpatialModel& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.input.SpatialModel)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  region_.MergeFrom(from.region_);
  obstacle_.MergeFrom(from.obstacle_);
}

void SpatialModel::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.input.SpatialModel)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void SpatialModel::CopyFrom(const SpatialModel& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.input.SpatialModel)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool SpatialModel::IsInitialized() const {
  if (!::PROTOBUF_NAMESPACE_ID::internal::AllAreInitialized(region_)) return false;
  if (!::PROTOBUF_NAMESPACE_ID::internal::AllAreInitialized(obstacle_)) return false;
  return true;
}

void SpatialModel::InternalSwap(SpatialModel* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  region_.InternalSwap(&other->region_);
  obstacle_.InternalSwap(&other->obstacle_);
}

::PROTOBUF_NAMESPACE_ID::Metadata SpatialModel::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace input
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::input::SpatialModel_Shape* Arena::CreateMaybeMessage< ::lm::input::SpatialModel_Shape >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::input::SpatialModel_Shape >(arena);
}
template<> PROTOBUF_NOINLINE ::lm::input::SpatialModel* Arena::CreateMaybeMessage< ::lm::input::SpatialModel >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::input::SpatialModel >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
