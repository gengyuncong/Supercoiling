// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/message/PingTarget.proto

#include "lm/message/PingTarget.pb.h"

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
namespace lm {
namespace message {
class PingTargetDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<PingTarget> _instance;
} _PingTarget_default_instance_;
}  // namespace message
}  // namespace lm
static void InitDefaultsscc_info_PingTarget_lm_2fmessage_2fPingTarget_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::message::_PingTarget_default_instance_;
    new (ptr) ::lm::message::PingTarget();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::message::PingTarget::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_PingTarget_lm_2fmessage_2fPingTarget_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_PingTarget_lm_2fmessage_2fPingTarget_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2fmessage_2fPingTarget_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2fmessage_2fPingTarget_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2fmessage_2fPingTarget_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2fmessage_2fPingTarget_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::message::PingTarget, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::message::PingTarget, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::message::PingTarget, id_),
  0,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 6, sizeof(::lm::message::PingTarget)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::message::_PingTarget_default_instance_),
};

const char descriptor_table_protodef_lm_2fmessage_2fPingTarget_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\033lm/message/PingTarget.proto\022\nlm.messag"
  "e\"\030\n\nPingTarget\022\n\n\002id\030\001 \002(\005"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2fmessage_2fPingTarget_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2fmessage_2fPingTarget_2eproto_sccs[1] = {
  &scc_info_PingTarget_lm_2fmessage_2fPingTarget_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2fmessage_2fPingTarget_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fmessage_2fPingTarget_2eproto = {
  false, false, descriptor_table_protodef_lm_2fmessage_2fPingTarget_2eproto, "lm/message/PingTarget.proto", 67,
  &descriptor_table_lm_2fmessage_2fPingTarget_2eproto_once, descriptor_table_lm_2fmessage_2fPingTarget_2eproto_sccs, descriptor_table_lm_2fmessage_2fPingTarget_2eproto_deps, 1, 0,
  schemas, file_default_instances, TableStruct_lm_2fmessage_2fPingTarget_2eproto::offsets,
  file_level_metadata_lm_2fmessage_2fPingTarget_2eproto, 1, file_level_enum_descriptors_lm_2fmessage_2fPingTarget_2eproto, file_level_service_descriptors_lm_2fmessage_2fPingTarget_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2fmessage_2fPingTarget_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2fmessage_2fPingTarget_2eproto)), true);
namespace lm {
namespace message {

// ===================================================================

void PingTarget::InitAsDefaultInstance() {
}
class PingTarget::_Internal {
 public:
  using HasBits = decltype(std::declval<PingTarget>()._has_bits_);
  static void set_has_id(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000001) ^ 0x00000001) != 0;
  }
};

PingTarget::PingTarget(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.message.PingTarget)
}
PingTarget::PingTarget(const PingTarget& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  id_ = from.id_;
  // @@protoc_insertion_point(copy_constructor:lm.message.PingTarget)
}

void PingTarget::SharedCtor() {
  id_ = 0;
}

PingTarget::~PingTarget() {
  // @@protoc_insertion_point(destructor:lm.message.PingTarget)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void PingTarget::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void PingTarget::ArenaDtor(void* object) {
  PingTarget* _this = reinterpret_cast< PingTarget* >(object);
  (void)_this;
}
void PingTarget::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void PingTarget::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const PingTarget& PingTarget::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_PingTarget_lm_2fmessage_2fPingTarget_2eproto.base);
  return *internal_default_instance();
}


void PingTarget::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.message.PingTarget)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  id_ = 0;
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* PingTarget::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // required int32 id = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          _Internal::set_has_id(&has_bits);
          id_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
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

::PROTOBUF_NAMESPACE_ID::uint8* PingTarget::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.message.PingTarget)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required int32 id = 1;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt32ToArray(1, this->_internal_id(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.message.PingTarget)
  return target;
}

size_t PingTarget::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.message.PingTarget)
  size_t total_size = 0;

  // required int32 id = 1;
  if (_internal_has_id()) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
        this->_internal_id());
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void PingTarget::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.message.PingTarget)
  GOOGLE_DCHECK_NE(&from, this);
  const PingTarget* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<PingTarget>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.message.PingTarget)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.message.PingTarget)
    MergeFrom(*source);
  }
}

void PingTarget::MergeFrom(const PingTarget& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.message.PingTarget)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  if (from._internal_has_id()) {
    _internal_set_id(from._internal_id());
  }
}

void PingTarget::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.message.PingTarget)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void PingTarget::CopyFrom(const PingTarget& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.message.PingTarget)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool PingTarget::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  return true;
}

void PingTarget::InternalSwap(PingTarget* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  swap(id_, other->id_);
}

::PROTOBUF_NAMESPACE_ID::Metadata PingTarget::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace message
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::message::PingTarget* Arena::CreateMaybeMessage< ::lm::message::PingTarget >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::message::PingTarget >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>