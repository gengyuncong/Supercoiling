// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/message/StartedOutputPerformanceSignaler.proto

#include "lm/message/StartedOutputPerformanceSignaler.pb.h"

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
extern PROTOBUF_INTERNAL_EXPORT_lm_2fmessage_2fEndpoint_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_Endpoint_lm_2fmessage_2fEndpoint_2eproto;
namespace lm {
namespace message {
class StartedOutputPerformanceSignalerDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<StartedOutputPerformanceSignaler> _instance;
} _StartedOutputPerformanceSignaler_default_instance_;
}  // namespace message
}  // namespace lm
static void InitDefaultsscc_info_StartedOutputPerformanceSignaler_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::message::_StartedOutputPerformanceSignaler_default_instance_;
    new (ptr) ::lm::message::StartedOutputPerformanceSignaler();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::message::StartedOutputPerformanceSignaler::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_StartedOutputPerformanceSignaler_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_StartedOutputPerformanceSignaler_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto}, {
      &scc_info_Endpoint_lm_2fmessage_2fEndpoint_2eproto.base,}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::message::StartedOutputPerformanceSignaler, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::message::StartedOutputPerformanceSignaler, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::message::StartedOutputPerformanceSignaler, address_),
  0,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 6, sizeof(::lm::message::StartedOutputPerformanceSignaler)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::message::_StartedOutputPerformanceSignaler_default_instance_),
};

const char descriptor_table_protodef_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n1lm/message/StartedOutputPerformanceSig"
  "naler.proto\022\nlm.message\032\031lm/message/Endp"
  "oint.proto\"I\n StartedOutputPerformanceSi"
  "gnaler\022%\n\007address\030\001 \002(\0132\024.lm.message.End"
  "point"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto_deps[1] = {
  &::descriptor_table_lm_2fmessage_2fEndpoint_2eproto,
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto_sccs[1] = {
  &scc_info_StartedOutputPerformanceSignaler_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto = {
  false, false, descriptor_table_protodef_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto, "lm/message/StartedOutputPerformanceSignaler.proto", 165,
  &descriptor_table_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto_once, descriptor_table_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto_sccs, descriptor_table_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto_deps, 1, 1,
  schemas, file_default_instances, TableStruct_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto::offsets,
  file_level_metadata_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto, 1, file_level_enum_descriptors_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto, file_level_service_descriptors_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto)), true);
namespace lm {
namespace message {

// ===================================================================

void StartedOutputPerformanceSignaler::InitAsDefaultInstance() {
  ::lm::message::_StartedOutputPerformanceSignaler_default_instance_._instance.get_mutable()->address_ = const_cast< ::lm::message::Endpoint*>(
      ::lm::message::Endpoint::internal_default_instance());
}
class StartedOutputPerformanceSignaler::_Internal {
 public:
  using HasBits = decltype(std::declval<StartedOutputPerformanceSignaler>()._has_bits_);
  static const ::lm::message::Endpoint& address(const StartedOutputPerformanceSignaler* msg);
  static void set_has_address(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000001) ^ 0x00000001) != 0;
  }
};

const ::lm::message::Endpoint&
StartedOutputPerformanceSignaler::_Internal::address(const StartedOutputPerformanceSignaler* msg) {
  return *msg->address_;
}
void StartedOutputPerformanceSignaler::clear_address() {
  if (address_ != nullptr) address_->Clear();
  _has_bits_[0] &= ~0x00000001u;
}
StartedOutputPerformanceSignaler::StartedOutputPerformanceSignaler(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.message.StartedOutputPerformanceSignaler)
}
StartedOutputPerformanceSignaler::StartedOutputPerformanceSignaler(const StartedOutputPerformanceSignaler& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  if (from._internal_has_address()) {
    address_ = new ::lm::message::Endpoint(*from.address_);
  } else {
    address_ = nullptr;
  }
  // @@protoc_insertion_point(copy_constructor:lm.message.StartedOutputPerformanceSignaler)
}

void StartedOutputPerformanceSignaler::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_StartedOutputPerformanceSignaler_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto.base);
  address_ = nullptr;
}

StartedOutputPerformanceSignaler::~StartedOutputPerformanceSignaler() {
  // @@protoc_insertion_point(destructor:lm.message.StartedOutputPerformanceSignaler)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void StartedOutputPerformanceSignaler::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
  if (this != internal_default_instance()) delete address_;
}

void StartedOutputPerformanceSignaler::ArenaDtor(void* object) {
  StartedOutputPerformanceSignaler* _this = reinterpret_cast< StartedOutputPerformanceSignaler* >(object);
  (void)_this;
}
void StartedOutputPerformanceSignaler::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void StartedOutputPerformanceSignaler::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const StartedOutputPerformanceSignaler& StartedOutputPerformanceSignaler::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_StartedOutputPerformanceSignaler_lm_2fmessage_2fStartedOutputPerformanceSignaler_2eproto.base);
  return *internal_default_instance();
}


void StartedOutputPerformanceSignaler::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.message.StartedOutputPerformanceSignaler)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000001u) {
    GOOGLE_DCHECK(address_ != nullptr);
    address_->Clear();
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* StartedOutputPerformanceSignaler::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // required .lm.message.Endpoint address = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10)) {
          ptr = ctx->ParseMessage(_internal_mutable_address(), ptr);
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

::PROTOBUF_NAMESPACE_ID::uint8* StartedOutputPerformanceSignaler::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.message.StartedOutputPerformanceSignaler)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required .lm.message.Endpoint address = 1;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        1, _Internal::address(this), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.message.StartedOutputPerformanceSignaler)
  return target;
}

size_t StartedOutputPerformanceSignaler::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.message.StartedOutputPerformanceSignaler)
  size_t total_size = 0;

  // required .lm.message.Endpoint address = 1;
  if (_internal_has_address()) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
        *address_);
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

void StartedOutputPerformanceSignaler::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.message.StartedOutputPerformanceSignaler)
  GOOGLE_DCHECK_NE(&from, this);
  const StartedOutputPerformanceSignaler* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<StartedOutputPerformanceSignaler>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.message.StartedOutputPerformanceSignaler)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.message.StartedOutputPerformanceSignaler)
    MergeFrom(*source);
  }
}

void StartedOutputPerformanceSignaler::MergeFrom(const StartedOutputPerformanceSignaler& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.message.StartedOutputPerformanceSignaler)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  if (from._internal_has_address()) {
    _internal_mutable_address()->::lm::message::Endpoint::MergeFrom(from._internal_address());
  }
}

void StartedOutputPerformanceSignaler::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.message.StartedOutputPerformanceSignaler)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void StartedOutputPerformanceSignaler::CopyFrom(const StartedOutputPerformanceSignaler& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.message.StartedOutputPerformanceSignaler)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool StartedOutputPerformanceSignaler::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  return true;
}

void StartedOutputPerformanceSignaler::InternalSwap(StartedOutputPerformanceSignaler* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  swap(address_, other->address_);
}

::PROTOBUF_NAMESPACE_ID::Metadata StartedOutputPerformanceSignaler::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace message
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::message::StartedOutputPerformanceSignaler* Arena::CreateMaybeMessage< ::lm::message::StartedOutputPerformanceSignaler >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::message::StartedOutputPerformanceSignaler >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
