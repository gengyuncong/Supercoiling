// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/message/ProcessWorkUnitOutput.proto

#include "lm/message/ProcessWorkUnitOutput.pb.h"

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
extern PROTOBUF_INTERNAL_EXPORT_lm_2fmessage_2fWorkUnitOutput_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<9> scc_info_WorkUnitOutput_lm_2fmessage_2fWorkUnitOutput_2eproto;
namespace lm {
namespace message {
class ProcessWorkUnitOutputDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<ProcessWorkUnitOutput> _instance;
} _ProcessWorkUnitOutput_default_instance_;
}  // namespace message
}  // namespace lm
static void InitDefaultsscc_info_ProcessWorkUnitOutput_lm_2fmessage_2fProcessWorkUnitOutput_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::message::_ProcessWorkUnitOutput_default_instance_;
    new (ptr) ::lm::message::ProcessWorkUnitOutput();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::message::ProcessWorkUnitOutput::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_ProcessWorkUnitOutput_lm_2fmessage_2fProcessWorkUnitOutput_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_ProcessWorkUnitOutput_lm_2fmessage_2fProcessWorkUnitOutput_2eproto}, {
      &scc_info_WorkUnitOutput_lm_2fmessage_2fWorkUnitOutput_2eproto.base,}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2fmessage_2fProcessWorkUnitOutput_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2fmessage_2fProcessWorkUnitOutput_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2fmessage_2fProcessWorkUnitOutput_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2fmessage_2fProcessWorkUnitOutput_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::message::ProcessWorkUnitOutput, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::message::ProcessWorkUnitOutput, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::message::ProcessWorkUnitOutput, work_unit_id_),
  PROTOBUF_FIELD_OFFSET(::lm::message::ProcessWorkUnitOutput, part_output_),
  0,
  ~0u,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 7, sizeof(::lm::message::ProcessWorkUnitOutput)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::message::_ProcessWorkUnitOutput_default_instance_),
};

const char descriptor_table_protodef_lm_2fmessage_2fProcessWorkUnitOutput_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n&lm/message/ProcessWorkUnitOutput.proto"
  "\022\nlm.message\032\037lm/message/WorkUnitOutput."
  "proto\"^\n\025ProcessWorkUnitOutput\022\024\n\014work_u"
  "nit_id\030\001 \002(\003\022/\n\013part_output\030\002 \003(\0132\032.lm.m"
  "essage.WorkUnitOutput"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2fmessage_2fProcessWorkUnitOutput_2eproto_deps[1] = {
  &::descriptor_table_lm_2fmessage_2fWorkUnitOutput_2eproto,
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2fmessage_2fProcessWorkUnitOutput_2eproto_sccs[1] = {
  &scc_info_ProcessWorkUnitOutput_lm_2fmessage_2fProcessWorkUnitOutput_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2fmessage_2fProcessWorkUnitOutput_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fmessage_2fProcessWorkUnitOutput_2eproto = {
  false, false, descriptor_table_protodef_lm_2fmessage_2fProcessWorkUnitOutput_2eproto, "lm/message/ProcessWorkUnitOutput.proto", 181,
  &descriptor_table_lm_2fmessage_2fProcessWorkUnitOutput_2eproto_once, descriptor_table_lm_2fmessage_2fProcessWorkUnitOutput_2eproto_sccs, descriptor_table_lm_2fmessage_2fProcessWorkUnitOutput_2eproto_deps, 1, 1,
  schemas, file_default_instances, TableStruct_lm_2fmessage_2fProcessWorkUnitOutput_2eproto::offsets,
  file_level_metadata_lm_2fmessage_2fProcessWorkUnitOutput_2eproto, 1, file_level_enum_descriptors_lm_2fmessage_2fProcessWorkUnitOutput_2eproto, file_level_service_descriptors_lm_2fmessage_2fProcessWorkUnitOutput_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2fmessage_2fProcessWorkUnitOutput_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2fmessage_2fProcessWorkUnitOutput_2eproto)), true);
namespace lm {
namespace message {

// ===================================================================

void ProcessWorkUnitOutput::InitAsDefaultInstance() {
}
class ProcessWorkUnitOutput::_Internal {
 public:
  using HasBits = decltype(std::declval<ProcessWorkUnitOutput>()._has_bits_);
  static void set_has_work_unit_id(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000001) ^ 0x00000001) != 0;
  }
};

void ProcessWorkUnitOutput::clear_part_output() {
  part_output_.Clear();
}
ProcessWorkUnitOutput::ProcessWorkUnitOutput(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  part_output_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.message.ProcessWorkUnitOutput)
}
ProcessWorkUnitOutput::ProcessWorkUnitOutput(const ProcessWorkUnitOutput& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_),
      part_output_(from.part_output_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  work_unit_id_ = from.work_unit_id_;
  // @@protoc_insertion_point(copy_constructor:lm.message.ProcessWorkUnitOutput)
}

void ProcessWorkUnitOutput::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_ProcessWorkUnitOutput_lm_2fmessage_2fProcessWorkUnitOutput_2eproto.base);
  work_unit_id_ = PROTOBUF_LONGLONG(0);
}

ProcessWorkUnitOutput::~ProcessWorkUnitOutput() {
  // @@protoc_insertion_point(destructor:lm.message.ProcessWorkUnitOutput)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void ProcessWorkUnitOutput::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void ProcessWorkUnitOutput::ArenaDtor(void* object) {
  ProcessWorkUnitOutput* _this = reinterpret_cast< ProcessWorkUnitOutput* >(object);
  (void)_this;
}
void ProcessWorkUnitOutput::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void ProcessWorkUnitOutput::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const ProcessWorkUnitOutput& ProcessWorkUnitOutput::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_ProcessWorkUnitOutput_lm_2fmessage_2fProcessWorkUnitOutput_2eproto.base);
  return *internal_default_instance();
}


void ProcessWorkUnitOutput::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.message.ProcessWorkUnitOutput)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  part_output_.Clear();
  work_unit_id_ = PROTOBUF_LONGLONG(0);
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* ProcessWorkUnitOutput::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // required int64 work_unit_id = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          _Internal::set_has_work_unit_id(&has_bits);
          work_unit_id_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated .lm.message.WorkUnitOutput part_output = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 18)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_part_output(), ptr);
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

::PROTOBUF_NAMESPACE_ID::uint8* ProcessWorkUnitOutput::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.message.ProcessWorkUnitOutput)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required int64 work_unit_id = 1;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt64ToArray(1, this->_internal_work_unit_id(), target);
  }

  // repeated .lm.message.WorkUnitOutput part_output = 2;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_part_output_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(2, this->_internal_part_output(i), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.message.ProcessWorkUnitOutput)
  return target;
}

size_t ProcessWorkUnitOutput::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.message.ProcessWorkUnitOutput)
  size_t total_size = 0;

  // required int64 work_unit_id = 1;
  if (_internal_has_work_unit_id()) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int64Size(
        this->_internal_work_unit_id());
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .lm.message.WorkUnitOutput part_output = 2;
  total_size += 1UL * this->_internal_part_output_size();
  for (const auto& msg : this->part_output_) {
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

void ProcessWorkUnitOutput::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.message.ProcessWorkUnitOutput)
  GOOGLE_DCHECK_NE(&from, this);
  const ProcessWorkUnitOutput* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<ProcessWorkUnitOutput>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.message.ProcessWorkUnitOutput)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.message.ProcessWorkUnitOutput)
    MergeFrom(*source);
  }
}

void ProcessWorkUnitOutput::MergeFrom(const ProcessWorkUnitOutput& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.message.ProcessWorkUnitOutput)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  part_output_.MergeFrom(from.part_output_);
  if (from._internal_has_work_unit_id()) {
    _internal_set_work_unit_id(from._internal_work_unit_id());
  }
}

void ProcessWorkUnitOutput::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.message.ProcessWorkUnitOutput)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void ProcessWorkUnitOutput::CopyFrom(const ProcessWorkUnitOutput& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.message.ProcessWorkUnitOutput)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool ProcessWorkUnitOutput::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  if (!::PROTOBUF_NAMESPACE_ID::internal::AllAreInitialized(part_output_)) return false;
  return true;
}

void ProcessWorkUnitOutput::InternalSwap(ProcessWorkUnitOutput* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  part_output_.InternalSwap(&other->part_output_);
  swap(work_unit_id_, other->work_unit_id_);
}

::PROTOBUF_NAMESPACE_ID::Metadata ProcessWorkUnitOutput::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace message
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::message::ProcessWorkUnitOutput* Arena::CreateMaybeMessage< ::lm::message::ProcessWorkUnitOutput >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::message::ProcessWorkUnitOutput >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>