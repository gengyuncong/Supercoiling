// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/io/ffpilot/FFPilotOutput.proto

#include "lm/io/ffpilot/FFPilotOutput.pb.h"

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
extern PROTOBUF_INTERNAL_EXPORT_lm_2fio_2fffpilot_2fFFPilotStageOutput_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_FFPilotStageOutput_lm_2fio_2fffpilot_2fFFPilotStageOutput_2eproto;
namespace lm {
namespace io {
namespace ffpilot {
class FFPilotOutputDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<FFPilotOutput> _instance;
} _FFPilotOutput_default_instance_;
}  // namespace ffpilot
}  // namespace io
}  // namespace lm
static void InitDefaultsscc_info_FFPilotOutput_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::io::ffpilot::_FFPilotOutput_default_instance_;
    new (ptr) ::lm::io::ffpilot::FFPilotOutput();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::io::ffpilot::FFPilotOutput::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_FFPilotOutput_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_FFPilotOutput_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto}, {
      &scc_info_FFPilotStageOutput_lm_2fio_2fffpilot_2fFFPilotStageOutput_2eproto.base,}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotOutput, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotOutput, stage_output_),
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, -1, sizeof(::lm::io::ffpilot::FFPilotOutput)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::io::ffpilot::_FFPilotOutput_default_instance_),
};

const char descriptor_table_protodef_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n!lm/io/ffpilot/FFPilotOutput.proto\022\rlm."
  "io.ffpilot\032&lm/io/ffpilot/FFPilotStageOu"
  "tput.proto\"H\n\rFFPilotOutput\0227\n\014stage_out"
  "put\030d \003(\0132!.lm.io.ffpilot.FFPilotStageOu"
  "tput"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto_deps[1] = {
  &::descriptor_table_lm_2fio_2fffpilot_2fFFPilotStageOutput_2eproto,
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto_sccs[1] = {
  &scc_info_FFPilotOutput_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto = {
  false, false, descriptor_table_protodef_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto, "lm/io/ffpilot/FFPilotOutput.proto", 164,
  &descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto_once, descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto_sccs, descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto_deps, 1, 1,
  schemas, file_default_instances, TableStruct_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto::offsets,
  file_level_metadata_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto, 1, file_level_enum_descriptors_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto, file_level_service_descriptors_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto)), true);
namespace lm {
namespace io {
namespace ffpilot {

// ===================================================================

void FFPilotOutput::InitAsDefaultInstance() {
}
class FFPilotOutput::_Internal {
 public:
};

void FFPilotOutput::clear_stage_output() {
  stage_output_.Clear();
}
FFPilotOutput::FFPilotOutput(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  stage_output_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.io.ffpilot.FFPilotOutput)
}
FFPilotOutput::FFPilotOutput(const FFPilotOutput& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      stage_output_(from.stage_output_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:lm.io.ffpilot.FFPilotOutput)
}

void FFPilotOutput::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_FFPilotOutput_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto.base);
}

FFPilotOutput::~FFPilotOutput() {
  // @@protoc_insertion_point(destructor:lm.io.ffpilot.FFPilotOutput)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void FFPilotOutput::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void FFPilotOutput::ArenaDtor(void* object) {
  FFPilotOutput* _this = reinterpret_cast< FFPilotOutput* >(object);
  (void)_this;
}
void FFPilotOutput::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void FFPilotOutput::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const FFPilotOutput& FFPilotOutput::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_FFPilotOutput_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto.base);
  return *internal_default_instance();
}


void FFPilotOutput::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.io.ffpilot.FFPilotOutput)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  stage_output_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* FFPilotOutput::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated .lm.io.ffpilot.FFPilotStageOutput stage_output = 100;
      case 100:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 34)) {
          ptr -= 2;
          do {
            ptr += 2;
            ptr = ctx->ParseMessage(_internal_add_stage_output(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<802>(ptr));
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

::PROTOBUF_NAMESPACE_ID::uint8* FFPilotOutput::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.io.ffpilot.FFPilotOutput)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated .lm.io.ffpilot.FFPilotStageOutput stage_output = 100;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_stage_output_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(100, this->_internal_stage_output(i), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.io.ffpilot.FFPilotOutput)
  return target;
}

size_t FFPilotOutput::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.io.ffpilot.FFPilotOutput)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .lm.io.ffpilot.FFPilotStageOutput stage_output = 100;
  total_size += 2UL * this->_internal_stage_output_size();
  for (const auto& msg : this->stage_output_) {
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

void FFPilotOutput::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.io.ffpilot.FFPilotOutput)
  GOOGLE_DCHECK_NE(&from, this);
  const FFPilotOutput* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<FFPilotOutput>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.io.ffpilot.FFPilotOutput)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.io.ffpilot.FFPilotOutput)
    MergeFrom(*source);
  }
}

void FFPilotOutput::MergeFrom(const FFPilotOutput& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.io.ffpilot.FFPilotOutput)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  stage_output_.MergeFrom(from.stage_output_);
}

void FFPilotOutput::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.io.ffpilot.FFPilotOutput)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void FFPilotOutput::CopyFrom(const FFPilotOutput& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.io.ffpilot.FFPilotOutput)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool FFPilotOutput::IsInitialized() const {
  if (!::PROTOBUF_NAMESPACE_ID::internal::AllAreInitialized(stage_output_)) return false;
  return true;
}

void FFPilotOutput::InternalSwap(FFPilotOutput* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  stage_output_.InternalSwap(&other->stage_output_);
}

::PROTOBUF_NAMESPACE_ID::Metadata FFPilotOutput::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace ffpilot
}  // namespace io
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::io::ffpilot::FFPilotOutput* Arena::CreateMaybeMessage< ::lm::io::ffpilot::FFPilotOutput >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::io::ffpilot::FFPilotOutput >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>