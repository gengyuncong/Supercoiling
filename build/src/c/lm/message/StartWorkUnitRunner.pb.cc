// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/message/StartWorkUnitRunner.proto

#include "lm/message/StartWorkUnitRunner.pb.h"

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
class StartWorkUnitRunnerDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<StartWorkUnitRunner> _instance;
} _StartWorkUnitRunner_default_instance_;
}  // namespace message
}  // namespace lm
static void InitDefaultsscc_info_StartWorkUnitRunner_lm_2fmessage_2fStartWorkUnitRunner_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::message::_StartWorkUnitRunner_default_instance_;
    new (ptr) ::lm::message::StartWorkUnitRunner();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::message::StartWorkUnitRunner::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_StartWorkUnitRunner_lm_2fmessage_2fStartWorkUnitRunner_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_StartWorkUnitRunner_lm_2fmessage_2fStartWorkUnitRunner_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2fmessage_2fStartWorkUnitRunner_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2fmessage_2fStartWorkUnitRunner_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2fmessage_2fStartWorkUnitRunner_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2fmessage_2fStartWorkUnitRunner_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::message::StartWorkUnitRunner, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::message::StartWorkUnitRunner, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::message::StartWorkUnitRunner, work_unit_runner_id_),
  PROTOBUF_FIELD_OFFSET(::lm::message::StartWorkUnitRunner, use_cpu_affinity_),
  PROTOBUF_FIELD_OFFSET(::lm::message::StartWorkUnitRunner, cpu_),
  PROTOBUF_FIELD_OFFSET(::lm::message::StartWorkUnitRunner, gpu_),
  PROTOBUF_FIELD_OFFSET(::lm::message::StartWorkUnitRunner, me_solver_),
  PROTOBUF_FIELD_OFFSET(::lm::message::StartWorkUnitRunner, diffusion_pde_solver_),
  2,
  3,
  ~0u,
  ~0u,
  0,
  1,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 11, sizeof(::lm::message::StartWorkUnitRunner)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::message::_StartWorkUnitRunner_default_instance_),
};

const char descriptor_table_protodef_lm_2fmessage_2fStartWorkUnitRunner_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n$lm/message/StartWorkUnitRunner.proto\022\n"
  "lm.message\"\227\001\n\023StartWorkUnitRunner\022\033\n\023wo"
  "rk_unit_runner_id\030\001 \002(\005\022\030\n\020use_cpu_affin"
  "ity\030\002 \001(\010\022\013\n\003cpu\030\003 \003(\005\022\013\n\003gpu\030\004 \003(\005\022\021\n\tm"
  "e_solver\030\005 \002(\t\022\034\n\024diffusion_pde_solver\030\013"
  " \002(\t"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2fmessage_2fStartWorkUnitRunner_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2fmessage_2fStartWorkUnitRunner_2eproto_sccs[1] = {
  &scc_info_StartWorkUnitRunner_lm_2fmessage_2fStartWorkUnitRunner_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2fmessage_2fStartWorkUnitRunner_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fmessage_2fStartWorkUnitRunner_2eproto = {
  false, false, descriptor_table_protodef_lm_2fmessage_2fStartWorkUnitRunner_2eproto, "lm/message/StartWorkUnitRunner.proto", 204,
  &descriptor_table_lm_2fmessage_2fStartWorkUnitRunner_2eproto_once, descriptor_table_lm_2fmessage_2fStartWorkUnitRunner_2eproto_sccs, descriptor_table_lm_2fmessage_2fStartWorkUnitRunner_2eproto_deps, 1, 0,
  schemas, file_default_instances, TableStruct_lm_2fmessage_2fStartWorkUnitRunner_2eproto::offsets,
  file_level_metadata_lm_2fmessage_2fStartWorkUnitRunner_2eproto, 1, file_level_enum_descriptors_lm_2fmessage_2fStartWorkUnitRunner_2eproto, file_level_service_descriptors_lm_2fmessage_2fStartWorkUnitRunner_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2fmessage_2fStartWorkUnitRunner_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2fmessage_2fStartWorkUnitRunner_2eproto)), true);
namespace lm {
namespace message {

// ===================================================================

void StartWorkUnitRunner::InitAsDefaultInstance() {
}
class StartWorkUnitRunner::_Internal {
 public:
  using HasBits = decltype(std::declval<StartWorkUnitRunner>()._has_bits_);
  static void set_has_work_unit_runner_id(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_use_cpu_affinity(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
  static void set_has_me_solver(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_diffusion_pde_solver(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000007) ^ 0x00000007) != 0;
  }
};

StartWorkUnitRunner::StartWorkUnitRunner(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  cpu_(arena),
  gpu_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.message.StartWorkUnitRunner)
}
StartWorkUnitRunner::StartWorkUnitRunner(const StartWorkUnitRunner& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_),
      cpu_(from.cpu_),
      gpu_(from.gpu_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  me_solver_.UnsafeSetDefault(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited());
  if (from._internal_has_me_solver()) {
    me_solver_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), from._internal_me_solver(),
      GetArena());
  }
  diffusion_pde_solver_.UnsafeSetDefault(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited());
  if (from._internal_has_diffusion_pde_solver()) {
    diffusion_pde_solver_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), from._internal_diffusion_pde_solver(),
      GetArena());
  }
  ::memcpy(&work_unit_runner_id_, &from.work_unit_runner_id_,
    static_cast<size_t>(reinterpret_cast<char*>(&use_cpu_affinity_) -
    reinterpret_cast<char*>(&work_unit_runner_id_)) + sizeof(use_cpu_affinity_));
  // @@protoc_insertion_point(copy_constructor:lm.message.StartWorkUnitRunner)
}

void StartWorkUnitRunner::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_StartWorkUnitRunner_lm_2fmessage_2fStartWorkUnitRunner_2eproto.base);
  me_solver_.UnsafeSetDefault(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited());
  diffusion_pde_solver_.UnsafeSetDefault(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited());
  ::memset(&work_unit_runner_id_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&use_cpu_affinity_) -
      reinterpret_cast<char*>(&work_unit_runner_id_)) + sizeof(use_cpu_affinity_));
}

StartWorkUnitRunner::~StartWorkUnitRunner() {
  // @@protoc_insertion_point(destructor:lm.message.StartWorkUnitRunner)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void StartWorkUnitRunner::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
  me_solver_.DestroyNoArena(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited());
  diffusion_pde_solver_.DestroyNoArena(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited());
}

void StartWorkUnitRunner::ArenaDtor(void* object) {
  StartWorkUnitRunner* _this = reinterpret_cast< StartWorkUnitRunner* >(object);
  (void)_this;
}
void StartWorkUnitRunner::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void StartWorkUnitRunner::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const StartWorkUnitRunner& StartWorkUnitRunner::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_StartWorkUnitRunner_lm_2fmessage_2fStartWorkUnitRunner_2eproto.base);
  return *internal_default_instance();
}


void StartWorkUnitRunner::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.message.StartWorkUnitRunner)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cpu_.Clear();
  gpu_.Clear();
  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      me_solver_.ClearNonDefaultToEmpty();
    }
    if (cached_has_bits & 0x00000002u) {
      diffusion_pde_solver_.ClearNonDefaultToEmpty();
    }
  }
  if (cached_has_bits & 0x0000000cu) {
    ::memset(&work_unit_runner_id_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&use_cpu_affinity_) -
        reinterpret_cast<char*>(&work_unit_runner_id_)) + sizeof(use_cpu_affinity_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* StartWorkUnitRunner::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // required int32 work_unit_runner_id = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          _Internal::set_has_work_unit_runner_id(&has_bits);
          work_unit_runner_id_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional bool use_cpu_affinity = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 16)) {
          _Internal::set_has_use_cpu_affinity(&has_bits);
          use_cpu_affinity_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated int32 cpu = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 24)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_cpu(::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr));
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<24>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 26) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedInt32Parser(_internal_mutable_cpu(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated int32 gpu = 4;
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 32)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_gpu(::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr));
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<32>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 34) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedInt32Parser(_internal_mutable_gpu(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // required string me_solver = 5;
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 42)) {
          auto str = _internal_mutable_me_solver();
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::InlineGreedyStringParser(str, ptr, ctx);
          #ifndef NDEBUG
          ::PROTOBUF_NAMESPACE_ID::internal::VerifyUTF8(str, "lm.message.StartWorkUnitRunner.me_solver");
          #endif  // !NDEBUG
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // required string diffusion_pde_solver = 11;
      case 11:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 90)) {
          auto str = _internal_mutable_diffusion_pde_solver();
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::InlineGreedyStringParser(str, ptr, ctx);
          #ifndef NDEBUG
          ::PROTOBUF_NAMESPACE_ID::internal::VerifyUTF8(str, "lm.message.StartWorkUnitRunner.diffusion_pde_solver");
          #endif  // !NDEBUG
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

::PROTOBUF_NAMESPACE_ID::uint8* StartWorkUnitRunner::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.message.StartWorkUnitRunner)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required int32 work_unit_runner_id = 1;
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt32ToArray(1, this->_internal_work_unit_runner_id(), target);
  }

  // optional bool use_cpu_affinity = 2;
  if (cached_has_bits & 0x00000008u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteBoolToArray(2, this->_internal_use_cpu_affinity(), target);
  }

  // repeated int32 cpu = 3;
  for (int i = 0, n = this->_internal_cpu_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt32ToArray(3, this->_internal_cpu(i), target);
  }

  // repeated int32 gpu = 4;
  for (int i = 0, n = this->_internal_gpu_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt32ToArray(4, this->_internal_gpu(i), target);
  }

  // required string me_solver = 5;
  if (cached_has_bits & 0x00000001u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_me_solver().data(), static_cast<int>(this->_internal_me_solver().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "lm.message.StartWorkUnitRunner.me_solver");
    target = stream->WriteStringMaybeAliased(
        5, this->_internal_me_solver(), target);
  }

  // required string diffusion_pde_solver = 11;
  if (cached_has_bits & 0x00000002u) {
    ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::VerifyUTF8StringNamedField(
      this->_internal_diffusion_pde_solver().data(), static_cast<int>(this->_internal_diffusion_pde_solver().length()),
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::SERIALIZE,
      "lm.message.StartWorkUnitRunner.diffusion_pde_solver");
    target = stream->WriteStringMaybeAliased(
        11, this->_internal_diffusion_pde_solver(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.message.StartWorkUnitRunner)
  return target;
}

size_t StartWorkUnitRunner::RequiredFieldsByteSizeFallback() const {
// @@protoc_insertion_point(required_fields_byte_size_fallback_start:lm.message.StartWorkUnitRunner)
  size_t total_size = 0;

  if (_internal_has_me_solver()) {
    // required string me_solver = 5;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
        this->_internal_me_solver());
  }

  if (_internal_has_diffusion_pde_solver()) {
    // required string diffusion_pde_solver = 11;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
        this->_internal_diffusion_pde_solver());
  }

  if (_internal_has_work_unit_runner_id()) {
    // required int32 work_unit_runner_id = 1;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
        this->_internal_work_unit_runner_id());
  }

  return total_size;
}
size_t StartWorkUnitRunner::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.message.StartWorkUnitRunner)
  size_t total_size = 0;

  if (((_has_bits_[0] & 0x00000007) ^ 0x00000007) == 0) {  // All required fields are present.
    // required string me_solver = 5;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
        this->_internal_me_solver());

    // required string diffusion_pde_solver = 11;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::StringSize(
        this->_internal_diffusion_pde_solver());

    // required int32 work_unit_runner_id = 1;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
        this->_internal_work_unit_runner_id());

  } else {
    total_size += RequiredFieldsByteSizeFallback();
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated int32 cpu = 3;
  {
    size_t data_size = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      Int32Size(this->cpu_);
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_cpu_size());
    total_size += data_size;
  }

  // repeated int32 gpu = 4;
  {
    size_t data_size = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      Int32Size(this->gpu_);
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_gpu_size());
    total_size += data_size;
  }

  // optional bool use_cpu_affinity = 2;
  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000008u) {
    total_size += 1 + 1;
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void StartWorkUnitRunner::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.message.StartWorkUnitRunner)
  GOOGLE_DCHECK_NE(&from, this);
  const StartWorkUnitRunner* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<StartWorkUnitRunner>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.message.StartWorkUnitRunner)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.message.StartWorkUnitRunner)
    MergeFrom(*source);
  }
}

void StartWorkUnitRunner::MergeFrom(const StartWorkUnitRunner& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.message.StartWorkUnitRunner)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cpu_.MergeFrom(from.cpu_);
  gpu_.MergeFrom(from.gpu_);
  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x0000000fu) {
    if (cached_has_bits & 0x00000001u) {
      _internal_set_me_solver(from._internal_me_solver());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_set_diffusion_pde_solver(from._internal_diffusion_pde_solver());
    }
    if (cached_has_bits & 0x00000004u) {
      work_unit_runner_id_ = from.work_unit_runner_id_;
    }
    if (cached_has_bits & 0x00000008u) {
      use_cpu_affinity_ = from.use_cpu_affinity_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void StartWorkUnitRunner::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.message.StartWorkUnitRunner)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void StartWorkUnitRunner::CopyFrom(const StartWorkUnitRunner& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.message.StartWorkUnitRunner)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool StartWorkUnitRunner::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  return true;
}

void StartWorkUnitRunner::InternalSwap(StartWorkUnitRunner* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  cpu_.InternalSwap(&other->cpu_);
  gpu_.InternalSwap(&other->gpu_);
  me_solver_.Swap(&other->me_solver_, &::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
  diffusion_pde_solver_.Swap(&other->diffusion_pde_solver_, &::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(StartWorkUnitRunner, use_cpu_affinity_)
      + sizeof(StartWorkUnitRunner::use_cpu_affinity_)
      - PROTOBUF_FIELD_OFFSET(StartWorkUnitRunner, work_unit_runner_id_)>(
          reinterpret_cast<char*>(&work_unit_runner_id_),
          reinterpret_cast<char*>(&other->work_unit_runner_id_));
}

::PROTOBUF_NAMESPACE_ID::Metadata StartWorkUnitRunner::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace message
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::message::StartWorkUnitRunner* Arena::CreateMaybeMessage< ::lm::message::StartWorkUnitRunner >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::message::StartWorkUnitRunner >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>