// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/io/ffpilot/FFPilotStage.proto

#include "lm/io/ffpilot/FFPilotStage.pb.h"

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
extern PROTOBUF_INTERNAL_EXPORT_lm_2fio_2fffpilot_2fFFPilotPhase_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<2> scc_info_FFPilotPhase_lm_2fio_2fffpilot_2fFFPilotPhase_2eproto;
extern PROTOBUF_INTERNAL_EXPORT_lm_2fio_2fffpilot_2fFFPilotStageOutput_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_FFPilotStageOutput_lm_2fio_2fffpilot_2fFFPilotStageOutput_2eproto;
extern PROTOBUF_INTERNAL_EXPORT_lm_2ftypes_2fOrderParameters_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_OrderParameters_lm_2ftypes_2fOrderParameters_2eproto;
extern PROTOBUF_INTERNAL_EXPORT_lm_2ftypes_2fTilings_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_Tiling_lm_2ftypes_2fTilings_2eproto;
namespace lm {
namespace io {
namespace ffpilot {
class FFPilotStageDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<FFPilotStage> _instance;
} _FFPilotStage_default_instance_;
}  // namespace ffpilot
}  // namespace io
}  // namespace lm
static void InitDefaultsscc_info_FFPilotStage_lm_2fio_2fffpilot_2fFFPilotStage_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::io::ffpilot::_FFPilotStage_default_instance_;
    new (ptr) ::lm::io::ffpilot::FFPilotStage();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::io::ffpilot::FFPilotStage::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<4> scc_info_FFPilotStage_lm_2fio_2fffpilot_2fFFPilotStage_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 4, 0, InitDefaultsscc_info_FFPilotStage_lm_2fio_2fffpilot_2fFFPilotStage_2eproto}, {
      &scc_info_Tiling_lm_2ftypes_2fTilings_2eproto.base,
      &scc_info_OrderParameters_lm_2ftypes_2fOrderParameters_2eproto.base,
      &scc_info_FFPilotPhase_lm_2fio_2fffpilot_2fFFPilotPhase_2eproto.base,
      &scc_info_FFPilotStageOutput_lm_2fio_2fffpilot_2fFFPilotStageOutput_2eproto.base,}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2fio_2fffpilot_2fFFPilotStage_2eproto[1];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2fio_2fffpilot_2fFFPilotStage_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2fio_2fffpilot_2fFFPilotStage_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2fio_2fffpilot_2fFFPilotStage_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, basin_id_),
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, tiling_id_),
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, replicate_id_),
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, tiling_),
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, order_parameter_),
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, is_pilot_stage_),
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, needs_pilot_stage_),
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, ffpilot_phases_),
  PROTOBUF_FIELD_OFFSET(::lm::io::ffpilot::FFPilotStage, pilot_stage_output_),
  3,
  4,
  5,
  0,
  1,
  6,
  7,
  ~0u,
  2,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 14, sizeof(::lm::io::ffpilot::FFPilotStage)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::io::ffpilot::_FFPilotStage_default_instance_),
};

const char descriptor_table_protodef_lm_2fio_2fffpilot_2fFFPilotStage_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n lm/io/ffpilot/FFPilotStage.proto\022\rlm.i"
  "o.ffpilot\032 lm/io/ffpilot/FFPilotPhase.pr"
  "oto\032&lm/io/ffpilot/FFPilotStageOutput.pr"
  "oto\032\036lm/types/OrderParameters.proto\032\026lm/"
  "types/Tilings.proto\"\331\002\n\014FFPilotStage\022\020\n\010"
  "basin_id\030\001 \001(\003\022\021\n\ttiling_id\030\002 \001(\004\022\027\n\014rep"
  "licate_id\030\004 \001(\004:\0010\022 \n\006tiling\030\003 \001(\0132\020.lm."
  "types.Tiling\0222\n\017order_parameter\030\005 \001(\0132\031."
  "lm.types.OrderParameters\022\035\n\016is_pilot_sta"
  "ge\030\025 \001(\010:\005false\022 \n\021needs_pilot_stage\030\026 \001"
  "(\010:\005false\0224\n\016ffpilot_phases\030\311\001 \003(\0132\033.lm."
  "io.ffpilot.FFPilotPhase\022>\n\022pilot_stage_o"
  "utput\030\255\002 \001(\0132!.lm.io.ffpilot.FFPilotStag"
  "eOutput"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto_deps[4] = {
  &::descriptor_table_lm_2fio_2fffpilot_2fFFPilotPhase_2eproto,
  &::descriptor_table_lm_2fio_2fffpilot_2fFFPilotStageOutput_2eproto,
  &::descriptor_table_lm_2ftypes_2fOrderParameters_2eproto,
  &::descriptor_table_lm_2ftypes_2fTilings_2eproto,
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto_sccs[1] = {
  &scc_info_FFPilotStage_lm_2fio_2fffpilot_2fFFPilotStage_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto = {
  false, false, descriptor_table_protodef_lm_2fio_2fffpilot_2fFFPilotStage_2eproto, "lm/io/ffpilot/FFPilotStage.proto", 527,
  &descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto_once, descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto_sccs, descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto_deps, 1, 4,
  schemas, file_default_instances, TableStruct_lm_2fio_2fffpilot_2fFFPilotStage_2eproto::offsets,
  file_level_metadata_lm_2fio_2fffpilot_2fFFPilotStage_2eproto, 1, file_level_enum_descriptors_lm_2fio_2fffpilot_2fFFPilotStage_2eproto, file_level_service_descriptors_lm_2fio_2fffpilot_2fFFPilotStage_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2fio_2fffpilot_2fFFPilotStage_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2fio_2fffpilot_2fFFPilotStage_2eproto)), true);
namespace lm {
namespace io {
namespace ffpilot {

// ===================================================================

void FFPilotStage::InitAsDefaultInstance() {
  ::lm::io::ffpilot::_FFPilotStage_default_instance_._instance.get_mutable()->tiling_ = const_cast< ::lm::types::Tiling*>(
      ::lm::types::Tiling::internal_default_instance());
  ::lm::io::ffpilot::_FFPilotStage_default_instance_._instance.get_mutable()->order_parameter_ = const_cast< ::lm::types::OrderParameters*>(
      ::lm::types::OrderParameters::internal_default_instance());
  ::lm::io::ffpilot::_FFPilotStage_default_instance_._instance.get_mutable()->pilot_stage_output_ = const_cast< ::lm::io::ffpilot::FFPilotStageOutput*>(
      ::lm::io::ffpilot::FFPilotStageOutput::internal_default_instance());
}
class FFPilotStage::_Internal {
 public:
  using HasBits = decltype(std::declval<FFPilotStage>()._has_bits_);
  static void set_has_basin_id(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
  static void set_has_tiling_id(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
  static void set_has_replicate_id(HasBits* has_bits) {
    (*has_bits)[0] |= 32u;
  }
  static const ::lm::types::Tiling& tiling(const FFPilotStage* msg);
  static void set_has_tiling(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static const ::lm::types::OrderParameters& order_parameter(const FFPilotStage* msg);
  static void set_has_order_parameter(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static void set_has_is_pilot_stage(HasBits* has_bits) {
    (*has_bits)[0] |= 64u;
  }
  static void set_has_needs_pilot_stage(HasBits* has_bits) {
    (*has_bits)[0] |= 128u;
  }
  static const ::lm::io::ffpilot::FFPilotStageOutput& pilot_stage_output(const FFPilotStage* msg);
  static void set_has_pilot_stage_output(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
};

const ::lm::types::Tiling&
FFPilotStage::_Internal::tiling(const FFPilotStage* msg) {
  return *msg->tiling_;
}
const ::lm::types::OrderParameters&
FFPilotStage::_Internal::order_parameter(const FFPilotStage* msg) {
  return *msg->order_parameter_;
}
const ::lm::io::ffpilot::FFPilotStageOutput&
FFPilotStage::_Internal::pilot_stage_output(const FFPilotStage* msg) {
  return *msg->pilot_stage_output_;
}
void FFPilotStage::clear_tiling() {
  if (tiling_ != nullptr) tiling_->Clear();
  _has_bits_[0] &= ~0x00000001u;
}
void FFPilotStage::clear_order_parameter() {
  if (order_parameter_ != nullptr) order_parameter_->Clear();
  _has_bits_[0] &= ~0x00000002u;
}
void FFPilotStage::clear_ffpilot_phases() {
  ffpilot_phases_.Clear();
}
void FFPilotStage::clear_pilot_stage_output() {
  if (pilot_stage_output_ != nullptr) pilot_stage_output_->Clear();
  _has_bits_[0] &= ~0x00000004u;
}
FFPilotStage::FFPilotStage(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  ffpilot_phases_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.io.ffpilot.FFPilotStage)
}
FFPilotStage::FFPilotStage(const FFPilotStage& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_),
      ffpilot_phases_(from.ffpilot_phases_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  if (from._internal_has_tiling()) {
    tiling_ = new ::lm::types::Tiling(*from.tiling_);
  } else {
    tiling_ = nullptr;
  }
  if (from._internal_has_order_parameter()) {
    order_parameter_ = new ::lm::types::OrderParameters(*from.order_parameter_);
  } else {
    order_parameter_ = nullptr;
  }
  if (from._internal_has_pilot_stage_output()) {
    pilot_stage_output_ = new ::lm::io::ffpilot::FFPilotStageOutput(*from.pilot_stage_output_);
  } else {
    pilot_stage_output_ = nullptr;
  }
  ::memcpy(&basin_id_, &from.basin_id_,
    static_cast<size_t>(reinterpret_cast<char*>(&needs_pilot_stage_) -
    reinterpret_cast<char*>(&basin_id_)) + sizeof(needs_pilot_stage_));
  // @@protoc_insertion_point(copy_constructor:lm.io.ffpilot.FFPilotStage)
}

void FFPilotStage::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_FFPilotStage_lm_2fio_2fffpilot_2fFFPilotStage_2eproto.base);
  ::memset(&tiling_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&needs_pilot_stage_) -
      reinterpret_cast<char*>(&tiling_)) + sizeof(needs_pilot_stage_));
}

FFPilotStage::~FFPilotStage() {
  // @@protoc_insertion_point(destructor:lm.io.ffpilot.FFPilotStage)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void FFPilotStage::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
  if (this != internal_default_instance()) delete tiling_;
  if (this != internal_default_instance()) delete order_parameter_;
  if (this != internal_default_instance()) delete pilot_stage_output_;
}

void FFPilotStage::ArenaDtor(void* object) {
  FFPilotStage* _this = reinterpret_cast< FFPilotStage* >(object);
  (void)_this;
}
void FFPilotStage::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void FFPilotStage::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const FFPilotStage& FFPilotStage::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_FFPilotStage_lm_2fio_2fffpilot_2fFFPilotStage_2eproto.base);
  return *internal_default_instance();
}


void FFPilotStage::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.io.ffpilot.FFPilotStage)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  ffpilot_phases_.Clear();
  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000007u) {
    if (cached_has_bits & 0x00000001u) {
      GOOGLE_DCHECK(tiling_ != nullptr);
      tiling_->Clear();
    }
    if (cached_has_bits & 0x00000002u) {
      GOOGLE_DCHECK(order_parameter_ != nullptr);
      order_parameter_->Clear();
    }
    if (cached_has_bits & 0x00000004u) {
      GOOGLE_DCHECK(pilot_stage_output_ != nullptr);
      pilot_stage_output_->Clear();
    }
  }
  if (cached_has_bits & 0x000000f8u) {
    ::memset(&basin_id_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&needs_pilot_stage_) -
        reinterpret_cast<char*>(&basin_id_)) + sizeof(needs_pilot_stage_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* FFPilotStage::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // optional int64 basin_id = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          _Internal::set_has_basin_id(&has_bits);
          basin_id_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional uint64 tiling_id = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 16)) {
          _Internal::set_has_tiling_id(&has_bits);
          tiling_id_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional .lm.types.Tiling tiling = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 26)) {
          ptr = ctx->ParseMessage(_internal_mutable_tiling(), ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional uint64 replicate_id = 4 [default = 0];
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 32)) {
          _Internal::set_has_replicate_id(&has_bits);
          replicate_id_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional .lm.types.OrderParameters order_parameter = 5;
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 42)) {
          ptr = ctx->ParseMessage(_internal_mutable_order_parameter(), ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional bool is_pilot_stage = 21 [default = false];
      case 21:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 168)) {
          _Internal::set_has_is_pilot_stage(&has_bits);
          is_pilot_stage_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional bool needs_pilot_stage = 22 [default = false];
      case 22:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 176)) {
          _Internal::set_has_needs_pilot_stage(&has_bits);
          needs_pilot_stage_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated .lm.io.ffpilot.FFPilotPhase ffpilot_phases = 201;
      case 201:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 74)) {
          ptr -= 2;
          do {
            ptr += 2;
            ptr = ctx->ParseMessage(_internal_add_ffpilot_phases(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<1610>(ptr));
        } else goto handle_unusual;
        continue;
      // optional .lm.io.ffpilot.FFPilotStageOutput pilot_stage_output = 301;
      case 301:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 106)) {
          ptr = ctx->ParseMessage(_internal_mutable_pilot_stage_output(), ptr);
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

::PROTOBUF_NAMESPACE_ID::uint8* FFPilotStage::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.io.ffpilot.FFPilotStage)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // optional int64 basin_id = 1;
  if (cached_has_bits & 0x00000008u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt64ToArray(1, this->_internal_basin_id(), target);
  }

  // optional uint64 tiling_id = 2;
  if (cached_has_bits & 0x00000010u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt64ToArray(2, this->_internal_tiling_id(), target);
  }

  // optional .lm.types.Tiling tiling = 3;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        3, _Internal::tiling(this), target, stream);
  }

  // optional uint64 replicate_id = 4 [default = 0];
  if (cached_has_bits & 0x00000020u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt64ToArray(4, this->_internal_replicate_id(), target);
  }

  // optional .lm.types.OrderParameters order_parameter = 5;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        5, _Internal::order_parameter(this), target, stream);
  }

  // optional bool is_pilot_stage = 21 [default = false];
  if (cached_has_bits & 0x00000040u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteBoolToArray(21, this->_internal_is_pilot_stage(), target);
  }

  // optional bool needs_pilot_stage = 22 [default = false];
  if (cached_has_bits & 0x00000080u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteBoolToArray(22, this->_internal_needs_pilot_stage(), target);
  }

  // repeated .lm.io.ffpilot.FFPilotPhase ffpilot_phases = 201;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_ffpilot_phases_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(201, this->_internal_ffpilot_phases(i), target, stream);
  }

  // optional .lm.io.ffpilot.FFPilotStageOutput pilot_stage_output = 301;
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(
        301, _Internal::pilot_stage_output(this), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.io.ffpilot.FFPilotStage)
  return target;
}

size_t FFPilotStage::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.io.ffpilot.FFPilotStage)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .lm.io.ffpilot.FFPilotPhase ffpilot_phases = 201;
  total_size += 2UL * this->_internal_ffpilot_phases_size();
  for (const auto& msg : this->ffpilot_phases_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x000000ffu) {
    // optional .lm.types.Tiling tiling = 3;
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *tiling_);
    }

    // optional .lm.types.OrderParameters order_parameter = 5;
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *order_parameter_);
    }

    // optional .lm.io.ffpilot.FFPilotStageOutput pilot_stage_output = 301;
    if (cached_has_bits & 0x00000004u) {
      total_size += 2 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(
          *pilot_stage_output_);
    }

    // optional int64 basin_id = 1;
    if (cached_has_bits & 0x00000008u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int64Size(
          this->_internal_basin_id());
    }

    // optional uint64 tiling_id = 2;
    if (cached_has_bits & 0x00000010u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt64Size(
          this->_internal_tiling_id());
    }

    // optional uint64 replicate_id = 4 [default = 0];
    if (cached_has_bits & 0x00000020u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt64Size(
          this->_internal_replicate_id());
    }

    // optional bool is_pilot_stage = 21 [default = false];
    if (cached_has_bits & 0x00000040u) {
      total_size += 2 + 1;
    }

    // optional bool needs_pilot_stage = 22 [default = false];
    if (cached_has_bits & 0x00000080u) {
      total_size += 2 + 1;
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

void FFPilotStage::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.io.ffpilot.FFPilotStage)
  GOOGLE_DCHECK_NE(&from, this);
  const FFPilotStage* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<FFPilotStage>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.io.ffpilot.FFPilotStage)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.io.ffpilot.FFPilotStage)
    MergeFrom(*source);
  }
}

void FFPilotStage::MergeFrom(const FFPilotStage& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.io.ffpilot.FFPilotStage)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  ffpilot_phases_.MergeFrom(from.ffpilot_phases_);
  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x000000ffu) {
    if (cached_has_bits & 0x00000001u) {
      _internal_mutable_tiling()->::lm::types::Tiling::MergeFrom(from._internal_tiling());
    }
    if (cached_has_bits & 0x00000002u) {
      _internal_mutable_order_parameter()->::lm::types::OrderParameters::MergeFrom(from._internal_order_parameter());
    }
    if (cached_has_bits & 0x00000004u) {
      _internal_mutable_pilot_stage_output()->::lm::io::ffpilot::FFPilotStageOutput::MergeFrom(from._internal_pilot_stage_output());
    }
    if (cached_has_bits & 0x00000008u) {
      basin_id_ = from.basin_id_;
    }
    if (cached_has_bits & 0x00000010u) {
      tiling_id_ = from.tiling_id_;
    }
    if (cached_has_bits & 0x00000020u) {
      replicate_id_ = from.replicate_id_;
    }
    if (cached_has_bits & 0x00000040u) {
      is_pilot_stage_ = from.is_pilot_stage_;
    }
    if (cached_has_bits & 0x00000080u) {
      needs_pilot_stage_ = from.needs_pilot_stage_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void FFPilotStage::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.io.ffpilot.FFPilotStage)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void FFPilotStage::CopyFrom(const FFPilotStage& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.io.ffpilot.FFPilotStage)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool FFPilotStage::IsInitialized() const {
  if (_internal_has_tiling()) {
    if (!tiling_->IsInitialized()) return false;
  }
  if (_internal_has_order_parameter()) {
    if (!order_parameter_->IsInitialized()) return false;
  }
  if (_internal_has_pilot_stage_output()) {
    if (!pilot_stage_output_->IsInitialized()) return false;
  }
  return true;
}

void FFPilotStage::InternalSwap(FFPilotStage* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  ffpilot_phases_.InternalSwap(&other->ffpilot_phases_);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(FFPilotStage, needs_pilot_stage_)
      + sizeof(FFPilotStage::needs_pilot_stage_)
      - PROTOBUF_FIELD_OFFSET(FFPilotStage, tiling_)>(
          reinterpret_cast<char*>(&tiling_),
          reinterpret_cast<char*>(&other->tiling_));
}

::PROTOBUF_NAMESPACE_ID::Metadata FFPilotStage::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace ffpilot
}  // namespace io
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::io::ffpilot::FFPilotStage* Arena::CreateMaybeMessage< ::lm::io::ffpilot::FFPilotStage >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::io::ffpilot::FFPilotStage >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
