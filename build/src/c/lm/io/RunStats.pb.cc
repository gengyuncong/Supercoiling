// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/io/RunStats.proto

#include "lm/io/RunStats.pb.h"

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
extern PROTOBUF_INTERNAL_EXPORT_lm_2fio_2fRunStats_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_RunStats_TrajectoryStat_lm_2fio_2fRunStats_2eproto;
namespace lm {
namespace io {
class RunStats_TrajectoryStatDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<RunStats_TrajectoryStat> _instance;
} _RunStats_TrajectoryStat_default_instance_;
class RunStatsDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<RunStats> _instance;
} _RunStats_default_instance_;
}  // namespace io
}  // namespace lm
static void InitDefaultsscc_info_RunStats_lm_2fio_2fRunStats_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::io::_RunStats_default_instance_;
    new (ptr) ::lm::io::RunStats();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::io::RunStats::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_RunStats_lm_2fio_2fRunStats_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_RunStats_lm_2fio_2fRunStats_2eproto}, {
      &scc_info_RunStats_TrajectoryStat_lm_2fio_2fRunStats_2eproto.base,}};

static void InitDefaultsscc_info_RunStats_TrajectoryStat_lm_2fio_2fRunStats_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::lm::io::_RunStats_TrajectoryStat_default_instance_;
    new (ptr) ::lm::io::RunStats_TrajectoryStat();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::lm::io::RunStats_TrajectoryStat::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_RunStats_TrajectoryStat_lm_2fio_2fRunStats_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_RunStats_TrajectoryStat_lm_2fio_2fRunStats_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_lm_2fio_2fRunStats_2eproto[2];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_lm_2fio_2fRunStats_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_lm_2fio_2fRunStats_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_lm_2fio_2fRunStats_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::lm::io::RunStats_TrajectoryStat, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::io::RunStats_TrajectoryStat, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::io::RunStats_TrajectoryStat, trajectory_id_),
  PROTOBUF_FIELD_OFFSET(::lm::io::RunStats_TrajectoryStat, run_time_),
  0,
  1,
  PROTOBUF_FIELD_OFFSET(::lm::io::RunStats, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::lm::io::RunStats, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::lm::io::RunStats, trajectory_stat_),
  PROTOBUF_FIELD_OFFSET(::lm::io::RunStats, total_runt_time_),
  PROTOBUF_FIELD_OFFSET(::lm::io::RunStats, count_trajectories_),
  ~0u,
  0,
  1,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 7, sizeof(::lm::io::RunStats_TrajectoryStat)},
  { 9, 17, sizeof(::lm::io::RunStats)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::io::_RunStats_TrajectoryStat_default_instance_),
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::lm::io::_RunStats_default_instance_),
};

const char descriptor_table_protodef_lm_2fio_2fRunStats_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\024lm/io/RunStats.proto\022\005lm.io\"\277\001\n\010RunSta"
  "ts\0227\n\017trajectory_stat\030\001 \003(\0132\036.lm.io.RunS"
  "tats.TrajectoryStat\022\032\n\017total_runt_time\030\002"
  " \001(\001:\0010\022\035\n\022count_trajectories\030\003 \001(\004:\0010\032\?"
  "\n\016TrajectoryStat\022\030\n\rtrajectory_id\030\001 \001(\004:"
  "\0010\022\023\n\010run_time\030\002 \001(\001:\0010"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_lm_2fio_2fRunStats_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_lm_2fio_2fRunStats_2eproto_sccs[2] = {
  &scc_info_RunStats_lm_2fio_2fRunStats_2eproto.base,
  &scc_info_RunStats_TrajectoryStat_lm_2fio_2fRunStats_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_lm_2fio_2fRunStats_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fio_2fRunStats_2eproto = {
  false, false, descriptor_table_protodef_lm_2fio_2fRunStats_2eproto, "lm/io/RunStats.proto", 223,
  &descriptor_table_lm_2fio_2fRunStats_2eproto_once, descriptor_table_lm_2fio_2fRunStats_2eproto_sccs, descriptor_table_lm_2fio_2fRunStats_2eproto_deps, 2, 0,
  schemas, file_default_instances, TableStruct_lm_2fio_2fRunStats_2eproto::offsets,
  file_level_metadata_lm_2fio_2fRunStats_2eproto, 2, file_level_enum_descriptors_lm_2fio_2fRunStats_2eproto, file_level_service_descriptors_lm_2fio_2fRunStats_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_lm_2fio_2fRunStats_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_lm_2fio_2fRunStats_2eproto)), true);
namespace lm {
namespace io {

// ===================================================================

void RunStats_TrajectoryStat::InitAsDefaultInstance() {
}
class RunStats_TrajectoryStat::_Internal {
 public:
  using HasBits = decltype(std::declval<RunStats_TrajectoryStat>()._has_bits_);
  static void set_has_trajectory_id(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_run_time(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
};

RunStats_TrajectoryStat::RunStats_TrajectoryStat(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.io.RunStats.TrajectoryStat)
}
RunStats_TrajectoryStat::RunStats_TrajectoryStat(const RunStats_TrajectoryStat& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&trajectory_id_, &from.trajectory_id_,
    static_cast<size_t>(reinterpret_cast<char*>(&run_time_) -
    reinterpret_cast<char*>(&trajectory_id_)) + sizeof(run_time_));
  // @@protoc_insertion_point(copy_constructor:lm.io.RunStats.TrajectoryStat)
}

void RunStats_TrajectoryStat::SharedCtor() {
  ::memset(&trajectory_id_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&run_time_) -
      reinterpret_cast<char*>(&trajectory_id_)) + sizeof(run_time_));
}

RunStats_TrajectoryStat::~RunStats_TrajectoryStat() {
  // @@protoc_insertion_point(destructor:lm.io.RunStats.TrajectoryStat)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void RunStats_TrajectoryStat::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void RunStats_TrajectoryStat::ArenaDtor(void* object) {
  RunStats_TrajectoryStat* _this = reinterpret_cast< RunStats_TrajectoryStat* >(object);
  (void)_this;
}
void RunStats_TrajectoryStat::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void RunStats_TrajectoryStat::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const RunStats_TrajectoryStat& RunStats_TrajectoryStat::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_RunStats_TrajectoryStat_lm_2fio_2fRunStats_2eproto.base);
  return *internal_default_instance();
}


void RunStats_TrajectoryStat::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.io.RunStats.TrajectoryStat)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    ::memset(&trajectory_id_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&run_time_) -
        reinterpret_cast<char*>(&trajectory_id_)) + sizeof(run_time_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* RunStats_TrajectoryStat::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // optional uint64 trajectory_id = 1 [default = 0];
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          _Internal::set_has_trajectory_id(&has_bits);
          trajectory_id_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional double run_time = 2 [default = 0];
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 17)) {
          _Internal::set_has_run_time(&has_bits);
          run_time_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
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

::PROTOBUF_NAMESPACE_ID::uint8* RunStats_TrajectoryStat::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.io.RunStats.TrajectoryStat)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // optional uint64 trajectory_id = 1 [default = 0];
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt64ToArray(1, this->_internal_trajectory_id(), target);
  }

  // optional double run_time = 2 [default = 0];
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(2, this->_internal_run_time(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.io.RunStats.TrajectoryStat)
  return target;
}

size_t RunStats_TrajectoryStat::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.io.RunStats.TrajectoryStat)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    // optional uint64 trajectory_id = 1 [default = 0];
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt64Size(
          this->_internal_trajectory_id());
    }

    // optional double run_time = 2 [default = 0];
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 + 8;
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

void RunStats_TrajectoryStat::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.io.RunStats.TrajectoryStat)
  GOOGLE_DCHECK_NE(&from, this);
  const RunStats_TrajectoryStat* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<RunStats_TrajectoryStat>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.io.RunStats.TrajectoryStat)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.io.RunStats.TrajectoryStat)
    MergeFrom(*source);
  }
}

void RunStats_TrajectoryStat::MergeFrom(const RunStats_TrajectoryStat& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.io.RunStats.TrajectoryStat)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      trajectory_id_ = from.trajectory_id_;
    }
    if (cached_has_bits & 0x00000002u) {
      run_time_ = from.run_time_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void RunStats_TrajectoryStat::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.io.RunStats.TrajectoryStat)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void RunStats_TrajectoryStat::CopyFrom(const RunStats_TrajectoryStat& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.io.RunStats.TrajectoryStat)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool RunStats_TrajectoryStat::IsInitialized() const {
  return true;
}

void RunStats_TrajectoryStat::InternalSwap(RunStats_TrajectoryStat* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(RunStats_TrajectoryStat, run_time_)
      + sizeof(RunStats_TrajectoryStat::run_time_)
      - PROTOBUF_FIELD_OFFSET(RunStats_TrajectoryStat, trajectory_id_)>(
          reinterpret_cast<char*>(&trajectory_id_),
          reinterpret_cast<char*>(&other->trajectory_id_));
}

::PROTOBUF_NAMESPACE_ID::Metadata RunStats_TrajectoryStat::GetMetadata() const {
  return GetMetadataStatic();
}


// ===================================================================

void RunStats::InitAsDefaultInstance() {
}
class RunStats::_Internal {
 public:
  using HasBits = decltype(std::declval<RunStats>()._has_bits_);
  static void set_has_total_runt_time(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_count_trajectories(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
};

RunStats::RunStats(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  trajectory_stat_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:lm.io.RunStats)
}
RunStats::RunStats(const RunStats& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_),
      trajectory_stat_(from.trajectory_stat_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&total_runt_time_, &from.total_runt_time_,
    static_cast<size_t>(reinterpret_cast<char*>(&count_trajectories_) -
    reinterpret_cast<char*>(&total_runt_time_)) + sizeof(count_trajectories_));
  // @@protoc_insertion_point(copy_constructor:lm.io.RunStats)
}

void RunStats::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_RunStats_lm_2fio_2fRunStats_2eproto.base);
  ::memset(&total_runt_time_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&count_trajectories_) -
      reinterpret_cast<char*>(&total_runt_time_)) + sizeof(count_trajectories_));
}

RunStats::~RunStats() {
  // @@protoc_insertion_point(destructor:lm.io.RunStats)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void RunStats::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void RunStats::ArenaDtor(void* object) {
  RunStats* _this = reinterpret_cast< RunStats* >(object);
  (void)_this;
}
void RunStats::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void RunStats::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const RunStats& RunStats::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_RunStats_lm_2fio_2fRunStats_2eproto.base);
  return *internal_default_instance();
}


void RunStats::Clear() {
// @@protoc_insertion_point(message_clear_start:lm.io.RunStats)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  trajectory_stat_.Clear();
  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    ::memset(&total_runt_time_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&count_trajectories_) -
        reinterpret_cast<char*>(&total_runt_time_)) + sizeof(count_trajectories_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* RunStats::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated .lm.io.RunStats.TrajectoryStat trajectory_stat = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_trajectory_stat(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<10>(ptr));
        } else goto handle_unusual;
        continue;
      // optional double total_runt_time = 2 [default = 0];
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 17)) {
          _Internal::set_has_total_runt_time(&has_bits);
          total_runt_time_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else goto handle_unusual;
        continue;
      // optional uint64 count_trajectories = 3 [default = 0];
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 24)) {
          _Internal::set_has_count_trajectories(&has_bits);
          count_trajectories_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
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

::PROTOBUF_NAMESPACE_ID::uint8* RunStats::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:lm.io.RunStats)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated .lm.io.RunStats.TrajectoryStat trajectory_stat = 1;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_trajectory_stat_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(1, this->_internal_trajectory_stat(i), target, stream);
  }

  cached_has_bits = _has_bits_[0];
  // optional double total_runt_time = 2 [default = 0];
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(2, this->_internal_total_runt_time(), target);
  }

  // optional uint64 count_trajectories = 3 [default = 0];
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt64ToArray(3, this->_internal_count_trajectories(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:lm.io.RunStats)
  return target;
}

size_t RunStats::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:lm.io.RunStats)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .lm.io.RunStats.TrajectoryStat trajectory_stat = 1;
  total_size += 1UL * this->_internal_trajectory_stat_size();
  for (const auto& msg : this->trajectory_stat_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    // optional double total_runt_time = 2 [default = 0];
    if (cached_has_bits & 0x00000001u) {
      total_size += 1 + 8;
    }

    // optional uint64 count_trajectories = 3 [default = 0];
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::UInt64Size(
          this->_internal_count_trajectories());
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

void RunStats::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:lm.io.RunStats)
  GOOGLE_DCHECK_NE(&from, this);
  const RunStats* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<RunStats>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:lm.io.RunStats)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:lm.io.RunStats)
    MergeFrom(*source);
  }
}

void RunStats::MergeFrom(const RunStats& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:lm.io.RunStats)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  trajectory_stat_.MergeFrom(from.trajectory_stat_);
  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      total_runt_time_ = from.total_runt_time_;
    }
    if (cached_has_bits & 0x00000002u) {
      count_trajectories_ = from.count_trajectories_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void RunStats::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:lm.io.RunStats)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void RunStats::CopyFrom(const RunStats& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:lm.io.RunStats)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool RunStats::IsInitialized() const {
  return true;
}

void RunStats::InternalSwap(RunStats* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  trajectory_stat_.InternalSwap(&other->trajectory_stat_);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(RunStats, count_trajectories_)
      + sizeof(RunStats::count_trajectories_)
      - PROTOBUF_FIELD_OFFSET(RunStats, total_runt_time_)>(
          reinterpret_cast<char*>(&total_runt_time_),
          reinterpret_cast<char*>(&other->total_runt_time_));
}

::PROTOBUF_NAMESPACE_ID::Metadata RunStats::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace io
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::lm::io::RunStats_TrajectoryStat* Arena::CreateMaybeMessage< ::lm::io::RunStats_TrajectoryStat >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::io::RunStats_TrajectoryStat >(arena);
}
template<> PROTOBUF_NOINLINE ::lm::io::RunStats* Arena::CreateMaybeMessage< ::lm::io::RunStats >(Arena* arena) {
  return Arena::CreateMessageInternal< ::lm::io::RunStats >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>