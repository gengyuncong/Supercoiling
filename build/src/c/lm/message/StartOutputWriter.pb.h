// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/message/StartOutputWriter.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2fmessage_2fStartOutputWriter_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2fmessage_2fStartOutputWriter_2eproto

#include <limits>
#include <string>

#include <google/protobuf/port_def.inc>
#if PROTOBUF_VERSION < 3013000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers. Please update
#error your headers.
#endif
#if 3013000 < PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers. Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/port_undef.inc>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/arena.h>
#include <google/protobuf/arenastring.h>
#include <google/protobuf/generated_message_table_driven.h>
#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/inlined_string_field.h>
#include <google/protobuf/metadata_lite.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>  // IWYU pragma: export
#include <google/protobuf/extension_set.h>  // IWYU pragma: export
#include <google/protobuf/unknown_field_set.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_lm_2fmessage_2fStartOutputWriter_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2fmessage_2fStartOutputWriter_2eproto {
  static const ::PROTOBUF_NAMESPACE_ID::internal::ParseTableField entries[]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::AuxiliaryParseTableField aux[]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::ParseTable schema[1]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::FieldMetadata field_metadata[];
  static const ::PROTOBUF_NAMESPACE_ID::internal::SerializationTable serialization_table[];
  static const ::PROTOBUF_NAMESPACE_ID::uint32 offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fmessage_2fStartOutputWriter_2eproto;
namespace lm {
namespace message {
class StartOutputWriter;
class StartOutputWriterDefaultTypeInternal;
extern StartOutputWriterDefaultTypeInternal _StartOutputWriter_default_instance_;
}  // namespace message
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::message::StartOutputWriter* Arena::CreateMaybeMessage<::lm::message::StartOutputWriter>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace message {

// ===================================================================

class StartOutputWriter PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.message.StartOutputWriter) */ {
 public:
  inline StartOutputWriter() : StartOutputWriter(nullptr) {}
  virtual ~StartOutputWriter();

  StartOutputWriter(const StartOutputWriter& from);
  StartOutputWriter(StartOutputWriter&& from) noexcept
    : StartOutputWriter() {
    *this = ::std::move(from);
  }

  inline StartOutputWriter& operator=(const StartOutputWriter& from) {
    CopyFrom(from);
    return *this;
  }
  inline StartOutputWriter& operator=(StartOutputWriter&& from) noexcept {
    if (GetArena() == from.GetArena()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }

  inline const ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet& unknown_fields() const {
    return _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance);
  }
  inline ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet* mutable_unknown_fields() {
    return _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
  }

  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* descriptor() {
    return GetDescriptor();
  }
  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* GetDescriptor() {
    return GetMetadataStatic().descriptor;
  }
  static const ::PROTOBUF_NAMESPACE_ID::Reflection* GetReflection() {
    return GetMetadataStatic().reflection;
  }
  static const StartOutputWriter& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const StartOutputWriter* internal_default_instance() {
    return reinterpret_cast<const StartOutputWriter*>(
               &_StartOutputWriter_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(StartOutputWriter& a, StartOutputWriter& b) {
    a.Swap(&b);
  }
  inline void Swap(StartOutputWriter* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(StartOutputWriter* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline StartOutputWriter* New() const final {
    return CreateMaybeMessage<StartOutputWriter>(nullptr);
  }

  StartOutputWriter* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<StartOutputWriter>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const StartOutputWriter& from);
  void MergeFrom(const StartOutputWriter& from);
  PROTOBUF_ATTRIBUTE_REINITIALIZES void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  const char* _InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) final;
  ::PROTOBUF_NAMESPACE_ID::uint8* _InternalSerialize(
      ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  inline void SharedCtor();
  inline void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(StartOutputWriter* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.message.StartOutputWriter";
  }
  protected:
  explicit StartOutputWriter(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2fmessage_2fStartOutputWriter_2eproto);
    return ::descriptor_table_lm_2fmessage_2fStartOutputWriter_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kOutputFilenameFieldNumber = 3,
    kOutputWriterClassFieldNumber = 4,
    kUseCpuAffinityFieldNumber = 1,
    kCpuFieldNumber = 2,
  };
  // required string output_filename = 3;
  bool has_output_filename() const;
  private:
  bool _internal_has_output_filename() const;
  public:
  void clear_output_filename();
  const std::string& output_filename() const;
  void set_output_filename(const std::string& value);
  void set_output_filename(std::string&& value);
  void set_output_filename(const char* value);
  void set_output_filename(const char* value, size_t size);
  std::string* mutable_output_filename();
  std::string* release_output_filename();
  void set_allocated_output_filename(std::string* output_filename);
  private:
  const std::string& _internal_output_filename() const;
  void _internal_set_output_filename(const std::string& value);
  std::string* _internal_mutable_output_filename();
  public:

  // required string output_writer_class = 4;
  bool has_output_writer_class() const;
  private:
  bool _internal_has_output_writer_class() const;
  public:
  void clear_output_writer_class();
  const std::string& output_writer_class() const;
  void set_output_writer_class(const std::string& value);
  void set_output_writer_class(std::string&& value);
  void set_output_writer_class(const char* value);
  void set_output_writer_class(const char* value, size_t size);
  std::string* mutable_output_writer_class();
  std::string* release_output_writer_class();
  void set_allocated_output_writer_class(std::string* output_writer_class);
  private:
  const std::string& _internal_output_writer_class() const;
  void _internal_set_output_writer_class(const std::string& value);
  std::string* _internal_mutable_output_writer_class();
  public:

  // optional bool use_cpu_affinity = 1 [default = false];
  bool has_use_cpu_affinity() const;
  private:
  bool _internal_has_use_cpu_affinity() const;
  public:
  void clear_use_cpu_affinity();
  bool use_cpu_affinity() const;
  void set_use_cpu_affinity(bool value);
  private:
  bool _internal_use_cpu_affinity() const;
  void _internal_set_use_cpu_affinity(bool value);
  public:

  // optional int32 cpu = 2 [default = 0];
  bool has_cpu() const;
  private:
  bool _internal_has_cpu() const;
  public:
  void clear_cpu();
  ::PROTOBUF_NAMESPACE_ID::int32 cpu() const;
  void set_cpu(::PROTOBUF_NAMESPACE_ID::int32 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::int32 _internal_cpu() const;
  void _internal_set_cpu(::PROTOBUF_NAMESPACE_ID::int32 value);
  public:

  // @@protoc_insertion_point(class_scope:lm.message.StartOutputWriter)
 private:
  class _Internal;

  // helper for ByteSizeLong()
  size_t RequiredFieldsByteSizeFallback() const;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr output_filename_;
  ::PROTOBUF_NAMESPACE_ID::internal::ArenaStringPtr output_writer_class_;
  bool use_cpu_affinity_;
  ::PROTOBUF_NAMESPACE_ID::int32 cpu_;
  friend struct ::TableStruct_lm_2fmessage_2fStartOutputWriter_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// StartOutputWriter

// optional bool use_cpu_affinity = 1 [default = false];
inline bool StartOutputWriter::_internal_has_use_cpu_affinity() const {
  bool value = (_has_bits_[0] & 0x00000004u) != 0;
  return value;
}
inline bool StartOutputWriter::has_use_cpu_affinity() const {
  return _internal_has_use_cpu_affinity();
}
inline void StartOutputWriter::clear_use_cpu_affinity() {
  use_cpu_affinity_ = false;
  _has_bits_[0] &= ~0x00000004u;
}
inline bool StartOutputWriter::_internal_use_cpu_affinity() const {
  return use_cpu_affinity_;
}
inline bool StartOutputWriter::use_cpu_affinity() const {
  // @@protoc_insertion_point(field_get:lm.message.StartOutputWriter.use_cpu_affinity)
  return _internal_use_cpu_affinity();
}
inline void StartOutputWriter::_internal_set_use_cpu_affinity(bool value) {
  _has_bits_[0] |= 0x00000004u;
  use_cpu_affinity_ = value;
}
inline void StartOutputWriter::set_use_cpu_affinity(bool value) {
  _internal_set_use_cpu_affinity(value);
  // @@protoc_insertion_point(field_set:lm.message.StartOutputWriter.use_cpu_affinity)
}

// optional int32 cpu = 2 [default = 0];
inline bool StartOutputWriter::_internal_has_cpu() const {
  bool value = (_has_bits_[0] & 0x00000008u) != 0;
  return value;
}
inline bool StartOutputWriter::has_cpu() const {
  return _internal_has_cpu();
}
inline void StartOutputWriter::clear_cpu() {
  cpu_ = 0;
  _has_bits_[0] &= ~0x00000008u;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 StartOutputWriter::_internal_cpu() const {
  return cpu_;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 StartOutputWriter::cpu() const {
  // @@protoc_insertion_point(field_get:lm.message.StartOutputWriter.cpu)
  return _internal_cpu();
}
inline void StartOutputWriter::_internal_set_cpu(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _has_bits_[0] |= 0x00000008u;
  cpu_ = value;
}
inline void StartOutputWriter::set_cpu(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_set_cpu(value);
  // @@protoc_insertion_point(field_set:lm.message.StartOutputWriter.cpu)
}

// required string output_filename = 3;
inline bool StartOutputWriter::_internal_has_output_filename() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool StartOutputWriter::has_output_filename() const {
  return _internal_has_output_filename();
}
inline void StartOutputWriter::clear_output_filename() {
  output_filename_.ClearToEmpty(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
  _has_bits_[0] &= ~0x00000001u;
}
inline const std::string& StartOutputWriter::output_filename() const {
  // @@protoc_insertion_point(field_get:lm.message.StartOutputWriter.output_filename)
  return _internal_output_filename();
}
inline void StartOutputWriter::set_output_filename(const std::string& value) {
  _internal_set_output_filename(value);
  // @@protoc_insertion_point(field_set:lm.message.StartOutputWriter.output_filename)
}
inline std::string* StartOutputWriter::mutable_output_filename() {
  // @@protoc_insertion_point(field_mutable:lm.message.StartOutputWriter.output_filename)
  return _internal_mutable_output_filename();
}
inline const std::string& StartOutputWriter::_internal_output_filename() const {
  return output_filename_.Get();
}
inline void StartOutputWriter::_internal_set_output_filename(const std::string& value) {
  _has_bits_[0] |= 0x00000001u;
  output_filename_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), value, GetArena());
}
inline void StartOutputWriter::set_output_filename(std::string&& value) {
  _has_bits_[0] |= 0x00000001u;
  output_filename_.Set(
    &::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), ::std::move(value), GetArena());
  // @@protoc_insertion_point(field_set_rvalue:lm.message.StartOutputWriter.output_filename)
}
inline void StartOutputWriter::set_output_filename(const char* value) {
  GOOGLE_DCHECK(value != nullptr);
  _has_bits_[0] |= 0x00000001u;
  output_filename_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), ::std::string(value),
              GetArena());
  // @@protoc_insertion_point(field_set_char:lm.message.StartOutputWriter.output_filename)
}
inline void StartOutputWriter::set_output_filename(const char* value,
    size_t size) {
  _has_bits_[0] |= 0x00000001u;
  output_filename_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), ::std::string(
      reinterpret_cast<const char*>(value), size), GetArena());
  // @@protoc_insertion_point(field_set_pointer:lm.message.StartOutputWriter.output_filename)
}
inline std::string* StartOutputWriter::_internal_mutable_output_filename() {
  _has_bits_[0] |= 0x00000001u;
  return output_filename_.Mutable(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
}
inline std::string* StartOutputWriter::release_output_filename() {
  // @@protoc_insertion_point(field_release:lm.message.StartOutputWriter.output_filename)
  if (!_internal_has_output_filename()) {
    return nullptr;
  }
  _has_bits_[0] &= ~0x00000001u;
  return output_filename_.ReleaseNonDefault(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
}
inline void StartOutputWriter::set_allocated_output_filename(std::string* output_filename) {
  if (output_filename != nullptr) {
    _has_bits_[0] |= 0x00000001u;
  } else {
    _has_bits_[0] &= ~0x00000001u;
  }
  output_filename_.SetAllocated(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), output_filename,
      GetArena());
  // @@protoc_insertion_point(field_set_allocated:lm.message.StartOutputWriter.output_filename)
}

// required string output_writer_class = 4;
inline bool StartOutputWriter::_internal_has_output_writer_class() const {
  bool value = (_has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool StartOutputWriter::has_output_writer_class() const {
  return _internal_has_output_writer_class();
}
inline void StartOutputWriter::clear_output_writer_class() {
  output_writer_class_.ClearToEmpty(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
  _has_bits_[0] &= ~0x00000002u;
}
inline const std::string& StartOutputWriter::output_writer_class() const {
  // @@protoc_insertion_point(field_get:lm.message.StartOutputWriter.output_writer_class)
  return _internal_output_writer_class();
}
inline void StartOutputWriter::set_output_writer_class(const std::string& value) {
  _internal_set_output_writer_class(value);
  // @@protoc_insertion_point(field_set:lm.message.StartOutputWriter.output_writer_class)
}
inline std::string* StartOutputWriter::mutable_output_writer_class() {
  // @@protoc_insertion_point(field_mutable:lm.message.StartOutputWriter.output_writer_class)
  return _internal_mutable_output_writer_class();
}
inline const std::string& StartOutputWriter::_internal_output_writer_class() const {
  return output_writer_class_.Get();
}
inline void StartOutputWriter::_internal_set_output_writer_class(const std::string& value) {
  _has_bits_[0] |= 0x00000002u;
  output_writer_class_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), value, GetArena());
}
inline void StartOutputWriter::set_output_writer_class(std::string&& value) {
  _has_bits_[0] |= 0x00000002u;
  output_writer_class_.Set(
    &::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), ::std::move(value), GetArena());
  // @@protoc_insertion_point(field_set_rvalue:lm.message.StartOutputWriter.output_writer_class)
}
inline void StartOutputWriter::set_output_writer_class(const char* value) {
  GOOGLE_DCHECK(value != nullptr);
  _has_bits_[0] |= 0x00000002u;
  output_writer_class_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), ::std::string(value),
              GetArena());
  // @@protoc_insertion_point(field_set_char:lm.message.StartOutputWriter.output_writer_class)
}
inline void StartOutputWriter::set_output_writer_class(const char* value,
    size_t size) {
  _has_bits_[0] |= 0x00000002u;
  output_writer_class_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), ::std::string(
      reinterpret_cast<const char*>(value), size), GetArena());
  // @@protoc_insertion_point(field_set_pointer:lm.message.StartOutputWriter.output_writer_class)
}
inline std::string* StartOutputWriter::_internal_mutable_output_writer_class() {
  _has_bits_[0] |= 0x00000002u;
  return output_writer_class_.Mutable(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
}
inline std::string* StartOutputWriter::release_output_writer_class() {
  // @@protoc_insertion_point(field_release:lm.message.StartOutputWriter.output_writer_class)
  if (!_internal_has_output_writer_class()) {
    return nullptr;
  }
  _has_bits_[0] &= ~0x00000002u;
  return output_writer_class_.ReleaseNonDefault(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
}
inline void StartOutputWriter::set_allocated_output_writer_class(std::string* output_writer_class) {
  if (output_writer_class != nullptr) {
    _has_bits_[0] |= 0x00000002u;
  } else {
    _has_bits_[0] &= ~0x00000002u;
  }
  output_writer_class_.SetAllocated(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), output_writer_class,
      GetArena());
  // @@protoc_insertion_point(field_set_allocated:lm.message.StartOutputWriter.output_writer_class)
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace message
}  // namespace lm

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2fmessage_2fStartOutputWriter_2eproto
