// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: lm/io/ffpilot/FFPilotOutput.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto

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
#include "lm/io/ffpilot/FFPilotStageOutput.pb.h"
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto {
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
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto;
namespace lm {
namespace io {
namespace ffpilot {
class FFPilotOutput;
class FFPilotOutputDefaultTypeInternal;
extern FFPilotOutputDefaultTypeInternal _FFPilotOutput_default_instance_;
}  // namespace ffpilot
}  // namespace io
}  // namespace lm
PROTOBUF_NAMESPACE_OPEN
template<> ::lm::io::ffpilot::FFPilotOutput* Arena::CreateMaybeMessage<::lm::io::ffpilot::FFPilotOutput>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace lm {
namespace io {
namespace ffpilot {

// ===================================================================

class FFPilotOutput PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:lm.io.ffpilot.FFPilotOutput) */ {
 public:
  inline FFPilotOutput() : FFPilotOutput(nullptr) {}
  virtual ~FFPilotOutput();

  FFPilotOutput(const FFPilotOutput& from);
  FFPilotOutput(FFPilotOutput&& from) noexcept
    : FFPilotOutput() {
    *this = ::std::move(from);
  }

  inline FFPilotOutput& operator=(const FFPilotOutput& from) {
    CopyFrom(from);
    return *this;
  }
  inline FFPilotOutput& operator=(FFPilotOutput&& from) noexcept {
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
  static const FFPilotOutput& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const FFPilotOutput* internal_default_instance() {
    return reinterpret_cast<const FFPilotOutput*>(
               &_FFPilotOutput_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(FFPilotOutput& a, FFPilotOutput& b) {
    a.Swap(&b);
  }
  inline void Swap(FFPilotOutput* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(FFPilotOutput* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline FFPilotOutput* New() const final {
    return CreateMaybeMessage<FFPilotOutput>(nullptr);
  }

  FFPilotOutput* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<FFPilotOutput>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const FFPilotOutput& from);
  void MergeFrom(const FFPilotOutput& from);
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
  void InternalSwap(FFPilotOutput* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "lm.io.ffpilot.FFPilotOutput";
  }
  protected:
  explicit FFPilotOutput(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto);
    return ::descriptor_table_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kStageOutputFieldNumber = 100,
  };
  // repeated .lm.io.ffpilot.FFPilotStageOutput stage_output = 100;
  int stage_output_size() const;
  private:
  int _internal_stage_output_size() const;
  public:
  void clear_stage_output();
  ::lm::io::ffpilot::FFPilotStageOutput* mutable_stage_output(int index);
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotStageOutput >*
      mutable_stage_output();
  private:
  const ::lm::io::ffpilot::FFPilotStageOutput& _internal_stage_output(int index) const;
  ::lm::io::ffpilot::FFPilotStageOutput* _internal_add_stage_output();
  public:
  const ::lm::io::ffpilot::FFPilotStageOutput& stage_output(int index) const;
  ::lm::io::ffpilot::FFPilotStageOutput* add_stage_output();
  const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotStageOutput >&
      stage_output() const;

  // @@protoc_insertion_point(class_scope:lm.io.ffpilot.FFPilotOutput)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotStageOutput > stage_output_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  friend struct ::TableStruct_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// FFPilotOutput

// repeated .lm.io.ffpilot.FFPilotStageOutput stage_output = 100;
inline int FFPilotOutput::_internal_stage_output_size() const {
  return stage_output_.size();
}
inline int FFPilotOutput::stage_output_size() const {
  return _internal_stage_output_size();
}
inline ::lm::io::ffpilot::FFPilotStageOutput* FFPilotOutput::mutable_stage_output(int index) {
  // @@protoc_insertion_point(field_mutable:lm.io.ffpilot.FFPilotOutput.stage_output)
  return stage_output_.Mutable(index);
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotStageOutput >*
FFPilotOutput::mutable_stage_output() {
  // @@protoc_insertion_point(field_mutable_list:lm.io.ffpilot.FFPilotOutput.stage_output)
  return &stage_output_;
}
inline const ::lm::io::ffpilot::FFPilotStageOutput& FFPilotOutput::_internal_stage_output(int index) const {
  return stage_output_.Get(index);
}
inline const ::lm::io::ffpilot::FFPilotStageOutput& FFPilotOutput::stage_output(int index) const {
  // @@protoc_insertion_point(field_get:lm.io.ffpilot.FFPilotOutput.stage_output)
  return _internal_stage_output(index);
}
inline ::lm::io::ffpilot::FFPilotStageOutput* FFPilotOutput::_internal_add_stage_output() {
  return stage_output_.Add();
}
inline ::lm::io::ffpilot::FFPilotStageOutput* FFPilotOutput::add_stage_output() {
  // @@protoc_insertion_point(field_add:lm.io.ffpilot.FFPilotOutput.stage_output)
  return _internal_add_stage_output();
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::lm::io::ffpilot::FFPilotStageOutput >&
FFPilotOutput::stage_output() const {
  // @@protoc_insertion_point(field_list:lm.io.ffpilot.FFPilotOutput.stage_output)
  return stage_output_;
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__

// @@protoc_insertion_point(namespace_scope)

}  // namespace ffpilot
}  // namespace io
}  // namespace lm

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_lm_2fio_2fffpilot_2fFFPilotOutput_2eproto