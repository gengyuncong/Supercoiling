// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: robertslab/pbuf/NDArray.proto

#include "robertslab/pbuf/NDArray.pb.h"

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
namespace robertslab {
namespace pbuf {
class NDArrayDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<NDArray> _instance;
} _NDArray_default_instance_;
}  // namespace pbuf
}  // namespace robertslab
static void InitDefaultsscc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::robertslab::pbuf::_NDArray_default_instance_;
    new (ptr) ::robertslab::pbuf::NDArray();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::robertslab::pbuf::NDArray::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_robertslab_2fpbuf_2fNDArray_2eproto[1];
static const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* file_level_enum_descriptors_robertslab_2fpbuf_2fNDArray_2eproto[3];
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_robertslab_2fpbuf_2fNDArray_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_robertslab_2fpbuf_2fNDArray_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::robertslab::pbuf::NDArray, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::robertslab::pbuf::NDArray, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::robertslab::pbuf::NDArray, array_order_),
  PROTOBUF_FIELD_OFFSET(::robertslab::pbuf::NDArray, byte_order_),
  PROTOBUF_FIELD_OFFSET(::robertslab::pbuf::NDArray, data_type_),
  PROTOBUF_FIELD_OFFSET(::robertslab::pbuf::NDArray, shape_),
  PROTOBUF_FIELD_OFFSET(::robertslab::pbuf::NDArray, data_),
  PROTOBUF_FIELD_OFFSET(::robertslab::pbuf::NDArray, compressed_deflate_),
  PROTOBUF_FIELD_OFFSET(::robertslab::pbuf::NDArray, compressed_snappy_),
  1,
  2,
  3,
  ~0u,
  0,
  4,
  5,
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 12, sizeof(::robertslab::pbuf::NDArray)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::robertslab::pbuf::_NDArray_default_instance_),
};

const char descriptor_table_protodef_robertslab_2fpbuf_2fNDArray_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\035robertslab/pbuf/NDArray.proto\022\017roberts"
  "lab.pbuf\"\210\005\n\007NDArray\022C\n\013array_order\030\001 \001("
  "\0162#.robertslab.pbuf.NDArray.ArrayOrder:\t"
  "ROW_MAJOR\022K\n\nbyte_order\030\002 \001(\0162\".robertsl"
  "ab.pbuf.NDArray.ByteOrder:\023LITTLE_ENDIAN"
  "_ORDER\0224\n\tdata_type\030\003 \002(\0162!.robertslab.p"
  "buf.NDArray.DataType\022\r\n\005shape\030\004 \003(\r\022\014\n\004d"
  "ata\030\005 \002(\014\022!\n\022compressed_deflate\030\006 \001(\010:\005f"
  "alse\022 \n\021compressed_snappy\030\007 \001(\010:\005false\"="
  "\n\nArrayOrder\022\r\n\tROW_MAJOR\020\000\022\020\n\014COLUMN_MA"
  "JOR\020\001\022\016\n\nIMPL_ORDER\020\002\":\n\tByteOrder\022\027\n\023LI"
  "TTLE_ENDIAN_ORDER\020\000\022\024\n\020BIG_ENDIAN_ORDER\020"
  "\001\"\327\001\n\010DataType\022\010\n\004int8\020\000\022\t\n\005int16\020\001\022\t\n\005i"
  "nt32\020\002\022\t\n\005int64\020\003\022\t\n\005uint8\020\004\022\n\n\006uint16\020\005"
  "\022\n\n\006uint32\020\006\022\n\n\006uint64\020\007\022\013\n\007float16\020\010\022\013\n"
  "\007float32\020\t\022\013\n\007float64\020\n\022\r\n\tcomplex64\020\013\022\016"
  "\n\ncomplex128\020\014\022\006\n\002S8\020\r\022\007\n\003S16\020\016\022\007\n\003S32\020\017"
  "\022\007\n\003S64\020\020\022\010\n\004S128\020\021"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto_sccs[1] = {
  &scc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto = {
  false, false, descriptor_table_protodef_robertslab_2fpbuf_2fNDArray_2eproto, "robertslab/pbuf/NDArray.proto", 699,
  &descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto_once, descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto_sccs, descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto_deps, 1, 0,
  schemas, file_default_instances, TableStruct_robertslab_2fpbuf_2fNDArray_2eproto::offsets,
  file_level_metadata_robertslab_2fpbuf_2fNDArray_2eproto, 1, file_level_enum_descriptors_robertslab_2fpbuf_2fNDArray_2eproto, file_level_service_descriptors_robertslab_2fpbuf_2fNDArray_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_robertslab_2fpbuf_2fNDArray_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto)), true);
namespace robertslab {
namespace pbuf {
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* NDArray_ArrayOrder_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto);
  return file_level_enum_descriptors_robertslab_2fpbuf_2fNDArray_2eproto[0];
}
bool NDArray_ArrayOrder_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
    case 2:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)
constexpr NDArray_ArrayOrder NDArray::ROW_MAJOR;
constexpr NDArray_ArrayOrder NDArray::COLUMN_MAJOR;
constexpr NDArray_ArrayOrder NDArray::IMPL_ORDER;
constexpr NDArray_ArrayOrder NDArray::ArrayOrder_MIN;
constexpr NDArray_ArrayOrder NDArray::ArrayOrder_MAX;
constexpr int NDArray::ArrayOrder_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* NDArray_ByteOrder_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto);
  return file_level_enum_descriptors_robertslab_2fpbuf_2fNDArray_2eproto[1];
}
bool NDArray_ByteOrder_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)
constexpr NDArray_ByteOrder NDArray::LITTLE_ENDIAN_ORDER;
constexpr NDArray_ByteOrder NDArray::BIG_ENDIAN_ORDER;
constexpr NDArray_ByteOrder NDArray::ByteOrder_MIN;
constexpr NDArray_ByteOrder NDArray::ByteOrder_MAX;
constexpr int NDArray::ByteOrder_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)
const ::PROTOBUF_NAMESPACE_ID::EnumDescriptor* NDArray_DataType_descriptor() {
  ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&descriptor_table_robertslab_2fpbuf_2fNDArray_2eproto);
  return file_level_enum_descriptors_robertslab_2fpbuf_2fNDArray_2eproto[2];
}
bool NDArray_DataType_IsValid(int value) {
  switch (value) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
    case 17:
      return true;
    default:
      return false;
  }
}

#if (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)
constexpr NDArray_DataType NDArray::int8;
constexpr NDArray_DataType NDArray::int16;
constexpr NDArray_DataType NDArray::int32;
constexpr NDArray_DataType NDArray::int64;
constexpr NDArray_DataType NDArray::uint8;
constexpr NDArray_DataType NDArray::uint16;
constexpr NDArray_DataType NDArray::uint32;
constexpr NDArray_DataType NDArray::uint64;
constexpr NDArray_DataType NDArray::float16;
constexpr NDArray_DataType NDArray::float32;
constexpr NDArray_DataType NDArray::float64;
constexpr NDArray_DataType NDArray::complex64;
constexpr NDArray_DataType NDArray::complex128;
constexpr NDArray_DataType NDArray::S8;
constexpr NDArray_DataType NDArray::S16;
constexpr NDArray_DataType NDArray::S32;
constexpr NDArray_DataType NDArray::S64;
constexpr NDArray_DataType NDArray::S128;
constexpr NDArray_DataType NDArray::DataType_MIN;
constexpr NDArray_DataType NDArray::DataType_MAX;
constexpr int NDArray::DataType_ARRAYSIZE;
#endif  // (__cplusplus < 201703) && (!defined(_MSC_VER) || _MSC_VER >= 1900)

// ===================================================================

void NDArray::InitAsDefaultInstance() {
}
class NDArray::_Internal {
 public:
  using HasBits = decltype(std::declval<NDArray>()._has_bits_);
  static void set_has_array_order(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static void set_has_byte_order(HasBits* has_bits) {
    (*has_bits)[0] |= 4u;
  }
  static void set_has_data_type(HasBits* has_bits) {
    (*has_bits)[0] |= 8u;
  }
  static void set_has_data(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static void set_has_compressed_deflate(HasBits* has_bits) {
    (*has_bits)[0] |= 16u;
  }
  static void set_has_compressed_snappy(HasBits* has_bits) {
    (*has_bits)[0] |= 32u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000009) ^ 0x00000009) != 0;
  }
};

NDArray::NDArray(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  shape_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:robertslab.pbuf.NDArray)
}
NDArray::NDArray(const NDArray& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_),
      shape_(from.shape_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  data_.UnsafeSetDefault(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited());
  if (from._internal_has_data()) {
    data_.Set(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), from._internal_data(),
      GetArena());
  }
  ::memcpy(&array_order_, &from.array_order_,
    static_cast<size_t>(reinterpret_cast<char*>(&compressed_snappy_) -
    reinterpret_cast<char*>(&array_order_)) + sizeof(compressed_snappy_));
  // @@protoc_insertion_point(copy_constructor:robertslab.pbuf.NDArray)
}

void NDArray::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto.base);
  data_.UnsafeSetDefault(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited());
  ::memset(&array_order_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&compressed_snappy_) -
      reinterpret_cast<char*>(&array_order_)) + sizeof(compressed_snappy_));
}

NDArray::~NDArray() {
  // @@protoc_insertion_point(destructor:robertslab.pbuf.NDArray)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void NDArray::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
  data_.DestroyNoArena(&::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited());
}

void NDArray::ArenaDtor(void* object) {
  NDArray* _this = reinterpret_cast< NDArray* >(object);
  (void)_this;
}
void NDArray::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void NDArray::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const NDArray& NDArray::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_NDArray_robertslab_2fpbuf_2fNDArray_2eproto.base);
  return *internal_default_instance();
}


void NDArray::Clear() {
// @@protoc_insertion_point(message_clear_start:robertslab.pbuf.NDArray)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  shape_.Clear();
  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000001u) {
    data_.ClearNonDefaultToEmpty();
  }
  if (cached_has_bits & 0x0000003eu) {
    ::memset(&array_order_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&compressed_snappy_) -
        reinterpret_cast<char*>(&array_order_)) + sizeof(compressed_snappy_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* NDArray::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // optional .robertslab.pbuf.NDArray.ArrayOrder array_order = 1 [default = ROW_MAJOR];
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          ::PROTOBUF_NAMESPACE_ID::uint64 val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::robertslab::pbuf::NDArray_ArrayOrder_IsValid(val))) {
            _internal_set_array_order(static_cast<::robertslab::pbuf::NDArray_ArrayOrder>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(1, val, mutable_unknown_fields());
          }
        } else goto handle_unusual;
        continue;
      // optional .robertslab.pbuf.NDArray.ByteOrder byte_order = 2 [default = LITTLE_ENDIAN_ORDER];
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 16)) {
          ::PROTOBUF_NAMESPACE_ID::uint64 val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::robertslab::pbuf::NDArray_ByteOrder_IsValid(val))) {
            _internal_set_byte_order(static_cast<::robertslab::pbuf::NDArray_ByteOrder>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(2, val, mutable_unknown_fields());
          }
        } else goto handle_unusual;
        continue;
      // required .robertslab.pbuf.NDArray.DataType data_type = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 24)) {
          ::PROTOBUF_NAMESPACE_ID::uint64 val = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
          if (PROTOBUF_PREDICT_TRUE(::robertslab::pbuf::NDArray_DataType_IsValid(val))) {
            _internal_set_data_type(static_cast<::robertslab::pbuf::NDArray_DataType>(val));
          } else {
            ::PROTOBUF_NAMESPACE_ID::internal::WriteVarint(3, val, mutable_unknown_fields());
          }
        } else goto handle_unusual;
        continue;
      // repeated uint32 shape = 4;
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 32)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_shape(::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr));
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<32>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 34) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedUInt32Parser(_internal_mutable_shape(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // required bytes data = 5;
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 42)) {
          auto str = _internal_mutable_data();
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::InlineGreedyStringParser(str, ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional bool compressed_deflate = 6 [default = false];
      case 6:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 48)) {
          _Internal::set_has_compressed_deflate(&has_bits);
          compressed_deflate_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // optional bool compressed_snappy = 7 [default = false];
      case 7:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 56)) {
          _Internal::set_has_compressed_snappy(&has_bits);
          compressed_snappy_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
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

::PROTOBUF_NAMESPACE_ID::uint8* NDArray::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:robertslab.pbuf.NDArray)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // optional .robertslab.pbuf.NDArray.ArrayOrder array_order = 1 [default = ROW_MAJOR];
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteEnumToArray(
      1, this->_internal_array_order(), target);
  }

  // optional .robertslab.pbuf.NDArray.ByteOrder byte_order = 2 [default = LITTLE_ENDIAN_ORDER];
  if (cached_has_bits & 0x00000004u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteEnumToArray(
      2, this->_internal_byte_order(), target);
  }

  // required .robertslab.pbuf.NDArray.DataType data_type = 3;
  if (cached_has_bits & 0x00000008u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteEnumToArray(
      3, this->_internal_data_type(), target);
  }

  // repeated uint32 shape = 4;
  for (int i = 0, n = this->_internal_shape_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteUInt32ToArray(4, this->_internal_shape(i), target);
  }

  // required bytes data = 5;
  if (cached_has_bits & 0x00000001u) {
    target = stream->WriteBytesMaybeAliased(
        5, this->_internal_data(), target);
  }

  // optional bool compressed_deflate = 6 [default = false];
  if (cached_has_bits & 0x00000010u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteBoolToArray(6, this->_internal_compressed_deflate(), target);
  }

  // optional bool compressed_snappy = 7 [default = false];
  if (cached_has_bits & 0x00000020u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteBoolToArray(7, this->_internal_compressed_snappy(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:robertslab.pbuf.NDArray)
  return target;
}

size_t NDArray::RequiredFieldsByteSizeFallback() const {
// @@protoc_insertion_point(required_fields_byte_size_fallback_start:robertslab.pbuf.NDArray)
  size_t total_size = 0;

  if (_internal_has_data()) {
    // required bytes data = 5;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::BytesSize(
        this->_internal_data());
  }

  if (_internal_has_data_type()) {
    // required .robertslab.pbuf.NDArray.DataType data_type = 3;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::EnumSize(this->_internal_data_type());
  }

  return total_size;
}
size_t NDArray::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:robertslab.pbuf.NDArray)
  size_t total_size = 0;

  if (((_has_bits_[0] & 0x00000009) ^ 0x00000009) == 0) {  // All required fields are present.
    // required bytes data = 5;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::BytesSize(
        this->_internal_data());

    // required .robertslab.pbuf.NDArray.DataType data_type = 3;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::EnumSize(this->_internal_data_type());

  } else {
    total_size += RequiredFieldsByteSizeFallback();
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated uint32 shape = 4;
  {
    size_t data_size = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      UInt32Size(this->shape_);
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_shape_size());
    total_size += data_size;
  }

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000006u) {
    // optional .robertslab.pbuf.NDArray.ArrayOrder array_order = 1 [default = ROW_MAJOR];
    if (cached_has_bits & 0x00000002u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::EnumSize(this->_internal_array_order());
    }

    // optional .robertslab.pbuf.NDArray.ByteOrder byte_order = 2 [default = LITTLE_ENDIAN_ORDER];
    if (cached_has_bits & 0x00000004u) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::EnumSize(this->_internal_byte_order());
    }

  }
  if (cached_has_bits & 0x00000030u) {
    // optional bool compressed_deflate = 6 [default = false];
    if (cached_has_bits & 0x00000010u) {
      total_size += 1 + 1;
    }

    // optional bool compressed_snappy = 7 [default = false];
    if (cached_has_bits & 0x00000020u) {
      total_size += 1 + 1;
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

void NDArray::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:robertslab.pbuf.NDArray)
  GOOGLE_DCHECK_NE(&from, this);
  const NDArray* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<NDArray>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:robertslab.pbuf.NDArray)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:robertslab.pbuf.NDArray)
    MergeFrom(*source);
  }
}

void NDArray::MergeFrom(const NDArray& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:robertslab.pbuf.NDArray)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  shape_.MergeFrom(from.shape_);
  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x0000003fu) {
    if (cached_has_bits & 0x00000001u) {
      _internal_set_data(from._internal_data());
    }
    if (cached_has_bits & 0x00000002u) {
      array_order_ = from.array_order_;
    }
    if (cached_has_bits & 0x00000004u) {
      byte_order_ = from.byte_order_;
    }
    if (cached_has_bits & 0x00000008u) {
      data_type_ = from.data_type_;
    }
    if (cached_has_bits & 0x00000010u) {
      compressed_deflate_ = from.compressed_deflate_;
    }
    if (cached_has_bits & 0x00000020u) {
      compressed_snappy_ = from.compressed_snappy_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void NDArray::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:robertslab.pbuf.NDArray)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void NDArray::CopyFrom(const NDArray& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:robertslab.pbuf.NDArray)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool NDArray::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  return true;
}

void NDArray::InternalSwap(NDArray* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  shape_.InternalSwap(&other->shape_);
  data_.Swap(&other->data_, &::PROTOBUF_NAMESPACE_ID::internal::GetEmptyStringAlreadyInited(), GetArena());
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(NDArray, compressed_snappy_)
      + sizeof(NDArray::compressed_snappy_)
      - PROTOBUF_FIELD_OFFSET(NDArray, array_order_)>(
          reinterpret_cast<char*>(&array_order_),
          reinterpret_cast<char*>(&other->array_order_));
}

::PROTOBUF_NAMESPACE_ID::Metadata NDArray::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace pbuf
}  // namespace robertslab
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::robertslab::pbuf::NDArray* Arena::CreateMaybeMessage< ::robertslab::pbuf::NDArray >(Arena* arena) {
  return Arena::CreateMessageInternal< ::robertslab::pbuf::NDArray >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>