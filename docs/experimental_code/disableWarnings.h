/* Haven't been able to get this to work as a header, but it does have a lot of good compiler-specific
 * info about various annoying warnings.
 *
 * Based on code found at https://eigen.tuxfamily.org/dox/DisableStupidWarnings_8h_source.html */

#ifndef LMES_WARNINGS_DISABLED
#define LMES_WARNINGS_DISABLED

#ifdef _MSC_VER
    /** 4100 - unreferenced formal parameter (occurred e.g. in aligned_allocator::destroy(pointer p))
    4101 - unreferenced local variable
    4127 - conditional expression is constant
    4181 - qualifier applied to reference type ignored
    4211 - nonstandard extension used : redefined extern to static
    4244 - 'argument' : conversion from 'type1' to 'type2', possible loss of data
    4273 - QtAlignedMalloc, inconsistent DLL linkage
    4324 - structure was padded due to declspec(align())
    4503 - decorated name length exceeded, name was truncated
    4512 - assignment operator could not be generated
    4522 - 'class' : multiple assignment operators specified
    4700 - uninitialized local variable 'xyz' used
    4714 - function marked as __forceinline not inlined
    4717 - 'function' : recursive on all control paths, function will cause runtime stack overflow
    4800 - 'type' : forcing value to bool 'true' or 'false' (performance warning) */
//    #pragma warning( push )
//    #pragma warning( disable : 4100 4101 4127 4181 4211 4244 4273 4324 4503 4512 4522 4700 4714 4717 4800)

#elif defined __INTEL_COMPILER
    /** 858  - type qualifier on return type is meaningless
           Compiler generates this warning when it finds a type qualifier applied to the return type of a function.
           C++ allows the type-qualifier to be const or volatile , but either qualifier applied to a function return type is meaningless,
           because functions can only return rvalues and the type qualifiers apply only to lvalues.
    2196 - routine is both "inline" and "noinline" ("noinline" assumed)
           ICC 12 generates this warning even without any inline keyword, when defining class methods 'inline' i.e. inside of class body
           typedef that may be a reference type.
    279  - controlling expression is constant
           ICC 12 generates this warning on assert(constant_expression_depending_on_template_params) and frankly this is a legitimate use case.
    1684 - conversion from pointer to same-sized integral type (potential portability problem)
    2259 - non-pointer conversion from "Eigen::Index={ptrdiff_t={long}}" to "int" may lose significant bits */
    #pragma warning push
    #pragma warning disable 858 2196

#elif defined __clang__
    /** -Wconstant-logical-operand - warning: use of logical && with constant operand; switch to bitwise & or remove constant
        this is really a stupid warning as it warns on compile-time expressions involving enums */
//    #pragma clang diagnostic push
//    #pragma clang diagnostic ignored "-Wconstant-logical-operand"

#elif defined __GNUC__ && __GNUC__>=6
//    #pragma GCC diagnostic push
//    #pragma GCC diagnostic ignored "-Wignored-attributes"

#endif

#if defined __NVCC__
//  // Disable the "statement is unreachable" message
//  #pragma diag_suppress code_is_unreachable
//  // Disable the "dynamic initialization in unreachable code" message
//  #pragma diag_suppress initialization_not_reachable
//  // Disable the "invalid error number" message that we get with older versions of nvcc
//  #pragma diag_suppress 1222
//  // Disable the "calling a __host__ function from a __host__ __device__ function is not allowed" messages (yes, there are many of them and they seem to change with every version of the compiler)
//  #pragma diag_suppress 2527
//  #pragma diag_suppress 2529
//  #pragma diag_suppress 2651
//  #pragma diag_suppress 2653
//  #pragma diag_suppress 2668
//  #pragma diag_suppress 2669
//  #pragma diag_suppress 2670
//  #pragma diag_suppress 2671
//  #pragma diag_suppress 2735
//  #pragma diag_suppress 2737
#endif

/** reenable warnings */
//#if defined __clang__
//#pragma clang diagnostic pop
//#endif /* defined __clang__ */

#endif /* LMES_WARNINGS_DISABLED */
