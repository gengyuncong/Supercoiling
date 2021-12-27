/*
 * Copyright 2012-2019 Johns Hopkins University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Developed by: Roberts Group
 *               Johns Hopkins University
 *               http://biophysics.jhu.edu/roberts/
 *
 * Author(s): Elijah Roberts, Max Klein
 */

#ifndef LM_LIMIT_LIMITCHECKMACROS
#define LM_LIMIT_LIMITCHECKMACROS

// non-standard compliant check macro, similar to the inequality checking macros from the avx library
//define check_limit_MIN_EXCLUSIVE(val, limitVal) __extension__ ({ (val < limitVal); })

// specializations of check_limit with regards to stoppingCondition and includeEndpoint for the basic min/max limits
#define check_limit_MIN_EXCLUSIVE(val, limitVal, checkBool) checkBool = (val <  limitVal);
#define check_limit_MIN_INCLUSIVE( val, limitVal, checkBool) checkBool = (val <= limitVal);
#define check_limit_MAX_EXCLUSIVE(val, limitVal, checkBool) checkBool = (val >  limitVal);
#define check_limit_MAX_INCLUSIVE( val, limitVal, checkBool) checkBool = (val >= limitVal);

// specializations of check_limit with regards to stoppingCondition and includeEndpoint for second degree limits (limits that depend on both previous and present value)
#define check_limit_CROSSED_EXCLUSIVE(prevVal, val, limitVal, checkBool) checkBool = ((prevVal >= limitVal && val <  limitVal) || (prevVal <= limitVal && val >  limitVal));
#define check_limit_CROSSED_INCLUSIVE( prevVal, val, limitVal, checkBool) checkBool = ((prevVal >  limitVal && val <= limitVal) || (prevVal <  limitVal && val >= limitVal));
#define check_limit_CROSSED_DECREASING_EXCLUSIVE(prevVal, val, limitVal, checkBool) checkBool = (prevVal >= limitVal && val <  limitVal);
#define check_limit_CROSSED_DECREASING_INCLUSIVE( prevVal, val, limitVal, checkBool) checkBool = (prevVal >  limitVal && val <= limitVal);
#define check_limit_CROSSED_INCREASING_EXCLUSIVE(prevVal, val, limitVal, checkBool) checkBool = (prevVal <= limitVal && val >  limitVal);
#define check_limit_CROSSED_INCREASING_INCLUSIVE( prevVal, val, limitVal, checkBool) checkBool = (prevVal <  limitVal && val >= limitVal);

#ifdef OPT_AVX
#include <immintrin.h>
// AVX versions
// specializations of check_limit with regards to stoppingCondition and includeEndpoint for the basic min/max limits
#define check_simple_limit_avx(valueArr, valueID, limitValue, tmpBool, checkBool, opCode) \
    tmpBool   = _mm256_cmp_pd(_mm256_load_pd(&valueArr[valueID*DOUBLES_PER_AVX]), limitValue, opCode); \
    checkBool = _mm256_movemask_pd(tmpBool);

#define check_limit_avx_MIN_EXCLUSIVE(valueArr, valueID, limitValue, tmpBool, checkBool) check_simple_limit_avx(valueArr, valueID, limitValue, tmpBool, checkBool, _CMP_LT_OQ)
#define check_limit_avx_MIN_INCLUSIVE( valueArr, valueID, limitValue, tmpBool, checkBool) check_simple_limit_avx(valueArr, valueID, limitValue, tmpBool, checkBool, _CMP_LE_OQ)
#define check_limit_avx_MAX_EXCLUSIVE(valueArr, valueID, limitValue, tmpBool, checkBool) check_simple_limit_avx(valueArr, valueID, limitValue, tmpBool, checkBool, _CMP_GT_OQ)
#define check_limit_avx_MAX_INCLUSIVE( valueArr, valueID, limitValue, tmpBool, checkBool) check_simple_limit_avx(valueArr, valueID, limitValue, tmpBool, checkBool, _CMP_GE_OQ)

// specializations of check_limit with regards to stoppingCondition and includeEndpoint for second degree limits (limits that depend on both previous and present value)
#define check_second_degree_limit_avx(previousValueArr, valueArr, valueID, limitValue, previousTmpBool, tmpBool, checkBool, previousOpCode, opCode) \
    previousTmpBool = _mm256_cmp_pd(_mm256_load_pd(&previousValueArr[valueID*DOUBLES_PER_AVX]), limitValue, previousOpCode); \
    tmpBool         = _mm256_cmp_pd(_mm256_load_pd(&valueArr[valueID*DOUBLES_PER_AVX]), limitValue, opCode); \
    checkBool       = _mm256_movemask_pd(previousTmpBool)&_mm256_movemask_pd(tmpBool);

#define check_limit_avx_DECREASING_EXCLUSIVE(previousValueArr, valueArr, valueID, limitValue, previousTmpBool, tmpBool, checkBool) check_second_degree_limit_avx(previousValueArr, valueArr, valueID, limitValue, previousTmpBool, tmpBool, checkBool, _CMP_GE_OQ, _CMP_LT_OQ)
#define check_limit_avx_DECREASING_INCLUSIVE( previousValueArr, valueArr, valueID, limitValue, previousTmpBool, tmpBool, checkBool) check_second_degree_limit_avx(previousValueArr, valueArr, valueID, limitValue, previousTmpBool, tmpBool, checkBool, _CMP_GT_OQ, _CMP_LE_OQ)
#define check_limit_avx_INCREASING_EXCLUSIVE(previousValueArr, valueArr, valueID, limitValue, previousTmpBool, tmpBool, checkBool) check_second_degree_limit_avx(previousValueArr, valueArr, valueID, limitValue, previousTmpBool, tmpBool, checkBool, _CMP_LE_OQ, _CMP_GT_OQ)
#define check_limit_avx_INCREASING_INCLUSIVE( previousValueArr, valueArr, valueID, limitValue, previousTmpBool, tmpBool, checkBool) check_second_degree_limit_avx(previousValueArr, valueArr, valueID, limitValue, previousTmpBool, tmpBool, checkBool, _CMP_LT_OQ, _CMP_GE_OQ)
#endif /* OPT_AVX */

#endif /* LM_LIMIT_LIMITCHECKMACROS */
