/*
 * Copyright 2016 Johns Hopkins University
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

#include "robertslab/Types.h"

template<> const char* printf_format_string<int8_t>()   {return "% 2d";}
template<> const char* printf_format_string<int32_t>()  {return "% 2d";}
template<> const char* printf_format_string<int64_t>()  {return "% 2lld";}
template<> const char* printf_format_string<uint8_t>()  {return "%u";}
template<> const char* printf_format_string<uint32_t>() {return "%u";}
template<> const char* printf_format_string<uint64_t>() {return "%llu";}
template<> const char* printf_format_string<float>()    {return "%0.3e";}
template<> const char* printf_format_string<double>()   {return "%0.3e";}
