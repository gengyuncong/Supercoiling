/*
 * Copyright 2019 Johns Hopkins University
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
 * Author(s): Elijah Roberts
 */

#include <cstdio>
#include <stdexcept>
#include <string>
#include <sys/stat.h>

#include "lm/Exceptions.h"
#include "lm/io/File.h"

using std::string;

namespace lm {
namespace io {

File::File(const string filename)
:filename(filename),fp(NULL),line(NULL),lineSize(0)
{
    filestatsRet=stat(filename.c_str(), &filestats);
}

File::~File()
{
    close();
}

bool File::isFile()
{
    return S_ISREG(filestats.st_mode);

}

bool File::isDir()
{
    return S_ISDIR(filestats.st_mode);
}

bool File::exists()
{
    return (filestatsRet==0);
}

string File::getFilename()
{
    return filename;
}

int64_t File::getSize()
{
    return filestats.st_size;
}

void File::openTruncate()
{
    if (fp == NULL)
    {
        fp = fopen(filename.c_str(), "w");
        if (fp == NULL) throw lm::IOException(string("Could not open file for writing truncate: ")+filename);
    }
}

void File::openAppend()
{
    if (fp == NULL)
    {
        fp = fopen(filename.c_str(), "a");
        if (fp == NULL) throw lm::IOException(string("Could not open file for writing append: ")+filename);
    }
}

void File::openRead()
{
    if (fp == NULL)
    {
        fp = fopen(filename.c_str(), "r");
        if (fp == NULL) throw lm::IOException(string("Could not open file for reading: ")+filename);
    }
}

size_t File::read(void* buffer, size_t count)
{
    if (fp != NULL)
    {
        size_t ret=fread(buffer, sizeof(unsigned char), count, fp);
        if (count != ret && ferror(fp) != 0) throw lm::IOException(string("Error reading from file: ")+filename);
        return ret;
    }
    return 0;
}

void File::readFully(void* vbuffer, size_t length)
{
    unsigned char* buffer = (unsigned char*)vbuffer;
    size_t bytesRead;
    do
    {
        bytesRead=read(buffer, length);
        buffer += bytesRead;
        length -= bytesRead;
    }
    while (length > 0 && bytesRead != 0);
}

std::string File::readLine()
{
    if (fp != NULL)
    {
        ssize_t ret=getline(&line, &lineSize, fp);
        if (ret < 0) THROW_EXCEPTION(lm::RuntimeException, string("Error reading from file: ")+filename);

        // Remove any trailing newline and carriage returns.
        while (ret > 0 && (line[ret-1] == '\n' || line[ret-1] == '\r'))
        {
            line[ret-1] = 0;
            ret--;
        }

        return string(line);
    }
    THROW_EXCEPTION(lm::RuntimeException, "File not open.");
}

void File::skip(uint64_t length)
{
    if (fp != NULL)
    {
        if (fseek(fp, length, SEEK_CUR) != 0) throw lm::IOException(string("Error seeking in file: ")+filename);
    }
}

bool File::isEof()
{
    if (fp != NULL)
    {
        return (ftell(fp) == getSize());
    }
    return true;
}

size_t File::write(void* buffer, size_t count)
{
    if (fp != NULL)
    {
        size_t ret=fwrite(buffer, sizeof(unsigned char), count, fp);
        if (count != ret) throw lm::IOException(string("Error writing to file: ")+filename);
        return ret;
    }
    return 0;
}

void File::flush()
{
    if (fp != NULL)
    {
        int ret=fflush(fp);
        if (ret != 0) throw lm::IOException(string("Could not flush file: ")+filename);
    }
}

void File::close()
{
    if (fp != NULL)
    {
        int ret=fclose(fp);
        fp = NULL;
        if (ret != 0) throw lm::IOException(string("Could not close file: ")+filename);
    }

    if (line != NULL)
    {
        free(line);
        line = NULL;
        lineSize = 0;
    }
}

}
}
