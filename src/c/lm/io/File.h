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


#ifndef LM_IO_FILE_H
#define LM_IO_FILE_H

#include <cstdio>
#include <map>
#include <string>
#include <sys/stat.h>


using std::string;

namespace lm {
namespace io {

class File
{
public:
    File(const string filename);
    virtual ~File();
    virtual bool isFile();
    virtual bool isDir();
    virtual bool exists();
    virtual string getFilename();
    virtual int64_t getSize();
    virtual void openTruncate();
    virtual void openAppend();
    virtual void openRead();
    virtual size_t read(void* buffer, size_t length);
    virtual void readFully(void* buffer, size_t length);
    virtual string readLine();
    virtual void skip(uint64_t length);
    virtual bool isEof();
    virtual size_t write(void* buffer, size_t length);
    virtual void flush();
    virtual void close();

protected:
    string filename;
    struct stat filestats;
    int filestatsRet;
    FILE* fp;
    char* line;
    size_t lineSize;
};

}
}

#endif
