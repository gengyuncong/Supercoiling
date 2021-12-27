/*
 * University of Illinois Open Source License
 * Copyright 2012-2016 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to
 * do so, subject to the following conditions:
 *
 * - Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimers.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimers in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the names of the Roberts Group, Johns Hopkins University,
 * nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts, Max Klein
 */

#ifndef LM_MESSAGE_COMMUNICATOR_H
#define LM_MESSAGE_COMMUNICATOR_H

#include <string>
#include <google/protobuf/message.h>

#include "lm/message/Endpoint.pb.h"
#include "lm/message/Message.pb.h"

namespace lm {
namespace message {

class Communicator
{
public:
    static bool initializeDefaultSubclass();
    static void finalizeDefaultSubclass(bool abort=false);
    static Communicator* createObjectOfDefaultSubclass(bool isSupervisor=false);
    static std::string printableAddress(const lm::message::Endpoint& endpoint);

public:
    Communicator();
    virtual ~Communicator();

    // accessors
    virtual std::string getHostname() const;
    virtual Endpoint getSourceAddress() const;
    virtual Endpoint getHypervisorAddress() const=0;

    // send and receive messages
    virtual void sendMessage(Endpoint dest, lm::message::Message* msg, int sleepMilliseconds=0) const=0;
    virtual void receiveMessage(lm::message::Message* msg, int sleepMilliseconds=0) const=0;

protected:
    virtual bool initializeClass()=0;
    virtual void finalizeClass(bool abort)=0;
    virtual Endpoint constructObject(bool isSupervisor)=0;

protected:
    Endpoint sourceAddress;
};

}
}
#endif // LM_MESSAGE_COMMUNICATOR_H
