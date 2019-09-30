/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef ALPS_MPSPSCAN_SCHEDULER_HPP
#define ALPS_MPSPSCAN_SCHEDULER_HPP

#include "dmrg/utils/time_stopper.h"
#include "libpscan/options.hpp"
#include <boost/filesystem.hpp>

enum TaskStatusFlag {
    TaskNotStarted, TaskRunning, TaskHalted, TaskFinished
};

struct TaskDescriptor {
    TaskStatusFlag status;
    boost::filesystem::path in;
    boost::filesystem::path out;
};

class Scheduler {
public:
    Scheduler(const Options&);
    void run();
    
private:
    void parse_job_file(const boost::filesystem::path&);
    void checkpoint_status() const;
    
    boost::filesystem::path outfilepath;
    boost::filesystem::path infilepath;
    
    time_stopper stop_callback;
    bool write_xml;
    std::vector<TaskDescriptor> tasks;
};

#endif