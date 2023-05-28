#   Copyright 2023 Zuse Institute Berlin
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 16:33:53 2017

@author: gongotar
"""

import datetime
import math
import utils
import matlab.engine

class Scheduler:
    """Scheduler containing timeline of checkpoints"""
    
    
    # the time interval of the scheduler
    SCHEDULE_INTERVAL = datetime.timedelta(seconds = 400)
    # Jobs natural interval precision based on job requests
    DO_PRECISION_INTERVAL = 5
    # Estimated scheduling duration for a job in milliseconds
    SCHEDULE_DURATION = 2000
    # Time interval in milliseconds for which a checkpoint request can be matched to a planed checkpoint
    REQUEST_TIME_TRESHOLD = 500
    
    def __init__(self, system):
        # empty timeline at the beginning
        self.timeline = []
        # Set the system running jobs
        self.system = system
        # Initialize the set of jobs last scheduled checkpoints
        self.jobs_last_cp = {}
        # Initialize the list of all jobs being scheduled
        self.jobs = []
        # Initialize the list of all job requests by time and result (approved/rejected)
        self.job_requests = {}
        # Connect to the matlab engine
        self.eng = matlab.engine.start_matlab()
        # mark the scheduler start time to detect the delay when receiving a new job
        self.scheduler_start_time = datetime.datetime.now()
        # define the next schedule milstone: until when do we want to schedule checkpoints
        self.next_schedule_milestone = self.scheduler_start_time + Scheduler.SCHEDULE_INTERVAL
    
    def get_scheduler_time(self):
        """Returns the passed time from starting scheduler in milliseconds"""
        passed_duration = datetime.datetime.now() - self.scheduler_start_time
        return passed_duration.days * 86400000 + \
            passed_duration.seconds * 1000 + \
            int(passed_duration.microseconds / 1000)
    
        
    def extend_scheduler_milestone(self):
        """Extends the milestone of the scheduler to schedule further checkpoints"""
        
        self.next_schedule_milestone = self.next_schedule_milestone + Scheduler.SCHEDULE_INTERVAL
        # TODO: add next checkpoints until the next milestone
    
    def compute_optimal_interval(self, job):
        """Computes the optimum checkpointing interval of the given job considering"""
        
        MTTI = job['MTTI']
        restart_time = job['restart_time']
        
        write_duration = job['cp_size'] * job['mpi_size'] * 1000 / self.system['write_bandwidth']
        
        optimum_interval = MTTI;
        if MTTI > write_duration / 2:
            var = write_duration / (2 * MTTI);
            optimum_interval = math.sqrt(2 * MTTI * write_duration) * (1 + math.sqrt(var)/3 + var/9) - write_duration;
        
        prev_step = optimum_interval - optimum_interval % job['natural_intv']
        next_step = prev_step + job['natural_intv']
        
        prev_expected_lifetime = utils.expected_lifetime(optimum_interval, prev_step, write_duration, MTTI, restart_time)
        next_expected_lifetime = utils.expected_lifetime(optimum_interval, next_step, write_duration, MTTI, restart_time)
        
        if prev_expected_lifetime > next_expected_lifetime:
            return int(next_step)
        else:
            return int(prev_step)


    def find_first_conflict(self, timeline):
        """Finds the first conflict of checkpoints in the given timeline"""
        conflict = []
        timeline.sort(key=lambda cp: cp['time'], reverse=False)
        max_end = timeline[0]['time'] + timeline[0]['duration']
        
        prev_end = 0
        next_start = 0
        
        for cp in timeline:    
            if cp['time'] < max_end:
                conflict.append(cp)
                if cp['time'] + cp['duration'] > max_end:
                    max_end = cp['time'] + cp['duration']
            elif len(conflict) == 1:
                prev_end = conflict[0]['time'] + conflict[0]['duration']
                max_end = cp['time'] + cp['duration']
                conflict[0] = cp
            else:
                next_start = cp['time']
                break
                
        window = [prev_end, next_start]
        
        return window, (conflict if len(conflict) > 1 else None)
    
    def optimize_conflict(self, eng, conflict, window):
        """Optimizes the given conflict by moving checkpoints limited to the given window"""

        movements = utils.solve_optimization_problem(eng, conflict, window)
        offsets = {}
        
        for i, mov in enumerate(movements):
            int_mov = int(mov[0])
            conflict[i]['time'] = conflict[i]['time'] + int_mov
            if conflict[i]['owner']['id'] not in offsets or offsets[conflict[i]['owner']['id']]['index'] < conflict[i]['index'] + 1:
                offsets[conflict[i]['owner']['id']] = {conflict[i]['index'] + 1: int_mov}
                
        return offsets
                
    def update_timeline(self, timeline, offsets):
        """Updates the given timeline according to the given offsets (movement steps for each job after an index)"""
        # TODO: implement
        for cp in timeline:
            if cp['owner']['id'] in offsets:
                index = list(offsets[cp['owner']['id']].keys())[0]
                if cp['index'] >= index:
                    cp['time'] = cp['time'] + offsets[cp['owner']['id']][index]
        
    
    def optimize_timeline(self, timeline):
        """Optimizes the timeline of checkpoints and constructs a schedule"""
        
        window, conflict = self.find_first_conflict(timeline)
        while conflict is not None:
            offsets = self.optimize_conflict(self.eng, conflict, window)
            self.update_timeline(timeline, offsets)
            window, conflict = self.find_first_conflict(timeline)
        
        timeline.sort(key=lambda cp: cp['time'], reverse=False)
        #print('++++++++++++++++++')
        #print(timeline)
    def add_new_job_to_scheduler(self, job):
        """Adds the new job to the checkpoint plan of scheduler (updates the scheduler)"""
        
        # Fetch the job offset corresponding to the scheduler start time
        job_start_offset = self.get_scheduler_time() + self.SCHEDULE_DURATION
        
        if utils.contains(self.jobs, lambda j: j['id'] == job['id']):
            self.finalize(self, job['id'], False)
        
        # Add the job to all jobs on the scheduler
        job['offset'] = job_start_offset
        self.jobs.append(job)
        self.job_requests[job['id']] = []
        
        # Compute optimal interval of the given job in milliseconds
        MTTI = utils.fetch_job_MTTI(job, self.system)
        restart_time = utils.fetch_restart_time(job, self.system)
        job['MTTI'] = MTTI
        job['restart_time'] = restart_time
        optimal_interval = self.compute_optimal_interval(job)
        job['opt_intv'] = optimal_interval

        # Compute the checkpointing interval in the schedule
        job_cp_start = job_start_offset + optimal_interval

        if self.next_schedule_milestone <= datetime.datetime.now():
            sched_intv = int(Scheduler.SCHEDULE_INTERVAL.days * 86400000 + \
                Scheduler.SCHEDULE_INTERVAL.seconds * 1000 +\
                Scheduler.SCHEDULE_INTERVAL.microseconds / 1000)
            scheduler_cp_end = job_start_offset + sched_intv
            self.next_schedule_milestone = datetime.datetime.now() + Scheduler.SCHEDULE_INTERVAL
        else:
            cp_duration = self.next_schedule_milestone - self.scheduler_start_time
            scheduler_cp_end = int(cp_duration.days * 86400000 + \
                cp_duration.seconds * 1000 +\
                cp_duration.microseconds / 1000)
        
        # Checkpoint write duration assuming availablity of all write bandwidth
        write_duration = int(job['cp_size'] * job['mpi_size'] * 1000 / self.system['write_bandwidth'])
        
        # Generate checkpoints and append them to the timeline of checkpoints
        index = 0
        cp_times = range(job_cp_start, scheduler_cp_end, optimal_interval + write_duration)
        for ii in cp_times:
            checkpoint = {'owner': job, 'index': index, 'time': ii, 'size': job['cp_size'] * job['mpi_size'], 'duration': write_duration}
            self.timeline.append(checkpoint)
            index+=1
            if index == len(cp_times):
                self.jobs_last_cp[job['id']] = checkpoint
        
        # Optimize the timeline considering the new job
        self.optimize_timeline(self.timeline)
        
        return self.get_scheduler_time() - job_start_offset

            
    def need_checkpoint(self, cp_request):
        """Checks the scheduler and returns 1 if checkpoint request is accepted, 0 if not"""
        # TODO: only consider requests that are rejected in the precision process, accepted requests take longer
        # Precise the job natural interval after a defined number of requests based on frequency of jobs requests
        
        result = 0
        request_time = self.get_scheduler_time()
        
        for cp in self.timeline:
            if cp['time'] - request_time < self.REQUEST_TIME_TRESHOLD and cp['owner']['id'] == cp_request['job']:
                result = 1
                break
        
        # Note the request time and answer for later precision process
        self.job_requests[cp_request['job']].append({self.get_scheduler_time(): result})

        return 1
    
    def finalize_job(self, job_id, reschedule):
        """Finalizes a job and removes all job entries from scheduler"""
        
        self.timeline = [cp for cp in self.timeline if cp['owner']['id'] != job_id]
        del self.jobs_last_cp[job_id]
        del self.job_requests[job_id]
        for i, o in enumerate(self.jobs):
            if o['id'] == job_id:
                del self.jobs[i]
                break
            
        if reschedule and len(self.timeline) > 1:
            # Optimizes the schedule again considering the removed job
            self.optimize_timeline(self.timeline)
    
    
    
