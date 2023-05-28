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
Created on Wed Aug 16 14:15:26 2017

@author: gongotar
"""

#!/usr/bin/env python

import socket
import scheduler as sch
import threading as th
import traceback
from time import sleep

    
def init_system():
    """Initializes the attributes of the system running jobs and returns it"""
    
    system = {
        'write_bandwidth': 2 * 1024 * 1024 * 1024,  # 2 GB/s
        'read_bandwidth': 3 * 1024 * 1024 * 1024,   # 3 GB/s
        'hosts': {
            'ubuntu': 
                {'core1': 'core_1_info',
                 'core2': 'core_2_info',
                 'RAM': 'RAM_info'}
            }
        }
                
    return system

# Listening server ip and port
TCP_IP = '130.73.144.90' #'192.168.56.101'
TCP_PORT = 3030
BUFFER_SIZE = 512  # Normally 1024, but we want fast response

system = init_system()
scheduler = sch.Scheduler(system)
lock = th.Lock()

def need_checkpoint(first_mpi, checkpoint_request_count):
    """Checks the scheduler and returns 1 if checkpoint request is accepted, 0 if not"""
    
    cp_request = {'job': first_mpi, 'cp_req_count': checkpoint_request_count}
    result = scheduler.need_checkpoint(cp_request)
    return result

def add_new_job_to_scheduler(first_mpi, mpis, checkpoint_size, natural_interval, mpi_size):
    """Adds the new job to the checkpoint plan of scheduler (updates the scheduler)"""


    job = {'id': first_mpi, 'cp_size': checkpoint_size, 'natural_intv': natural_interval, 'mpis': mpis, 'mpi_size': mpi_size}

    elapsed_time = scheduler.add_new_job_to_scheduler(job)
    return elapsed_time

def finalize_job(job_id):
    scheduler.finalize_job(job_id, True)

def fetch_info_from_data(data):
    """Fetches the information from data by parsing the received data"""
    
    # Splits data by ':'
    splitted = data.decode().split('#')
    
    # Type of the message would be the first item (checkpoint request, init request, ...)
    msg_type = splitted[0]
    
    # The second item would be the process id of the first mpi rank of the job
    first_mpi = splitted[1]
    
    if msg_type == 'cpr':
        # The next item would be the checkpoint request count
        checkpoint_request_count = int(splitted[4])
        return msg_type, first_mpi, checkpoint_request_count
    elif msg_type == 'init':
        mpi_host_cpus = [x.strip() for x in splitted[2].split(',')]
        mpis = []
        for mpi_host_cpu in mpi_host_cpus:
            splitted_mpi_host_cpu = mpi_host_cpu.split(':')
            mpis.append({
                'mpi': splitted_mpi_host_cpu[0],
                'host': splitted_mpi_host_cpu[1], 
                'core': splitted_mpi_host_cpu[2]
                 })

        # The next item would be the expected checkpoint size
        checkpoint_size = int(splitted[3])
        # The next item would be the job natural interval
        natural_interval = int(splitted[4])
        # The next item would be the job mpi size
        mpi_size = int(splitted[5])
        return msg_type, first_mpi, mpis, checkpoint_size, natural_interval, mpi_size
    elif msg_type == 'fin':
        return msg_type, first_mpi
        

def process_received_connection(conn, addr):
    """Processes the received connection as a seprated thread"""
    
    print('Connection address:', addr)
    first_mpi = None
    try:
        # Receives data from client
        data = conn.recv(BUFFER_SIZE)
        while data:
            # Parses the data and fetches the needed information from it
            msg = fetch_info_from_data(data)
            msg_type = msg[0]
            first_mpi = msg[1]
            # get lock
            lock.acquire()
            # If message type is checkpoint request then check if its accepted
            if msg_type == 'cpr':
                # Passes the fetched information to the scheduler and returns the result
                checkpoint_request_count = msg[2]
                needCP = need_checkpoint(first_mpi, checkpoint_request_count)
                msg = str(needCP)
                #  Release lock
                lock.release()
            # Else if the message type is intro then add the new job to the scheduler
            elif msg_type == 'init':
                # Passes the fetched information to the scheduler for initialization
                mpis = msg[2]
                checkpoint_size = msg[3]
                natural_interval = msg[4]
                mpi_size = msg[5]
                elapsed_time = add_new_job_to_scheduler(first_mpi, mpis, checkpoint_size, natural_interval, mpi_size)
                to_wait = scheduler.SCHEDULE_DURATION - elapsed_time
                # Release lock and wait to stay synced with scheduler
                lock.release()
                sleep(float(to_wait)/1000.0)
                msg = 'done'
                
            elif msg_type == 'fin':
                finalize_job(first_mpi)
                # Release lock
                lock.release()
                msg = 'done'
            
            # Sends the result (permission to write checkpoint) back to the client
            conn.send(msg.encode('ascii'))  # echo
            data = conn.recv(BUFFER_SIZE)
        
    except Exception as e:
        print('released in  exception: ', first_mpi)
        if lock.locked():
            lock.release()
        traceback.print_exc()
        
    # Closes the connection when the process which is writing checkpoints ends
    conn.close()
    print('Connection ', addr, ' closed')

    
def main():
    """Main function"""

    # Creates a socket on the given ip and port and listens to it
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    s.bind((TCP_IP, TCP_PORT))
    print('listening for connections ...')
    s.listen(10)
    conn = None
    try:
        # Accepts connections from clients and processes them in a separated thread
        while 1:
            conn, addr = s.accept()
            clinet_thread = th.Thread(target=process_received_connection, args=(conn, addr))
            clinet_thread.start()
    except Exception as e:
        print('Exception occured')
        s.close()
        traceback.print_exc()
    except KeyboardInterrupt:
        if conn:
            conn.close()
        s.close()
        print('socket closed')
main()
