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
Created on Tue Aug 22 13:43:06 2017

@author: gongotar
"""
import datetime
import numpy as np
import utils
import math

def scope_test():
    def do_local():
        spam = "local spam"

    def do_nonlocal():
        nonlocal spam
        spam = "nonlocal spam"

    def do_global():
        global spam
        spam = "global spam"

    spam = "test spam"
    do_local()
    print("After local assignment:", spam)
    do_nonlocal()
    print("After nonlocal assignment:", spam)
    do_global()
    print("After global assignment:", spam)


scope_test()
print("In global scope:", spam)

class Dog:

    kind = []         # class variable shared by all instances

    def __init__(self, name):
        self.name = name    # instance variable unique to each instance
        
        
d = Dog('Fido')
e = Dog('Buddy')

d.kind.append('sss')
print(e.kind)

num = 34
#print(num.__class__.__doc__)

print(datetime.datetime.now())

system = {
    'write_bandwidth': 2, 
    'read_bandwidth': 3, 
    'hosts': {
        'ubuntu': 
            {'core1': 'core_1_info',
             'core2': 'core_2_info',
             'RAM': 'RAM_info'}
        }
    }
         
A = np.array([1, 4, 3, 7, 7])
B = np.array([4, 6, 4, 7, 3])

print(sum(np.multiply(A,B)))

print('necy')

def optim_intv(job, w):
    MTTI = job['MTTI']
    
    write_duration = w
    optimum_interval = MTTI;
    if MTTI > write_duration / 2:
        var = write_duration / (2 * MTTI);
        optimum_interval = math.sqrt(2 * MTTI * write_duration) * (1 + math.sqrt(var)/3 + var/9) - write_duration;
    job['opt_intv'] = optimum_interval
            
job1 = {'id': '12233', 'cp_size': 4, 'natural_intv': 1, 'mpis': 2, 'mpi_size': 2, 'MTTI': 60.0, 'restart_time': 2.0}
optim_intv(job1, 3)
job2 = {'id': '4442', 'cp_size': 4, 'natural_intv': 1, 'mpis': 2, 'mpi_size': 2, 'MTTI': 60.0, 'restart_time': 2.0}
optim_intv(job2, 2)
job3 = {'id': '354343', 'cp_size': 2, 'natural_intv': 1, 'mpis': 2, 'mpi_size': 2, 'MTTI': 60.0, 'restart_time': 2.0}
optim_intv(job3, 3)

conflict = [{'owner': job1, 'index': 0, 'time': 33.0, 'size': job1['cp_size'] * job1['mpi_size'], 'duration': 3.0},
             {'owner': job2, 'index': 0, 'time': 35.0, 'size': job2['cp_size'] * job2['mpi_size'], 'duration': 2.0},
             {'owner': job3, 'index': 0, 'time': 31.0, 'size': job3['cp_size'] * job3['mpi_size'], 'duration': 3.0},
             {'owner': job3, 'index': 1, 'time': 36.0, 'size': job3['cp_size'] * job3['mpi_size'], 'duration': 3.0}]

window = [30, 43]

#if eng is None:
 #   eng = matlab.engine.start_matlab()
utils.solve_optimization_problem(eng, conflict, window)
