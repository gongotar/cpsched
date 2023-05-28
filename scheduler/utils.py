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
Created on Thu Aug 24 12:21:58 2017

@author: gongotar
"""
import math
import numpy as np
#import matlab.engine
#from scipy.optimize import minimize

def contains(list, filter):
    """Checks if the given list contains a member passing the given filter"""
    for x in list:
        if filter(x):
            return True
    return False


def fetch_job_MTTI(job, system):
    """Computes the Meantime to failure of the given job in millis"""
    
    # TODO: Compute the MTTI of the given job
    MTTI = 4 * 3600 * 1000
    return MTTI

def fetch_restart_time(job, system):
    """Computes the restart time of the job on the given system in case of failure in millis"""
    
    # TODO: Compute the restart time of the given job on the given system
    R = 300*1000
    return R

def expected_lifetime(computation_time, cp_interval, write_duration, MTTI, restart_time):
        """Computes the expected lifetime of computing "computation_time" according to Daly's approach
        
            Considering the given checkpointing interval "cp_interval",
            write duration and also mean time to failure of the system and restart time
            of the system"""

        lifetime = MTTI * math.exp(restart_time / MTTI) \
            * (math.exp((cp_interval + write_duration) / MTTI) - 1) * computation_time / cp_interval
            
        return lifetime
 
'''
def d_expected_lifetime(cp_interval, write_duration, MTTI):
    """Computes the derivative of expected lifetime"""

    dlifetime = (cp_interval - MTTI) * math.exp((cp_interval + write_duration) / MTTI) + MTTI
        
    return dlifetime

def d2_expected_lifetime(cp_interval, write_duration, MTTI):
    """Computes the second derivative of expected lifetime"""

    d2lifetime = (cp_interval / MTTI) * math.exp((cp_interval + write_duration) / MTTI)
        
    return d2lifetime

    
def no_overlaps_constraint(x, conflict):
    K = 100000 
    for i in range(0, len(x)):
        ti = conflict[i]['time']
        oi = conflict[i]['owner']['opt_intv']
        wi = conflict[i]['duration']
        xi = x[i] - oi
        for j in range(i + 1, len(x)):
            tj = conflict[j]['time']
            wj = conflict[j]['duration']
            oj = conflict[j]['owner']['opt_intv']
            xj = x[j] - oj
            Aij = ti + xi - tj - xj < 0
            rij = ti + xi + wi - tj - xj - K * (1 - Aij)
            if rij > 0:
                print('1: ', rij)
                return -1
            Aji = tj + xj - ti - xi < 0
            rji = tj + xj + wj - ti - xi - K * (1 - Aji)
            if rji > 0:
                print('2: ', rji, ' j ', j,' i ', i, ' A ', Aji, ' xj ', xj, ' xi ', xi, ' tj ', tj, ' ti' , ti)
                return -1
    
    print(x, 'answer accepted')
    return 1

def fill_constraints_bounds(conflict, window):
    """Fills the values of G and h according to cvxopt library"""

    num = len(conflict)
    bnds = []
    cons = []#{'type': 'ineq', 'fun': no_overlaps_constraint, 'args': (conflict,)}
    
    for i in range(num):
        ti = conflict[i]['time']
        oi = conflict[i]['owner']['opt_intv']
        wi = conflict[i]['duration']
        #print('min ', ti - si + window[0], ' max ', ti - si + window[1] - wi)
        bnds.append((max(oi - ti + window[0], 1), oi - ti + window[1] - wi))
        

        for j in range(i + 1, num):
            tj = conflict[j]['time']
            wj = conflict[j]['duration']
            oj = conflict[j]['owner']['opt_intv']
            K = 10000.0 #sys.maxsize

            cons.append({'type': 'ineq', 'fun': lambda x: -(ti + x[i] - oi + wi - tj - x[j] + oj - K * (1.0 - (ti + x[i] - oi - tj - x[j] + oj < 0)))})
            cons.append({'type': 'ineq', 'fun': lambda x: -(tj + x[j] - oj + wj - ti - x[i] + oi - K * (1.0 - (tj + x[j] - oj - ti - x[i] + oi < 0)))})

    positive_pointer = 0
    negative_pointer = 1
    
    negative_value = -1
    positive_value = 1

    
    for ii in range(rows):
        
        if ii == rows / 2:
            negative_value, positive_value = positive_value, negative_value
            positive_pointer = 0
            negative_pointer = 1
            
        if ii > rows / 2:
            h[ii] = C[positive_pointer][negative_pointer]
        else:
            h[ii] = C[negative_pointer][positive_pointer]
        
        G[ii][positive_pointer] = positive_value
        G[ii][negative_pointer] = negative_value
        
        
        if negative_pointer < cols - 1:
            negative_pointer += 1
        else:
            positive_pointer += 1
            negative_pointer = positive_pointer + 1

    return cons, bnds

            
def lifetimes(computation_times, cp_intervals, write_durations, MTTIs, restart_times):
    """Computes the expected lifetimes considering the given array of parameters"""
    lifetimes = np.zeros(len(cp_intervals))
    for ii in range(len(cp_intervals)):
        lifetimes[ii] = expected_lifetime(computation_times[ii], cp_intervals[ii], write_durations[ii], MTTIs[ii], restart_times[ii])
    
    return lifetimes

def obj_func(x, conflict):
    optx = np.zeros(len(conflict))
    wx = np.zeros(len(conflict))
    mx = np.zeros(len(conflict))
    rx = np.zeros(len(conflict))
    
    for i in range(len(conflict)):
        optx[i] = conflict[i]['owner']['opt_intv']
        wx[i] = conflict[i]['duration']
        mx[i] = conflict[i]['owner']['MTTI']
        rx[i] = conflict[i]['owner']['restart_time']
        
    return sum(lifetimes(optx, x, wx, mx, rx))

def d_obj_func(x, conflict):
    d_sum = 0
    
    for i in range(len(conflict)):
        d_sum += d_expected_lifetime(x[i], conflict[i]['duration'], conflict[i]['owner']['MTTI'])
        
    return d_sum
 '''
 
def evaluate(solution, conflict):
    '''Evaluates the given solution against the given conflict'''
    # TODO:
    pass

def checkParallelWrite(solution, conflict):
    '''Checks the solution and conflict, returns the resulting parallel writes'''
    #TODO:    
    pass

def relax(X, conflict, priorities):
    '''Changes the given X based on natural interval of jobs and priorities'''
    #TODO:
    pass

def applyMovements(movements, conflict):
    '''Appleis the given movements set on the given conflict'''
    pass

def getConflicts(checkpoints):
    '''returns conflicts of the given checkpoints set'''
    pass

def chooseMostPrio(conflict):
    '''Chooses the checkpoint with most priority to stay at the optimal place'''
    '''The checkpoint whose movement causes most loss is the most prio checkpoint'''
    pass

def moveAllButPrio(conflict, prio):
    '''Moves all checkpoints in the conflict one step further but the given prio checkpoint which stays where it is'''
    pass

def fix(X, conflict, priorities, solution):
    '''Fixes the found solution in doubles considering natural intervals of jobs'''
    # TODO:
#    relaxedSolution = relax(X, conflict, priorities)
#    
#    parallelWriteSets = checkParallelWrite(relaxedSolution, conflict)
#    if(parallelWriteSets is None):
#        rate = evaluate(relaxedSolution, conflict)
#        if rate < solution['rate']:
#            solution['rate'] = rate
#            solution['x'] = relaxedSolution
#    else:
#        for parallelWrites in parallelWriteSets:
#            pass

    loss = np.zeros(len(conflict))
    movements = np.zeros(len(conflict))
    for i, c in enumerate(conflict):
        x = X[i]
        cp_interval1 = c['owner']['opt_intv'] + x - (x % c['owner']['natural_intv'])
        cp_interval2 = cp_interval1 + c['owner']['natural_intv']
        l1 = expected_lifetime(c['owner']['opt_intv'], cp_interval1, c['duration'], c['owner']['MTTI'], c['owner']['restart_time'])
        l2 = expected_lifetime(c['owner']['opt_intv'], cp_interval2, c['duration'], c['owner']['MTTI'], c['owner']['restart_time'])
        l_opt = expected_lifetime(c['owner']['opt_intv'], c['owner']['opt_intv'] + x, c['duration'], c['owner']['MTTI'], c['owner']['restart_time'])
        loss[i] = abs(l1 - l2)
        if l1 > l2:
            movements[i] = cp_interval2 - c['owner']['opt_intv']
            loss[i] = l1 - l_opt
        else:
            movements[i] = cp_interval1 - c['owner']['opt_intv']
            loss[i] = l2 - l_opt
         
    applyMovements(movements, conflict)
    conflicts = getConflicts(conflict)
    
    while conflicts is not None:
        for conflict in conflicts:
            mostPrio = chooseMostPrio(conflict)
            moveAllButPrio(conflict, mostPrio)
        conflict = getConflicts(conflict)
        

def solve_optimization_problem(eng, conflict, window):
    """Solves the given optimization problem"""
 
    optx = np.zeros(len(conflict))
    mx = np.zeros(len(conflict))
    rx = np.zeros(len(conflict))

    for i in range(len(conflict)):
        optx[i] = conflict[i]['owner']['opt_intv']
        mx[i] = conflict[i]['owner']['MTTI']
        rx[i] = conflict[i]['owner']['restart_time']


    eng.workspace['p_optx'] = optx.tolist()
    eng.workspace['p_mx'] = mx.tolist()
    eng.workspace['p_rx'] = rx.tolist()
    eng.workspace['p_conflict'] = conflict
    eng.workspace['p_window'] = window
    
    eng.optimization(nargout = 0)
    
    x = eng.eval('x')
    print('fval: ', eng.eval('fval'))
    
    print('x: ', x)
    
    solution = {'rate': 1000000, 'x': None}
    fix(x, conflict, np.zeros(len(conflict)), solution)
    
    #eng.constraints(matlab.double([0, 1, -1]))
    #print(eng.eval('msg'))
    return solution['x']
    '''
    
    optx = np.zeros(len(conflict))
    wx = np.zeros(len(conflict))
    mx = np.zeros(len(conflict))
    rx = np.zeros(len(conflict))
    
    for i in range(len(conflict)):
        optx[i] = conflict[i]['owner']['opt_intv']
        wx[i] = conflict[i]['duration']
        mx[i] = conflict[i]['owner']['MTTI']
        rx[i] = conflict[i]['owner']['restart_time']
        
    fun = lambda x: sum(lifetimes(optx, x, wx, mx, rx)) 

    optx = np.zeros(len(conflict))
    
    for i in range(len(conflict)):
        optx[i] = conflict[i]['owner']['opt_intv']
        
    cons, bnds = fill_constraints_bounds(conflict, window)
    res = minimize(obj_func, optx, args =(conflict,), method='SLSQP', bounds=bnds,
               constraints=cons)
    
    print(optx)
    print(res)

    m = 1
    n = 2
    optx = []
    wx = []
    mx = []
    rx = []
    G, h =  fill_G_h(conflict, window)
    
    def F(x = None, z = None):
        if x is None: return m, matrix(1.0, (n,1))
        if 0 in x: return None
        f = sum(expected_lifetime(optx, x, wx, mx, rx))
        Df = sum(d_expected_lifetime(x, wx, mx)).T
        if z is None: return f, Df
        H = spdiag(z[0] * sum(d2_expected_lifetime(x, wx, mx)))
        return f, Df, H
    
    sol = solvers.cp(F, G, h)
    print(sol['x'])


        m = 1
    n = 2
    optx = []
    wx = []
    mx = []
    rx = []
    G, h =  fill_G_h(no_overlap, window)
    
    def F(x):
        if x is None: return m, matrix(1.0, (n,1))
        else:
            return None
        if 0 in x: return None
        f = sum(expected_lifetime(optx, x, wx, mx, rx))
        Df = sum(d_expected_lifetime(x, wx, mx)).T
        if z is None: return f, Df
        
        H = spdiag(z[0] * sum(d2_expected_lifetime(x, wx, mx)))
        return f, Df, H
    
    sol = solvers.cp(F, G, h)
    print(sol['x'])
    '''
