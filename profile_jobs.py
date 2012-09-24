import itertools
import numpy as np
import re
import sys
import time
import os
from collections import namedtuple

Event = namedtuple('Event', 'timestamp type step job')
Proc  = namedtuple('Proc', 'step start stop job')

step_mapping = {
    'Run bowtie on genome'              : 'Run Bowtie on genome',
    'Parse genome Bowtie output'        : 'Run Bowtie on genome',
    'Run bowtie on transcriptome'       : 'Run Bowtie on transcriptome',
    'Parse transcriptome Bowtie output' : 'Run Bowtie on transcriptome',
    'Run blat on unmapped reads'        : 'Run BLAT',
    'Parse blat output'                 : 'Run BLAT',
    'Run mdust on unmapped reads'       : 'Run BLAT',
    }

def parse_log_file(filename, job_name):

    time_re = "(\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2})"
    pat = re.compile(time_re + ".*RUM\.Workflow.*(START|FINISH)\s+(.*)")
    
    with open(filename) as f:
        for line in f:
            m = pat.match(line)
            if (m is not None):
                (tm, type, step) = m.groups()
                t = time.strptime(tm, "%Y/%m/%d %H:%M:%S")
                e = Event(t, type, step, job_name)
                yield e

def build_timings(events):
    stack = []
    timings = []
    for e in events:
        if e.type == 'START':
            stack.append(e)
        elif e.type == 'FINISH':
            prev = stack.pop()
            yield Proc(e.step, prev.timestamp, e.timestamp, e.job)

def append_filenames(dest, dirname, filenames):
    pat = re.compile("rum_(postproc|\d+).*log$")
    for filename in filenames:
        if pat.match(filename):
            dest.append(dirname + "/" + filename)

def load_events_for_job(path, job_name):
    log_files = []
    os.path.walk(path, append_filenames, log_files)
    event_iters = [parse_log_file(f, job_name) for f in log_files]
    return itertools.chain(*event_iters)

def group_times_by_step(procs):

    times_for_step = {}

    for p in procs:
        if p.step not in times_for_step:
            times_for_step[p.step] = []
        start = time.mktime(p.start)
        stop  = time.mktime(p.stop)
            
        times_for_step[p.step].append(stop - start)

    return times_for_step

def summarize_steps(steps):
    for (step, times) in steps:
        chunks = len(times)
        total  = sum(times)
        med    = np.median(times)

def rename_steps(step_to_times):
    result = {}
    for step in step_to_times:
        times = np.array(step_to_times[step])

        new_step = step_mapping[step] if step in step_mapping else step
        print "New step is " + new_step
        if new_step in result:
            print "Adding " + str(result[new_step]) + " to " + str(times)
            result[new_step] += times
        else:
            result[new_step] = times
    return result

def main():

    dtype=[('step', 'S100')]

    jobs = []

    for arg in sys.argv[1:]:
        
        (job_name, dir) = arg.split("=")
        jobs.append((job_name, dir))

        dtype.extend([
                ('%s_chunks'   % job_name, int),
                ('%s_total_time'  % job_name, int),
                ('%s_median_time' % job_name, int)])

    times = {}

    step_idx = {}

    times_for_job = {}

    seen_steps = set()
    steps      = []

    for (job_name, dir) in jobs:
        procs = list(build_timings(load_events_for_job(dir, job_name)))
        grouped = group_times_by_step(procs)
        renamed = rename_steps(grouped)
        times_for_job[job_name] = renamed
        for p in procs:
            step = step_mapping[p.step] if p.step in step_mapping else p.step
            if step not in seen_steps:
                steps.append(step)
                seen_steps.add(step)
    
    table = []

    for step in steps:
        row = (step,)
        for (job_name, dir) in jobs:
            if step in times_for_job[job_name]:
                times = times_for_job[job_name][step]
                row += (len(times), 
                        sum(times),
                        np.median(times))
            else:
                row += (0, 0, 0)

        table.append(row)

    table=np.array(table, dtype=dtype)
    print table

main()
