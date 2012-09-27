import itertools
import numpy as np
import re
import sys
import time
import os
from lxml.html import builder as E
from lxml.html.builder import TR, TD, TH
import lxml


from collections import namedtuple
from numpy.lib.recfunctions import append_fields

seconds_per_hour = 60.0 * 60.0

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
    'Sort junctions (all, bed) by location' : 'Sort junctions',
    'Sort junctions (all, rum) by location' : 'Sort junctions',
    'Sort junctions (high-quality, bed) by location' : 'Sort junctions',
    'Finish mapping stats'      : 'Other post-processing',
    'Merge SAM headers'          : 'Other post-processing',
    'Sort junctions'             : 'Other post-processing',
    'Get inferred internal exons'             : 'Other post-processing',
    'Merge quants'             : 'Other post-processing',
    'Sort junctions'             : 'Other post-processing',
    'Remove duplicates from NU' : 'Other processing',
    'Make unmapped reads file for blat' : 'Other processing',
    
    }

def parse_log_file(filename, job_name):

    """Parse the specified log file and return an iterator of Events
    from the file."""

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
            if prev.step != e.step:
                raise Exception(
                    """I have a FINISH event for the START event of a
                    different step""")
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

def new_step_name(step):
    if step in step_mapping:
        return step_mapping[step]
    return step

def procs_to_array(procs):

    steps      = []
    seen_steps = set()
    for p in procs:
        step = new_step_name(p.step)
        if step not in seen_steps:
            seen_steps.add(step)
            steps.append(step)

    times_for_step = group_times_by_step(procs)
    times_for_step = rename_steps(times_for_step)

    table = []
    for s in steps:
        times = np.array(times_for_step[s])
        table.append((s, len(times), sum(times), max(times)))
    
    return np.array(table, dtype=[
            ('step', 'S100'),
            ('chunks', int),
            ('total',  float),
            ('median', float)])

def rename_steps(step_to_times):
    result = {}
    for step in step_to_times:
        times = np.array(step_to_times[step])

        new_step = step_mapping[step] if step in step_mapping else step
        if new_step in result:
            result[new_step] += times
        else:
            result[new_step] = times
    return result

def main():

    jobs = []

    for arg in sys.argv[1:]:
        
        (job_name, dir) = arg.split("=")
        jobs.append((job_name, dir))

    stats = {}

    for (job_name, dir) in jobs:
        procs = list(build_timings(load_events_for_job(dir, job_name)))
        job_stats = procs_to_array(procs)
        if job_name not in stats:
            stats[job_name] = []
        stats[job_name].append(job_stats)

    tables = []

    seen_job_name = set()
    for (job_name, path) in jobs:
        if job_name not in seen_job_name:
            seen_job_name.add(job_name)
            copies = stats[job_name]
            merged = merge_copies(stats[job_name])
            tables.append((job_name, merged))

    table = make_final_table(tables)
    print_table('rum_profile/index.html', table, [x[0] for x in tables])

def merge_copies(copies):

    result = []

    steps = copies[0]['step']
    print "--- Merging %d copies ---" % len(copies)
    for step in range(len(copies[0])):
        steps   = set([c[step]['step']   for c in copies])
        chunks  = set([c[step]['chunks'] for c in copies])
        totals  =     [c[step]['total']  for c in copies]
        medians =     [c[step]['median'] for c in copies]

        if (len(steps) != 1):
            raise Exception("Different steps for copies of same job")
        if (len(chunks) != 1):
            raise Exception("Different number of chunks for copies of same job")

        row = (list(steps)[0],
               list(chunks)[0],
               np.median(totals),
               np.median(medians))
        result.append(row)
    return np.array(result, dtype=copies[0].dtype)

def make_final_table(jobs):

    stats = [j[1] for j in jobs]

    columns = []
    did_step = False
    steps = [(step,) for step in jobs[0][1]['step']]
    result = np.array(steps, dtype=[('step', 'S100')])

    for (name, table) in jobs:
        
        result = append_fields(
            result,
            ['%s_chunks' % name,
             '%s_cpu'  % name,
             '%s_wc' % name],
            [table['chunks'],
             table['total'],
             table['median']
             ])

    for (name, table) in jobs:

        cpu_times  = result['%s_cpu' % name]
        wc_times = result['%s_wc' % name]

        cpu_intensity  = (cpu_times  - min(cpu_times))  / (max(cpu_times)  - min(cpu_times))
        wc_intensity = (wc_times - min(wc_times)) / (max(wc_times) - min(wc_times))
                
        cpu_intensity = (255 - (cpu_intensity * 255))
        wc_intensity = (255 - (wc_intensity * 255))

        result = append_fields(
            result,
            ['%s_cpu_pct'        % name,
             '%s_wc_pct'       % name,
             '%s_cpu_intensity'  % name,
             '%s_wc_intensity' % name,
             ],
            [100 * cpu_times / sum(cpu_times),
             100 * wc_times / sum(wc_times),
             cpu_intensity,
             wc_intensity
             ])

    return result
#    return np.array(list(result) + [footer], dtype=result.dtype)

def calc_gain(table, job_names):
    baseline_cpu_secs = table['%s_cpu' % job_names[0]]
    baseline_wc_secs = table['%s_wc' % job_names[0]]

    for name in job_names[1:]:
        cpu_gain = baseline_cpu_secs  - table['%s_cpu' % name]
        median_gain  = baseline_wc_secs - table['%s_wc' % name]

        cpu_pct_gain = cpu_gain / sum(baseline_cpu_secs)
        wc_pct_gain  = median_gain  / sum(baseline_wc_secs)

        table = append_fields(table, 
                      ['%s_cpu_gain' % name,
                       '%s_cpu_pct_gain' % name,
                       '%s_wc_gain' % name,
                       '%s_wc_pct_gain' % name,
                       ],
                      [cpu_gain, 100. * cpu_pct_gain,
                       median_gain,  100. * wc_pct_gain]
                      )
    return table

def td_hours(seconds):
    return TD('%.2f' % (seconds / seconds_per_hour))

def td_percent(pct, bgcolor=None):
    if bgcolor is None:
        return TD('%.2f%%' % pct)
    else:
        return TD('%.2f%%' % pct, bgcolor=bgcolor)

def print_table(filename, table, job_names):

    table = calc_gain(table, job_names)

    metrics = ['cpu', 'wc']

    baseline = job_names[0]

    with open(filename, 'w') as out:

        headers = []

        # Build the header row
        for j in job_names:
            headers.append('%s chunks' % j)
            for m in metrics:
                headers.extend(['%s %s hours' % (j, m), '(%)'])

            # All jobs except the baseline get "hours gained" and
            # "percent hours gained" columns
            if j != baseline:
                for m in metrics:
                    headers.extend(['%s hours gained' % (m), '(%)'])

        header_row = TR(TH('Step'))        
        for h in headers:
            header_row.append(TH(h))

        data_rows = []

        for row in table:
            
            tr = TR(TD(str(row['step'])))

            for j in job_names:

                tr.append(TD(str(row['%s_chunks' % j])))

                for m in metrics:
                    intensity  = row['%s_%s_intensity' % (j, m)]
                    bgcolor  = '#ff%02x%02x' % (intensity, intensity)
                    tr.append(td_hours(row['%s_%s'  % (j, m)]))
                    tr.append(td_percent(row['%s_%s_pct' % (j, m)], bgcolor))

                if j != baseline:
                    for m in metrics:
                        tr.append(td_hours(row['%s_%s_gain' % (j, m)]))
                        tr.append(td_percent(row['%s_%s_pct_gain' % (j, m)]))

            data_rows.append(tr)

        summary = [TD('Totals')]
        for j in job_names:
            summary.append(TD(''))
            for m in metrics:
                summary.extend([
                        td_hours(sum(table['%s_%s' % (j, m)])),
                        td_percent(sum(table['%s_%s_pct' % (j, m)]))])
            
            if j != baseline:
                for m in metrics:
                    summary.extend([
                            td_hours(sum(table['%s_%s_gain' % (j, m)])),
                            td_percent(sum(table['%s_%s_pct_gain' % (j, m)]))])

        rows = [header_row]
        rows.extend(data_rows)
        rows.append(TR(*summary))

        html = E.HTML(
            E.HEAD(
                E.LINK(rel='stylesheet', type='text/css', href='profile.css'),
                E.TITLE('RUM Job Profile')),
            E.BODY(
                E.H1("RUM Job Profile"),
                E.TABLE(*rows)))
                
        out.write(lxml.html.tostring(html))

main()
