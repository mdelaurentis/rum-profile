import itertools
import numpy as np
import re
import sys
import time
import os

from collections import namedtuple
from numpy.lib.recfunctions import append_fields

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
        chunks = len(times)
        total  = sum(times)
        med    = np.median(times)
        table.append((s, chunks, total, med))
    
    return np.array(table, dtype=[
            ('step', 'S100'),
            ('chunks', int),
            ('total',  float),
            ('median', float)])


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
        if new_step in result:
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
                ('%s_chunks'      % job_name, int),
                ('%s_total_time'  % job_name, int),
                ('%s_median_time' % job_name, int)])

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
    print steps
    result = np.array(steps, dtype=[('step', 'S100')])

    for (name, table) in jobs:
        print "Name is " + name
        
        result = append_fields(
            result,
            ['%s_chunks' % name,
             '%s_total'  % name,
             '%s_median' % name],
            [table['chunks'],
             table['total'],
             table['median']
             ])

    for (name, table) in jobs:

        totals  = result['%s_total' % name]
        medians = result['%s_median' % name]

        total_intensity  = (totals  - min(totals))  / (max(totals)  - min(totals))
        median_intensity = (medians - min(medians)) / (max(medians) - min(medians))
                
        total_intensity = (255 - (total_intensity * 255))
        median_intensity = (255 - (median_intensity * 255))

        print "Total intensity is " + str(total_intensity)

        result = append_fields(
            result,
            ['%s_total_pct'        % name,
             '%s_median_pct'       % name,
             '%s_total_intensity'  % name,
             '%s_median_intensity' % name,
             ],
            [100 * result['%s_total'  % name] / sum(result['%s_total'  % name]),
             100 * result['%s_median' % name] / sum(result['%s_median' % name]),
             total_intensity,
             median_intensity
             ])


    return result
#    return np.array(list(result) + [footer], dtype=result.dtype)

def calc_gain(table, job_names):
    baseline_total_secs = table['%s_total' % job_names[0]]
    baseline_median_secs = table['%s_median' % job_names[0]]

    for name in job_names[1:]:
        stacked_gain = baseline_total_secs  - table['%s_total' % name]
        median_gain  = baseline_median_secs - table['%s_median' % name]

        stacked_pct_gain = stacked_gain / sum(baseline_total_secs)
        median_pct_gain  = median_gain  / sum(baseline_median_secs)

        table = append_fields(table, 
                      ['%s_stacked_gain' % name,
                       '%s_stacked_pct_gain' % name,
                       '%s_median_gain' % name,
                       '%s_median_pct_gain' % name,
                       ],
                      [stacked_gain, 100. * stacked_pct_gain,
                       median_gain,  100. * median_pct_gain]
                      )
    return table

def print_table(filename, table, job_names):

    table = calc_gain(table, job_names)
    print table.dtype
    with open(filename, 'w') as out:
        out.write("""

<html>
  <head>
    <link rel="stylesheet" type="text/css" href="profile.css"></link>
  </head>
  <body>
    <table>
      <thead>

        <tr>
          <th>Step</th>
""")

        is_first = True

        for j in job_names:
            out.write("<th>%s chunks</th>" % j)
            out.write("<th>%s total</th>"  % j)
            out.write("<th>(%)</th>")
            out.write("<th>%s wallclock</th>" % j)
            out.write("<th>(%)</th>")
            if not is_first:
                out.write('<th>Stacked gain</th>')
                out.write('<th>(%)</th>')
                out.write('<th>Median gain</th>')
                out.write('<th>(%)</th>')
            is_first = False
        
        out.write("""
        </tr>
      </thead>
      <tbody>""")

        for row in table:
            out.write("""
        <tr>
          <td>%s</td>""" % row['step'])
            
            is_first = True

            for j in job_names:

                total_intensity = row['%s_total_intensity' % j]
                median_intensity = row['%s_median_intensity' % j]

                total_color  = '#ff%02x%02x' % (total_intensity, total_intensity)
                median_color = '#ff%02x%02x' % (median_intensity, median_intensity)
                
                out.write("<td>%d</td>" % row['%s_chunks' % j])
                out.write("<td>%d</td>" % row['%s_total'  % j])
                out.write("<td bgcolor='%s'>%.2f%%</td>" % (total_color, row['%s_total_pct'  % j]))
                out.write("<td>%d</td>" % row['%s_median' % j])
                out.write("<td bgcolor='%s'>%.2f%%</td>" % (median_color, row['%s_median_pct' % j]))

                if not is_first:
                    out.write('<td>%d</td>' % row['%s_stacked_gain' % j])
                    out.write('<td>%.2f</td>' % row['%s_stacked_pct_gain' % j])
                    out.write('<td>%d</td>' % row['%s_median_gain' % j])
                    out.write('<td>%.2f</td>' % row['%s_median_pct_gain' % j])
                is_first = False
            out.write("""
        </tr>""")

        out.write("<tr><td>Total</td>")

        is_first = True
        for j in job_names:

            total_intensity = row['%s_total_intensity' % j]
            median_intensity = row['%s_median_intensity' % j]

            total_color  = '#ff%02x%02x' % (total_intensity, total_intensity)
            median_color = '#ff%02x%02x' % (median_intensity, median_intensity)
            
            out.write("<td></td>")
            out.write("<td>%d</td>" % sum(table['%s_total' % j]))
            out.write("<td>%.2f%%</td>" % sum(table['%s_total_pct' % j]))
            out.write("<td>%d</td>" % sum(table['%s_median' % j]))
            out.write("<td>%.2f%%</td>" % sum(table['%s_median_pct' % j]))
            
            if not is_first:
                out.write('<td>%d</td>' % sum(table['%s_stacked_gain' % j]))
                out.write('<td>%.2f</td>' % sum(table['%s_stacked_pct_gain' % j]))
                out.write('<td>%d</td>' % sum(table['%s_median_gain' % j]))
                out.write('<td>%.2f</td>' % sum(table['%s_median_pct_gain' % j]))
            is_first = False
        
        out.write("""
        </tr>
      </tbody>
    </table>
  </body>
</html>
""")


main()
