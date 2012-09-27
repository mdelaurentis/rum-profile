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

    metrics = ['cpu', 'wc']
    for metric in metrics:

        filename = 'rum_profile/%s.html' % metric
        print_table(filename, table, [x[0] for x in tables], metric)

    print_help_page()

def print_help_page():
    with open('rum_profile/help.html', 'w') as f:
        contents = E.DIV(
            E.P("""
This shows the amount of time spent on each step for one or more RUM
jobs. You should be able to use the output to gain some insight into
the performance of a single job, or to compare the running times for
two or more jobs in order to understand the effects of changes made to
the RUM code, to the input data, or to the system on which RUM was
run."""),

            E.P("""
The times for a single job are are compiled by parsing the log files
after the job is finished. For each step, we consider the duration to
be from the time the script was actually invoked to the time the
script exited. This means that any latency caused by the cluster's
scheduling system will not be reflected in the times."""),

            E.P("""
"CPU time" for a step measures the total time for the step across all
chunks. If you have grouped multiple jobs together, the CPU time
reported here is the median value for all the jobs. This is intended
to show the total amount of computing time the job would take for a
typical invocation. The total CPU time for the job is the sum of all
these median values.
"""),

            E.P("""
"Wallclock time" for a step is the time of the chunk that took the
longest time to complete the step. We use the maximum value in order
model the worst-case scenario, where one of the chunks takes much
longer than the other chunks, and becomes the limiting factor for the
running time of the whole job. If you have grouped jobs together,
wallclock time is the median value for all jobs. This is intended to
show the maximum duration of a typical invocation of the job."""),

            E.P("""
The times for all steps are highlighted in varying shades of yelloq in
order to indicate each step's impact on the total running time,
relative to the other steps. The step with the longest time is pure
yellow, the step with the shortest time is pure white, and the colors
for the other steps are scaled linearly according to the running
time. This is intended to allow you to quickly identify hot spots, or
steps that took a very long time compared to other steps, and which
might benefit from optimization or parameter tuning.
"""),

            E.P("""
If you have run two or more job groups, we use the first group as a
baseline and compare the running time of all the other jobs to the
baseline. The "improvement" for a step shows the number of hours saved
relative to the baseline. The percent improvement is the improvement
divided by the total running time for the baseline job. This is
intended to show the impact that improving the running time of one
step has on the running time of the entire job. For example suppose
the baseline job took 100 hours, 30 of which were spent on the "Run
BLAT" step, and that the "Run BLAT" step in the comparison job took
only 20 hours. The improvement is 30 - 20 = 10 hours, and the percent
improvement is (30 - 20) / 100 = 10%.  If a step in the new job is
slower, the improvement and percent improvement will be negative."""),

            E.P("""
The improvement for each step is colored green or red according to the
degree to which that step improved or degraded the performance
compared to the baseline."""),


            E.P(E.STRONG("Note:"),
                """
These numbers do not include the "preprocessing" phase at all. Prior
to RUM 2.0.3, it is difficult to determine from the log files exactly
when pre-processing begins and ends. RUM 2.0.3 and greater will
clearly include these times in the log file, so future versions of the
profiling program will be able to incorporate times for
preprocessing.""")
            )
        


        help_page = template('help',
                             contents)
        f.write(lxml.html.tostring(help_page))


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

def td_hours(seconds, bgcolor=None):
    if bgcolor is None:
        return TD('%.2f' % (seconds / seconds_per_hour),
                  CLASS='numeric')
    else:
        return TD('%.2f' % (seconds / seconds_per_hour),
                  CLASS='numeric',
                  bgcolor=bgcolor)

def td_percent(pct, bgcolor=None):
    if bgcolor is None:
        return TD('%.2f%%' % pct,
                  CLASS='numeric')
    else:
        return TD('%.2f%%' % pct, bgcolor=bgcolor,
                  CLASS='numeric')

metric_names = {
    'wc' : 'Wallclock',
    'cpu' : 'CPU'
}

def gain_pct_to_bgcolor(pct):
    intensity = 255 * (pct / 100.)
    if (intensity < 0):
        intensity = 255 + intensity
        return '#ff%02x%02x' % (intensity, intensity)
    else:
        intensity = 255 - intensity
        return '#%02xff%02x' % (intensity, intensity)
        
def template(name, contents):

    cpu_class  = 'active' if name == 'cpu'  else ''
    wc_class   = 'active' if name == 'wc'   else ''
    help_class = 'active' if name == 'help' else ''

    return E.HTML(
        E.HEAD(
            E.LINK(rel='stylesheet', type='text/css', href='css/bootstrap.css'),
            E.LINK(rel='stylesheet', type='text/css', href='profile.css'),
            E.SCRIPT(src='js/bootstrap.min.js'),
            E.TITLE('RUM Job Profile')),

        E.BODY(

            E.DIV(
                E.DIV(
                    E.DIV(
                        E.A(E.SPAN(CLASS='icon-bar'),
                            E.SPAN(CLASS='icon-bar'),
                            E.SPAN(CLASS='icon-bar'),
                            CLASS='btn btn-navbar'),
                        E.A('RUM Profile', CLASS='brand', href='#'),
                        E.DIV(
                            E.UL(
                                E.LI(
                                    E.A('CPU time', href='cpu.html'),
                                    CLASS=cpu_class),
                                E.LI(
                                    E.A('Wallclock time', href='wc.html'),
                                    CLASS=wc_class),
                                E.LI(
                                    E.A('Help', href='help.html'),
                                    CLASS=help_class),
                                
                                CLASS='nav'),
                            CLASS='nav-collapse collapse'),
                        CLASS='container'),
                    CLASS='navbar-inner'),
                CLASS='navbar navbar-inverse navbar-fixed-top'),
            E.BR(),
            E.BR(),
            E.BR(),
            E.DIV(contents,
                  CLASS='container')))
    
def print_table(filename, table, job_names, metric):

    table = calc_gain(table, job_names)

    baseline = job_names[0]

    with open(filename, 'w') as out:

        top_headers = TR(TH(''))
        headers    =  TR(TH('Steps')) 

        # Build the header row
        for j in job_names:

            colspan = 1

            headers.append(TH('chunks'))
            colspan += 2
            headers.append(TH('hours', colspan='2'))

            # All jobs except the baseline get "hours gained" and
            # "percent hours gained" columns
            if j != baseline:
                colspan += 2
                headers.append(TH('improvement', colspan='2'))

            top_headers.append(TH(j, colspan=str(colspan)))

        data_rows = []

        for row in table:
            
            tr = TR(TD(str(row['step'])))

            for j in job_names:

                tr.append(TD(str(row['%s_chunks' % j]),
                             CLASS='numeric'))

                intensity  = row['%s_%s_intensity' % (j, metric)]
                bgcolor  = '#ffff%02x' % intensity
                tr.append(td_hours(row['%s_%s'  % (j, metric)], bgcolor))
                tr.append(td_percent(row['%s_%s_pct' % (j, metric)], bgcolor))

                if j != baseline:
                    hours = row['%s_%s_gain' % (j, metric)]
                    pct   = row['%s_%s_pct_gain' % (j, metric)]

                    intensity = 255 * (pct / 100.)

                    bgcolor = gain_pct_to_bgcolor(pct)
                    tr.append(td_hours(hours, bgcolor=bgcolor))
                    tr.append(td_percent(pct, bgcolor=bgcolor))

            data_rows.append(tr)

        summary = [TD('Totals')]
        for j in job_names:
            summary.append(TD(''))
            summary.extend([
                    td_hours(sum(table['%s_%s' % (j, metric)])),
                    td_percent(sum(table['%s_%s_pct' % (j, metric)]))])
            
            if j != baseline:
                hours   = sum(table['%s_%s_gain' % (j, metric)])
                pct     = sum(table['%s_%s_pct_gain' % (j, metric)])
                bgcolor = gain_pct_to_bgcolor(pct)
                summary.extend([
                        td_hours(hours, bgcolor),
                        td_percent(pct, bgcolor)])

        rows = [top_headers, headers]
        rows.extend(data_rows)
        rows.append(TR(*summary))

        html = template(metric, E.DIV(
                E.TABLE(*rows)))
                
        out.write(lxml.html.tostring(html))

main()
