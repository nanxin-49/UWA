"""Manual D1b resource prototype transport and Linux job supervisor."""
import csv
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import signal
import statistics
import subprocess
import sys
import tarfile
import time

from d1_package import sha, verify


GIB = 1024**3
REMOTE = Path('/home/misc/cwx')
MATLAB = Path('/home/misc/MATLAB/R2025b/bin/matlab')
SCOPE = 'D1B_RESOURCE_PROTOTYPE_ONLY'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def write(path, body):
    temp = path.with_suffix('.pending')
    temp.write_text(json.dumps(body, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    temp.replace(path)


def preflight(root):
    manifest = read(root/'manifest.json')
    cfg = read(root/'config.json')
    assert manifest['scope'] == cfg['scope'] == SCOPE and cfg['full_d1b_enabled'] is False
    verify(root, manifest['files'])
    assert cfg['prototype_half_widths_m']==[6,12,24] and cfg['max_nodes']==1536
    assert cfg['pool_workers']==[0,2,4] and cfg['repeats']==2
    assert cfg['max_wall_s']==1800 and cfg['min_available_GiB']==18
    assert cfg['rss_cap_GiB']==12 and cfg['reserve_fraction']==.3
    return cfg, sha(root/'manifest.json')


def owned_root(root):
    assert root.parent == REMOTE and root.resolve() == root and root.is_dir()
    assert re.fullmatch(r'd1b_v2_[A-Za-z0-9_-]{1,64}',root.name)


def available():
    with Path('/proc/meminfo').open() as handle:
        return int(next(row.split()[1] for row in handle if row.startswith('MemAvailable:')))*1024


def pressure(name):
    lines = Path('/proc/pressure',name).read_text().splitlines()
    full = next(row for row in lines if row.startswith('full '))
    return float(next(token[6:] for token in full.split() if token.startswith('avg10=')))


def proc(pid):
    folder = Path('/proc',str(pid))
    line = (folder/'stat').read_text()
    data = line[line.rindex(')')+2:].split()
    return dict(pid=pid,state=data[0],ppid=int(data[1]),pgid=int(data[2]),sid=int(data[3]),
                start=int(data[19]),cpu_ticks=int(data[11])+int(data[12]),threads=int(data[17]))


def owned_snapshot(leader,tracked):
    records = {}
    for name in Path('/proc').iterdir():
        if name.name.isdecimal():
            try:
                row = proc(int(name.name))
                records[row['pid']] = row
            except (ProcessLookupError, FileNotFoundError, PermissionError):
                continue
    active = {pid: row for pid,row in records.items() if row['sid']==leader or tracked.get(pid)==row['start']}
    while True:
        extra = {pid: row for pid,row in records.items() if row['ppid'] in active and pid not in active}
        if not extra:
            break
        active.update(extra)
    for row in active.values():
        pid = row['pid']
        tracked[pid] = row['start']
        try:
            lines = Path('/proc',str(pid),'status').read_text().splitlines()
            row['rss_bytes'] = int(next((line.split()[1] for line in lines if line.startswith('VmRSS:')),0))*1024
        except (ProcessLookupError, FileNotFoundError):
            row['rss_bytes'] = 0
    return active


def signal_owned(active, sig):
    for p in active.values():
        try:
            if proc(p['pid'])['start']==p['start'] and p['state']!='Z':
                os.kill(p['pid'],sig)
        except (ProcessLookupError, FileNotFoundError):
            continue


def cleanup_owned(child, tracked):
    for sig,delay in [(signal.SIGTERM,15),(signal.SIGKILL,5)]:
        deadline = time.monotonic()+delay
        while True:
            child.poll()
            active = {pid:p for pid,p in owned_snapshot(child.pid,tracked).items() if p['state']!='Z'}
            if not active:
                return child.wait()
            signal_owned(active,sig)
            if time.monotonic()>=deadline:
                break
            time.sleep(.25)
    raise RuntimeError('Task-owned processes still present after cleanup')


def summarize(out, stage_resources):
    modes = ['serial','process_2','process_4']
    summary = {}
    for workload in ['L06','L12','L24','receivers']:
        rows = {}
        for mode in modes:
            ids = [f'{mode}_r{repeat}_{workload}' for repeat in (1,2)]
            samples = [read(out/mode/f'{task}.json') for task in ids]
            assert all(r['status']=='PASS_PROTOTYPE_EQUIVALENCE' and r['failed_checks']==0 for r in samples)
            if workload=='receivers':
                keys = ['PS_setup_s','PS_receive_s','FullFK_setup_s','FullFK_receive_s']
                values = {key: [r[key] for r in samples] for key in keys}
            else:
                keys = ['base_assembly_s','near_replacement_s','dense_solve_s','receivers_s','boundary_residual_s']
                values = {key: [r['timings_s'][key] for r in samples] for key in keys}
                values.update({key: [r[key] for r in samples] for key in ['bie_wall_s','independent_audit_s']})
            resources = [stage_resources.get(task,{}) for task in ids]
            assert all(r.get('samples',0)>0 for r in resources)
            timing = {key: dict(samples_s=v,median_s=statistics.median(v)) for key,v in values.items()}
            rows[mode] = dict(timings=timing,resources=resources)
        for mode in modes:
            for key,val in rows[mode]['timings'].items():
                val['speedup_vs_serial'] = rows['serial']['timings'][key]['median_s']/max(val['median_s'],1e-9)
        summary[workload] = rows
    pools = {mode:read(out/mode/'pool.json') for mode in modes}
    # Screening estimates from the largest pilot, with the receiver count held fixed.
    forecasts = {}
    for mode in modes:
        times = summary['L24'][mode]['timings']
        high = {key:max(row['samples_s']) for key,row in times.items()}
        rows = []
        for name,n in [('C0',16000),('C1',19200),('C2',19200),('C3',22400),('C4',21760)]:
            ratio = n/1536
            components = dict(assembly_s=high['base_assembly_s']*ratio**2,
                              near_s=high['near_replacement_s']*ratio,
                              dense_solve_s=high['dense_solve_s']*ratio**3,
                              boundary_s=high['boundary_residual_s']*ratio**2,
                              independent_19_s=high['independent_audit_s']*ratio,
                              receive_s=high['receivers_s']*ratio)
            rows.append(dict(layer=name,nodes=n,matrix_bytes=16*n*n,screening_s=components,
                             screening_total_s=sum(components.values())))
        forecasts[mode] = rows
    write(out/'resource_summary.json',dict(status='REVIEW_REQUIRED',workloads=summary,pools=pools,
          timing_extrapolation_only=forecasts,
          rss_sampling_period_s=.25,full_domain_admitted=False))


def seal_stages(root,checkpoint,completed):
    out = root/'output'; folder = out/'stage_seals'; folder.mkdir(exist_ok=True)
    for stage in checkpoint.get('completed',[]):
        if stage in completed:
            continue
        if stage=='SURFACES':
            files = sorted(p for p in (out/'surfaces').rglob('*') if p.is_file())
        elif stage=='PCT_CHECKOUT':
            files = [out/'pct.json']
        else:
            hit = re.fullmatch(r'(serial|process_[24])_r[12]_(L06|L12|L24|receivers)',stage)
            assert hit, 'Unexpected stage identity'
            files = [out/hit[1]/f'{stage}.json',out/hit[1]/'pool.json']
        assert files
        write(folder/f'{stage}.json',dict(stage=stage,manifest_sha256=checkpoint['manifest_sha256'],
              completed_at_observed_elapsed_s=checkpoint['elapsed_s'],
              files=[dict(path=p.relative_to(root).as_posix(),sha256=sha(p)) for p in files]))
        completed.add(stage)


def run(root):
    import fcntl
    out = root/'output'; t0 = time.monotonic()
    status = dict(status='FAILED_BEFORE_MATLAB',matlab_exit_code=None,supervisor_exit_code=78)
    child = None; tracked = {}; stage_resources = {}; registered_groups = {}; peak = 0; baseline = None
    stopping = []; completed = set(); last_cpu = {}
    signal.signal(signal.SIGTERM,lambda *_: stopping.append('SUPERVISOR_TERM'))
    signal.signal(signal.SIGINT,lambda *_: stopping.append('SUPERVISOR_INT'))
    try:
        with (REMOTE/'d1b_v2_serial.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            cfg,binding = preflight(root)
            baseline = available(); reserve = math.ceil(baseline*cfg['reserve_fraction']); cap = cfg['rss_cap_GiB']*GIB
            affinity = len(os.sched_getaffinity(0))
            disk = shutil.disk_usage(out).free
            load1 = os.getloadavg()[0]
            memory,io = pressure('memory'),pressure('io')
            assert MATLAB.is_file() and os.access(MATLAB,os.X_OK), 'R2025b executable unavailable'
            assert baseline>=cfg['min_available_GiB']*GIB and cap<=baseline-reserve, 'RAM preflight failed'
            assert disk>=5*GIB and os.statvfs(out).f_favail>=5000, 'Disk/inode preflight failed'
            assert affinity>=5 and load1<=21 and memory<=5 and io<=10, 'CPU/load/PSI preflight failed'
            write(out/'host.json',dict(cpu_affinity_count=affinity,initial_available_bytes=baseline,
                  initial_load1=load1,initial_memory_full_avg10=memory,initial_io_full_avg10=io,disk_free_bytes=disk,
                  reserve_bytes=reserve,owned_rss_cap_bytes=cap,cpu_ticks_per_second=os.sysconf('SC_CLK_TCK'),
                  matlab_executable=str(MATLAB),manifest_sha256=binding))
            env = dict(os.environ,D1B_ROOT=str(root),D1B_OUTPUT=str(out),D1B_MANIFEST_SHA=binding,
                       OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
            expr = ("addpath(fullfile(getenv('D1B_ROOT'),'code')); "
                    "run_d1b_resource_prototype(getenv('D1B_ROOT'),getenv('D1B_OUTPUT'),getenv('D1B_MANIFEST_SHA'))")
            with (out/'matlab.log').open('xb') as log, (out/'resources.csv').open('x',newline='') as stream:
                child = subprocess.Popen([str(MATLAB),'-singleCompThread','-batch',expr],cwd=root,
                                         stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                                         env=env,start_new_session=True)
                tracked[child.pid] = proc(child.pid)['start']
                writer = csv.writer(stream)
                writer.writerow(['elapsed_s','stage','total_RSS_bytes','group_RSS_bytes','main_registered_workers_RSS_bytes',
                                 'available_bytes','load1','memory_full_avg10','io_full_avg10','owned_processes_json'])
                bad_pressure = 0; missing = 0; min_available = baseline
                while child.poll() is None:
                    elapsed = time.monotonic()-t0
                    cp_path = out/'checkpoint.json'
                    cp = read(cp_path) if cp_path.is_file() else {}
                    if cp:
                        assert cp['manifest_sha256']==binding
                        seal_stages(root,cp,completed)
                    pool = out/'active_pool.json'
                    worker_ids = read(pool).get('worker_pid',[]) if pool.is_file() else []
                    active = owned_snapshot(child.pid,tracked)
                    now = available(); min_available = min(min_available,now)
                    stage = cp.get('active_stage','MATLAB_STARTUP')
                    main_pid = cp.get('main_pid',child.pid)
                    group_rss = sum(row['rss_bytes'] for row in active.values() if row['pgid']==child.pid)
                    total = sum(row['rss_bytes'] for row in active.values())
                    registered_rss = sum(active[pid]['rss_bytes'] for pid in {main_pid,*worker_ids} if pid in active)
                    peak = max(total,peak)
                    pressure_m,pressure_io = pressure('memory'),pressure('io')
                    writer.writerow([round(elapsed,3),stage,total,group_rss,registered_rss,now,os.getloadavg()[0],
                                     pressure_m,pressure_io,json.dumps(list(active.values()),separators=(',',':'))]); stream.flush()
                    item = stage_resources.setdefault(stage,dict(samples=0,peak_owned_RSS_bytes=0,peak_group_RSS_bytes=0,
                                                                  peak_main_worker_RSS_bytes=0,min_available_bytes=baseline,
                                                                  sampled_peak_CPU_cores=0))
                    item['samples'] += 1
                    item['peak_owned_RSS_bytes'] = max(item['peak_owned_RSS_bytes'],total)
                    item['peak_group_RSS_bytes'] = max(item['peak_group_RSS_bytes'],group_rss)
                    item['peak_main_worker_RSS_bytes'] = max(item['peak_main_worker_RSS_bytes'],registered_rss)
                    item['min_available_bytes'] = min(item['min_available_bytes'],now)
                    last = last_cpu.get(stage)
                    cpu = {pid:(row['start'],row['cpu_ticks']) for pid,row in active.items()}
                    if last:
                        ticks = sum(value[1]-last[1][pid][1] for pid,value in cpu.items()
                                    if pid in last[1] and value[0]==last[1][pid][0])
                        cores = ticks/os.sysconf('SC_CLK_TCK')/max(elapsed-last[0],.001)
                        item['sampled_peak_CPU_cores'] = max(item['sampled_peak_CPU_cores'],cores)
                    last_cpu[stage] = (elapsed,cpu)
                    bad_pressure = bad_pressure+1 if pressure_m>10 or pressure_io>30 else 0
                    if main_pid in active and all(pid in active for pid in worker_ids):
                        missing = 0
                        if worker_ids:
                            registered_groups[str(len(worker_ids))] = sorted(worker_ids)
                    else:
                        missing += 1
                    checks = [(now<reserve,'AVAILABLE_RAM_RESERVE'),(total>cap,'OWNED_RSS_CAP'),
                              (elapsed>cfg['max_wall_s'],'WALL_TIME'),(bad_pressure>=40,'HOST_PRESSURE'),
                              (missing>=8,'WORKER_OR_MAIN_NOT_OWNED')]
                    reasons = stopping+[label for fail,label in checks if fail]
                    if reasons:
                        status['stop_reason'] = reasons[0]
                        break
                    time.sleep(.25)
                status['matlab_exit_code'] = cleanup_owned(child,tracked)
                status.update(min_available_bytes=min_available,peak_owned_RSS_bytes=peak,
                              registered_workers_seen=registered_groups)
                if 'stop_reason' in status:
                    status.update(status='RESOURCE_STOPPED',supervisor_exit_code=78)
                    return 78
                status['status'] = 'MATLAB_OR_AUDIT_FAILED'
                assert status['matlab_exit_code']==0, 'Nonzero MATLAB exit'
                cp = read(out/'checkpoint.json')
                assert cp['status']=='PROTOTYPE_NUMERICS_PASS' and cp['failed_checks']==0
                assert cp['manifest_sha256']==binding and len(cp['completed'])==26
                seal_stages(root,cp,completed)
                assert len(completed)==26
                assert set(registered_groups)=={'2','4'}, 'Missing live worker ownership evidence'
                assert all(registered_groups[str(n)]==sorted(read(out/f'process_{n}'/'pool.json')['worker_pid']) for n in [2,4])
                summarize(out,stage_resources)
                assert preflight(root)[1]==binding
                status.update(status='PROTOTYPE_COMPLETE_REVIEW_REQUIRED',supervisor_exit_code=0)
                return 0
    except Exception as exc:
        status['error'] = f'{type(exc).__name__}: {exc}'
        if child is not None and status['matlab_exit_code'] is None:
            try:
                status['matlab_exit_code'] = cleanup_owned(child,tracked)
            except Exception as cleanup_error:
                status['cleanup_error'] = str(cleanup_error)
        return 78
    finally:
        status.update(elapsed_s=time.monotonic()-t0,initial_available_bytes=baseline,
                      peak_owned_RSS_bytes=peak,next_allowed_step='STOP; user review; no complete D1b solve')
        write(out/'stage_resources.json',stage_resources)
        archive_files = [p for p in sorted(out.rglob('*')) if p.is_file() and not p.is_symlink()
                         and 'pct_jobs_2' not in p.parts and 'pct_jobs_4' not in p.parts]
        status['stage_files_sha256'] = [dict(path=p.relative_to(root).as_posix(),sha256=sha(p)) for p in archive_files]
        write(out/'exit.json',status)


def start(root):
    preflight(root)
    out = root/'output'; out.mkdir(mode=0o700)
    with (root/'supervisor.log').open('xb') as log:
        job = subprocess.Popen([sys.executable,'-B',str(root/'code/d1b_remote.py'),'run',str(root)],
                               stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    write(root/'started.json',dict(supervisor=proc(job.pid),started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    print(f'STARTED {root.name} supervisor_pid={job.pid}')


def status(root):
    out = root/'output'
    if (out/'exit.json').is_file():
        final = read(out/'exit.json'); final.pop('stage_files_sha256',None)
        print(json.dumps(final,indent=2)); return
    if not (root/'started.json').is_file():
        print('UPLOADED_OR_LAUNCH_INCOMPLETE'); return
    original = read(root/'started.json')['supervisor']
    try:
        current = proc(original['pid']); alive = current['start']==original['start'] and current['state']!='Z'
    except (FileNotFoundError,ProcessLookupError):
        alive = False
    cp_path = out/'checkpoint.json'
    cp = read(cp_path) if cp_path.is_file() else {}
    print(json.dumps(dict(status='RUNNING' if alive else 'SUPERVISOR_EXIT_MISSING',
                          active_stage=cp.get('active_stage','STARTUP'),completed=cp.get('completed',[]),
                          error=cp.get('error'),supervisor_pid=original['pid']),indent=2))


def fetch(root):
    out = root/'output'; exit_info = read(out/'exit.json')
    assert 'cleanup_error' not in exit_info, 'Process cleanup unresolved'
    verify(root,exit_info['stage_files_sha256'])
    preflight(root)
    started = read(root/'started.json')['supervisor']
    try:
        current = proc(started['pid'])
        assert current['start']!=started['start'] or current['state']=='Z', 'Supervisor is still running'
    except (FileNotFoundError,ProcessLookupError):
        pass
    result = root/'retrieval'; result.mkdir(exist_ok=True)
    archive = result/f'result_{time.time_ns()}.tar'
    files = [root/r['path'] for r in exit_info['stage_files_sha256']]
    files += [out/'exit.json',root/'supervisor.log',root/'started.json',root/'manifest.json',root/'config.json']
    manifest = dict(files=[dict(path=p.relative_to(root).as_posix(),sha256=sha(p)) for p in files])
    index = archive.with_suffix('.json'); write(index,manifest)
    with tarfile.open(archive,'x') as tar:
        for p in files:
            tar.add(p,arcname=p.relative_to(root).as_posix(),recursive=False)
        tar.add(index,arcname='retrieval_manifest.json',recursive=False)
    print(f'ARCHIVE_PATH={archive}\nARCHIVE_SHA256={sha(archive)}')


def extract(archive,dest,expected=None):
    assert not dest.exists()
    if expected:
        assert sha(archive)==expected
    with tarfile.open(archive,'r') as tar:
        members = tar.getmembers()
        names = set()
        for item in members:
            name = item.name; relative = PurePosixPath(name)
            assert item.isfile() and not relative.is_absolute() and '..' not in relative.parts
            assert re.fullmatch(r'[A-Za-z0-9_./-]+',name) and relative.as_posix()==name
            assert name.lower() not in names
            names.add(name.lower())
        dest.mkdir(parents=True)
        for item in members:
            path = dest/item.name; path.parent.mkdir(parents=True,exist_ok=True)
            with tar.extractfile(item) as src, path.open('xb') as dst:
                shutil.copyfileobj(src,dst)
    manifest = read(dest/'retrieval_manifest.json')
    assert names=={r['path'].lower() for r in manifest['files']}|{'retrieval_manifest.json'}
    verify(dest,manifest['files'])
    print('FETCH_VERIFIED='+str(dest))


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Prototype requires Python assertions enabled')
    op = sys.argv[1]
    path = Path(sys.argv[2]).absolute()
    if op=='verify':
        print('INPUT_MANIFEST_SHA256='+preflight(path)[1])
    elif op=='extract':
        extract(path,Path(sys.argv[3]).absolute(),sys.argv[4])
    else:
        owned_root(path)
        if op=='run':
            sys.exit(run(path))
        {'start':start,'status':status,'fetch':fetch}[op](path)
