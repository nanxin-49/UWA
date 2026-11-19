"""Manual scale pilot supervisor. Import/dry-run never connects or launches."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import time
import d1b_remote as base
from d1_package import sha, verify

GIB=1024**3
REMOTE=Path('/home/misc/cwx')
SCOPE='D1B_SCALE_PILOT_ONLY'
read,write=base.read,base.write


def preflight(root):
    cfg=read(root/'config.json'); manifest=read(root/'manifest.json')
    assert manifest['scope']==cfg['scope']==SCOPE and cfg['full_d1b_enabled'] is False
    verify(root,manifest['files'])
    assert cfg['nodes']==[4096,8192] and cfg['half_widths_m']==[64,128] and cfg['control_nodes']==1536
    assert cfg['workers']==2 and cfg['threads']==1 and not cfg['automatic_retry'] and not cfg['increase_workers']
    assert cfg['max_wall_s']==1800 and cfg['action_cutoff_s']==1770
    assert cfg['min_available_GiB']==18 and cfg['rss_cap_GiB']==12 and cfg['reserve_fraction']==.3
    assert cfg['min_disk_GiB']==20 and cfg['output_cap_GiB']==8
    assert cfg['source']==dict(frequency_hz=4000,c_mps=1500,rho_kgpm3=1025,v0_mps=.001,total_width_m=.1,z_source_m=1490)
    assert cfg['PM']==dict(period_m=640,U_PM_mps=6,seed=260001,N=160,low_frequency_enabled=False,plateau_m=220,support_m=240)
    assert cfg['BIE']==dict(ppw=12,panel_order=8,self_order=96,h_m=-2,window_plateau_ratio=.85)
    assert cfg['receiver_count']==2279 and cfg['receiver_depth_m']==6 and cfg['receiver_half_width_m']==200
    assert cfg['gates']==dict(matrix=1e-11,rhs=1e-11,density=1e-10,receiver=1e-11,midpoint_equivalence=1e-10,independent_equivalence=1e-10,linear=1e-10,boundary=1e-4)
    archive=root/'input_package.tar'
    if archive.exists():
        assert sha(archive)==(root/'input_package.tar.sha256').read_text().split()[0]
        with tarfile.open(archive) as tar:
            members=tar.getmembers(); expected={r['path']:r['sha256'] for r in manifest['files']}
            expected['manifest.json']=sha(root/'manifest.json')
            assert len(members)==len(expected) and {m.name for m in members}==set(expected)
            for m in members:
                assert m.isfile() and hashlib.sha256(tar.extractfile(m).read()).hexdigest()==expected[m.name]
    return cfg,sha(root/'manifest.json')


def owned_root(root):
    assert root.parent==REMOTE and root.resolve()==root and root.is_dir()
    assert re.fullmatch(r'd1b_scale_[A-Za-z0-9_-]{1,64}',root.name)


def stop_reason(cfg,available,total,elapsed,baseline,disk,output_bytes):
    for bad,label in [(available<math.ceil(baseline*.3),'AVAILABLE_RAM_RESERVE'),
                      (total>cfg['rss_cap_GiB']*GIB,'OWNED_RSS_CAP'),
                      (elapsed>=cfg['action_cutoff_s'],'WALL_ACTION_CUTOFF'),
                      (disk<2*GIB,'RUNNING_DISK_RESERVE'),
                      (output_bytes>cfg['output_cap_GiB']*GIB,'OUTPUT_CAP')]:
        if bad: return label
    return None


def cpu_delta(previous,active,elapsed,hz):
    if previous is None: return 0.0
    then,old=previous
    ticks=sum(max(0,p['cpu_ticks']-old[pid][1]) for pid,p in active.items()
              if pid in old and old[pid][0]==p['start'])
    return ticks/hz/max(elapsed-then,.001)


def verify_case(out,name,binding):
    folder=out/name; seal=read(folder/'seal.json'); audit=read(folder/'audit.json')
    assert seal['status']=='AUDITED_PASS' and seal['manifest_sha256']==binding
    assert audit['status']=='PASS_INDEPENDENT_SAVED_AUDIT' and audit['failed_checks']==0
    assert audit['linear_relative']<1e-10 and audit['midpoint_relative']<1e-4 and audit['independent_relative']<1e-4
    assert seal['data_sha256']==audit['data_sha256']==sha(folder/'data.mat')
    assert seal['audit_sha256']==sha(folder/'audit.json')
    assert seal['independent_sha256']==sha(folder/'independent.mat')
    expected = {'control_serial': (1536,0), 'control_process2': (1536,2), 'N4096': (4096,2), 'N8192': (8192,2)}
    assert (audit['nodes'],audit['workers'])==expected[name] and seal['nodes']==audit['nodes']
    return audit,seal


def summarize(out,resources,binding):
    sealed={name:verify_case(out,name,binding) for name in ['control_serial','control_process2','N4096','N8192']}
    cases={name:value[0] for name,value in sealed.items()}
    times={name:a['timings_s'] for name,a in cases.items()}
    components=['base_assembly_s','near_replacement_s','dense_solve_s','boundary_residual_s','receivers_s']
    rows=[]
    for n in [4096,8192]:
        name=f'N{n}'; a=cases[name]; t=times[name]
        total=sum(t[k] for k in components)+a['independent_audit_s']
        wall=sealed[name][1]['total_case_s']
        phase_set=['base_assembly','near_replacement','dense_solve','boundary_residual','independent_GL112']
        rss = {p:resources[f'{name}.{p}']['peak_owned_RSS_bytes'] for p in phase_set}
        rows.append(dict(nodes=n,case_wall_s=wall,numerical_components_s=total,dense_solve_s=t['dense_solve_s'],
                         dense_fraction=t['dense_solve_s']/total,
                         dense_fraction_of_case_wall=t['dense_solve_s']/wall,
                         **{f'{p}_sampled_peak_RSS_bytes':v for p,v in rss.items()},
                         ideal_2_to_4_worker_saving_upper_bound=(t['base_assembly_s']+t['boundary_residual_s']+t['receivers_s']+a['independent_audit_s'])/(2*total)))
    forecasts=[]
    exponents={'base_assembly_s':2,'near_replacement_s':1,'dense_solve_s':3,'boundary_residual_s':2,'receivers_s':1}
    for name,n in [('N4096',4096),('N8192',8192)]:
        values={k:times[name][k]*(22400/n)**p for k,p in exponents.items()}
        values['independent_audit_s']=cases[name]['independent_audit_s']*(22400/n)
        forecasts.append(dict(source_nodes=n,target_nodes=22400,components_s=values,total_s=sum(values.values())))
    ratio=8192/4096
    observed_orders={k:math.log(times['N8192'][k]/times['N4096'][k])/math.log(ratio) for k in components}
    # Peak RSS is an observation. A runtime+N^2 extrapolation remains an estimate.
    case_peak=lambda name:max(r['peak_owned_RSS_bytes'] for key,r in resources.items() if key.startswith(name+'.'))
    peak4,peak8=case_peak('N4096'),case_peak('N8192')
    coeff=max(0,(peak8-peak4)/(8192**2-4096**2)); intercept=max(0,peak8-coeff*8192**2)
    memory_estimate=intercept+coeff*22400**2
    result=dict(status='PILOT_MEASURED_REVIEW_REQUIRED',cases=rows,observed_two_size_orders=observed_orders,
                serial_control_speedup=read(out/'control_serial/seal.json')['total_case_s']/read(out/'control_process2/seal.json')['total_case_s'],
                forecasts_22400_screening_only=forecasts,screening_22400_RSS_bytes=memory_estimate,
                screening_22400_RSS_with_50_percent_margin_bytes=1.5*memory_estimate,
                full_22400_admitted=False,more_workers_admitted=False,
                limitations=['one repetition; not a timing confidence interval','different finite BIE windows; no target-field qualification',
                             'large serial equivalence is all base rows, not a duplicate full serial solve',
                             'communication bytes measured; isolated network time is not identifiable',
                             'RSS sampling may miss short peaks; shared pages counted per process'])
    write(out/'resource_summary.json',result)
    with (out/'timing_summary.csv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    lines=['# D1b-v2 4096/8192 BIE试点实测报告','',
           '两档保存数据独立审计通过。完整22400节点仍未准入；本报告不授权启动。','',
           '|节点|整档含审计与写盘/s|数值阶段合计/s|直接解/s|直接解占整档|直接解占数值阶段|理想2→4 worker节省|','|---|---:|---:|---:|---:|---:|---:|']
    for r in rows: lines.append(f"|{r['nodes']}|{r['case_wall_s']:.3f}|{r['numerical_components_s']:.3f}|{r['dense_solve_s']:.3f}|{r['dense_fraction_of_case_wall']:.2%}|{r['dense_fraction']:.2%}|{r['ideal_2_to_4_worker_saving_upper_bound']:.2%}|")
    lines+=['','上限假定组装/边界/接收/独立审计耗时减半、无新增通信，是乐观界而非4-worker实测。近项替换和稠密直接解仍在单线程客户端；本次没有4-worker试验。',
            '',f'22400阶段时间的两档锚点外推为 {min(r["total_s"] for r in forecasts)/60:.2f}–{max(r["total_s"] for r in forecasts)/60:.2f} min（非置信区间）。',
            f'按两档总RSS拟合的22400筛查约 {memory_estimate/GIB:.2f} GiB，加50%余量为 {1.5*memory_estimate/GIB:.2f} GiB；单矩阵为 {16*22400**2/GIB:.3f} GiB。',
            '完整运行需另审预算、现场可用RAM、LU/临时副本峰值、C0–C4及C2/非局域敏感性。不能将这个试点作为完整PM声场资格。',
            '', '回收后用 render_d1b_scale_saved 从已关闭的 resources.csv 生成 FIG/PNG。阶段峰值见 stage_resources.json，通信见 phase_events.jsonl。',
            '', '一致性范围：1536节点完整串行/2-worker矩阵、RHS、密度、接收、中点和独立残差；4096/8192逐 worker 核全部基础矩阵行、固定密度接收/中点各12行及全部19个GL112点。近项替换及直接解在客户端串行完成。']
    (out/'measured_report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def run(root):
    import fcntl
    out=root/'output'; t0=time.monotonic(); child=None; tracked={}; resources={}; stopping=[]
    status=dict(status='FAILED_BEFORE_MATLAB',matlab_exit_code=None,supervisor_exit_code=78)
    signal.signal(signal.SIGTERM,lambda *_:stopping.append('SUPERVISOR_TERM'))
    signal.signal(signal.SIGINT,lambda *_:stopping.append('SUPERVISOR_INT'))
    try:
        with (REMOTE/'d1b_v2_serial.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            cfg,binding=preflight(root); baseline=base.available(); reserve=math.ceil(.3*baseline)
            assert base.MATLAB.is_file() and os.access(base.MATLAB,os.X_OK)
            assert baseline>=18*GIB and 12*GIB<=baseline-reserve
            assert shutil.disk_usage(root).free>=20*GIB and os.statvfs(root).f_favail>=5000
            assert len(os.sched_getaffinity(0))>=5 and os.getloadavg()[0]<=21
            assert base.pressure('memory')<=5 and base.pressure('io')<=10
            write(out/'host.json',dict(initial_available_bytes=baseline,reserve_bytes=reserve,task_RSS_cap_bytes=12*GIB,
                                      hz=os.sysconf('SC_CLK_TCK'),MATLAB=str(base.MATLAB),manifest_sha256=binding))
            env=dict(os.environ,D1B_SCALE_ROOT=str(root),D1B_SCALE_OUTPUT=str(out),D1B_SCALE_MANIFEST_SHA=binding,
                     OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
            expr="addpath(fullfile(getenv('D1B_SCALE_ROOT'),'code')); run_d1b_scale_pilot(getenv('D1B_SCALE_ROOT'),getenv('D1B_SCALE_OUTPUT'),getenv('D1B_SCALE_MANIFEST_SHA'),'execute')"
            previous=None; bad_pressure=0; missing=0; minimum=baseline; peak=0; registered_seen=set(); output_bytes=0; last_disk_sample=-10
            with (out/'matlab.log').open('xb') as log, (out/'resources.csv').open('x',newline='') as stream:
                child=subprocess.Popen([str(base.MATLAB),'-singleCompThread','-batch',expr],cwd=root,env=env,
                                       stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                status['status']='MATLAB_OR_AUDIT_FAILED'
                tracked[child.pid]=base.proc(child.pid)['start']
                w=csv.writer(stream); w.writerow(['elapsed_s','stage','total_RSS_bytes','group_RSS_bytes','main_registered_workers_RSS_bytes',
                    'available_bytes','reserve_bytes','CPU_cores','main_worker_CPU_cores','load1','memory_full_avg10','io_full_avg10','owned_processes_json'])
                while child.poll() is None:
                    elapsed=time.monotonic()-t0
                    cp=read(out/'checkpoint.json') if (out/'checkpoint.json').exists() else {}
                    ph=read(out/'phase.json') if (out/'phase.json').exists() else {}
                    stage=ph.get('case_id','startup')+'.'+ph.get('phase',cp.get('active_stage','MATLAB_STARTUP'))
                    active=base.owned_snapshot(child.pid,tracked); pool=read(out/'active_pool.json') if (out/'active_pool.json').exists() else {}
                    ids=pool.get('worker_pid',[]); ids=[ids] if isinstance(ids,(int,float)) else ids
                    assert len(ids) in [0,2] and len(set(ids))==len(ids)
                    main=cp.get('main_pid',child.pid); selected={main,*ids}
                    now=base.available(); minimum=min(minimum,now)
                    total=sum(p['rss_bytes'] for p in active.values()); peak=max(peak,total)
                    group=sum(p['rss_bytes'] for p in active.values() if p['pgid']==child.pid)
                    registered=sum(p['rss_bytes'] for pid,p in active.items() if pid in selected)
                    interval=elapsed-previous[0] if previous else 0
                    cpu=cpu_delta(previous,active,elapsed,os.sysconf('SC_CLK_TCK'))
                    registered_cpu=cpu_delta(previous,{pid:p for pid,p in active.items() if pid in selected},elapsed,os.sysconf('SC_CLK_TCK'))
                    previous=(elapsed,{pid:(p['start'],p['cpu_ticks']) for pid,p in active.items()})
                    pm,pi=base.pressure('memory'),base.pressure('io')
                    w.writerow([elapsed,stage,total,group,registered,now,reserve,cpu,registered_cpu,os.getloadavg()[0],pm,pi,json.dumps(list(active.values()),separators=(',',':'))]); stream.flush()
                    r=resources.setdefault(stage,dict(samples=0,peak_owned_RSS_bytes=0,peak_main_worker_RSS_bytes=0,min_available_bytes=baseline,CPU_core_seconds=0,peak_CPU_cores=0))
                    r['samples']+=1; r['peak_owned_RSS_bytes']=max(r['peak_owned_RSS_bytes'],total)
                    r['peak_main_worker_RSS_bytes']=max(r['peak_main_worker_RSS_bytes'],registered); r['min_available_bytes']=min(r['min_available_bytes'],now)
                    r['CPU_core_seconds']+=cpu*interval; r['peak_CPU_cores']=max(r['peak_CPU_cores'],cpu)
                    bad_pressure=bad_pressure+1 if pm>10 or pi>30 else 0
                    if main in active and all(pid in active for pid in ids):
                        missing=0; registered_seen.update(ids)
                    else: missing+=1
                    if elapsed-last_disk_sample>=2:
                        output_bytes=sum(p.stat().st_size for p in out.rglob('*') if p.is_file() and not p.is_symlink()); last_disk_sample=elapsed
                    reason=stop_reason(cfg,now,total,elapsed,baseline,shutil.disk_usage(root).free,output_bytes)
                    if bad_pressure>=40: reason='HOST_PRESSURE'
                    if missing>=8: reason='MAIN_OR_WORKER_NOT_OWNED'
                    if stopping: reason=stopping[0]
                    if reason:
                        status.update(status='RESOURCE_STOPPED',stop_reason=reason)
                        break
                    time.sleep(.25)
                status['matlab_exit_code']=base.cleanup_owned(child,tracked)
                status.update(peak_owned_RSS_bytes=peak,min_available_bytes=minimum,registered_workers_seen=sorted(registered_seen))
            assert 'stop_reason' not in status, 'Resource/ownership stop'
            assert status['matlab_exit_code']==0, 'MATLAB nonzero exit'
            cp=read(out/'checkpoint.json')
            assert cp['status']=='PILOT_NUMERICS_PASS_REVIEW_REQUIRED' and cp['failed_checks']==0
            assert cp['completed']==['CONTROL_EQUIVALENCE_PASS','N4096_INDEPENDENT_PASS','N8192_INDEPENDENT_PASS']
            pool_info=read(out/'pool.json')
            assert pool_info['profile']=='Processes' and pool_info['compute_threads']==[1,1]
            assert len(registered_seen)==2 and registered_seen==set(pool_info['worker_pid'])
            assert read(out/'pct.json')['checkout']==1
            summarize(out,resources,binding); assert preflight(root)[1]==binding
            assert time.monotonic()-t0<cfg['max_wall_s'], 'Closeout exceeded total budget'
            status.update(status='SCALE_PILOT_COMPLETE_REVIEW_REQUIRED',supervisor_exit_code=0)
    except Exception as exc:
        status['error']=f'{type(exc).__name__}: {exc}'
        if child is not None and status['matlab_exit_code'] is None:
            try: status['matlab_exit_code']=base.cleanup_owned(child,tracked)
            except Exception as err: status['cleanup_error']=str(err)
    finally:
        write(out/'stage_resources.json',resources)
        status.update(elapsed_s=time.monotonic()-t0,next_allowed_step='STOP; manual review; never 22400 or automatic retry')
        files=[p for p in sorted(out.rglob('*')) if p.is_file() and not p.is_symlink() and 'pct_jobs_2' not in p.parts]
        status['stage_files_sha256']=[dict(path=p.relative_to(root).as_posix(),sha256=sha(p)) for p in files]
        write(out/'exit.json',status)
    return status['supervisor_exit_code']


def start(root):
    preflight(root); (root/'output').mkdir(mode=0o700)
    with (root/'supervisor.log').open('xb') as log:
        job=subprocess.Popen([sys.executable,'-B',str(root/'code/d1b_scale_remote.py'),'run',str(root)],
                             stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    write(root/'started.json',dict(supervisor=base.proc(job.pid),started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    print(f'STARTED {root.name} supervisor_pid={job.pid}')


def fetch(root):
    final=read(root/'output/exit.json'); assert 'cleanup_error' not in final
    verify(root,final['stage_files_sha256']); preflight(root)
    original=read(root/'started.json')['supervisor']
    try:
        p=base.proc(original['pid']); assert p['start']!=original['start'] or p['state']=='Z','Supervisor still running'
    except (FileNotFoundError,ProcessLookupError): pass
    result=root/'retrieval'; result.mkdir(exist_ok=True); archive=result/f'result_{time.time_ns()}.tar'
    files=[root/r['path'] for r in final['stage_files_sha256']]+[root/'output/exit.json',root/'started.json',root/'supervisor.log',root/'manifest.json',root/'config.json']
    index=archive.with_suffix('.json'); write(index,dict(files=[dict(path=p.relative_to(root).as_posix(),sha256=sha(p)) for p in files]))
    with tarfile.open(archive,'x') as tar:
        for p in files: tar.add(p,arcname=p.relative_to(root).as_posix(),recursive=False)
        tar.add(index,arcname='retrieval_manifest.json',recursive=False)
    print(f'ARCHIVE_PATH={archive}\nARCHIVE_SHA256={sha(archive)}')


if __name__=='__main__':
    if not __debug__: raise RuntimeError('Assertions must remain enabled')
    op=sys.argv[1]; path=Path(sys.argv[2]).absolute()
    if op=='verify': print('INPUT_MANIFEST_SHA256='+preflight(path)[1])
    elif op=='extract': base.extract(path,Path(sys.argv[3]).absolute(),sys.argv[4])
    else:
        owned_root(path)
        if op=='run': sys.exit(run(path))
        {'start':start,'status':base.status,'fetch':fetch}[op](path)
