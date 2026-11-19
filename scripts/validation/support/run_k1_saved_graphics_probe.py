"""One approved <=180 s saved-array graphics probe; no solver on path."""
import json, math, os, subprocess, sys, time
from pathlib import Path

root=Path(sys.argv[1]).resolve(); assert root.name.startswith('k1_saved_graphics_diag_') and root.parent==Path('/home/misc/cwx')
sys.path.insert(0,str(root/'code'))
import k1_pm_supervisor as sup
import k1_pm_preflight as pre
start=time.monotonic(); cfg=pre.read(root/'config.json'); cfg['_root']=str(root)
pre.verify_rows(root,pre.read(root/'diagnostic_manifest.json'))
baseline=max(119021486080,sup.available()); child=None; ident=None
result=dict(status='LINUX_SAVED_GRAPHICS_FAILED_STOPPED',failed_checks=1,new_acoustic_solves=0,old_used_wall_s=666,RAM_baseline_bytes=baseline)
try:
    snap=sup.admission(cfg,baseline); prefs=root/'prefs'; prefs.mkdir()
    cmd=['/home/misc/MATLAB/R2025b/bin/matlab','-singleCompThread','-batch',f"addpath('{root}/code'); probe_k1_saved_graphics('{root}');"]
    sup.write(root/'command.json',dict(argv=cmd,admission=snap))
    with (root/'stdout.log').open('wb') as out,(root/'stderr.log').open('wb') as err:
        child=subprocess.Popen(cmd,cwd=root,env=dict(os.environ,MATLAB_PREFDIR=str(prefs),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),stdout=out,stderr=err,start_new_session=True)
        ident=dict(pid=child.pid,process_group=child.pid,start_ticks=sup.start_ticks(child.pid)); assert sup.owned_alive(ident)
        sup.write(root/'identity.json',ident); bad=0; stop=None; peak=0; low=snap['MemAvailable_bytes']
        with (root/'resources.csv').open('w',buffering=1) as f:
            f.write('elapsed_s,MemAvailable_bytes,owned_group_RSS_bytes,memory_full_avg10,io_full_avg10\n')
            while child.poll() is None:
                mem=sup.available(); rss=sup.group_rss(child.pid); mp=sup.pressure('memory'); ip=sup.pressure('io')
                peak=max(peak,rss); low=min(low,mem); bad=bad+1 if mp>10 or ip>30 else 0
                f.write(f'{time.monotonic()-start:.3f},{mem},{rss},{mp},{ip}\n')
                if mem<.3*baseline: stop='RAM_STOP'
                elif rss>cfg['task_RSS_cap_bytes']: stop='RSS_STOP'
                elif bad>=10: stop='PRESSURE_STOP'
                elif time.monotonic()-start>=180: stop='WALL_180S_STOP'
                if stop: sup.stop_owned(child,ident); break
                time.sleep(1)
        code=child.wait()
    result.update(matlab_exit_code=code,stop_reason=stop,peak_RSS_bytes=peak,minimum_available_bytes=low)
    assert code==0 and stop is None
    receipt=pre.read(root/'plot_receipt.json'); assert receipt['failed_checks']==0 and receipt['new_acoustic_solves']==0
    result.update(status='LINUX_SAVED_GRAPHICS_ORIGINAL_GATES_PASS',failed_checks=0)
except Exception as exc:
    result['error']=repr(exc)
    if child is not None and child.poll() is None:
        try: sup.stop_owned(child,ident)
        except Exception as error: result['stop_refused_or_failed']=repr(error)
finally:
    elapsed=time.monotonic()-start
    result.update(elapsed_s=elapsed,charged_wall_s=math.ceil(elapsed),K1_used_wall_s_after=666+math.ceil(elapsed),all_owned_children_reaped=child is None or child.poll() is not None)
    sup.write(root/'receipt.json',result)
print(json.dumps(result))
sys.exit(0 if result['failed_checks']==0 else 1)
