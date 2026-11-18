"""Read-only K1 host resource inventory. No MATLAB, network, or job mutation."""
import datetime
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time


def command(argv):
    proc = subprocess.run(argv, capture_output=True, text=True, timeout=10, env=dict(os.environ, LC_ALL='C'))
    return dict(exit_code=proc.returncode, stdout=proc.stdout.strip(), stderr=proc.stderr.strip())


def read(path):
    p = Path(path)
    return p.read_text().strip() if p.is_file() else None


def cpu_ticks():
    values = [int(x) for x in Path('/proc/stat').read_text().splitlines()[0].split()[1:9]]
    return sum(values), values[3] + values[4]


def main():
    root = Path('/home/misc/cwx')
    assert root.is_dir()
    a = cpu_ticks()
    time.sleep(1)
    b = cpu_ticks()
    memory = {}
    for line in Path('/proc/meminfo').read_text().splitlines():
        key, value = line.split(':', 1)
        if key in ['MemTotal', 'MemAvailable', 'SwapTotal', 'SwapFree', 'Dirty', 'Writeback', 'Committed_AS', 'CommitLimit']:
            memory[key + '_bytes'] = int(value.split()[0])*1024
    proc = command(['ps', '-eo', 'pid,ppid,pcpu,pmem,rss,comm', '--sort=-pcpu'])
    lines = proc['stdout'].splitlines()[1:]
    processes = []
    matlab = []
    for line in lines:
        item = line.split(None, 5)
        if len(item) != 6:
            continue
        row = dict(pid=int(item[0]), ppid=int(item[1]), CPU_percent=float(item[2]), RAM_percent=float(item[3]), RSS_bytes=int(item[4])*1024, comm=item[5])
        if len(processes) < 12:
            processes.append(row)
        if item[5].lower() in ['matlab', 'matlab_helper']:
            matlab.append(row)
    vfs = os.statvfs(root)
    executable = Path('/home/misc/MATLAB/R2025b/bin/matlab')
    limits = {}
    for name in ['RLIMIT_NOFILE', 'RLIMIT_NPROC', 'RLIMIT_AS', 'RLIMIT_STACK', 'RLIMIT_CPU']:
        limits[name] = resource.getrlimit(getattr(resource, name))
    lscpu = command(['lscpu', '--json'])
    result = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  scope='read-only host resources; no solver or MATLAB launched',
                  lscpu=json.loads(lscpu['stdout']) if lscpu['exit_code'] == 0 else lscpu,
                  affinity_logical_CPUs=sorted(os.sched_getaffinity(0)), load_average=os.getloadavg(),
                  CPU_busy_fraction_1s=1-(b[1]-a[1])/max(1,b[0]-a[0]), memory=memory,
                  pressure={k:read('/proc/pressure/'+k) for k in ['cpu','memory','io']},
                  top_processes_resources_only=processes, matlab_processes_resources_only=matlab,
                  process_scope='No arguments, open files, credentials or other-user sessions inspected',
                  disk=dict(path=str(root),free_bytes=vfs.f_bavail*vfs.f_frsize,total_bytes=vfs.f_blocks*vfs.f_frsize,
                            free_inodes=vfs.f_favail,total_inodes=vfs.f_files),
                  limits=limits,
                  cgroups={p:read('/sys/fs/cgroup/'+p) for p in ['cpu.max','cpuset.cpus.effective','memory.max','memory.current','memory.swap.max']},
                  matlab=dict(executable_path=str(executable), exists=executable.is_file(), is_executable=os.access(executable,os.X_OK),
                              PCT_directory_exists=Path('/home/misc/MATLAB/R2025b/toolbox/parallel').is_dir(),
                              license_environment_present={k:bool(os.getenv(k)) for k in ['MLM_LICENSE_FILE','LM_LICENSE_FILE']}),
                  required_tools={k:bool(shutil.which(k)) for k in ['python3','flock','setsid','nohup','tar','sha256sum','timeout']})
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
