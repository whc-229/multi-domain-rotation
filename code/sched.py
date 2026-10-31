"""带优先级的调度器（第三批 hi、第四批 lo 共用；取代 t9_sched.py 的作用）。
作业文件格式同 run_jobs.py：每行 "<完成标志文件>\t<超时秒>\t<命令>"；完成标志已存在就跳过，中断后重跑同一条命令即可续跑。
- 全机上限：正在运行的作业进程（python3 code/<JOBS 里的脚本>）总数 < cap 才启动；
- hi：从启动到全部作业派出，在 /tmp/mdrr_sched_wait/ 下挂一个 "<自己的 PID>" 文件；
- lo：只要有活着的 hi 在等名额就不启动；另外 lo 自己同时运行的作业数不超过 lo_cap（默认 = cap），长作业不会把核全占住；
- "检查名额 + 启动"在全局文件锁里做（与 t9_sched 同一把锁），多个调度器同时运行也不会超过上限；某个作业的命令已经有进程在跑就跳过。
停调度器时连子进程一起停：调度器把每个作业放进自己的进程组（start_new_session），日志里记下 PID = 进程组号，
需要停时 kill -- -<PGID>（只停记录过的组），再停调度器本身。
作业命令里写 ana 环境 python 的绝对路径（/home/hk/miniconda3/envs/ana/bin/python3）；启动时先检查它能 import numpy / scipy / highspy。
用法：nohup python3 code/sched.py <日志目录> <hi|lo> <cap> [lo_cap] -- <作业文件>... > <调度日志> 2>&1 &"""
import os, sys, subprocess, time, fcntl

JOBS = (b"t6_job.py", b"t6_rerun.py", b"t8_job.py", b"t9_job.py", b"t10_variants.py", b"t11_job.py", b"t12_job.py",
        b"t13_job.py", b"t13_rhmip.py", b"t14_gphh.py", b"t15_job.py", b"t16_job.py",
        b"t17_job.py", b"t18_job.py", b"t12b_es_job.py")
LOCK = "/tmp/mdrr_idle_sched.lock"
WAIT = "/tmp/mdrr_sched_wait"
ENV = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")


def procs_args():
    """正在运行的作业进程的参数表（python3 之后的部分）；只读 /proc，不动任何进程。"""
    out = []
    for p in os.listdir("/proc"):
        if not p.isdigit(): continue
        try: a = open(f"/proc/{p}/cmdline", "rb").read().split(b"\0")
        except OSError: continue
        if len(a) >= 2 and a[0].endswith(b"python3") and any(a[1].endswith(x) for x in JOBS):
            out.append(b" ".join(x for x in a[1:] if x).decode())
    return out


def hi_waiting():
    """是否有活着的 hi 调度器在等名额（PID 文件对应的进程还在）。"""
    if not os.path.isdir(WAIT): return False
    for f in os.listdir(WAIT):
        if f.isdigit() and os.path.exists(f"/proc/{f}"): return True
    return False


def main():
    """逐个作业：等到有名额再启动（后台运行）；全部启动后等它们结束。"""
    a = sys.argv[1:]; sep = a.index("--")
    logdir, prio, cap = a[0], a[1], int(a[2]); lo_cap = int(a[3]) if sep > 3 else cap; files = a[sep + 1:]
    assert prio in ("hi", "lo")
    # 环境检查（10-30 晚审核）：作业命令里的 python 必须能 import numpy / scipy / highspy，否则不启动（会话重启后 PATH 可能换成别的解释器）
    for jf in files:
        for line in open(jf):
            py = line.rstrip("\n").split("\t")[2].split("&& ")[-1].split()[0]
            r = subprocess.run(["bash", "-c", f"{py} -c 'import numpy, scipy, highspy'"], env=ENV, capture_output=True)
            if r.returncode != 0: sys.exit(f"作业文件 {jf} 里的解释器 {py} 不能 import numpy / scipy / highspy，调度器不启动")
            break
    os.makedirs(logdir, exist_ok=True); os.makedirs(WAIT, exist_ok=True)
    me = os.path.join(WAIT, str(os.getpid()))
    if prio == "hi": open(me, "w").close()   # hi 从启动到全部作业派出都挂着标记（10-31 修正：原来只在等名额时挂，一波作业同时结束时 lo 会在两次派出之间插进来）
    procs = []; k = 0
    for jf in files:
        for line in open(jf):
            done, tlim, cmd = line.rstrip("\n").split("\t")
            k += 1
            if os.path.exists(done): continue
            tail = cmd.split("&& ")[-1].split(None, 1)[1].strip()     # 去掉解释器（可能是绝对路径），与 /proc 里 python3 之后的部分比
            lk = open(LOCK, "a")
            while True:
                fcntl.flock(lk, fcntl.LOCK_EX)
                run = procs_args(); mine = sum(p.poll() is None for p in procs)
                if os.path.exists(done) or tail in run: break
                if len(run) < cap and (prio == "hi" or (not hi_waiting() and mine < lo_cap)): break
                fcntl.flock(lk, fcntl.LOCK_UN); time.sleep(30)
            if os.path.exists(done) or tail in procs_args():
                fcntl.flock(lk, fcntl.LOCK_UN); lk.close(); continue
            log = open(os.path.join(logdir, f"{k:05d}.log"), "w"); log.write(cmd + "\n"); log.flush()
            p = subprocess.Popen(["/usr/bin/time", "-v", "timeout", tlim, "bash", "-c", cmd], stdout=log, stderr=subprocess.STDOUT, env=ENV,
                                 start_new_session=True)
            procs.append(p)
            print(time.strftime("%m-%d %H:%M:%S"), "启动", k, "PGID", p.pid, cmd.split("&& ")[-1], flush=True)
            time.sleep(3)      # 等子进程出现在 /proc 里再放锁、数下一次
            fcntl.flock(lk, fcntl.LOCK_UN); lk.close()
    if os.path.exists(me): os.remove(me)
    for p in procs: p.wait()
    print(time.strftime("%m-%d %H:%M:%S"), "全部结束", flush=True)


if __name__ == "__main__":
    main()
