"""空闲核调度器（T9、T6b 用）：只用空闲的核，不挤占 T8。
按作业文件顺序逐个启动作业（格式同 run_jobs.py："<完成标志文件>\t<超时秒>\t<命令>"），
启动前要求机器上正在运行的 t6_job / t6_rerun / t8_job / t9_job 作业总数 < 上限（默认 12），否则每 60 秒再看一次。
"检查名额 + 启动"在全局文件锁里做，多个调度器同时运行也不会超过上限；某个作业的命令已经有进程在跑（例如上一个调度器启动的）就跳过。
T8 的调度器并行数固定，T9 只占用它们空出来的名额；T8 跑完后 T9 自然用上更多核。
完成标志已存在就跳过，所以中断后重跑同一条命令即可续跑。
用法：nohup python3 code/t9_sched.py <日志目录> <上限> <作业文件>... > <调度日志> 2>&1 &"""
import os, sys, subprocess, time, fcntl

JOBS = (b"t6_job.py", b"t6_rerun.py", b"t8_job.py", b"t9_job.py")
LOCK = "/tmp/mdrr_idle_sched.lock"
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


def running():
    """数正在运行的作业进程。"""
    return len(procs_args())


def main():
    """逐个作业：等到有空闲名额再启动，后台运行；全部启动后等它们结束。"""
    logdir, cap, files = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    os.makedirs(logdir, exist_ok=True)
    procs = []
    k = 0
    for jf in files:
        for line in open(jf):
            done, tlim, cmd = line.rstrip("\n").split("\t")
            k += 1
            if os.path.exists(done): continue
            tail = cmd.split("&& ")[-1].replace("python3 ", "", 1).strip()
            lk = open(LOCK, "a")
            while True:
                fcntl.flock(lk, fcntl.LOCK_EX)
                if os.path.exists(done) or tail in procs_args() or running() < cap: break
                fcntl.flock(lk, fcntl.LOCK_UN); time.sleep(60)
            if os.path.exists(done) or tail in procs_args():
                fcntl.flock(lk, fcntl.LOCK_UN); lk.close(); continue
            log = open(os.path.join(logdir, f"{k:05d}.log"), "w"); log.write(cmd + "\n"); log.flush()
            p = subprocess.Popen(["/usr/bin/time", "-v", "timeout", tlim, "bash", "-c", cmd], stdout=log, stderr=subprocess.STDOUT, env=ENV)
            procs.append(p)
            print(time.strftime("%m-%d %H:%M:%S"), "启动", k, "PID", p.pid, cmd.split("&& ")[-1], flush=True)
            time.sleep(20)     # 等子进程出现在 /proc 里再放锁、数下一次
            fcntl.flock(lk, fcntl.LOCK_UN); lk.close()
    for p in procs: p.wait()
    print(time.strftime("%m-%d %H:%M:%S"), "全部结束", flush=True)


if __name__ == "__main__":
    main()
