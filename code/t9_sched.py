"""T9 的调度器：只用空闲的核，不挤占 T8。
按作业文件顺序逐个启动 T9 作业（格式同 run_jobs.py："<完成标志文件>\t<超时秒>\t<命令>"），
启动前要求机器上正在运行的 t6_job / t8_job / t9_job 作业总数 < 上限（默认 12），否则每 60 秒再看一次。
T8 的调度器并行数固定，T9 只占用它们空出来的名额；T8 跑完后 T9 自然用上更多核。
完成标志已存在就跳过，所以中断后重跑同一条命令即可续跑。
用法：nohup python3 code/t9_sched.py <日志目录> <上限> <作业文件>... > <调度日志> 2>&1 &"""
import os, sys, subprocess, time

ENV = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")


def running():
    """数正在运行的作业进程（按 python3 code/t6_job.py 等命令行计数；只读 /proc，不动任何进程）。"""
    n = 0
    for p in os.listdir("/proc"):
        if not p.isdigit(): continue
        try: a = open(f"/proc/{p}/cmdline", "rb").read().split(b"\0")
        except OSError: continue
        if len(a) >= 2 and a[0].endswith(b"python3") and any(a[1].endswith(x) for x in (b"t6_job.py", b"t8_job.py", b"t9_job.py")): n += 1
    return n


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
            while running() >= cap: time.sleep(60)
            log = open(os.path.join(logdir, f"{k:05d}.log"), "w"); log.write(cmd + "\n"); log.flush()
            p = subprocess.Popen(["/usr/bin/time", "-v", "timeout", tlim, "bash", "-c", cmd], stdout=log, stderr=subprocess.STDOUT, env=ENV)
            procs.append(p)
            print(time.strftime("%m-%d %H:%M:%S"), "启动", k, "PID", p.pid, cmd.split("&& ")[-1], flush=True)
            time.sleep(20)     # 等子进程出现在 /proc 里再数下一次
    for p in procs: p.wait()
    print(time.strftime("%m-%d %H:%M:%S"), "全部结束", flush=True)


if __name__ == "__main__":
    main()
