"""通用作业调度：作业文件每行 "<完成标志文件>\t<超时秒>\t<命令>"。
- 完成标志文件已存在就跳过，所以中断后重跑同一条命令即可续跑；
- 每个作业由它自己写自己的结果文件（单写入者），调度器只写 <日志目录>/<序号>.log；
- 每个作业单独限时（timeout），超时记为未完成，不重试；
- 用 /usr/bin/time -v 记录峰值内存，便于判断并行度。
用法：nohup python3 code/run_jobs.py <作业文件> <并行数> <日志目录> > <调度日志> 2>&1 &"""
import os, sys, subprocess, time, json
from concurrent.futures import ThreadPoolExecutor, as_completed

ENV = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")


def run_one(a):
    """运行一个作业；返回 (序号, 返回码, 用时)。"""
    n, done, tlim, cmd, logdir = a
    if os.path.exists(done): return n, "skip", 0.0
    t0 = time.time()
    log = os.path.join(logdir, f"{n:05d}.log")
    with open(log, "w") as f:
        f.write(cmd + "\n"); f.flush()
        r = subprocess.run(["/usr/bin/time", "-v", "timeout", str(tlim), "bash", "-c", cmd], stdout=f, stderr=subprocess.STDOUT, env=ENV)
    return n, r.returncode, time.time() - t0


def main():
    jobfile, P, logdir = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    os.makedirs(logdir, exist_ok=True)
    jobs = []
    for n, line in enumerate(open(jobfile)):
        line = line.rstrip("\n")
        if not line or line.startswith("#"): continue
        done, tlim, cmd = line.split("\t", 2)
        jobs.append((n, done, int(tlim), cmd, logdir))
    todo = [j for j in jobs if not os.path.exists(j[1])]
    print(json.dumps(dict(event="start", total=len(jobs), todo=len(todo), workers=P, pid=os.getpid(), t=time.strftime("%F %T"))), flush=True)
    with ThreadPoolExecutor(P) as ex:
        for fut in as_completed([ex.submit(run_one, j) for j in todo]):
            n, rc, sec = fut.result()
            print(json.dumps(dict(job=n, rc=rc, sec=round(sec, 1), t=time.strftime("%F %T"))), flush=True)
    print(json.dumps(dict(event="end", t=time.strftime("%F %T"))), flush=True)


if __name__ == "__main__":
    main()
