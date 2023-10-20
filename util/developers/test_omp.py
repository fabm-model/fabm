import multiprocessing
import subprocess
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--omp_places", default=None, choices=("cores",))
    parser.add_argument("--plot", help="file to save scaling plot to")
    args, remaining_args = parser.parse_known_args()

    env = dict(os.environ)
    if args.omp_places:
        env.update(OMP_PLACES="cores")
    result = []
    for ncpus in range(1, 1 + multiprocessing.cpu_count()):
        print(f"Testing with {ncpus} cores...", end="", flush=True)
        p = subprocess.run(
            remaining_args,
            env=env | {"OMP_NUM_THREADS": str(ncpus)},
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding="ascii",
        )
        last_line = p.stdout.rsplit("\n", 2)[1].strip()
        assert last_line.startswith("Wall time:")
        nsec = float(last_line[10:-2].strip())
        print(f" {nsec} s")
        result.append([ncpus, nsec])

    if args.plot:
        from matplotlib import pyplot
        import numpy as np

        fig, ax = pyplot.subplots()
        result = np.array(result)
        ax.plot(result[:, 0], 100 * result[:, 1] / result[0, 1])
        ax.set_xlabel("number of cores")
        ax.set_ylabel("wall time relative to serial run (%)")
        ax.grid()
        fig.savefig(args.plot, dpi=150)
