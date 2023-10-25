import multiprocessing
import subprocess
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--omp_places", default=None, choices=("cores",))
    parser.add_argument("--omp_schedule", default=None)
    parser.add_argument("--plot", help="file to save scaling plot to")
    args, remaining_args = parser.parse_known_args()

    env = dict(os.environ)
    if args.omp_schedule:
        env.update(OMP_SCHEDULE=args.omp_schedule)
    if args.omp_places:
        env.update(OMP_PLACES=args.omp_places)
    result = []
    ref_ranges = None
    for ncpus in range(1, 1 + multiprocessing.cpu_count()):
        print(f"Testing with {ncpus} cores...", end="", flush=True)
        p = subprocess.run(
            remaining_args,
            env=env | {"OMP_NUM_THREADS": str(ncpus)},
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding="ascii",
        )
        nsec = None
        read_ranges = False
        ranges = {}
        for line in p.stdout.split("\n"):
            if line.startswith("Wall time:"):
                nsec = float(line[10:-2].strip())
            if read_ranges:
                read_ranges = line.startswith("  ")
                if read_ranges:
                    name, minval, maxval = line.strip().split(" ")
                    ranges[name] = (float(minval), float(maxval))
            else:
                read_ranges = line.startswith("Final variable ranges:") or ()
        print(f" {nsec} s")
        if ref_ranges is None:
            ref_ranges = ranges
        else:
            for name, (minref, maxref) in ref_ranges.items():
                minval, maxval = ranges[name]
                if minval != minref or maxval != maxref:
                    print(
                        f"  {name} mismatch: {minval} - {maxval} vs {minref} - {maxref}"
                    )
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
        ax.set_ylim(0.0, 100.0)
        fig.savefig(args.plot, dpi=150)
