import timeit
from functools import partial
from src.cluster import profile_cluster, PROPERTY_COMBINATIONS
from src.io.datasets import covid19_repertoire
from src.io.output import path_in_results


REPEAT = 1
NUMBER = 1


def main():
    data = covid19_repertoire()
    benchmark(data)


def time_function(func, *args):
    return min(timeit.Timer(partial(func, *args)).repeat(REPEAT, NUMBER))


def benchmark(data):
    results = f"""
For each combination of properties, a clustering is made {NUMBER} times and timed.
This process is repeated {REPEAT} times and the fastest time is recorded.

The following times are for clustering a dataset of {len(data)} TCR's {NUMBER} times

"""
    for combo in PROPERTY_COMBINATIONS:
        time = time_function(profile_cluster, data, combo)
        name = ' + '.join(combo)
        results += f'{name}: {time:.2f}s\n'
    with open(path_in_results(f'benchmarks/profile_method_{REPEAT}x{NUMBER}.txt'), 'w') as f:
        f.write(results)


if __name__ == '__main__':
    main()

