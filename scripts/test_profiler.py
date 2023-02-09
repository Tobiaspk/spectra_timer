from time import sleep
from spectra.utils.profiler import profile

@profile
def test():
    for i in range(10):
        sleep(0.1)
    s = 0
    for i in range(100000):
        s += 1
    return s

if __name__ == '__main__':
    test()