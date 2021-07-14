from group import *


def run(upto=6):
    hcn = HyperCubeGraph(2)
    hcn.explore()
    for dim in range(2, upto):
        hcn = hcn.extend()
        

if __name__ == "__main__":
    run()
