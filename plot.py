#!/usr/env/python3

import matplotlib.pyplot as plt
import pandas as pds


def main():

    bucket_sizes=(20,40,60,80,100)
    datas = [pds.read_table("test_{}.txt".format(bucket_size),
                           header=None, names=["nppc","ull{}".format(bucket_size), "phare{}".format(bucket_size)],
                           index_col=0, sep=" ") for bucket_size in bucket_sizes]
    fig,ax = plt.subplots()
    for nppc,data in zip(bucket_sizes,datas):
        data.plot(ax=ax)
    fig.savefig("perf.png")


if __name__ == "__main__":
    main()
