#!/usr/env/python3

import matplotlib.pyplot as plt
import pandas as pds


def main():

    bucket_sizes=(20,40,60,80,100)
    datas = [pds.read_table("test_{}.txt".format(bucket_size),
                           header=None, names=["nppc","ull{}".format(bucket_size), "phare{}".format(bucket_size)],
                           index_col=0, sep=" ") for bucket_size in bucket_sizes]
    fig,(ax1,ax2) = plt.subplots(ncols=2, figsize=(10,5))
    for nppc,data in zip(bucket_sizes,datas):
        for col in data.columns:
            if "phare" in col:
                continue
            else:
                data[col].plot(ax=ax1, label=f"{nppc}")

    for nppc,data in zip(bucket_sizes, datas):
        for col in data.columns:
            if "phare" in col:
                print(col)
                data[col].plot(ax=ax2, label=f"{nppc}", ls="-")
    for ax in (ax1,ax2):
        ax.legend(ncol=4)
        for b in bucket_sizes:
            ax.axvline(b, ls="--", color="k", alpha=0.2)
    ax1.set_title("ull")
    ax2.set_title("phare")

    fig.savefig("perf.png")


if __name__ == "__main__":
    main()
