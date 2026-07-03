import numpy as np

for name in [
    # [旧·服务器路径]
    # "/mnt/sda/gws_1020251255/data/processed/all_proteins_sci_site_scores_mean.npy",
    # "/mnt/sda/gws_1020251255/data/processed/all_proteins_sci_site_scores_top20.npy",
    # "/mnt/sda/gws_1020251255/data/processed/all_proteins_sci_site_scores_top50.npy"
    # [新·本地路径]
    "D:/文件/工作室/website/data/processed/all_proteins_sci_site_scores_mean.npy",
    "D:/文件/工作室/website/data/processed/all_proteins_sci_site_scores_top20.npy",
    "D:/文件/工作室/website/data/processed/all_proteins_sci_site_scores_top50.npy"
]:
    x=np.load(name)

    print("\n",name)
    print("shape:",x.shape)
    print("mean :",x.mean())
    print("std  :",x.std())
    print("min  :",x.min())
    print("max  :",x.max())
    print("p1   :",np.percentile(x,1))
    print("p5   :",np.percentile(x,5))
    print("p50  :",np.percentile(x,50))
    print("p95  :",np.percentile(x,95))
    print("p99  :",np.percentile(x,99))