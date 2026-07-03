import numpy as np
import pandas as pd

def check_data():
    # [旧·服务器路径]
    # mat = np.memmap(
    #     "/mnt/sda/gws_1020251255/data/processed/all_proteins_sci_site_matrices.dat",
    #     dtype=np.float32, mode="r", shape=(95142, 33, 33)
    # )
    # df = pd.read_csv("/mnt/sda/gws_1020251255/data/processed/all_proteins_site_metadata.csv")
    # [新·本地路径]
    mat = np.memmap(
        "D:/文件/工作室/website/data/processed/all_proteins_sci_site_matrices.dat",
        dtype=np.float32,
        mode="r",
        shape=(95142, 33, 33)
    )
    print("Matrix shape:", mat.shape)

    df = pd.read_csv("D:/文件/工作室/website/data/processed/all_proteins_site_metadata.csv")
    print("Metadata rows:", len(df))
    print("Metadata columns:", df.columns.tolist())

if __name__ == "__main__":
    check_data()