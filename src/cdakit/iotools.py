import subprocess

import pandas as pd


def to_format_table(df: pd.DataFrame, index_label="index"):
    csv_str = df.to_csv(
        None, sep=" ", float_format="%.6f", na_rep="NaN", index_label=index_label
    )
    fmt_proc = subprocess.Popen(
        ['column', '-t'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    table_str, _ = fmt_proc.communicate(csv_str)
    return table_str


def read_format_table(ftable):
    df = pd.read_table(ftable, sep=r"\s+", index_col=None)
    if "index" in df.columns:
        df = df.set_index("index")
    return df

