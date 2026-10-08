import pandas as pd
from pathlib import Path
from tqdm import tqdm

base_path = Path("/usr/local/storage/data_microarray")
out_path = base_path / "qc"

## miRNA
data = []
for gse in (base_path / "260927_normalized_miRNA").glob("*"):
    tmp = pd.read_csv(gse, sep="\t", index_col=0)
    tmp.set_index("GeneName" if gse.name in ["GSE137580_norm_data.tsv", "GSE137582_norm_data.tsv", "GSE205119_norm_data.tsv", "comb_GSE137580_GSE205119_norm_data.tsv"] else "Symbol", drop=True, inplace=True)
    tmp.index.rename(None, inplace=True)
    tmp = tmp.loc[:, tmp.columns.str.startswith("GSM")]
    data.append(tmp)
data_comb = data[0]
for dat in tqdm(data[1:]):
    data_comb = data_comb.merge(dat, how="outer", left_index=True, right_index=True)
data_comb.to_csv(out_path / "miRNA_samples_norm.tsv", sep="\t")

data = []
for gse in (base_path / "260927_normalized_mRNA").glob("*"):
    tmp = pd.read_csv(gse, sep="\t", index_col=0)
    tmp.set_index("Symbol", drop=True, inplace=True)
    tmp.index.rename(None, inplace=True)
    tmp = tmp.loc[:, tmp.columns.str.startswith("GSM")]
    data.append(tmp)
data_comb = data[0]
for dat in tqdm(data[1:]):
    data_comb = data_comb.merge(dat, how="outer", left_index=True, right_index=True)
data_comb.to_csv(out_path / "mRNA_samples_norm.tsv", sep="\t")