from pathlib import Path
import subprocess
from collections import defaultdict

import pandas as pd
import numpy as np

def main():
    base_path = Path("/usr/local/storage/data_microarray")

    base_path_data_E = base_path / "261004_normalized_mewis"
    base_path_data_D = base_path / "261004_dea_mewis"
    out_E_view = base_path / "E_view_miRNA_mewis"
    out_D_view = base_path / "D_view_miRNA_mewis"

    out_E_view.mkdir(parents=True, exist_ok=True)
    out_D_view.mkdir(parents=True, exist_ok=True)

    problems = defaultdict(list)

    exemptions = ["GSE89858"]

    for file in base_path_data_E.glob("*.tsv"):
        print()
        out_E_sample_concept = base_path  / "fuzzy_concepts_miRNA_mewis" / f"{file.name.split(".")[0]}_E_sample"
        out_E_sample_result = base_path / "fuzzy_res_miRNA_mewis" / f"{file.name.split(".")[0]}_E_sample"

        out_E_sample_concept.mkdir(parents=True, exist_ok=True)
        out_E_sample_result.mkdir(parents=True, exist_ok=True)

        # convert raw files to expected format (shift in first line of files, which work well in pandas)
        tmp = pd.read_csv(file, sep="\t", index_col=0)
        tmp.set_index("GeneName" if any([x in file.name for x in exemptions]) else "Symbol", drop=True, inplace=True)
        tmp.index.rename(None, inplace=True)
        tmp = tmp.loc[:, tmp.columns.str.startswith("GSM")]
        print("saving tmp file")
        tmp.to_csv(out_E_view / file.name, sep="\t")

        print("running {} on samples".format(file.name))

        try:
            subprocess.run([
                    "python", "main_concepts.py",
                    "--mtx", str(out_E_view / file.name),
                    "--config", "/home/f/flor/fuzzy_concept_E_sample_wise.json",
                    "--threads", "8",
                    "--output", str(out_E_sample_concept)
                ], check=True
            )

            subprocess.run([
                "python", "main_fuzzify.py",
                "--mtx", str(out_E_view / file.name),
                "--concept", str(out_E_sample_concept / "concepts_detailed.json"),
                "--config", "/home/f/flor/fuzzy_results_E_sample_wise.json",
                "--threads", "8",
                "--output", str(out_E_sample_result)
            ], check=True)
        except subprocess.CalledProcessError as e:
            problems[file.name].append(e)


    for file in base_path_data_D.glob("*.tsv"):
        print(f"working on {file}")
        out_D_logfc_concept = base_path / "fuzzy_concepts_miRNA_mewis" / f"{file.name}_D_logfc"
        out_D_padj_concept = base_path / "fuzzy_concepts_miRNA_mewis" / f"{file.name}_D_padj"
        out_D_logfc_result = base_path / "fuzzy_res_miRNA_mewis" / f"{file.name}_D_logfc"
        out_D_padj_result = base_path / "fuzzy_res_miRNA_mewis" / f"{file.name}_D_padj"

        out_D_logfc_concept.mkdir(parents=True, exist_ok=True)
        out_D_padj_concept.mkdir(parents=True, exist_ok=True)
        out_D_logfc_result.mkdir(parents=True, exist_ok=True)
        out_D_padj_result.mkdir(parents=True, exist_ok=True)

        data_path_logfc = out_D_view / f"{".".join(file.name.split(".")[:-1])}_logfc.tsv"
        data_path_padj = out_D_view / f"{".".join(file.name.split(".")[:-1])}_pval.tsv"

        dea = pd.read_csv(file, sep="\t")
        dea['-log10(p.adj)'] = -np.log10(dea['adj.P.Val'].values)
        dea_pval = pd.DataFrame(
            dea["-log10(p.adj)"].values, 
            columns=["-log10(p.adj)"], 
            index = dea["GeneName" if any([x in file.name for x in exemptions]) else "Symbol"].values
            )
        dea_logfc = pd.DataFrame(
            dea["logFC"].values,
            columns=["logFC"], 
            index = dea["GeneName" if any([x in file.name for x in exemptions]) else "Symbol"].values
            )

        dea_pval.to_csv(data_path_padj, sep="\t")
        dea_logfc.to_csv(data_path_logfc, sep="\t")
        print("running {} for logfc".format(file.name))

        try:
            subprocess.run([
                    "python", "main_concepts.py",
                    "--mtx", str(data_path_logfc),
                    "--config", "/home/f/flor/fuzzy_concept_D_logfc.json",
                    "--threads", "8",
                    "--output", str(out_D_logfc_concept)
                ], check=True
            )

            subprocess.run([
                "python", "main_fuzzify.py",
                "--mtx", str(data_path_logfc),
                "--concept", str(out_D_logfc_concept / "concepts_detailed.json"),
                "--config", "/home/f/flor/fuzzy_results_D.json",
                "--threads", "8",
                "--output", str(out_D_logfc_result)
            ], check=True)
        except subprocess.CalledProcessError as e:
            problems[file.name].append(e)

        print("running {} for p.adj".format(file.name))

        try:
            subprocess.run([
                    "python", "main_concepts.py",
                    "--mtx", str(data_path_padj),
                    "--config", "/home/f/flor/fuzzy_concept_D_padj.json",
                    "--threads", "8",
                    "--output", str(out_D_padj_concept)
                ], check=True
            )

            subprocess.run([
                "python", "main_fuzzify.py",
                "--mtx", str(data_path_padj),
                "--concept", str(out_D_padj_concept / "concepts_detailed.json"),
                "--config", "/home/f/flor/fuzzy_results_D.json",
                "--threads", "8",
                "--output", str(out_D_padj_result)
            ], check=True)
        except subprocess.CalledProcessError as e:
            problems[file.name].append(e)

    if problems:
        print(problems)
    else:
        print("Done")

if __name__ == "__main__":
    main()