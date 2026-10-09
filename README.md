# NEAP_Project
Hand-in for master seminar (Master Bioinformatics)

## Setup:

**R:**

Please execute this on the cluster to ensure usage of checked compiliers and compatible R version:
```
mkdir -p ~/.R
cp ~/.R/Makevars ~/.R/Makevars.gcc14.backup

GCC13="$(command -v gcc-13)"
GXX13="$(command -v g++-13)"
GFORTRAN13="$(command -v gfortran-13)"

GFORTRAN_LIBDIR="$(
  dirname "$("$GFORTRAN13" -print-file-name=libgfortran.so)"
)"

cat > ~/.R/Makevars <<EOF
CC = $GCC13
CXX = $GXX13
CXX11 = $GXX13
CXX14 = $GXX13
CXX17 = $GXX13
CXX20 = $GXX13

FC = $GFORTRAN13

CXX11STD = -std=gnu++11
CXX14STD = -std=gnu++14
CXX17STD = -std=gnu++17
CXX20STD = -std=gnu++20

FLIBS = -L$GFORTRAN_LIBDIR -lgfortran -lquadmath -lm
EOF

export LIBRARY_PATH="$GFORTRAN_LIBDIR${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="$GFORTRAN_LIBDIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export USE_BUNDLED_LIBUV=1
export R_LIBS_USER="$HOME/R/gcc13/%p-library/%v"

mkdir -p "$HOME/R/gcc13/x86_64-suse-linux-gnu-library/4.3"
```
Then start an R session and install the required packages:
```
install.packages(c("stringr", "tibble", "BiocManager", "readxl"))
BiocManager::install(c("affy", "affyPLM", "oligo", "limma", GEOquery", "ArrayExpress", "AgiMicroRna", "illuminaio", "multiMiR"))
BiocManager::install(c(
  "hgu133plus2.db",             
  "hgu133a2.db",                
  "hgu133a.db",                  
  "rat2302.db",    
  "htmg430a.db",
  "mouse4302.db",                   
  "hugene10sttranscriptcluster.db",
  "hugene20sttranscriptcluster.db", 
  "hugene21sttranscriptcluster.db",
  "mta10transcriptcluster.db",
  "mogene10sttranscriptcluster.db", 
  "ragene10sttranscriptcluster.db",
  "illuminaHumanv4.db",          
  "illuminaHumanv3.db",          
  "illuminaMousev2.db"     
), ask = FALSE, update = FALSE)
```
Furthermore ensure that https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16384 and https://www.affymetrix.com/api/downloads/lf/ivt/miRNA-4_0/ are downloaded and the paths are set correctly in `annotate_data` (microarray_functions.R)

**python:**

Create an environment with `python>3.12`. Afterwards, install the required packages:

```
pip install pandas seaborn tqdm ipykernel openpyxl
```

**Further installation notes:**

Files starting wiht `FuzzSys` are expected to be run from within a local clone of this repo: https://github.com/zimmerlab/FuzzSys

**File descriptions**
1. `analyse_views_miRNA`: Determine significantly changed (logfc D-view cat --, -, +, ++ and padj significant) miRNAs in available views
2. `analyse_views_mRNA`: See `analyse_views_miRNA`. Same procedure, only for mRNA
3. `analysis_PRISM_datasets`: visualization of atherosclerosis related datasets
4. `create_cond_matrix`: Recovery of test and reference conditions after DEA
5. `FuzzSys_miRNA`: Runs FuzzSys on miRNA datasets. Creates D-views (logfc and padj) and E-views(normalized)
6. `FuzzSys_mRNA`: See `FuzzSys_miRNA`. Same procedure, only for mRNA
7. `FuzzSys_miRNA_mewis`: See `FuzzSys_miRNA`. Same procedure, only for datasets of Mr. Mewis
8. `fuzzy_qc`: combines all miRNA and mRNA datasets for later Fuzzyfication to determine deviations of samples.
8. `microarray_functions`: Functions for execution of mircoarray preprocessing pipeline
9. `microarray_views_miRNA`: Preprocessing of miRNA microarrays and subsequent DEA
10. `microarray_views_mRNA`: See `microarray_views_miRNA`. Same procedure, only for mRNA
11. `microarray_views_miRNA_mewis`: See `microarray_views_miRNA`. Same procedure, only for datasets of Mr. Mewis
12. `miRNA_targets`: Target gene discovery for dysregulated mouse and human miRNAs
13. `view_analysis`: Analysis of miRNA and mRNA dysregulation in atherosclerotic contexts in ECs
14. `comb_code`: code executed to combine GSE137582 with GSE205119

The files are expected to be executed in the following order:
1. `microarray_views_miRNA` + `microarray_views_mRNA`
2. `FuzzSys_miRNA` + `FuzzSys_mRNA`
3. `analyse_views_miRNA` + `analyse_views_mRNA`
4. `create_cond_matrix`
5. `view_analysis` + `miRNA_targets`

# Documentation microarray_functions

**Main Functions:**

download_microarray_data(accession, out_dir, force=F)

- Accession: Accession number to download
- Out_dir: output directory for the download
- Force: If true, redownloads the file even when present

Download function for datasets from ArrayExpress and GEO. Checks for “GSE” and “E-MTAB” prefix


load_data(data_dir, metadata_file_path, sep=”.”, other.columns=”Detection”)

- Data_dir: directory with dataset
- Metadata_file_path: Path to metadata file
- Sep: separator used for separating names when parsing non-normalized Illumina files
- Other.columns: Prefix of other columns to parse for non-normalized Illumina files

Returns r object with raw data. Available microarray see Appendix 1.


background_correction(raw_data)

- Raw_data: R objects of class: AffyBatch, ExpressionSet, ExpressionFeatureSet, EList, EListRaw, RGList, uRNAList

Applies background correction to microarrays. Illumina microarrays of type Elist/EListRaw  limma::nec; uListRaw  AgiMicroRna::tgsNormalization; Agilent microarrays of type RGList/EListRaw/Elist  limma::backgroundCorrect + limma:: normalizeBetweenArrays; AffyBatch  affy::rma; ExpressionFeatureSet/GeneFeatureSet  oligo::rma


normalization(raw_data)

- Raw_data: R objects of class: AffyBatch, ExpressionSet, ExpressionFeatureSet, EList, EListRaw, RGList, uRNAList

Applies normalization to microarrays. Illumina microarrays of type Elist/EListRaw  limma::neqc; uListRaw  already corrected; Agilent microarrays of type RGList/EListRaw/Elist  limma::backgroundCorrect; AffyBatch  affy::rma; ExpressionFeatureSet/GeneFeatureSet  oligo::rma


annotate_data(data, install_missing=F)

- Data: R objects of class: AffyBatch, ExpressionSet, ExpressionFeatureSet, EList, EListRaw, RGList
- Install_missing: Whether to install missing annotation packages

Annotates microarray data with the appropriate annotation package. Available microarray see Appendix 1.


save_expression_matrix(data, save.dir, save.name)

- Data: R object returned by load_data, background_correction, normalization
- Save.dir: Directory to save expression values in
- Save.name: Name of file

Saves the expression values of the provided data as .tsv. Rownames are set as column “Symbol” from the gene information of the object; Columnames are set as the “sample” column of the sample level data.


clean_genes(data, symbol_col =”Symbol”, remove_controls=T, verbose=T)

- Data: R object of class EList, EListRaw, ExpressionSet, ExpressionFeatureSet. The object is assumed to have a column with the gene symbol saved in the respective feature-annotating matrix
- Symbol_col: column in feature-annotating matrix with gene symbol
- Remove_controls: remove features/genes which are probably controls
- Verbose: send updates through the terminal

This function cleans the features/genes of the microarray data of not annotated probes and control probes. Control probes are determined by checking the gene metadata for common registries of control probes, like “controlType” or “status”. Here control probes are determined through string matching. For Affymetrix microarray data the probe names are screened for probes starting with “AFFX-“. For Agilent microarrays the probe names are screened for “DarkCorner”, “SCorner", “NC[0-9]+_”. The control genes are removed if remove_controls is true. Probes with missing gene names are always removed.


collapse_duplicate_genes(data, symbol_col = “Symbol”, method, verbose=T)

- Data: R object of class EList, EListRaw, ExpressionSet, ExpressionFeatureSet. The object is assumed to have a column with the gene symbol saved in the respective feature-annotating matrix
- Symbol_col: column in feature-annotating matrix with gene symbol
- Method: one of “mean”, “median”, “max_mean”. Collapses duplicate probes and returns either the mean, median or max_mean of the probes. In case of max_mean the probe with the largest mean across samples is taken.
- Verbose: send updates through the terminal

Collapses probes detecting the same gene. Available methods of collapsing the values are mean, median and max_mean


combine_chip_datasets(…, dataset_names=NULL)

- …: datasets to combine
- Dataset_names: vector of dataset names to create a unique dataset identifier for accession

Combines datasets with each other. Chips need to be identical to avoid misssense probe name collisions. Cross-Chip microarray comparison is not recommended


run_dea(data, group_col, group_levels=NULL, contrast_str=NULL, covariate_cols=NULL, p_adjust_method=”BH”, save.view=F, save.dir=NULL)

- Data: processed microarray data; expects one of ExpressionSet, ExpressionFeatureSet, EListRaw, EList, uRNAList
- Group_col: column containing the contrast conditions
- Group_levels: specify the levels of the group_col. If NULL then the levels are set automatically
- Contrast_str: contrast string of form “<cond2>-<cond1>”. Cond1 and cond2 must be in group_col. If contrast_str=NULL and group_col contains exactly two conditions, then these are taken as contrast_str.
- P_adjust_method: method for multiple testing correction of DEA. Defaults to Benjamin-Hochberg correction.
- Save.view: If true, saves the DEA table
- Save.dir: Output directory for the dea. File is saved in the output directory as “<accession_nr>_<contrast_str>.tsv”. contrast_str is adapted to a valid file name string when saving. If save.view=TRUE, then required.

Runs DEA on normalized microarray data using limma.


run_pca(data, n_top=1000, center=T, scale.=F)

- Data: R Object of class ExpressionSet, ExpressionFeatureSet, EList or EListRaw
- N_top: number of most highly variable features to use for PCA
- Center: center the features
- Scale: scale the features by the standard deviation

Calculates the PCA representation of the given dataset




**Helper Functions**

read_metadata(metadata_file_path, sheet=1)
- Metadata_file_path: Path to excel file containing metadata
- Sheet: sheet number of excel file which contains the metadata

Read in excel file containing metadata.

Assumptions: The unique identifier is saved in the 4th column of the excel.


match_metadata_to_files(files, meta)
- Files: list of file paths
- Meta: metadata dataframe

Match file names to metadata rows


unpack_archives(dataset_dir)

- Dataset_dir: directory 

Unpack downloaded archives (.zip, .tar and .gz files)


Unpack_geo_seq_archives(dataset_dir, remove_archives=F, force=F, max_depth=10, quiet=F)

- Dataset_dir: directory with files downloaded from GEO
- Remove_archives: Delete archives from GEO after download
- Force: Overwrites previously unpacked data
- Max_depth: maximal depth of iterative unpacking
- Quiet: If false, no messages are printed

Unpack datasets downloaded from GEO. Iteratively checks subdirectories for further compressed datasets. Is sensitive to 10x sequencing data.


create_agilent_targets(data_dir)
- Data_dir: data directory to save the targets.txt in

Formats the targets.txt files for Agilent microarrays. The file is saved as data_dir / targets.txt


make_annotation_df(meta_df, file_paths)


make_targets(meta_df, file_paths, 	data_dir, col=”FileName”)
- Meta_df: metadata dataframe
- File_paths: list of paths to the agilent paths
- Data_dir: Data directory of the files
- Col: Name of column containing the file paths

Returns a dataframe with filenames matched to the respective samples


Read_chp_files(files,  meta)
- Files: list of files to read in
- Meta: Metadata for the respective files

Read-in .chp saved files


Get_supplier(data)

- Data: R object of one of the following classes: AffyBatch, ExpressionSet, ExpressionFeatureSet, EList, EListRaw, RGList, uRNAList

Returns the chip types if saved in the metadata. Assumes that the sample metadata contains “chip”


read_ilmn_matrix_safe(file_path, probe_id=NULL, expr=”AVG_signal”, sep=”_”, other.columns=”Detection”)
- File_path: path to .ilmn file
- Probe_id: Name of column containing the probe id. If NULL, then one of "ID_REF", "PROBE_ID", "ProbeID", "Probe" is taken, if any is present.
- Expr: column containing expression values
- Sep: separator of sample id and column function
- Other.columns: Prefix of other columns to parse

Returns an EListRaw function with expression values in E, gene information in genes, other information in other and sample level information in targets


Detection_pval(raw_data)
- Raw_data: raw_data returned by load_data

Determines the name of column containing pvalues for measurements


Filter_illumina(x, detection_name=”Detection PValue”)
- X: R Object of class EList or ElistRaw containing data from an Illumina chip
- Detection_name: name of matrix containing detection p-values.

Filters Illumina data for eligibly detected probes and removes not properly detected ones.


Pca_plot(pca_result, color_by, shape_by=NULL, pc_x=1, pc_y=2, point_size=1.5)
- Pca_result: calculated pca from run_pca
- Color_by: column of metadata to color the different samples
- Shape_by: column of metadata used for shapes
- Pca_x: PC component to use for x-axis
- Pca_y: PC component to use for y-axis
- Point_size: size of points

Plots the PCA


Urna_to_elist(x)
- X: dataset of class uRNAList

Converts uRNAList to Elist object

**Appendix 1: list of implemented microarray chips**

- Affymetrix Human Genome U133 Plus 2.0 Array
- Affymetrix Human Genome U133A 2.0 Array
- Affymetrix HT Human Genome U133A
- Affymetrix HT Mouse Genome 430A Array
- Affymetrix Rat Genome 230A
- Affymetrix Mouse Genome 430 2.0 Array
- Affymetrix Mouse Genome 430A 2.0 Array
- Affymetrix Mouse Transcriptome Array 1.0
- Affymetrix Human Gene 2.1 ST
- Affymetrix Human Gene 2.0 ST
- Affymetrix Human Gene 1.0 ST
- Affymetrix Mouse Gene 1.0 ST Array
- Affymetrix Mouse Exon-Junction Array
- Affymetrix Rat Gene 1.0 ST
- Affymetrix miRNA 3.0
- Affymetrix miRNA 4.0
- Illumina HumanHT-12 v4.0
- Illumina HumanHT-12 v3.0
- Illumina MouseRef-8 v2.0
- Illumina HumanRef-8 v2.0
- Illumina Human v2 MicroRNA
- Illumina Mouse v2 MicroRNA
- Agilent-084555 026652QM_RCUG_HomoSapiens
- Agilent-074809 SurePrint G3 Mouse GE v2
- Agilent-072363 SurePrint G3 Human GE v3
- Agilent-028005 SurePrint G3 Mouse GE
- Agilent miRNA microarray
- Agilent-019159 Rat miRNA
- Agilent-021827 Human miRNA Microarray G4470C
- Agilent-021827 Human miRNA Microarray v2
- Agilent-021827 Human miRNA Microarray v3
- Agilent-035430 mouse miRNA array
- Agilent-050340 Custom Rat miRNA Microarray
- Agilent-070155 Mouse miRNA Microarray
- Agilent-070156 Human miRNA
