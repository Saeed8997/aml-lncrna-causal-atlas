# Phase 0: Dataset Reconnaissance & Size/Format Audit
**Project:** Causal lncRNA–mRNA Regulatory Atlas of Acute Myeloid Leukemia Multi-Resolution Constraint-Based Discovery on Large-Scale Single-Cell Atlases with CRISPRi Functional Validation  
**Date:** October 14, 2026  
**Agent:** Computational Biology & MLOps Agent  
**Audit Policy:** Zero full-dataset downloads to local workstation. Strict metadata queries, Kaggle API reconnaissance, and HTTP/REST header checks.

---

## Executive Summary & Decision Matrix

| Dataset / Entity | Source Platform & Accession / UUID | Primary Format | Exact Size | In Kaggle Datasets? | Primary Recommendation |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **AML scAtlas (~748k cells)** | CELLxGENE `071b706a-7ea7-47a4-bddf-6457725839fc` | `.h5ad` | **6.40 GB** (6,398,033,807 B) | **No** (0 matches) | **Primary scRNA Atlas**: Direct cloud-mount to Kaggle. Load in `backed='r'` mode. |
| **AML scAtlas (Figshare)** | Figshare `27269946` | `.h5ad.gz` | **11.36 GB** (11,358,833,727 B) | **No** (0 matches) | Avoid on Kaggle (compressed size consumes 57% of `/kaggle/working`). Use CELLxGENE `.h5ad` instead. |
| **Zeng BoneMarrowMap** | CELLxGENE `f6c50495-3361-40ed-a819-fb9644396ed9` | `.h5ad` | **2.83 GB** (2,829,047,593 B) | **No** (0 matches) | Reference map (263,159 healthy hematopoiesis cells) for cell-state projection. |
| **Zeng 21 AML Cohorts** | Zenodo `18528442` (DOI: 10.5281/zenodo.18528442) | `.rds` (21x) + `.csv` | **33.5 GB** total (CSV: 863.4 MB) | **No** (0 matches) | Multi-cohort validation & benchmarking. Target individual `.rds` files per cohort. |
| **van Galen AML BM subset** | Kaggle: `thomasneff/aml-2019-healthy-bone-marrow-scrna-seq` | `.h5ad` + `.csv` | **1.72 GB** (`GSE116256_BM.h5ad`) | **Yes** | Useful for rapid sanity check, but contains only healthy BM cells from the 2019 study. |
| **lncRNA CRISPRi Screens** | NCBI GEO `GSE318789` / PMC `PMC13318546` | `.txt.gz` + `.xlsx` | **1.4 MB** (counts) + Supp. Tables | **No** (0 matches) | **Causal Ground Truth**: 7,996 lncRNAs screened in MOLM-13 across proliferation, differentiation, and venetoclax response. |

---

## 1. Kaggle Datasets Reconnaissance via API / MCP

We executed authenticated queries against the Kaggle Dataset registry for all designated keyphrases.

### Query Results Summary

| Search Query | Total Matches | Matching Dataset Slugs | Details & Formats |
| :--- | :---: | :--- | :--- |
| `"BoneMarrowMap"` | **0** | *None* | No existing public upload. |
| `"Zeng AML"` | **0** | *None* | No existing public upload. |
| `"scAML"` | **0** | *None* | No existing public upload. |
| `"AML scAtlas"` | **0** | *None* | No existing public upload. |
| `"748,679"` / `"748679"` | **0** | *None* | No existing public upload. |
| `"van Galen AML"` | **1** | `thomasneff/aml-2019-healthy-bone-marrow-scrna-seq` | Contains GSE116256 healthy bone marrow control cells. |
| `"van Galen"` | **1** | `thomasneff/aml-2019-healthy-bone-marrow-scrna-seq` | Same dataset as above. |
| `"GSE318789"` | **0** | *None* | No existing public upload. |
| `"AML CRISPRi lncRNA"` | **0** | *None* | No existing public upload. |
| `"AML CRISPRi"` | **0** | *None* | No existing public upload. |
| `"CRISPRi lncRNA"` | **0** | *None* | No existing public upload. |
| `"MOLM-13"` | **1** | `francescodeleonardis/graphpad-prism-file-manuscript-de-leonardis-et-al` | Unrelated GraphPad Prism files. |

---

## 2. Primary External Repositories (Metadata & API Inspection)

### 2.1. Zenodo Record `18528442` (Zeng et al., *Blood Cancer Discovery*, 2025)
- **Record Identifier:** `18528442`
- **DOI:** `10.5281/zenodo.18528442`
- **Associated Publication:** *Single-cell Transcriptional Atlas of Human Hematopoiesis Reveals Genetic and Hierarchy-Based Determinants of Aberrant AML Differentiation* (DOI: `10.1158/2643-3230.BCD-24-0342`)
- **Total Record Size:** **33.5 GB**
- **Biological Scale:** ~1.8 million cells (~1.2 million AML blasts + ~600,000 non-leukemic immune / healthy control cells) spanning **21 independent scRNA-seq studies**.
- **File Inventory & Structure:**
  - **Master Cell Metadata (`.csv`):** `scAML_21cohorts_Annotations_byCell.csv` — **863.4 MB**  
    Contains uniform cell-type classifications, BoneMarrowMap UMAP projections, pseudotime coordinates, and LSC signature scores for all 1.8M cells.
  - **Codebook / Dictionary (`.xlsx`):** `scAML_21cohorts_ColumnDescriptions.xlsx` — Metadata column definitions.
  - **Cohort Expression Objects (`.rds`):** **21 individual `.rds` files** (R Seurat objects with normalized/raw count matrices and BoneMarrowMap projections):
    - Major cohorts include: *Abbas et al.* (~2.1 GB), *Beneyto-Calabuig et al.* (~1.8 GB), *Triana et al.*, *van Galen et al.*, *Miles et al.*, *Benard et al.*, *Morita et al.* (individual sizes range from **~150 MB to ~2.5 GB** each, summing to **~32.6 GB**).

### 2.2. Figshare `27269946` vs. CELLxGENE Collection `071b706a-7ea7-47a4-bddf-6457725839fc` (AML scAtlas ~748k cells)
- **CELLxGENE Collection UUID:** `071b706a-7ea7-47a4-bddf-6457725839fc`
  - **Collection Name:** *Single-Cell Atlas of Acute Myeloid Leukemia (AML scAtlas)*
  - **Dataset Title:** `All Cells`
  - **Total Cells:** **748,679**
  - **File Format:** Native AnnData HDF5 (`.h5ad`)
  - **Exact File Size:** **6,398,033,807 bytes** (**5.96 GiB / 6.40 GB**)
  - **Direct Download URL:**  
    `https://datasets.cellxgene.cziscience.com/be1f3ced-6b8f-47b0-a3e2-7edc8d466abd.h5ad`
- **Figshare Article `27269946`:**
  - **Title:** *A Single-Cell Transcriptomic Atlas of Acute Myeloid Leukemia*
  - **File Name:** `AML_scAtlas_Jan2026.h5ad.gz`
  - **Exact File Size:** **11,358,833,727 bytes** (**10.58 GiB / 11.36 GB**)
  - **Direct Download URL:**  
    `https://ndownloader.figshare.com/files/63342394`

### 2.3. CELLxGENE Collection `f6c50495-3361-40ed-a819-fb9644396ed9` (Zeng BoneMarrowMap Reference)
- **Collection Name:** *A Balanced Bone Marrow Reference Map of Hematopoietic Development*
- **Dataset Title:** `A (Balanced) Bone Marrow Reference Map of Hematopoietic Development`
- **Cell Count:** **263,159** healthy hematopoietic cells
- **File Format:** `.h5ad`
- **Exact File Size:** **2,829,047,593 bytes** (**2.63 GiB / 2.83 GB**)
- **Direct Download URL:**  
  `https://datasets.cellxgene.cziscience.com/2cf4ad28-177f-4326-ab49-c4cb012c295d.h5ad`

### 2.4. NCBI GEO `GSE318789` & PMC `PMC13318546` (lncRNA Functional CRISPRi Screen)
- **Accession:** NCBI GEO `GSE318789` / BioProject `PRJNA1420621`
- **Publication:** Cozzi E, Neddermeyer A, Miliara S, et al. *Identification of long non-coding RNAs involved in leukemogenesis and venetoclax response in acute myeloid leukemia through functional CRISPR-dCas9 interference screens.* *Non-coding RNA Res.* 2026; 20: 29-43. (PMID: `42388349`, PMCID: `PMC13318546`, DOI: `10.1016/j.ncrna.2026.05.002`).
- **Screen Targets:** **7,996 lncRNAs** screened in MOLM-13 across Proliferation (58 hits), Differentiation (4 hits), and Venetoclax response (23 hits).
- **GEO Read Count Matrix URL:**  
  `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318789/suppl/GSE318789_MOLM13_All_Readcount_Data.txt.gz` (1.4 MB)
