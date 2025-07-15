---

# 🧬 DHS-to-cCRE Pipeline

This pipeline processes open chromatin data (e.g., DNase-seq or ATAC-seq) to identify high-confidence **candidate cis-regulatory elements (cCREs)** by integrating chromatin accessibility and histone modification signals (e.g., H3K27ac).

It performs peak calling, signal enrichment scoring, DHS filtering, normalization, and final cCRE annotation using R-based quantile thresholds.

---

## 📁 Working Directory Structure

```
├── bams/
│   ├── *out_.bam
│   ├── *out_.bam.bai
├── signal_files/
│   └── S3norm_rc_bedgraph/
│       └── *_S3.bw
├── H3K27ac_data/
│   ├── *.bw
│   ├── *.bed
```

---

## 🧩 Pipeline Overview

1. **Peak Calling**
   Uses MACS3 on BAM files to generate DHS regions.

2. **Signal Enrichment**
   Calculates QPoisson enrichment scores across DHSs.

3. **Filtering DHSs**
   Retains only enriched DHSs (length >150 bp and signal > threshold).

4. **BigWig Signal Averaging**
   Uses `bigWigAverageOverBed` to extract signal intensities from ATAC/H3K27ac data.

5. **Quantile Filtering (R)**
   Identifies high-signal DHSs based on quantile thresholds in a context-specific manner.

6. **Final Output**
   Annotated BED files containing candidate cCREs per context.

---

## 📦 Requirements

Install the following dependencies:

* **Bash + GNU Coreutils**
* [MACS3](https://github.com/macs3-project/MACS)
* [BEDOPS](https://bedops.readthedocs.io/)
* [bigWigAverageOverBed](https://hgdownload.soe.ucsc.edu/admin/exe/)
* `R` with the following packages:

  * `data.table`
  * `dplyr`
* [`GNU parallel`](https://www.gnu.org/software/parallel/)

You can install R dependencies via:

```r
install.packages(c("data.table", "dplyr"))
```

---

## 🚀 Quickstart

```bash
# Step 1: Set your working directory
export WORKING_DIRECTORY="/your/working/directory"
export METADATA_FILE="data/metadata.tsv"

# Step 2: Call peaks with MACS3
bash scripts/call_macs3.sh

# Step 3: Process QPoisson signal and DHSs
bash scripts/process_signal.sh

# Step 4: Compute normalized signals (bigWigAverageOverBed required)
bash scripts/compute_qpoisson.sh

# Step 5: Run quantile filtering in R
Rscript scripts/compute_quantiles.R
```

---

## 📄 Input Files

* `metadata.tsv`: Sample info with columns:

  ```
  sample_id    replicate    context    tissue
  ```

* `.bam` files: Raw aligned reads.

* `.bw` files: Signal tracks for ATAC-seq and H3K27ac.

---

## 📤 Output Files

* `output/dhs/`: All peaks and filtered DHSs per sample.
* `output/qpoisson/`: DHSs with signal scores.
* `output/bigwig_signal/`: Normalized average signals per site.
* `output/cCREs/`: Final cCRE BED files per context.

---

## 📊 Final Output Format

Final annotated BED file (`output/cCREs/`) contains:

```
chr    start    end    name    score    strand
```

Where `name` is the unique cCRE ID and `score` is the averaged normalized signal.

---

## 📌 Notes

* You can adjust the `minP`, `signal threshold`, or `quantile` cutoffs in the R scripts as needed.
* Designed for DNase-seq or ATAC-seq input data but flexible to other open chromatin assays.

---



