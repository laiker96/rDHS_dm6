Got it! Here's a **single `README.md` file**, fully written in Markdown, ready to drop into your GitHub repository:

```markdown
# rDHS Processing Pipeline

A shell-based pipeline to identify and quantify **regulatory DNase hypersensitive sites (rDHSs)** using ATAC-seq or ChIP-seq data. It combines MACS3 peak calling, enrichment filtering, and signal quantification in a reproducible and parallelized workflow.

---

## 📁 Directory Structure

```

WORKING\_DIRECTORY/
├── bams/                           # Input BAM files (<ID>\_out.bam)
├── macs\_peaks/                     # Output from MACS3 and filtering
│   └── <ID>/                       # Per-sample results
├── DHSs/                           # Final rDHSs and signal-annotated output
│   └── Processed-DHSs/
├── signal\_files/
│   └── S3norm\_rc\_bedgraph/        # BigWig files (<ID>\_S3.bw)
├── METADATA.csv                   # Metadata file (sample ID in column 5)
├── run\_pipeline.sh                # Master pipeline script

````

---

## ⚙️ Requirements

- Bash ≥ 4
- [MACS3](https://github.com/macs3-project/MACS)
- [BEDTools](https://bedtools.readthedocs.io/)
- [`bigWigAverageOverBed`](http://hgdownload.soe.ucsc.edu/admin/exe/)
- GNU `parallel`
- `awk`, `sort`, `paste`, `wc`
- Python 3 (for `filter-long-double.py`)

---

## 🚀 Usage

```bash
bash run_pipeline.sh /absolute/path/to/WORKING_DIRECTORY /absolute/path/to/METADATA.csv
````

* The **working directory** contains all input/output.
* The **metadata CSV** should have sample IDs in the **5th column** (no header assumed by default).

---

## 🧪 Pipeline Overview

### 1. Peak Calling with MACS3

Calls peaks per sample using `--shift -75 --extsize 150` to accommodate ATAC-seq data.

### 2. qPoisson Enrichment Scores

Computes enrichment scores using `macs3 bdgcmp -m qpois`.

### 3. Signal Processing

Converts qPoisson (-log10) values to p-values, filters by intersection with MACS3 peaks.

### 4. Iterative DHS Calling

Applies a custom Python script (`filter-long-double.py`) to filter by enrichment thresholds (from `1e-2` to `1e-325`), merges peaks, and selects rDHSs up to 400 bp.

Outputs:

* `<ID>.DHSs.bed`: rDHSs
* `<ID>.Excluded.bed`: filtered-out regions

### 5. Signal Quantification

Uses UCSC’s `bigWigAverageOverBed` to quantify the normalized signal over each rDHS.

Output:

* `DHSs/Processed-DHSs/output.<ID>`: BED file with per-rDHS signal

---

## 📝 Output Format

Each `output.<ID>` file includes:

| Column | Description             |
| ------ | ----------------------- |
| 1      | Chromosome              |
| 2      | Start                   |
| 3      | End                     |
| 4      | DHS ID (`<ID>-<index>`) |
| 5      | Score from MACS3        |
| 6      | Mean signal             |

---

## 🧹 Temporary Files

Temporary/intermediate files are kept within `macs_peaks/<ID>/tmp_files/` and `DHSs/`. Cleaned automatically after final output is saved.

---

## 🧾 Citation

If you use this pipeline, please cite:

* Zhang et al., *MACS: Model-based Analysis of ChIP-Seq*, Genome Biology (2008)
* Quinlan & Hall, *BEDTools*, Bioinformatics (2010)
* UCSC Genome Browser tools

---

## 📄 License

This repository is open-source under the [MIT License](LICENSE).

---

## 🙋 Author

Developed by **\[Your Name]**
📫 [your.email@domain.com](mailto:your.email@domain.com)
🔗 GitHub: [github.com/yourusername](https://github.com/yourusername)

```

Let me know if you want to:
- Add an example dataset section
- Include command-line arguments for flexibility
- Add a visual diagram of the pipeline

Or I can generate a matching `LICENSE` or `filter-long-double.py` stub as well.
```

