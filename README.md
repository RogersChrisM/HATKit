# HATKit
**HemTools Auxiliary Transcription Factor Kit**

HATKit is a modular Python toolkit that extends the 
[HemTools](https://github.com/YongchengYAO/HemTools) pipeline with 
downstream analysis and visualization utilities for transcription factor 
binding and regulatory element data.

## Current Modules

### `diffpeak_plots.py`
Generates volcano plots, MA plots, and replicate-correlation regression 
from HemTools diffPeak output (DESeq2 results). Supports peak and 
differential gene expression (DGE) input formats.

**Outputs:**
- Volcano plots (base and annotated)
- MA plots (full and IQR-filtered)
- Replicate correlation plots with Pearson and Spearman statistics
- Up/down-regulated peak and gene lists (TSV and BED)

**Dependencies:**
- Python 3.9+
- bedtools (must be on PATH)

```bash
pip install matplotlib numpy pandas requests adjustText scipy tqdm
```

**Usage:**
```bash
python diffpeak_plots.py -f results.tsv [options]
```
| Flag | Default | Description |
|---|---|---|
| `-f` / `--file` | required | diffPeak result TSV |
| `-r` / `--regions` | none | BED file(s) to highlight (up to 3) |
| `--dge` | off | Input is differential gene expression |
| `-lt` / `--lfcThresh` | 2.0 | Log fold-change threshold |
| `-pt` / `--pThresh` | 0.01 | FDR threshold |
| `--labelVolcano` | off | Auto-label top DE genes |
| `--labelSpecific` | none | Specific genes to label |
| `--genome` | hg19 | UCSC assembly (hg19 or hg38) |
| `--reverseFC` | off | Swap background/treat assignment |
| `--debug` | off | Verbose logging, keep temp files |

---

### `process_utils.py`
HPC/cluster utility for monitoring and sequencing job execution.

- `waitForProcess(pid)` — blocks until a given PID completes before 
  allowing downstream steps to proceed

---

## Roadmap

HATkit is under active development alongside dissertation work. 
Planned modules include additional TF binding analysis and 
regulatory element utilities built on HemTools pipeline outputs.

## License
MIT
