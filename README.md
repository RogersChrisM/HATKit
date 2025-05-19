# HATKit
**HemTools Auxiliary Transcription Factor Kit**

HATKit is a modular toolkit designed to extend HemTools with custom scripts and workflows for analyzing transcription factor binding and regulatory elements. It supports reformatting, peak processing, motif integration, and more.

## Installation
```bash
pip install -e .
```

## Usage
Example command-line tool:
```bash
hatkit convert --input foo.bed --output bar.bw
```

## Modules
- `format_conversion`: Convert between genomic file formats
- `peak_analysis`: Merge, compare, or score TF peak data
- `motif_utils`: Handle motif file parsing and matching
