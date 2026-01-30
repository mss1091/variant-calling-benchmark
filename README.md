# BLG348E Term Project: Comparative Analysis of Germline and Somatic Variant Calling Pipelines

## Authors

This comprehensive bioinformatics analysis was successfully conducted by **Mehmet Sait Sever** ([@mssever](https://github.com/mssever)) and my dear friend **Tunahan Geçit** ([@tunahangecit](https://github.com/tunahangecit)) as part of the BLG348E course at Istanbul Technical University.

## Project Overview

This project performs a comprehensive comparative analysis of variant calling pipelines for both **germline** (Phase 1) and **somatic** (Phase 2) variant detection. The analysis evaluates different combinations of workflows (COSAP vs Sarek), variant callers (HaplotypeCaller, DeepVariant, MuTect2, Strelka2), and alignment tools (BWA, Bowtie) to assess their performance, accuracy, and detection patterns.

> **Note**: This project was developed and tested on a **Linux system** (Ubuntu), which is recommended for bioinformatics pipelines due to better terminal support and compatibility with tools like Docker, Nextflow, and bcftools.

## Project Structure

```
project/
├── data/                          # Input data files
│   ├── reference/                 # Reference genome (hg38)
│   ├── phase1/                    # Germline analysis data (NA12878)
│   └── phase2/                    # Somatic analysis data (tumor-normal pairs)
│
├── bed_files/                     # BED files for filtering
│   ├── phase1/                    # Phase 1 BED files
│   └── phase2/                    # Phase 2 BED files
│
├── scripts/                       # Analysis scripts
│   ├── phase1/                    # Germline variant calling scripts
│   └── phase2/                    # Somatic variant calling scripts
│
├── commands/                      # Pipeline execution scripts
│   ├── phase1/                    # Phase 1 command scripts
│   └── phase2/                    # Phase 2 command scripts
│
├── outputs/                       # Raw pipeline outputs
│   ├── phase1/                    # Phase 1 unfiltered VCFs
│   └── phase2/                    # Phase 2 unfiltered VCFs
│
├── results/                       # Final analysis results
│   ├── phase1/                    # Phase 1 filtered VCFs, metrics, visualizations
│   └── phase2/                    # Phase 2 filtered VCFs, metrics, visualizations
│
└── setup_cosap.sh                 # Environment setup script
```

## Phase 1: Germline Variant Calling

### Pipelines Evaluated

1. **P1**: COSAP + HaplotypeCaller
2. **P2**: COSAP + DeepVariant
3. **P3**: Sarek + HaplotypeCaller
4. **P4**: Sarek + DeepVariant

### Workflow

1. **Variant Calling**: Run each pipeline on NA12878 exome data
2. **Filtering**: Apply ploidy fix, exome filtering, and PASS filter
3. **Metrics Calculation**: Compare against GIAB truth set
4. **Visualization**: Generate comprehensive visualizations

### Key Visualizations

- Filtering counts (variant reduction at each step)
- Performance metrics (TP, FP, FN, Precision, Recall, F1-Score)
- Similarity matrix (Jaccard similarity between pipelines)

## Phase 2: Somatic Variant Calling

### Pipelines Evaluated

1. **P1**: BWA + COSAP + MuTect2
2. **P2**: BWA + COSAP + Strelka2
3. **P3**: BWA + Sarek + MuTect2
4. **P4**: BWA + Sarek + Strelka2
5. **P5**: Bowtie + COSAP + MuTect2
6. **P6**: Bowtie + COSAP + Strelka2
7. **P7**: Bowtie + Sarek + MuTect2
8. **P8**: Bowtie + Sarek + Strelka2

### Workflow

1. **Variant Calling**: Run each pipeline on tumor-normal pairs
2. **Filtering**: Apply ploidy fix, exome filtering, and PASS filter
3. **Metrics Calculation**: Compare against high-confidence somatic truth set
4. **Visualization**: Generate comprehensive visualizations

### Key Visualizations

- Filtering counts (variant reduction at each step)
- Performance metrics (TP, FP, FN, Precision, Recall, F1-Score)
- Similarity matrix (Jaccard similarity between pipelines)

## Key Findings

### Phase 1 (Germline)
- **Workflow Equivalence**: COSAP and Sarek produce identical results when using the same variant caller (Jaccard = 1.00)
- **Caller Effect**: DeepVariant shows slight advantage over HaplotypeCaller (F1: 0.56 vs 0.55)

### Phase 2 (Somatic)
- **Both Mapper and Caller Matter**: Detection patterns vary significantly across mapper-caller combinations
- **Precision-Recall Trade-off**: Bowtie+Strelka2 achieves best F1-scores; BWA-based pipelines offer better recall

## Dependencies

### System Requirements
- **Linux** (Ubuntu 20.04+ recommended)
- **Python 3.9+**
- **Docker** (for containerized tools)

### Workflow Managers

**COSAP** (Comparative Sequencing Analysis Platform) is an easy yet comprehensive pipeline creation tool for NGS data developed by [MBaysanLab](https://github.com/MBaysanLab) at Istanbul Technical University. It provides reproducibility and allows comparison of different pipeline configurations.

- **GitHub**: [https://github.com/MBaysanLab/cosap](https://github.com/MBaysanLab/cosap)
- **Documentation**: [https://docs.cosap.bio](https://docs.cosap.bio)
- **Installation**: [https://docs.cosap.bio/fundamentals/getting-set-up](https://docs.cosap.bio/fundamentals/getting-set-up)

**Sarek** is an nf-core pipeline for variant calling:
- **Documentation**: [https://nf-co.re/sarek](https://nf-co.re/sarek)
- Requires **Nextflow** workflow manager

### Variant Calling Tools
- **GATK** (HaplotypeCaller, MuTect2)
- **DeepVariant**
- **Strelka2**
- **BWA** / **Bowtie** (alignment)
- **bcftools** (VCF manipulation)

## Contact

For questions or issues, please contact:
- **Mehmet Sait Sever**: severm21@itu.edu.tr
- **Tunahan Geçit**: gecitt21@itu.edu.tr
- Or the course instructors at Istanbul Technical University
