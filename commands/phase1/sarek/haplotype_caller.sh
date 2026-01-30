#!/bin/bash
# Pipeline 3: Sarek + HaplotypeCaller
# Commands to run Pipeline 3 from start to finish

# Step 1: Run Sarek Pipeline
# Note: Sarek will start from variant_calling step since we have BAM files
# Using Docker profile (Singularity not available)
cd /home/mssever/Desktop/blg348e/project && \
nextflow run nf-core/sarek \
    -r master \
    -profile docker \
    --input scripts/phase1/sarek/samplesheet_haplotype_caller.csv \
    --step variant_calling \
    --tools haplotypecaller \
    --genome GATK.GRCh38 \
    --fasta data/reference/Homo_sapiens_assembly38.fasta \
    --outdir outputs/phase1/sarek/haplotype_caller \
    --skip_tools fastqc,samtools,mosdepth \
    -c scripts/phase1/sarek/nextflow.config \
    2>&1 | tee pipeline_run_sarek_haplotype.log

# Step 2: Use the output VCF file
# Sarek outputs VCFs in: outputs/phase1/sarek/haplotype_caller/variant_calling/haplotypecaller/NA12878_exome-1/NA12878_exome-1.haplotypecaller.vcf.gz
cd /home/mssever/Desktop/blg348e/project && \
VCF_FILE="outputs/phase1/sarek/haplotype_caller/variant_calling/haplotypecaller/NA12878_exome-1/NA12878_exome-1.haplotypecaller.vcf.gz" && \
if [ ! -f "$VCF_FILE" ]; then
    # Try to find it if path is different
    VCF_FILE=$(find outputs/phase1/sarek/haplotype_caller -name "*haplotypecaller*.vcf.gz" | head -1)
    if [ -z "$VCF_FILE" ] || [ ! -f "$VCF_FILE" ]; then
        echo "ERROR: Could not find Sarek HaplotypeCaller output VCF file"
        echo "Please check outputs/phase1/sarek/haplotype_caller/ directory"
        exit 1
    fi
fi && \
echo "Using VCF file: $VCF_FILE" && \
# Decompress for processing
gunzip -c "$VCF_FILE" > results/phase1/filtered/sarek/haplotype_caller/sarek_caller_uncompressed.vcf && \
INPUT_VCF="results/phase1/filtered/sarek/haplotype_caller/sarek_caller_uncompressed.vcf"

# Step 3: Fix Ploidy
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools +fixploidy "$INPUT_VCF" -O v -o results/phase1/filtered/sarek/haplotype_caller/sarek_fixed_ploidy.vcf

# Step 4: Compress and Index Ploidy-Fixed VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase1/filtered/sarek/haplotype_caller/sarek_fixed_ploidy.vcf > results/phase1/filtered/sarek/haplotype_caller/sarek_fixed_ploidy.vcf.gz && \
tabix -p vcf results/phase1/filtered/sarek/haplotype_caller/sarek_fixed_ploidy.vcf.gz

# Step 5: Exome BED Filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools view -R bed_files/phase1/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
    results/phase1/filtered/sarek/haplotype_caller/sarek_fixed_ploidy.vcf.gz \
    -O v -o results/phase1/filtered/sarek/haplotype_caller/sarek_exome_filtered.vcf

# Step 6: Compress and Index Exome-Filtered VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase1/filtered/sarek/haplotype_caller/sarek_exome_filtered.vcf > results/phase1/filtered/sarek/haplotype_caller/sarek_exome_filtered.vcf.gz && \
tabix -p vcf results/phase1/filtered/sarek/haplotype_caller/sarek_exome_filtered.vcf.gz

# Step 7: High-Confidence BED Filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools view -R bed_files/phase1/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
    results/phase1/filtered/sarek/haplotype_caller/sarek_exome_filtered.vcf.gz \
    -O v -o results/phase1/filtered/sarek/haplotype_caller/sarek_hc_filtered.vcf

# Step 8: PASS Filtering (Filter out LowQual)
cd /home/mssever/Desktop/blg348e/project && \
(grep "^#" results/phase1/filtered/sarek/haplotype_caller/sarek_hc_filtered.vcf; \
 grep -v "^#" results/phase1/filtered/sarek/haplotype_caller/sarek_hc_filtered.vcf | awk '($7 != "LowQual")') > \
 results/phase1/filtered/sarek/haplotype_caller/sarek_final_filtered.vcf

# Step 9: Compress and Index Final VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase1/filtered/sarek/haplotype_caller/sarek_final_filtered.vcf > results/phase1/filtered/sarek/haplotype_caller/sarek_final_filtered.vcf.gz && \
tabix -p vcf results/phase1/filtered/sarek/haplotype_caller/sarek_final_filtered.vcf.gz

# Step 10: Index Truth VCF (if not already indexed)
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
tabix -p vcf data/phase1/truth_vcf/NA12878_exome_hc_filtered.vcf.gz

# Step 11: Calculate Metrics (bcftools isec)
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools isec -p results/phase1/metrics/sarek/haplotype_caller/ -c both \
    results/phase1/filtered/sarek/haplotype_caller/sarek_final_filtered.vcf.gz \
    data/phase1/truth_vcf/NA12878_exome_hc_filtered.vcf.gz

# Step 12: Calculate Precision, Recall, F1-Score
cd /home/mssever/Desktop/blg348e/project && python3 << 'EOF'
import subprocess
tp = int(subprocess.check_output(['grep', '-v', '^#', 'results/phase1/metrics/sarek/haplotype_caller/0002.vcf']).decode().count('\n'))
fp = int(subprocess.check_output(['grep', '-v', '^#', 'results/phase1/metrics/sarek/haplotype_caller/0000.vcf']).decode().count('\n'))
fn = int(subprocess.check_output(['grep', '-v', '^#', 'results/phase1/metrics/sarek/haplotype_caller/0001.vcf']).decode().count('\n'))
precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0
f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
print("=" * 60)
print("Pipeline 3: Sarek + HaplotypeCaller - Metrics")
print("=" * 60)
print(f"True Positives (TP):  {tp}")
print(f"False Positives (FP): {fp}")
print(f"False Negatives (FN): {fn}")
print(f"Total in our VCF:     {tp + fp}")
print(f"Total in truth VCF:   {tp + fn}")
print("")
print(f"Precision:  {precision:.4f} ({precision*100:.2f}%)")
print(f"Recall:     {recall:.4f} ({recall*100:.2f}%)")
print(f"F1-Score:   {f1:.4f} ({f1*100:.2f}%)")
print("=" * 60)
EOF
