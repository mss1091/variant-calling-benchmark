#!/bin/bash
# Pipeline 7: Sarek + MuTect2 (Bowtie)
# Commands to run Pipeline 7 from start to finish

# Step 1: Run Sarek Pipeline
# Note: Sarek will start from variant_calling step since we have BAM files
# Using Docker profile
cd /home/mssever/Desktop/blg348e/project && \
nextflow run nf-core/sarek \
    -r master \
    -profile docker \
    --input scripts/phase2/bowtie/sarek/samplesheet_mutect2.csv \
    --step variant_calling \
    --tools mutect2 \
    --genome GATK.GRCh38 \
    --fasta data/reference/Homo_sapiens_assembly38.fasta \
    --outdir outputs/phase2/bowtie/sarek/mutect2 \
    --skip_tools fastqc,samtools,mosdepth \
    --nucleotides_per_second 500 \
    -c scripts/phase2/bowtie/sarek/nextflow.config \
    2>&1 | tee pipeline_run_phase2_p7.log

# Step 2: Use the output VCF file
# Sarek outputs MuTect2 VCFs in: outputs/phase2/bowtie/sarek/mutect2/variant_calling/mutect2/tumor_vs_normal/tumor_vs_normal.mutect2.vcf.gz
cd /home/mssever/Desktop/blg348e/project && \
VCF_FILE=$(find outputs/phase2/bowtie/sarek/mutect2 -name "*mutect2*.vcf.gz" | grep -v "filtered" | head -1) && \
if [ -z "$VCF_FILE" ] || [ ! -f "$VCF_FILE" ]; then
    echo "ERROR: Could not find Sarek MuTect2 output VCF file"
    echo "Please check outputs/phase2/bowtie/sarek/mutect2/ directory"
    exit 1
fi && \
echo "Using VCF file: $VCF_FILE" && \
# Decompress for processing
gunzip -c "$VCF_FILE" > results/phase2/filtered/bowtie/sarek/mutect2/sarek_mutect2_uncompressed.vcf && \
INPUT_VCF="results/phase2/filtered/bowtie/sarek/mutect2/sarek_mutect2_uncompressed.vcf"

# Step 3: Fix Ploidy
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
mkdir -p results/phase2/filtered/bowtie/sarek/mutect2 && \
bcftools +fixploidy "$INPUT_VCF" -O v -o results/phase2/filtered/bowtie/sarek/mutect2/sarek_fixed_ploidy.vcf

# Step 4: Compress and Index Ploidy-Fixed VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase2/filtered/bowtie/sarek/mutect2/sarek_fixed_ploidy.vcf > results/phase2/filtered/bowtie/sarek/mutect2/sarek_fixed_ploidy.vcf.gz && \
tabix -p vcf results/phase2/filtered/bowtie/sarek/mutect2/sarek_fixed_ploidy.vcf.gz

# Step 5: Exome BED Filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools view -R bed_files/phase2/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz \
results/phase2/filtered/bowtie/sarek/mutect2/sarek_fixed_ploidy.vcf.gz \
-O v -o results/phase2/filtered/bowtie/sarek/mutect2/sarek_exome_filtered.vcf

# Step 6: Compress and Index Exome-Filtered VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase2/filtered/bowtie/sarek/mutect2/sarek_exome_filtered.vcf > results/phase2/filtered/bowtie/sarek/mutect2/sarek_exome_filtered.vcf.gz && \
tabix -p vcf results/phase2/filtered/bowtie/sarek/mutect2/sarek_exome_filtered.vcf.gz

# Step 7: High-Confidence BED Filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools view -R bed_files/phase2/High-Confidence_Regions_v1.2.bed \
results/phase2/filtered/bowtie/sarek/mutect2/sarek_exome_filtered.vcf.gz \
-O v -o results/phase2/filtered/bowtie/sarek/mutect2/sarek_hc_filtered.vcf

# Step 8: Compress and Index HC-Filtered VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase2/filtered/bowtie/sarek/mutect2/sarek_hc_filtered.vcf > results/phase2/filtered/bowtie/sarek/mutect2/sarek_hc_filtered.vcf.gz && \
tabix -p vcf results/phase2/filtered/bowtie/sarek/mutect2/sarek_hc_filtered.vcf.gz

# Step 9: PASS Filtering (Keep only PASS variants)
# Using bcftools view -f 'PASS,.' for robust filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools view -f 'PASS,.' results/phase2/filtered/bowtie/sarek/mutect2/sarek_hc_filtered.vcf.gz -O v -o results/phase2/filtered/bowtie/sarek/mutect2/sarek_final_filtered.vcf

# Step 10: Compress and Index Final Filtered VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase2/filtered/bowtie/sarek/mutect2/sarek_final_filtered.vcf > results/phase2/filtered/bowtie/sarek/mutect2/sarek_final_filtered.vcf.gz && \
tabix -p vcf results/phase2/filtered/bowtie/sarek/mutect2/sarek_final_filtered.vcf.gz

# Step 10: Ensure Truth VCF is indexed (exome-filtered for fair comparison)
# Note: Truth set is already in HC regions (filename indicates this)
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
if [ ! -f data/phase2/truth_vcf/high-confidence_sSNV_exome_filtered.vcf.gz.tbi ]; then
    tabix -p vcf data/phase2/truth_vcf/high-confidence_sSNV_exome_filtered.vcf.gz
fi

# Step 11: Calculate Metrics (bcftools isec)
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
mkdir -p results/phase2/metrics/bowtie/sarek/mutect2 && \
bcftools isec -p results/phase2/metrics/bowtie/sarek/mutect2/ -c both \
results/phase2/filtered/bowtie/sarek/mutect2/sarek_final_filtered.vcf.gz \
data/phase2/truth_vcf/high-confidence_sSNV_exome_filtered.vcf.gz

# Step 12: Calculate Precision, Recall, F1-Score
cd /home/mssever/Desktop/blg348e/project && python3 << 'EOF'
import subprocess
import os

metrics_dir = "results/phase2/metrics/bowtie/sarek/mutect2"

# Count variants (excluding header lines)
def count_variants(vcf_file):
    if not os.path.exists(vcf_file):
        return 0
    try:
        result = subprocess.run(['grep', '-v', '^#', vcf_file], capture_output=True, text=True)
        return len([line for line in result.stdout.strip().split('\n') if line.strip()])
    except:
        return 0

tp = count_variants(f"{metrics_dir}/0002.vcf")
fp = count_variants(f"{metrics_dir}/0000.vcf")
fn = count_variants(f"{metrics_dir}/0001.vcf")

precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0
f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

print("=" * 60)
print("Pipeline 7: Sarek + MuTect2 (Bowtie) - Metrics")
print("=" * 60)
print(f"True Positives (TP):  {tp}")
print(f"False Positives (FP): {fp}")
print(f"False Negatives (FN): {fn}")
print(f"Total in our VCF:     {tp + fp}")
print(f"Total in truth VCF:   {tp + fn}")
print("-" * 60)
print(f"Precision:            {precision:.4f} ({precision*100:.2f}%)")
print(f"Recall:               {recall:.4f} ({recall*100:.2f}%)")
print(f"F1-Score:             {f1:.4f}")
print("=" * 60)
EOF

echo ""
echo "Pipeline 7 (COSAP + MuTect2, Bowtie) completed!"
echo "Results saved to: results/phase2/filtered/bowtie/sarek/mutect2/"
echo "Metrics saved to: results/phase2/metrics/bowtie/sarek/mutect2/"
