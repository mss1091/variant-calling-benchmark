#!/bin/bash
# Pipeline 8: Sarek + Strelka2 (Bowtie)
# Commands to run Pipeline 8 from start to finish

# Step 1: Run Sarek Pipeline
# Note: Sarek will start from variant_calling step since we have BAM files
# Using Docker profile
cd /home/mssever/Desktop/blg348e/project && \
nextflow run nf-core/sarek \
    -r master \
    -profile docker \
    --input scripts/phase2/bowtie/sarek/samplesheet_strelka2.csv \
    --step variant_calling \
    --tools strelka \
    --genome GATK.GRCh38 \
    --fasta data/reference/Homo_sapiens_assembly38.fasta \
    --outdir outputs/phase2/bowtie/sarek/strelka2 \
    --skip_tools fastqc,samtools,mosdepth \
    --nucleotides_per_second 500 \
    -c scripts/phase2/bowtie/sarek/nextflow.config.strelka2 \
    2>&1 | tee pipeline_run_phase2_p8.log

# Step 2: Use the output VCF files (Strelka outputs SNVs and Indels separately)
# Merge them into one VCF for filtering
cd /home/mssever/Desktop/blg348e/project && \
SNV_VCF=$(find outputs/phase2/bowtie/sarek/strelka2 -name "*somatic_snvs.vcf.gz" | head -1) && \
INDEL_VCF=$(find outputs/phase2/bowtie/sarek/strelka2 -name "*somatic_indels.vcf.gz" | head -1) && \
if [ -z "$SNV_VCF" ] || [ ! -f "$SNV_VCF" ]; then
    echo "ERROR: Could not find Sarek Strelka SNV output VCF file"
    exit 1
fi && \
if [ -z "$INDEL_VCF" ] || [ ! -f "$INDEL_VCF" ]; then
    echo "ERROR: Could not find Sarek Strelka Indel output VCF file"
    exit 1
fi && \
echo "Using SNV VCF: $SNV_VCF" && \
echo "Using Indel VCF: $INDEL_VCF" && \
mkdir -p results/phase2/filtered/bowtie/sarek/strelka2 && \
# Merge SNVs and Indels
bcftools concat -a "$SNV_VCF" "$INDEL_VCF" -O v -o results/phase2/filtered/bowtie/sarek/strelka2/sarek_strelka2_merged.vcf && \
INPUT_VCF="results/phase2/filtered/bowtie/sarek/strelka2/sarek_strelka2_merged.vcf"

# Step 3: Fix Ploidy
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools +fixploidy "$INPUT_VCF" -O v -o results/phase2/filtered/bowtie/sarek/strelka2/sarek_fixed_ploidy.vcf

# Step 4: Compress and Index Ploidy-Fixed VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase2/filtered/bowtie/sarek/strelka2/sarek_fixed_ploidy.vcf > results/phase2/filtered/bowtie/sarek/strelka2/sarek_fixed_ploidy.vcf.gz && \
tabix -p vcf results/phase2/filtered/bowtie/sarek/strelka2/sarek_fixed_ploidy.vcf.gz

# Step 5: Exome BED Filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools view -R bed_files/phase2/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz \
results/phase2/filtered/bowtie/sarek/strelka2/sarek_fixed_ploidy.vcf.gz \
-O v -o results/phase2/filtered/bowtie/sarek/strelka2/sarek_exome_filtered.vcf

# Step 6: Compress and Index Exome-Filtered VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase2/filtered/bowtie/sarek/strelka2/sarek_exome_filtered.vcf > results/phase2/filtered/bowtie/sarek/strelka2/sarek_exome_filtered.vcf.gz && \
tabix -p vcf results/phase2/filtered/bowtie/sarek/strelka2/sarek_exome_filtered.vcf.gz

# Step 7: High-Confidence BED Filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools view -R bed_files/phase2/High-Confidence_Regions_v1.2.bed \
results/phase2/filtered/bowtie/sarek/strelka2/sarek_exome_filtered.vcf.gz \
-O v -o results/phase2/filtered/bowtie/sarek/strelka2/sarek_hc_filtered.vcf

# Step 8: Compress and Index HC-Filtered VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase2/filtered/bowtie/sarek/strelka2/sarek_hc_filtered.vcf > results/phase2/filtered/bowtie/sarek/strelka2/sarek_hc_filtered.vcf.gz && \
tabix -p vcf results/phase2/filtered/bowtie/sarek/strelka2/sarek_hc_filtered.vcf.gz

# Step 9: PASS Filtering (Keep only PASS variants)
# Using bcftools view -f 'PASS,.' for robust filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bcftools view -f 'PASS,.' results/phase2/filtered/bowtie/sarek/strelka2/sarek_hc_filtered.vcf.gz -O v -o results/phase2/filtered/bowtie/sarek/strelka2/sarek_final_filtered.vcf

# Step 10: Compress and Index Final Filtered VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
bgzip -c results/phase2/filtered/bowtie/sarek/strelka2/sarek_final_filtered.vcf > results/phase2/filtered/bowtie/sarek/strelka2/sarek_final_filtered.vcf.gz && \
tabix -p vcf results/phase2/filtered/bowtie/sarek/strelka2/sarek_final_filtered.vcf.gz

# Step 10: Ensure Truth VCF is indexed (exome-filtered for fair comparison)
# Note: Truth set is already in HC regions (filename indicates this)
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
if [ ! -f data/phase2/truth_vcf/high-confidence_sSNV_exome_filtered.vcf.gz.tbi ]; then
    tabix -p vcf data/phase2/truth_vcf/high-confidence_sSNV_exome_filtered.vcf.gz
fi

# Step 11: Calculate Metrics (bcftools isec)
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
mkdir -p results/phase2/metrics/bowtie/sarek/strelka2 && \
bcftools isec -p results/phase2/metrics/bowtie/sarek/strelka2/ -c both \
results/phase2/filtered/bowtie/sarek/strelka2/sarek_final_filtered.vcf.gz \
data/phase2/truth_vcf/high-confidence_sSNV_exome_filtered.vcf.gz

# Step 12: Calculate Precision, Recall, F1-Score
cd /home/mssever/Desktop/blg348e/project && python3 << 'EOF'
import subprocess
import os

metrics_dir = "results/phase2/metrics/bowtie/sarek/strelka2"

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
print("Pipeline 8: Sarek + Strelka2 (Bowtie) - Metrics")
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
echo "Pipeline 8 (Sarek + Strelka2, Bowtie) completed!"
echo "Results saved to: results/phase2/filtered/bowtie/sarek/strelka2/"
echo "Metrics saved to: results/phase2/metrics/bowtie/sarek/strelka2/"
