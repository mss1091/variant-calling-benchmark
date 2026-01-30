#!/bin/bash
# Pipeline 2: COSAP + DeepVariant (BWA)
# Commands to run Pipeline 2 from start to finish

# Step 1: Run COSAP Pipeline
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && python scripts/phase1/cosap/phase1_deep_variant.py 2>&1 | tee pipeline_run_deepvariant.log

# Step 2: Fix Ploidy
# Note: DeepVariant outputs caller.g.vcf (GVCF format)
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && \
if [ -f "outputs/phase1/cosap/deep_variant/VCF/deepvariant/caller.g.vcf" ]; then
    bcftools +fixploidy outputs/phase1/cosap/deep_variant/VCF/deepvariant/caller.g.vcf -O v -o results/phase1/filtered/deep_variant/deepvariant_fixed_ploidy.vcf
elif [ -f "outputs/phase1/cosap/deep_variant/VCF/deepvariant/caller.vcf" ]; then
    bcftools +fixploidy outputs/phase1/cosap/deep_variant/VCF/deepvariant/caller.vcf -O v -o results/phase1/filtered/deep_variant/deepvariant_fixed_ploidy.vcf
else
    echo "ERROR: Could not find DeepVariant output VCF file"
    exit 1
fi

# Step 3: Compress and Index Ploidy-Fixed VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && bgzip -c results/phase1/filtered/deep_variant/deepvariant_fixed_ploidy.vcf > results/phase1/filtered/deep_variant/deepvariant_fixed_ploidy.vcf.gz && tabix -p vcf results/phase1/filtered/deep_variant/deepvariant_fixed_ploidy.vcf.gz

# Step 4: Exome BED Filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && bcftools view -R bed_files/phase1/nexterarapidcapture_expandedexome_targetedregions.bed.gz results/phase1/filtered/deep_variant/deepvariant_fixed_ploidy.vcf.gz -O v -o results/phase1/filtered/deep_variant/deepvariant_exome_filtered.vcf

# Step 5: Compress and Index Exome-Filtered VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && bgzip -c results/phase1/filtered/deep_variant/deepvariant_exome_filtered.vcf > results/phase1/filtered/deep_variant/deepvariant_exome_filtered.vcf.gz && tabix -p vcf results/phase1/filtered/deep_variant/deepvariant_exome_filtered.vcf.gz

# Step 6: High-Confidence BED Filtering
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && bcftools view -R bed_files/phase1/HG001_GRCh38_1_22_v4.2.1_benchmark.bed results/phase1/filtered/deep_variant/deepvariant_exome_filtered.vcf.gz -O v -o results/phase1/filtered/deep_variant/deepvariant_hc_filtered.vcf

# Step 7: PASS Filtering
# Note: DeepVariant outputs GVCF format with reference blocks (<*> in ALT column) that need to be filtered out
# Filter: Keep PASS or "." in FILTER, exclude RefCall, and exclude reference blocks (<*>)
cd /home/mssever/Desktop/blg348e/project && (grep "^#" results/phase1/filtered/deep_variant/deepvariant_hc_filtered.vcf; grep -v "^#" results/phase1/filtered/deep_variant/deepvariant_hc_filtered.vcf | awk '($5 != "<*>" && $7 != "RefCall" && ($7 == "PASS" || $7 == "."))') > results/phase1/filtered/deep_variant/deepvariant_final_filtered.vcf

# Step 8: Compress and Index Final VCF
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && bgzip -c results/phase1/filtered/deep_variant/deepvariant_final_filtered.vcf > results/phase1/filtered/deep_variant/deepvariant_final_filtered.vcf.gz && tabix -p vcf results/phase1/filtered/deep_variant/deepvariant_final_filtered.vcf.gz

# Step 9: Index Truth VCF (if not already indexed)
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && tabix -p vcf data/phase1/truth_vcf/NA12878_exome_hc_filtered.vcf.gz

# Step 10: Calculate Metrics (bcftools isec)
cd /home/mssever/Desktop/blg348e/project && source setup_cosap.sh && bcftools isec -p results/phase1/metrics/deep_variant/ -c both results/phase1/filtered/deep_variant/deepvariant_final_filtered.vcf.gz data/phase1/truth_vcf/NA12878_exome_hc_filtered.vcf.gz

# Step 11: Calculate Precision, Recall, F1-Score
cd /home/mssever/Desktop/blg348e/project && python3 << 'EOF'
import subprocess
tp = int(subprocess.check_output(['grep', '-v', '^#', 'results/phase1/metrics/deep_variant/0002.vcf']).decode().count('\n'))
fp = int(subprocess.check_output(['grep', '-v', '^#', 'results/phase1/metrics/deep_variant/0000.vcf']).decode().count('\n'))
fn = int(subprocess.check_output(['grep', '-v', '^#', 'results/phase1/metrics/deep_variant/0001.vcf']).decode().count('\n'))
precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0
f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
print("=" * 60)
print("Pipeline 2: COSAP + DeepVariant - Metrics")
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
