#!/usr/bin/env python3
"""Phase 1 Pipeline Visualization Generator"""

import subprocess
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import seaborn as sns
    sns.set_style("whitegrid")
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "phase1"
DATA_DIR = PROJECT_ROOT / "data" / "phase1"
OUTPUT_DIR = PROJECT_ROOT / "results" / "phase1" / "visualizations"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PIPELINES = {
    "P1": {"name": "COSAP + HaplotypeCaller", "final_vcf": RESULTS_DIR / "filtered" / "haplotype_caller" / "caller_final_filtered.vcf.gz"},
    "P2": {"name": "COSAP + DeepVariant", "final_vcf": RESULTS_DIR / "filtered" / "deep_variant" / "deepvariant_final_filtered.vcf.gz"},
    "P3": {"name": "Sarek + HaplotypeCaller", "final_vcf": RESULTS_DIR / "filtered" / "sarek" / "haplotype_caller" / "sarek_final_filtered.vcf.gz"},
    "P4": {"name": "Sarek + DeepVariant", "final_vcf": RESULTS_DIR / "filtered" / "sarek" / "deep_variant" / "sarek_deepvariant_final_filtered.vcf.gz"}
}

TRUTH_VCF = DATA_DIR / "truth_vcf" / "NA12878_exome_hc_filtered.vcf.gz"
PIPELINE_COLORS = {"P1": "#6BAED6", "P2": "#FD8D3C", "P3": "#78C679", "P4": "#9E9AC8"}
PIPELINE_COLOR_LIST = [PIPELINE_COLORS["P1"], PIPELINE_COLORS["P2"], PIPELINE_COLORS["P3"], PIPELINE_COLORS["P4"]]


def count_variants(vcf_path):
    if not vcf_path.exists():
        return 0
    try:
        cmd = f"bcftools view -H {vcf_path} 2>/dev/null | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return int(result.stdout.strip())
    except:
        return 0


def read_vcf_variants(vcf_path):
    variants = set()
    if not vcf_path.exists():
        return variants
    try:
        cmd = f"zcat {vcf_path}" if str(vcf_path).endswith('.gz') else f"cat {vcf_path}"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
        for line in process.stdout:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                variants.add(f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4].split(',')[0]}")
        process.wait()
    except Exception as e:
        print(f"Error reading {vcf_path}: {e}")
    return variants


def get_metrics_from_files():
    metrics = {}
    metrics_paths = {
        "P1": RESULTS_DIR / "metrics" / "haplotype_caller",
        "P2": RESULTS_DIR / "metrics" / "deep_variant",
        "P3": RESULTS_DIR / "metrics" / "sarek" / "haplotype_caller",
        "P4": RESULTS_DIR / "metrics" / "sarek" / "deep_variant"
    }
    
    for pid, mdir in metrics_paths.items():
        tp = count_variants(mdir / "0002.vcf")
        fp = count_variants(mdir / "0000.vcf")
        fn = count_variants(mdir / "0001.vcf")
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0
        
        metrics[pid] = {"TP": tp, "FP": fp, "FN": fn, "Precision": precision, "Recall": recall, "F1": f1}
    return metrics


def visualization_1_filtering_counts():
    print("Creating Visualization 1: Variant counts after filtering...")
    
    filtering_steps = ["Raw VCF\n(After Calling)", "After Ploidy\nFix", "After Exome\nFilter", "After HC\nFilter", "Final PASS\nFilter"]
    pipeline_counts_data = {
        "P1": [32850, 32850, 1586, 1586, 1433],
        "P2": [32961, 32961, 1644, 1644, 1477],
        "P3": [31819, 31819, 1586, 1586, 1433],
        "P4": [32963, 32963, 1644, 1644, 1477]
    }
    
    fig, ax = plt.subplots(figsize=(14, 7))
    x = np.arange(len(filtering_steps))
    width = 0.18
    
    for i, pid in enumerate(["P1", "P2", "P3", "P4"]):
        counts = pipeline_counts_data[pid]
        positions = x + i*width
        bars = ax.bar(positions, counts, width, label=pid, color=PIPELINE_COLOR_LIST[i])
        for bar, count in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 200, 
                   f'{count:,}', ha='center', va='bottom', fontsize=7, rotation=90)
    
    ax.set_xlabel('Filtering Step', fontweight='bold')
    ax.set_ylabel('Variant Count', fontweight='bold')
    ax.set_title('Variant Counts After Each Filtering Step', fontweight='bold')
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels(filtering_steps)
    ax.legend(title='Pipeline', loc='upper right')
    ax.set_ylim(0, 38000)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "1_filtering_counts.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR / '1_filtering_counts.png'}")


def visualization_2_metrics():
    print("Creating Visualization 2: Performance metrics...")
    
    metrics = get_metrics_from_files()
    pipeline_labels = ["P1", "P2", "P3", "P4"]
    
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.axis('tight')
    ax.axis('off')
    
    table_data = []
    for pid in pipeline_labels:
        table_data.append([
            pid, f"{metrics[pid]['TP']}", f"{metrics[pid]['FP']}", f"{metrics[pid]['FN']}",
            f"{metrics[pid]['Precision']:.2f}", f"{metrics[pid]['Recall']:.2f}", f"{metrics[pid]['F1']:.2f}"
        ])
    
    table = ax.table(cellText=table_data,
                     colLabels=['Pipeline', 'TP', 'FP', 'FN', 'Precision', 'Recall', 'F1'],
                     cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.5)
    
    for i, pid in enumerate(pipeline_labels):
        table[(i+1, 0)].set_facecolor(PIPELINE_COLOR_LIST[i])
        table[(i+1, 0)].set_text_props(weight='bold', color='white')
        for j in range(1, 7):
            table[(i+1, j)].set_facecolor('#f0f0f0')
    
    ax.set_title('Performance Metrics Summary Table', fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "2_summary_table.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR / '2_summary_table.png'}")


def visualization_3_similarity_matrix():
    print("Creating Visualization 3: Pipeline similarity matrix...")
    
    pipeline_order = ["P1", "P3", "P2", "P4"]
    pipeline_variants = {}
    
    for pid in pipeline_order:
        variants = read_vcf_variants(PIPELINES[pid]["final_vcf"])
        pipeline_variants[pid] = variants
        print(f"  {pid}: {len(variants)} variants")
    
    n = len(pipeline_order)
    similarity_matrix = np.zeros((n, n))
    
    for i, pid1 in enumerate(pipeline_order):
        for j, pid2 in enumerate(pipeline_order):
            set1, set2 = pipeline_variants[pid1], pipeline_variants[pid2]
            intersection = len(set1 & set2)
            union = len(set1 | set2)
            similarity = intersection / union if union > 0 else 1.0
            similarity_matrix[i, j] = similarity
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    annot = [[f'{similarity_matrix[i, j]:.2f}' for j in range(n)] for i in range(n)]
    labels = [f"{pid}\n{PIPELINES[pid]['name']}" for pid in pipeline_order]
    sns.heatmap(similarity_matrix, annot=annot, fmt='', cmap='YlOrRd',
               vmin=0, vmax=1, xticklabels=labels, yticklabels=labels,
               cbar_kws={'label': 'Jaccard Similarity'}, ax=ax, square=True,
               linewidths=0.5, linecolor='white')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    ax.set_xlabel('Pipeline', fontweight='bold')
    ax.set_ylabel('Pipeline', fontweight='bold')
    ax.set_title('Pipeline Similarity Matrix (Jaccard Index)', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "3_similarity_matrix.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR / '3_similarity_matrix.png'}")


def main():
    print("=" * 60)
    print("Phase 1 Pipeline Visualization Generator")
    print("=" * 60)
    
    try:
        visualization_1_filtering_counts()
        visualization_2_metrics()
        visualization_3_similarity_matrix()
        print("\nAll visualizations created successfully!")
        print(f"Output directory: {OUTPUT_DIR}")
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
