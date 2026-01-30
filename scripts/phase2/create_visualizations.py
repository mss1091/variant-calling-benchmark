#!/usr/bin/env python3
"""Phase 2 Pipeline Visualization Generator"""

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
RESULTS_DIR = PROJECT_ROOT / "results" / "phase2"
DATA_DIR = PROJECT_ROOT / "data" / "phase2"
OUTPUT_DIR = PROJECT_ROOT / "results" / "phase2" / "visualizations"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PIPELINES = {
    "P1": {"name": "BWA + COSAP + MuTect2", "mapper": "BWA", "workflow": "COSAP", "caller": "MuTect2",
           "final_vcf": RESULTS_DIR / "filtered" / "bwa" / "cosap" / "mutect2" / "mutect2_final_filtered.vcf.gz"},
    "P2": {"name": "BWA + COSAP + Strelka2", "mapper": "BWA", "workflow": "COSAP", "caller": "Strelka2",
           "final_vcf": RESULTS_DIR / "filtered" / "bwa" / "cosap" / "strelka2" / "strelka2_final_filtered.vcf.gz"},
    "P3": {"name": "BWA + Sarek + MuTect2", "mapper": "BWA", "workflow": "Sarek", "caller": "MuTect2",
           "final_vcf": RESULTS_DIR / "filtered" / "bwa" / "sarek" / "mutect2" / "sarek_final_filtered.vcf.gz"},
    "P4": {"name": "BWA + Sarek + Strelka2", "mapper": "BWA", "workflow": "Sarek", "caller": "Strelka2",
           "final_vcf": RESULTS_DIR / "filtered" / "bwa" / "sarek" / "strelka2" / "sarek_final_filtered.vcf.gz"},
    "P5": {"name": "Bowtie + COSAP + MuTect2", "mapper": "Bowtie", "workflow": "COSAP", "caller": "MuTect2",
           "final_vcf": RESULTS_DIR / "filtered" / "bowtie" / "cosap" / "mutect2" / "mutect2_final_filtered.vcf.gz"},
    "P6": {"name": "Bowtie + COSAP + Strelka2", "mapper": "Bowtie", "workflow": "COSAP", "caller": "Strelka2",
           "final_vcf": RESULTS_DIR / "filtered" / "bowtie" / "cosap" / "strelka2" / "strelka2_final_filtered.vcf.gz"},
    "P7": {"name": "Bowtie + Sarek + MuTect2", "mapper": "Bowtie", "workflow": "Sarek", "caller": "MuTect2",
           "final_vcf": RESULTS_DIR / "filtered" / "bowtie" / "sarek" / "mutect2" / "sarek_final_filtered.vcf.gz"},
    "P8": {"name": "Bowtie + Sarek + Strelka2", "mapper": "Bowtie", "workflow": "Sarek", "caller": "Strelka2",
           "final_vcf": RESULTS_DIR / "filtered" / "bowtie" / "sarek" / "strelka2" / "sarek_final_filtered.vcf.gz"},
}

TRUTH_VCF = DATA_DIR / "truth_vcf" / "high-confidence_sSNV_exome_filtered.vcf.gz"
PIPELINE_COLORS = {
    "P1": "#6BAED6", "P2": "#FD8D3C", "P3": "#78C679", "P4": "#9E9AC8",
    "P5": "#969696", "P6": "#FFC000", "P7": "#08519C", "P8": "#54278F"
}
PIPELINE_COLOR_LIST = [PIPELINE_COLORS[f"P{i}"] for i in range(1, 9)]


def count_variants(vcf_path):
    if not vcf_path or not vcf_path.exists():
        return 0
    try:
        cmd = f"bcftools view -H {vcf_path} 2>/dev/null | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return int(result.stdout.strip())
    except:
        return 0


def read_vcf_variants(vcf_path):
    variants = set()
    if not vcf_path or not vcf_path.exists():
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
        "P1": RESULTS_DIR / "metrics" / "bwa" / "cosap" / "mutect2",
        "P2": RESULTS_DIR / "metrics" / "bwa" / "cosap" / "strelka2",
        "P3": RESULTS_DIR / "metrics" / "bwa" / "sarek" / "mutect2",
        "P4": RESULTS_DIR / "metrics" / "bwa" / "sarek" / "strelka2",
        "P5": RESULTS_DIR / "metrics" / "bowtie" / "cosap" / "mutect2",
        "P6": RESULTS_DIR / "metrics" / "bowtie" / "cosap" / "strelka2",
        "P7": RESULTS_DIR / "metrics" / "bowtie" / "sarek" / "mutect2",
        "P8": RESULTS_DIR / "metrics" / "bowtie" / "sarek" / "strelka2"
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
        "P1": [34914, 34914, 25016, 25016, 12020],
        "P2": [218356, 218356, 155541, 155541, 23368],
        "P3": [32193, 32193, 23368, 23368, 20466],
        "P4": [216856, 217745, 155817, 155817, 23389],
        "P5": [31261, 31261, 23048, 23048, 11580],
        "P6": [177274, 177274, 129305, 129305, 4867],
        "P7": [29598, 29598, 21914, 21914, 19707],
        "P8": [176918, 177769, 129597, 129597, 4884]
    }
    
    fig, ax = plt.subplots(figsize=(16, 8))
    x = np.arange(len(filtering_steps))
    width = 0.1
    
    for i, pid in enumerate([f"P{j}" for j in range(1, 9)]):
        counts = pipeline_counts_data[pid]
        positions = x + i*width
        bars = ax.bar(positions, counts, width, label=pid, color=PIPELINE_COLOR_LIST[i])
        for bar, count in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1500, 
                   f'{count:,}', ha='center', va='bottom', fontsize=6, rotation=90)
    
    ax.set_xlabel('Filtering Step', fontweight='bold')
    ax.set_ylabel('Variant Count', fontweight='bold')
    ax.set_title('Variant Counts After Each Filtering Step', fontweight='bold')
    ax.set_xticks(x + width * 3.5)
    ax.set_xticklabels(filtering_steps)
    ax.legend(title='Pipeline', ncol=2, loc='upper right')
    ax.set_ylim(0, 250000)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "1_filtering_counts.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR / '1_filtering_counts.png'}")


def visualization_2_metrics():
    print("Creating Visualization 2: Performance metrics...")
    
    metrics = get_metrics_from_files()
    pipeline_ids = [f"P{i}" for i in range(1, 9)]
    
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.axis('tight')
    ax.axis('off')
    
    table_data = []
    for pid in pipeline_ids:
        m = metrics[pid]
        table_data.append([
            pid, f"{m['TP']}", f"{m['FP']}", f"{m['FN']}",
            f"{m['Precision']:.2f}", f"{m['Recall']:.2f}", f"{m['F1']:.2f}"
        ])
    
    table = ax.table(cellText=table_data,
                     colLabels=['Pipeline', 'TP', 'FP', 'FN', 'Precision', 'Recall', 'F1'],
                     cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    for i, pid in enumerate(pipeline_ids):
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
    print("Creating Visualization 3: Pipeline Similarity Matrix...")
    
    pipeline_order = ["P6", "P8", "P2", "P4", "P1", "P5", "P3", "P7"]
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
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    annot = [[f'{similarity_matrix[i, j]:.2f}' for j in range(n)] for i in range(n)]
    labels = [f"{pid}\n{PIPELINES[pid]['workflow']} + {PIPELINES[pid]['caller']}" for pid in pipeline_order]
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
    print("Phase 2 Pipeline Visualization Generator")
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
