#!/usr/bin/env python3
"""Phase 2: Somatic Variant Calling with MuTect2 (BWA + COSAP)"""

from cosap.workflows import Pipeline, PipelineRunner, BamReader, VariantCaller
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", "..", ".."))

REFERENCE_PATH = os.path.join(PROJECT_ROOT, "data", "reference", "Homo_sapiens_assembly38.fasta")
TUMOR_BAM_PATH = os.path.join(PROJECT_ROOT, "data", "phase2", "bam", "bwa", "unprocessed_tumor_bwa.sorted.bam")
NORMAL_BAM_PATH = os.path.join(PROJECT_ROOT, "data", "phase2", "bam", "bwa", "unprocessed_normal_bwa.sorted.bam")
OUTPUT_BASE_DIR = os.path.join(PROJECT_ROOT, "outputs", "phase2", "bwa", "cosap", "mutect2")

def run_pipeline():
    pipeline = Pipeline()
    os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)
    
    tumor_sample = BamReader(TUMOR_BAM_PATH, name="tumor")
    normal_sample = BamReader(NORMAL_BAM_PATH, name="normal")
    
    variant_caller = VariantCaller(
        library="MuTect2",
        name="mutect2",
        tumor=tumor_sample,
        germline=normal_sample
    )
    
    tumor_sample.next_step = variant_caller
    normal_sample.next_step = variant_caller
    
    pipeline.add(tumor_sample)
    pipeline.add(normal_sample)
    pipeline.add(variant_caller)
    pipeline_config = pipeline.build(workdir=OUTPUT_BASE_DIR)
    
    print("Running COSAP Pipeline: BWA + MuTect2")
    print(f"Tumor: {TUMOR_BAM_PATH}")
    print(f"Normal: {NORMAL_BAM_PATH}")
    print(f"Output: {OUTPUT_BASE_DIR}")
    
    runner = PipelineRunner(device="cpu")
    runner.run_pipeline(pipeline_config)
    print("Pipeline completed")

if __name__ == "__main__":
    run_pipeline()
