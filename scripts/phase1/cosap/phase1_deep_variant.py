#!/usr/bin/env python3
"""Phase 1: Germline Variant Calling with DeepVariant (COSAP)"""

from cosap.workflows import Pipeline, PipelineRunner, BamReader, VariantCaller
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))

REFERENCE_PATH = os.path.join(PROJECT_ROOT, "data", "reference", "Homo_sapiens_assembly38.fasta")
INPUT_BAM_PATH = os.path.join(PROJECT_ROOT, "data", "phase1", "bam", "NA12878_exome.bam")
OUTPUT_BASE_DIR = os.path.join(PROJECT_ROOT, "outputs", "phase1", "cosap", "deep_variant")

def run_pipeline():
    pipeline = Pipeline()
    os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)
    
    germline_sample = BamReader(INPUT_BAM_PATH, name="NA12878")
    variant_caller = VariantCaller(
        library="DeepVariant",
        name="caller",
        germline=germline_sample,
        gvcf=False
    )
    
    germline_sample.next_step = variant_caller
    pipeline.add(germline_sample)
    pipeline.add(variant_caller)
    pipeline_config = pipeline.build(workdir=OUTPUT_BASE_DIR)
    
    print("Running COSAP Pipeline: BWA + DeepVariant")
    print(f"Input: {INPUT_BAM_PATH}")
    print(f"Output: {OUTPUT_BASE_DIR}")
    
    runner = PipelineRunner(device="cpu")
    runner.run_pipeline(pipeline_config)
    print("Pipeline completed")

if __name__ == "__main__":
    run_pipeline()
