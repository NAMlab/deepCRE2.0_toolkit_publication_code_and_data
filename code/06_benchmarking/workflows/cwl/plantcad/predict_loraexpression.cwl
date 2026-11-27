#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

# Metadata for the tool
label: "Tokenisation of PlantCaduceus input data"
doc: |
  Performs tokenisation for PlantCaduceus input files.

  This tool processes a TSV file and produces a Parquet file for PlantCaduceus predictions.

# Input parameters
inputs:
  data_file:
    type: File
    doc: "Input TSV file"
  model_lora:
    type: Directory
    doc: "Name of the trained LoRA-fine-tuned PlantCaduceus model"
  model_base:
    type: Directory
    doc: "Name of the trained base PlantCaduceus model"
  task_type:
    type: string
    doc: "The type of prediction, either 'regression', 'classification', or 'multi_label'"

# Output files
outputs:
  prediction:
    type: File
    outputBinding:
      glob: "*.csv"
    doc: "Generated CSV file with prediction results"

baseCommand: sh

arguments:
  - -c
  - |
    python /app/plantcad/src/lora_fine_tune.py predict \
      --data_dir $(inputs.data_file.path) \
      --output_file $(inputs.data_file.basename.replace(/\.[^/.]+$/, "")).csv \
      --checkpoint_dir $(inputs.model_lora.path) \
      --model_name $(inputs.model_base.path) \
      --task_type $(inputs.task_type)

# Standard output and error handling
stdout: predict.log
stderr: predict.error.log

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "plantcad-cwl"

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 32000
  - class: cwltool:CUDARequirement
    cudaComputeCapability: "8.0"
    cudaDeviceCountMin: 1 # Can spawn one prediction per CUDA device
    cudaDeviceCountMax: 1
    cudaVersionMin: "12.4"

# Define namespaces
$namespaces:
  cwltool: "http://commonwl.org/cwltool#"
