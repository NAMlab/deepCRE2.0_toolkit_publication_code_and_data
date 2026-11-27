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
  model:
    type: Directory
    doc: "Name of the trained PlantCaduceus model"
  sequence_length:
    type: int
    doc: "Expected sequence length for the model"

# Output files
outputs:
  tokenized:
    type: File
    outputBinding:
      glob: "*.parquet"
    doc: "Generated Parquet file with tokenized input data"

baseCommand: sh

arguments:
  - -c
  - |
    python /app/plantcad/src/lora_fine_tune.py tokenize \
      --data_dir $(inputs.data_file.path) \
      -output_path $(inputs.data_file.basename.replace(/\.[^/.]+$/, "")).parquet \
      --model_name $(inputs.model.path) \
      --sequence_length $(inputs.sequence_length)

# Standard output and error handling
stdout: tokenize.log
stderr: tokenize.error.log

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
