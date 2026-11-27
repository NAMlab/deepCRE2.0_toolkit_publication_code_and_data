#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

# Metadata for the tool
label: "DeepCis Prediction"
doc: |
  Performs deepcis prediction for chunked sequences.

  This tool processes a Fasta file and produces a CSV file with DeepCis predictions.

# Base command - using shell to chain two Python scripts
baseCommand: [python, /app/predict.py]

# Input parameters
inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
    doc: "Input FASTA file"
  output_file:
    type: string?
    inputBinding:
      position: 2
      prefix: --output_csv
    doc: "Output FASTA file, defaults to the input file name"
  model:
    type: File
    inputBinding:
      prefix: --model
    doc: "Path to the trained DeepCis model"
  length:
    type: int?
    default: 250
    inputBinding:
      prefix: --length
    doc: "Expected length of sequence chunks"
  overlap:
    type: int
    inputBinding:
      prefix: --overlap
    doc: "Overlap of sequence chunks (1-length)"
  threshold:
    type: float?
    inputBinding:
      prefix: --threshold
    doc: "Applied threshold to binarize output"

# Output files
outputs:
  predictions:
    type: File
    outputBinding:
      glob: "*.csv"
    doc: "Generated CSV file with prediction results"

# Standard output and error handling
stdout: predict_deepcis.log
stderr: predict_deepcis.error.log

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "deepcis-cwl"

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 32000

# Define namespaces
$namespaces:
  cwltool: "http://commonwl.org/cwltool#"
