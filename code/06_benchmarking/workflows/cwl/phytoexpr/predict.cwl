#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

# Metadata for the tool
label: "PhytoExpr Prediction"
doc: |
  Performs phytoexpr prediction for chunked sequences.

  This tool processes a Fasta file and produces a CSV file with PhytoExpr predictions.

# Base command - using shell to chain two Python scripts
baseCommand: [python, /app/predict.py]

# Input parameters
inputs:
  input_file:
    type: File
    inputBinding:
      prefix: --input-fasta
      position: 1
    doc: "Input FASTA file"
  models:
    type:
    - File
    - Directory
    inputBinding:
      prefix: --models
      position: 2
    doc: "Path to the trained PhytoExpr model"
  index_path:
    type: string
    inputBinding:
      prefix: --index-path
    default: "local.fai"
    doc: "Path to the (created) FASTA fai file (by default a local file)"

# Output files
outputs:
  predictions:
    type: File
    outputBinding:
      glob: "*.csv"
    doc: "Generated CSV file with prediction results"

# Standard output and error handling
stdout: predict_phytoexpr.log
stderr: predict_phytoexpr.error.log

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "phytoexpr-cwl"

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 32000

# Define namespaces
$namespaces:
  cwltool: "http://commonwl.org/cwltool#"
