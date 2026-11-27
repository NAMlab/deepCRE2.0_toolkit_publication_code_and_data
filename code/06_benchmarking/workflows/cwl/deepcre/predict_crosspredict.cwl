#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

# Metadata for the tool
label: "DeepCRE Cross Prediction"
doc: |
  Runs DeepCRE cross-species predictions using trained models.

  This tool executes the input_preparation.py script to prepare the input JSON configuration file
  and then runs the deepcre_crosspredict.py script to generate prediction results in CSV format.

# Base command - using shell to chain two Python scripts
baseCommand: sh

arguments:
  - -c
  - |
    python /app/ModelPipelineIO/scripts/common/replace_in_file.py \
      --input-file $(inputs.json_config.path) \
      --output-file $(inputs.json_path) \
      --placeholder $(inputs.pangenome_path_placeholder) '"$(inputs.prediction_models_placeholder)"' \
      --target $(inputs.pangenome_path.path) '"$(inputs.prediction_models.map(f => f.path).join('","'))"' && \
    python /app/deepcre/src/deepCRE/deepcre_crosspredict.py \
      --input $(inputs.json_path)

# Input parameters
inputs:
  pangenome_index_path:
    type: File
    doc: "Pan-genome index CSV file"
  pangenome_path:
    type: Directory
    doc: "Pan-genome root directory path"
  pangenome_path_placeholder:
    type: string
    doc: "Placeholder for the pangenome path in JSON config file"
  json_config:
    type: File
    doc: "DeepCRE JSON configuration file path"
  json_path:
    type: string
    default: "input.json"
    doc: "Output JSON configuration file path"
  output_dir:
    type: string
    default: "."
    doc: "Directory for prediction outputs"
  prediction_models:
    type: File[]
    doc: "Prediction models to use (comma-separated)"
  prediction_models_placeholder:
    type: string
    doc: "Placeholder for the prediction model paths in JSON config file"
  assays_dir:
    type: Directory
    doc: "Path to the assays directory containing input files"

# Output files
outputs:
  prediction_results:
    type:
    - File
    - File[]
    outputBinding:
      glob: "*.csv"
    doc: "Generated prediction result CSV files"

# Standard output and error handling
stdout: deepcre_crosspredict.log
stderr: deepcre_crosspredict.error.log

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "deepcre-cwl"

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 4096
    ramMax: 16384
  # - class: cwltool:CUDARequirement
  #   cudaComputeCapability: "3.5"
  #   cudaDeviceCountMin: 1 # Can spawn one prediction per CUDA device
  #   cudaDeviceCountMax: 1
  #   cudaVersionMin: "11.7"

# Define namespaces
$namespaces:
  cwltool: "http://commonwl.org/cwltool#"
