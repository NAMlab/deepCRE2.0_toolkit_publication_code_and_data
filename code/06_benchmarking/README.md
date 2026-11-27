# DeepCRE model benchmarking

This repository applies three published models for predicting gene expression based on DNA sequence to genomes of *Arabidopsis thaliana* and *Zea mays*.
Each model is integrated in workflows with standardized inputs and outputs using the Common Workflow Language (CWL, <https://www.commonwl.org/>).

This benchmarking was performed by Tobias Pfennig and is based on the yet unpublished ModelPipeline (<https://github.com/pfennigt/ModelPipeline>)

| Models | Publication |
| --- | --- |
| DeepCRE | <https://doi.org/10.1038/s41467-024-47744-0> |
| PlantCaduceus | <https://doi.org/10.1038/s41467-024-47744-0> |
| PhytoExpr | <https://doi.org/10.1073/pnas.2319811121> |

## Requirements

- ARCitect or ARCcommander (<https://nfdi4plants.github.io/nfdi4plants.knowledgebase/arcitect/installation/>)
- git
- conda or mamba
- Access to the following (non-public) GitHub repositories:
  - <https://github.com/pfennigt/ModelPipelineIO>
- NVIDIA graphics card ([Not required for all models](#models)):

| | CUDA | [Compute capability](<https://developer.nvidia.com/cuda/gpus>) (CC) | Notes |
| --- | --- | --- | --- |
| **Min** | 11.7 | 3.5 | [Cannot run all models](#models) |
| **Recommended** | 12.4 | 8.0 | |

## Installation

Download this repository from <https://git.nfdi4plants.org/szymanskilab/ModelPipeline_DeepCRETK> using the ARCitect.
Navigate into the created project folder (ModelPipeline_DeepCRETK).

If you have git setup with SSH access ignore this paragraph, otherwise open the file `.gitmodules`.
Replace every instance of `git@github.com:` with `https://github.com/`.
You will be prompted for you GitHub username and password, *but only a personal access token will be accepted as your password* ([GitHub: personal access tokens](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-personal-access-token-classic)).

Install the required submodules using

```bash
git submodule init
git submodule update
```

Furthermore, cwltool and either Docker or Apptainer (v1) / Singularity (v2/v3) are required to run this pipeline.

### Install CWLtool

Install cwltool using conda/mamba:

```bash
conda env create -f environment.yml
OR
conda create -n modelpipeline python==3.11 pip cwltool

conda activate mpl
```

### Install Apptainer (recommended)

Install apptainer following the instructions at <https://apptainer.org/docs/admin/main/installation.html>

Create the containers as flat SIF files by running (this may take a while)

```bash
bash setup-apptainer.sh
```

### Install Docker (alternative to Apptainer)

Install Docker following the instructions at <https://docs.docker.com/>.

The CWL runs are executed in Docker containers for reproducibility and avoiding package conflicts.
Create the relevant containers by running (this may take a while)

```bash
bash setup-docker.sh
```

#### Install NVIDIA container toolkit

When using Docker, install the nvidia-container-toolkit

```bash
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
  && curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

sudo apt-get update

export NVIDIA_CONTAINER_TOOLKIT_VERSION=1.19.1-1
  sudo apt-get install -y \
      nvidia-container-toolkit=${NVIDIA_CONTAINER_TOOLKIT_VERSION} \
      nvidia-container-toolkit-base=${NVIDIA_CONTAINER_TOOLKIT_VERSION} \
      libnvidia-container-tools=${NVIDIA_CONTAINER_TOOLKIT_VERSION} \
      libnvidia-container1=${NVIDIA_CONTAINER_TOOLKIT_VERSION}
```

Configure the container runtime by using the nvidia-ctk command:

```bash
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker
```

## Models

The following models are available to run in `runs/<run_name>`.

| Model | Run name | Description  | GPU |
| --- | --- | --- | --- |
| DeepCRE | `deepcre` | Per-gene expression | (CUDA >= 11.7, CC >= 3.5)¹ |
| PlantCaduceus | `plantcad_loraexpression` | Per-gene expression | CUDA >= 12.4, CC >= 8.0 |
| PhytoExpr | `phytoexpr` | Per-gene expression | (CUDA >= 11.7, CC >= 3.5)¹ |

¹ GPU not used by default.

## Re-run analysis

The manuscript Figure 2 and Supplementary Figure 4 are created in the jupyter notebook `runs/analysis-deepcreComparison/analysis.ipynb`.
To recreate the figures, navigate into `runs/analysis-deepcreComparison` and install the conda environment with

```bash
conda env create -f environment.yml
```

Then start the notebook by running

```bash
conda activate analysis-deepcreComparison

jupyter notebook
```

A Jupyter interface will open in your browser.
There you can start the notebook by clicking on `analysis.ipynb`.
Re-run the analysis by running each cell (e.g. by `Shift + Enter`), figures will be output into the `figures` folder.

## Re-run predictions

The `run` directories contain the computed predictions of each model, as they are shown in the manuscript, in the respective `dataset` folder.
To rerun the predictions, follow the instructions below.

### Running workflows

To run the CWL workflows, navigate to the respective folder in `runs`.
You have the option to use either Apptainer/Singularity (recommended) or Docker to run the workflows.

Aliases for the commands below have been implemented to ease usage.
Run

```bash
source setup-run-aliases.sh
```

> [!NOTE]
> You have to run `source setup-run-aliases.sh` each time you open a new shell

And then move into the `runs` directory of a model.
Then you have access to the following commands in your current shell:

```bash
# Run the current workflow with Docker
mpl-run

# Run the current workflow with Docker and parallelisation
mpl-run --parallel
mpl-run -p

# cwltool commands can also be used to adapt the run
mpl-run -p --parallel-max 5 # Maximum 5 parallel step executions (does not limit the number of CPUs per execution)

# Run the current workflow with Apptainer / Singularity
mpl-run --singularity
mpl-run -s

# Run the current workflow with Apptainer / Singularity and parallelisation
mpl-run --singularity --parallel
mpl-run -s -p
mpl-run -sp

# Delete the contents of the dataset and logs folders (asks for confirmation)
mpl-run reset
```

The predictions will be output into the runs' `dataset` folder.
One file is created each combination of model, prediction track, build and genotype.
The accompanying `index.csv` file lists all files in that folder whith their respective prediction.

The `logs` folder contains any logs, print statements, etc. produced by the workflow steps.
The overall status of a run is logged in `logs/run/run.error.log`

#### Running with Apptainer / Singularity (recommended)

To use apptainer for running, you first have to set the following environment variables from the projects root directory:

```bash
export SINGULARITY_CACHEDIR=../../workflows/cwl/containers
export SINGULARITY_TMPDIR=$SINGULARITY_CACHEDIR/tmp
export SINGULARITY_PULLDIR=$SINGULARITY_CACHEDIR/pull
export CWL_SINGULARITY_CACHE=$SINGULARITY_CACHEDIR
```

Then use the above commands with the flags `--singularity --nv`.
Without parallelisation:

```bash
cwltool --singularity --outdir dataset --log-dir logs --timestamps --disable-color --enable-ext run.cwl run.yml > logs/run/run.log 2> logs/run/run.error.log
```

With parallelisation:

```bash
cwltool --singularity --parallel --strict-memory-limit --outdir dataset --log-dir logs --timestamps --disable-color --enable-ext run.cwl run.yml  > logs/run/run.log 2> logs/run/run.error.log
```

#### Running with Docker

To use docker for running, use

```bash
cwltool --outdir dataset --log-dir logs --timestamps --disable-color --enable-ext run.cwl run.yml  > logs/run/run.log 2> logs/run/run.error.log
```

Parallelisation can be employed using

```bash
cwltool --parallel --strict-memory-limit --outdir dataset --log-dir logs --timestamps --disable-color --enable-ext run.cwl run.yml  > logs/run/run.log 2> logs/run/run.error.log
```
