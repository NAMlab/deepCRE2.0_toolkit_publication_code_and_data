# Dockerfile for ModelPipelineIO CWL workflows
# Uses micromamba for faster dependency resolution and smaller image size

FROM mambaorg/micromamba:2.8.1-cuda11.7.1-ubuntu22.04

USER root

# Basic tools
RUN apt-get update && apt-get install -y \
    build-essential  \
    git \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
# Optional: if you want to customize the ENTRYPOINT and have a conda
# environment activated, then do this:
# ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "my_entrypoint_program"]

# You can modify the CMD statement as needed....
CMD ["/bin/bash"]

# Create and activate the deepis environment
COPY --chown=$MAMBA_USER:$MAMBA_USER cwl/deepcis/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install the DeepCis package in development mode
# WORKDIR /app
# COPY --chown=$MAMBA_USER ./deepcis /app/deepcis
# WORKDIR /app/deepcis

# RUN pip install -e .

# Copy the prediction script
WORKDIR /app
COPY --chown=$MAMBA_USER ./cwl/deepcis/predict.py /app/predict.py

# Install the ModelPipelineIO package
WORKDIR /app
COPY --chown=$MAMBA_USER ./ModelPipelineIO /app/ModelPipelineIO
WORKDIR /app/ModelPipelineIO

RUN pip install .

# Set the working director to the cwl/deepcis folder
COPY --chown=$MAMBA_USER ./cwl/deepcis /app/cwl/deepcis
WORKDIR /app/cwl/deepcis

# Set the path to include the micromamba directory
ENV PATH=/opt/mamba/bin:$PATH
