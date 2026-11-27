# Dockerfile for ModelPipelineIO CWL workflows
# Uses micromamba for faster dependency resolution and smaller image size

FROM mambaorg/micromamba:2.8.1-cuda12.4.1-ubuntu22.04

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

# Create and activate the PlantCAD environment
COPY --chown=$MAMBA_USER:$MAMBA_USER cwl/plantcad/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a --yes

# Install additional requirements
COPY --chown=$MAMBA_USER:$MAMBA_USER cwl/plantcad/additional-requirements.txt /tmp/additional-requirements.txt
RUN micromamba run -n base pip install -r /tmp/additional-requirements.txt --no-build-isolation

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Copy the PlantCAD folder and set the working directory to it
COPY --chown=$MAMBA_USER ./plantcad /app/plantcad
WORKDIR /app/plantcad

# Set the path to include the micromamba directory
ENV PATH=/opt/mamba/bin:$PATH
