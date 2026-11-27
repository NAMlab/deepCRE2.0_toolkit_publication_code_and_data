#! /bin/bash

# Install using the workflows folder as reference folder
cd workflows  && \

# Install the DeepCRE docker container
docker build -t deepcre-cwl -f cwl/deepcre/Dockerfile .  && \

# Install the DeepCis docker container
docker build -t deepcis-cwl -f cwl/deepcis/Dockerfile .  && \

# Install the PlantCAD docker container
docker build -t plantcad-cwl -f cwl/plantcad/Dockerfile .  && \

# Install the PlantCAD docker container
docker build -t phytoexpr-cwl -f cwl/phytoexpr/Dockerfile .  && \

# Install the PanGenomeTools docker container
cd PanGenomeTools  && \
docker build -t pangenometools-cwl .  && \

# Install the ModelPipelineIO docker container
cd ../ModelPipelineIO  && \
docker build -t modelpipelineio-cwl .  && \

echo "Build finished successfully"
