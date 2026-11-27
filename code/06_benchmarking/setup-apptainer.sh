# Build PanGenomeTools container
cd workflows/PanGenomeTools  && \
apptainer build --force ../cwl/containers/pangenometools-cwl_latest.sif Container.def  && \

# Build ModelPipelineIO container
cd ../ModelPipelineIO  && \
apptainer build --force ../cwl/containers/modelpipelineio-cwl_latest.sif Container.def  && \

# Build DeepCis container
cd ..  && \
apptainer build --force --nv cwl/containers/deepcis-cwl_latest.sif cwl/deepcis/Container.def  && \

# Build DeepCRE container
apptainer build --force --nv cwl/containers/deepcre-cwl_latest.sif cwl/deepcre/Container.def  && \

# Build PlantCAD container
apptainer build --force --nv cwl/containers/plantcad-cwl_latest.sif cwl/plantcad/Container.def  && \

# Build PhytoExpr container
apptainer build --force --nv cwl/containers/phytoexpr-cwl_latest.sif cwl/phytoexpr/Container.def  && \

echo "Build finished successfully"
