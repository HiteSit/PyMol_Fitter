# Start with the latest micromamba image
FROM mambaorg/micromamba:2.0.8

# Set environment variables to ensure conda environment activation during build
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Copy environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install OpenMM, OpenFF and dependencies
# This installs into the base environment, which is the recommended approach
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes
    
# Set working directory for your project
WORKDIR /app

# Copy your application code (directory)
COPY --chown=$MAMBA_USER:$MAMBA_USER ../pymol_docking_plugin/CDPK_Utils.py /app/pymol_docking_src/
COPY --chown=$MAMBA_USER:$MAMBA_USER ../pymol_docking_plugin/Docking_Engine.py /app/pymol_docking_src/
COPY --chown=$MAMBA_USER:$MAMBA_USER ../pymol_docking_plugin/Protein_Minimization.py /app/pymol_docking_src/
COPY --chown=$MAMBA_USER:$MAMBA_USER ../pymol_docking_plugin/Protein_Preparation.py /app/pymol_docking_src/

# The default entrypoint in the micromamba image automatically 
# activates the base environment
CMD ["python"]