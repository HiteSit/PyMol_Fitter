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
    
# Install Flask and additional requirements
RUN pip install flask gunicorn requests numpy

# Set working directory for your project
WORKDIR /app

# Create the source directory
RUN mkdir -p /app/pymol_docking_src

# Copy your application code (directory)
COPY --chown=$MAMBA_USER:$MAMBA_USER ../pymol_docking_server/pymol_docking_src/CDPK_Utils.py /app/pymol_docking_src/
COPY --chown=$MAMBA_USER:$MAMBA_USER ../pymol_docking_server/pymol_docking_src/Docking_Engine.py /app/pymol_docking_src/
COPY --chown=$MAMBA_USER:$MAMBA_USER ../pymol_docking_server/pymol_docking_src/Protein_Minimization.py /app/pymol_docking_src/
COPY --chown=$MAMBA_USER:$MAMBA_USER ../pymol_docking_server/pymol_docking_src/Protein_Preparation.py /app/pymol_docking_src/

# Copy the Flask application
COPY --chown=$MAMBA_USER:$MAMBA_USER ../pymol_docking_server/app.py /app/

# Expose port for Flask
EXPOSE 5000

# The default entrypoint in the micromamba image automatically 
# activates the base environment
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "app:app"]