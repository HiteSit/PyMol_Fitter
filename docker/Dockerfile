# Start with the latest micromamba image
FROM mambaorg/micromamba:2.0.8

USER root
RUN apt-get update && apt-get install -y \
    libgl1-mesa-glx \
    libglu1-mesa \
    libxrender1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
USER mambauser

# Set environment variables to ensure conda environment activation during build
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Copy environment file
COPY --chown=mambauser:mambauser docker/environment.yml /tmp/environment.yml

# Install dependencies using conda/mamba
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Cleanup
RUN micromamba clean --all --yes && \
    # Remove package cache
    rm -rf /opt/conda/pkgs/* && \
    # Remove documentation
    rm -rf /opt/conda/share/doc/* && \
    rm -rf /opt/conda/share/man/* && \
    # Remove tests
    find /opt/conda -name "tests" -type d -exec rm -rf {} +; true && \
    find /opt/conda -name "test" -type d -exec rm -rf {} +; true && \
    # Remove unnecessary Python cache files
    find /opt/conda -name "__pycache__" -type d -exec rm -rf {} +; true && \
    find /opt/conda -name "*.pyc" -delete && \
    # Remove examples
    find /opt/conda -name "examples" -type d -exec rm -rf {} +; true && \
    # Remove unused files from large packages
    rm -rf /opt/conda/lib/python*/site-packages/rdkit/Docs && \
    rm -rf /opt/conda/lib/python*/site-packages/rdkit/Contrib && \
    rm -rf /opt/conda/lib/python*/site-packages/pymol/data/movies

# Set working directory for your project
WORKDIR /app

# Create the source directory
RUN mkdir -p /app/pymol_docking_src

# Copy the server computational code
COPY --chown=mambauser:mambauser ../pymol_docking_server/pymol_docking_src/*.py /app/pymol_docking_src/

# Copy the Flask application
COPY --chown=mambauser:mambauser ../pymol_docking_server/app.py /app/

# Create a data directory for output files
RUN mkdir -p /app/data
VOLUME /app/data

# Expose port for Flask
EXPOSE 5000

# The default entrypoint in the micromamba image automatically 
# activates the base environment
CMD ["gunicorn", "--timeout", "600", "--bind", "0.0.0.0:5000", "app:app"]