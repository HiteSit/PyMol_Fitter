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
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "app:app"] 