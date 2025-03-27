# Docker Configuration for PyMOL Docking

This directory contains the Docker configuration for the PyMOL Docking server.

## Contents

- `Dockerfile`: Configuration for building the Docker image
- `docker-compose.yml`: Compose file for easy deployment
- `environment.yml`: Conda environment specification for dependencies
- `data/`: Directory mounted as a volume for persistent data storage

## Building and Running

### Using Docker Compose (Recommended)

```bash
# Build and start the container
docker-compose up -d

# Check if the container is running
docker ps | grep pymol-docking-server

# Stop the container
docker-compose down
```

### Using Docker Directly

```bash
# Build the image
docker build -t pymol-docking-server -f Dockerfile ..

# Run the container
docker run -p 5000:5000 -v $(pwd)/data:/app/data pymol-docking-server
```

## Troubleshooting

### Common Issues

1. **"environment.yml not found" error**:
   - Make sure you're running the commands from the `docker` directory
   - Check that the `environment.yml` file exists

2. **Path-related issues**:
   - The Docker build context is set to the parent directory (`..`)
   - All paths in the Dockerfile are relative to the build context

3. **Windows-specific issues**:
   - On Windows, use forward slashes (`/`) in paths
   - Use WSL2 backend for Docker Desktop
   - For volume mounts, you may need to use absolute paths

### Testing the Container

To check if the container is running correctly:

```bash
# Health check
curl http://localhost:5000/health

# Expected response: {"status":"healthy"}
``` 