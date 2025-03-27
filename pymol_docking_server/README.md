# PyMOL Docking Server

This directory contains the server-side code for the PyMOL Docking Docker plugin.

## Contents

- `app.py`: Flask API serving the computational functions
- `pymol_docking_src/`: Directory containing computational modules
  - `Docking_Engine.py`: Main docking functionality
  - `Protein_Preparation.py`: Protein preparation utilities
  - `Protein_Minimization.py`: Minimization functions
  - `CDPK_Utils.py`: Utility functions
- `examples/`: Example files for testing
- `requirements.txt`: Python dependencies for the server

## How It Works

The server provides a Flask API with the following endpoints:

- `/health`: Health check endpoint
- `/dock`: Perform docking
- `/minimize_complex`: Perform complex minimization
- `/dock_and_minimize`: Perform docking followed by minimization

Each endpoint accepts specific JSON payloads and returns results as base64-encoded files.

## Deployment

This code is intended to be deployed in a Docker container using the configuration in the `docker/` directory.

To run the server directly (not recommended for production):

```bash
# Install dependencies
pip install -r requirements.txt

# Run the server
python app.py
```

## Development

To add new computational features:

1. Add your code to the appropriate module in `pymol_docking_src/`
2. Expose the functionality via a new or existing endpoint in `app.py`
3. Rebuild the Docker container if you're using Docker

## Testing

You can test the API endpoints using curl:

```bash
# Health check
curl http://localhost:5000/health

# For other endpoints, you need to send appropriate JSON payloads
# See the API documentation in app.py for details
```

The `examples/` directory contains sample files you can use for testing. 