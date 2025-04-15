# PyMol Fitter Test Suite

This directory contains comprehensive tests for the PyMol Fitter plugin, divided into two main categories:

1. **Scientific Tests**: Focus on testing the core scientific functionality
2. **Application Tests**: Test the API endpoints, client-server interactions, and Docker integration

## Test Structure

```
tests/
├── conftest.py                 # Global pytest fixtures
├── requirements-test.txt       # Test dependencies
├── data/                       # Test data directory
├── scientific_tests/           # Scientific functionality tests
│   ├── test_cdpk_utils.py      # Tests for CDPK utilities
│   ├── test_docking_engine.py  # Tests for docking engine
│   ├── test_protein_minimization.py  # Tests for protein minimization
│   └── test_protein_preparation.py   # Tests for protein preparation
└── application_tests/          # API and integration tests
    ├── conftest.py             # Flask-specific fixtures
    ├── test_flask_api.py       # Tests for Flask API endpoints
    └── test_client.py          # Tests for PyMOL client

```

## Installation

First, install testing dependencies:

```bash
pip install -r tests/requirements-test.txt
```

### Required Dependencies

These tests require several specialized chemical and computational biology libraries:

```bash
# Core test dependencies
pip install pytest pytest-cov

# Python libraries for API and client tests
pip install flask requests 

# Scientific libraries for the core functionality tests
pip install biopython biotite rdkit datamol
conda install -c conda-forge pdbfixer  # Not available via pip

# CDPL library may need to be installed separately
# See https://cdpl-team.github.io/CDPL/download.html for instructions
```

### Handling Import Errors

If you encounter import errors like `No module named 'pymol_fitter_src'`, you might need to fix import paths in:

1. `pymol_fitter_server/app.py`: Change imports from:
   ```python
   from pymol_fitter_src.Docking_Engine import ...
   ```
   To:
   ```python
   from pymol_fitter_server.pymol_fitter_src.Docking_Engine import ...
   ```

2. Or add the project root to your Python path:
   ```bash
   export PYTHONPATH=/path/to/PyMol_Fitter:$PYTHONPATH
   ```

## Running Tests

### Running All Tests

```bash
# From the project root directory
pytest tests/
```

### Running Only Scientific Tests

```bash
pytest tests/scientific_tests/
```

### Running Only Application Tests

```bash
pytest tests/application_tests/
```

### Running Specific Test Files

```bash
pytest tests/scientific_tests/test_docking_engine.py
```

### Running Tests with Coverage

```bash
pytest --cov=pymol_fitter_server --cov=pymol_fitter_plugin tests/
```

To generate an HTML coverage report:

```bash
pytest --cov=pymol_fitter_server --cov=pymol_fitter_plugin --cov-report=html tests/
```

## Test Categories and Markers

Tests are marked with categories to allow selective execution:

- **slow**: Tests that take a long time to run (e.g., actual docking)

Run only fast tests:

```bash
pytest -k "not slow" tests/
```

Run only slow tests:

```bash
pytest -k "slow" tests/
```

## Testing the Docker Environment

Some tests require the Docker container to be running. To run these tests:

1. Start the Docker container:
   ```bash
   cd docker
   docker-compose --profile cpu up -d
   ```

2. Run tests with Docker integration:
   ```bash
   pytest tests/application_tests/
   ```

You can set a custom Docker endpoint by setting an environment variable:

```bash
export PYMOL_FITTER_DOCKER_URL=http://localhost:5000
pytest tests/application_tests/
```

## Adding New Tests

### Adding Scientific Tests

1. Create a new test file in `tests/scientific_tests/`
2. Import the relevant modules from `pymol_fitter_server.pymol_fitter_src`
3. Use the fixtures from `conftest.py` as needed
4. Mark slow tests with `@pytest.mark.slow`

### Adding Application Tests

1. Create a new test file in `tests/application_tests/`
2. Use the Flask client fixture from `application_tests/conftest.py`
3. Mock external dependencies for faster tests

## Testing with Different Example Files

To add new example files for testing:

1. Place the files in the `tests/data/` directory
2. Create a new fixture in `conftest.py` to access these files

For example, to use test data in your tests:

```python
# In your test file
def test_something(data_dir):
    """Test using data from the data directory."""
    protein_file = data_dir / "8gcy.pdb"
    ligand_file = data_dir / "Ligand.sdf"
    
    # Use these files in your test
    assert protein_file.exists()
    # ...
```

The current test data directory contains:
- Protein files: 8gcy.pdb, LAC3.pdb
- Ligand files: Ligand.sdf, Lig_Min.sdf
- Crystal files: Crystal.sdf, 8gcy_Crystal.sdf

Use these existing files when writing new tests to ensure consistency.

## Troubleshooting

### Missing Dependencies

If you encounter errors about missing Python modules, make sure all dependencies listed in `requirements-test.txt` are installed.

### Docker-Related Test Failures

If Docker-related tests fail:

1. Check that the Docker container is running: `docker ps`
2. Check the Docker container logs: `docker logs pymol-fitter`
3. Verify the API is accessible: `curl http://localhost:5000/health`

### Scientific Module Import Errors

Some scientific modules may require specific libraries that are only available in the Docker environment. Consider running these tests within the Docker container.
