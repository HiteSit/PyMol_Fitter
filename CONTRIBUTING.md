# Contributing to PyMOL Docking Docker

Thank you for your interest in contributing to the PyMOL Docking Docker project! This guide will help you get started with contributing to this project.

## Development Setup

To set up the development environment:

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/pymol-docking-docker.git
   cd pymol-docking-docker
   ```

2. Install the server dependencies:
   ```bash
   cd pymol_docking_server
   pip install -r requirements.txt
   ```

3. Install the PyMOL plugin development requirements:
   ```bash
   pip install requests pathlib
   ```

## Project Structure

- `pymol_docking_plugin/`: Client-side PyMOL plugin
- `pymol_docking_server/`: Server-side Flask API and computational code
- `docker/`: Docker configuration files

## Adding Features

### Server-Side

1. To add new computational features:
   - Add your code to the appropriate module in `pymol_docking_server/pymol_docking_src/`
   - Update or extend the Flask API in `pymol_docking_server/app.py` to expose the new functionality

2. Update Docker configuration if needed:
   - If you add new dependencies, update `docker/environment.yml`
   - If you add new files, update the `docker/Dockerfile`

### Client-Side

1. To add new functionality to the PyMOL plugin:
   - Add new methods to `pymol_docking_plugin/client.py` to call the server API
   - Update the PyMOL commands in `pymol_docking_plugin/__init__.py`
   - Update the GUI in `pymol_docking_plugin/GUI.ui` if needed

## Testing

Please test your changes thoroughly:

1. Server functionality:
   - Start the Docker container with your changes
   - Test the API endpoints directly

2. PyMOL plugin:
   - Install the updated plugin in PyMOL
   - Test the new/modified functionality

## Submitting Changes

1. Create a branch for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Commit your changes with descriptive commit messages:
   ```bash
   git commit -m "Add feature X to handle Y"
   ```

3. Push your branch:
   ```bash
   git push origin feature/your-feature-name
   ```

4. Create a pull request detailing:
   - What changes you've made
   - Why they are valuable
   - How they have been tested

## Code Style

- Follow PEP 8 for Python code
- Add docstrings to all functions and classes
- Use type hints where appropriate

## Future Development

Here are some areas we're looking to improve:

- [ ] Add PoseBusters integration for pose validation
- [ ] Add manual conformer docking capabilities
- [ ] Improve error handling and recovery
- [ ] Enhance visualization options in PyMOL

We welcome contributions in these areas or any others that would improve the project! 