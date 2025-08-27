# How to Contribute to ITACaRT

First off, thank you for considering contributing to ITACaRT! Every contribution is welcome and greatly appreciated.

This document is a guide to help you get involved with the project.

## Reporting Bugs and Suggesting Features

The best way to report a bug or suggest a new feature is to open an **Issue** on our GitHub repository.

When opening an issue, please provide as much detail as possible, including:
* A clear and concise description of the problem or suggestion.
* Steps to reproduce the bug.
* The expected behavior and what actually happened.
* The version of `itacart` and Python you are using.

## Setting Up the Development Environment

To contribute code, you will need to set up your local environment.

1.  **Fork** the repository to your personal GitHub account.
2.  **Clone** your fork to your local machine:
    ```bash
    git clone [https://github.com/YOUR-USERNAME/itacart.git](https://github.com/YOUR-USERNAME/itacart.git)
    cd itacart
    ```
3.  **Create a Virtual Environment** (recommended):
    ```bash
    python -m venv venv
    # On Windows (Git Bash)
    source venv/Scripts/activate
    # On Linux/macOS
    # source venv/bin/activate
    ```
4.  **Install Dependencies** for development:
    ```bash
    python -m pip install --upgrade pip
    # Installs the package in editable mode and dev dependencies
    pip install -e .[dev]
    ```
    *(Note: We will need to define a `[project.optional-dependencies]` table in `pyproject.toml` for the `[dev]` dependencies, including `pytest`, `black`, `isort`, etc.)*

5.  **Set up the pre-commit hooks** to ensure automatic code formatting:
    ```bash
    python -m pre_commit install
    ```

## Pull Request Workflow

1.  Create a new **branch** from `main` (or `develop`) for your changes:
    ```bash
    git checkout -b feature/my-new-feature
    ```
2.  Make your code changes.
3.  **Write tests** for your new functionality in the `tests/` directory.
4.  **Run tests** locally to ensure everything is working correctly:
    ```bash
    pytest
    ```
5.  **Commit** your changes. Thanks to `pre-commit`, your code will be automatically formatted.
6.  **Push** your branch to your fork:
    ```bash
    git push origin feature/my-new-feature
    ```
7.  Open a **Pull Request (PR)** to the `ITACaRT/itacart` repository. In the PR description, clearly explain what your change does and why it's needed.

## Coding Standards

* We use **black** for code formatting.
* We use **isort** for sorting imports.
* Please follow the **Google Style** for docstrings.

Thank you again for your contribution!
