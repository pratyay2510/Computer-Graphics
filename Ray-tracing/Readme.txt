# Setting Up and Running the Project with SCons

This guide provides step-by-step instructions on how to set up a virtual environment, install SCons, build the project, and run the grading script.

## Prerequisites
- Python (version 3.x recommended)
- Pip package manager

## Installation and Setup

### Step 1: Create a Virtual Python Environment
Before installing SCons, it is recommended to create an isolated Python environment.

Run the following command in the project directory:

"python -m venv env_name"

This will create a virtual environment named `venv`.

### Step 2: Activate the Virtual Environment
To activate the virtual environment:

- **On Windows (cmd or PowerShell):**
" env_name/Scripts/activate"

- **On macOS/Linux:**
" source env_name/bin/activate


### Step 3: Install SCons
Once the virtual environment is activated, install `scons` using `pip`:

"pip install scons"


### Step 4: Build the Project Using SCons
Navigate to the project directory and run:

"scons"


This will start the build process as per the `SConstruct` configuration.

### Step 5: Run the Grading Script
After building the project, you can run the grading script to evaluate the submission:

"python grading-script.py a"

