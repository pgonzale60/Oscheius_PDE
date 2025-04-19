# Oscheius_PDE Research Workflow

## Data Management Protocol

### Raw Data
1. All raw data should be stored in `data/raw/` 
2. Raw data files should never be modified
3. Document the source and acquisition method for each raw data file
4. Use descriptive file names with dates (YYYY-MM-DD format)

### Data Processing
1. All data processing scripts should be in the `analysis/` directory
2. Number scripts in order of execution (e.g., `01_data_processing.R`, `02_statistical_analysis.R`)
3. Processed data goes in `data/processed/`
4. Each processing step should be documented in the script headers

### Analysis
1. Store analysis scripts in the `analysis/` directory
2. Store results in the `results/` directory
3. Generate figures in the `figures/` directory
4. Use consistent naming conventions for all output files

## Reproducibility Guidelines

### R Environment
- Use renv for package management
- Use `renv::snapshot()` after adding new packages
- Use `renv::restore()` when cloning the repository

### Documentation
- Document all code with clear comments
- Include session information at the end of each script
- Maintain a lab notebook or research log for key decisions

### Version Control
- Commit regularly with meaningful commit messages
- Use branches for experimental features
- Never commit large data files to Git

## Publication Preparation
1. Create publication-ready figures in the `figures/` directory
2. Document statistical methods thoroughly
3. Ensure all data processing steps are reproducible
