# Contributing to the Oscheius_PDE Project

This document outlines the best practices for contributing to this genomics research project.

## Code Style

All R code should follow the [tidyverse style guide](https://style.tidyverse.org/), which includes:
- Use spaces, not tabs
- 2 spaces for indentation
- Max 80 characters per line
- Use snake_case for variable and function names
- Use descriptive variable names

## Git Workflow

1. Create feature branches for new analyses or features
2. Use meaningful commit messages that describe what changes were made
3. Pull and merge from the main branch regularly
4. Create pull requests for review before merging into the main branch

## Documentation

- Document all functions with roxygen2-style comments
- Include examples in function documentation
- Document data processing steps and decisions
- Update README.md when adding new features or analyses

## Data Management

- Never commit raw data to the repository if files are large
- Keep a separate document describing where raw data is stored and how to access it
- Document all data processing steps clearly
- Include metadata for all datasets

## Package Dependencies

- All R package dependencies must be managed through renv
- After adding a new package dependency:
  ```r
  renv::snapshot()
  ```
- When pulling changes with new dependencies:
  ```r
  renv::restore()
  ```

## Testing

- Create test files in the `tests/` directory
- Test all functions with representative data
- Use testthat for testing R functions
- Run tests before committing changes

## Code Review

All code should be reviewed by at least one other team member before being merged into the main branch. Code reviews should focus on:
- Correctness
- Clarity
- Reproducibility
- Performance
- Documentation

## Questions

If you have any questions about contributing, please contact Pablo Gonzalez de la Rosa (pgonzale60@gmail.com).
