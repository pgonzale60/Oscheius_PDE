# Helper functions for tests

#' Create a mock test environment
#'
#' @param reset Whether to reset the existing test environment
#' @return A list with paths to mock files and test data
create_test_environment <- function(reset = FALSE) {
  # Create a test directory in tempdir
  test_dir <- file.path(tempdir(), "oscheius_test")
  if (reset && dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Return the test directory
  return(test_dir)
}

#' Mock the renv functions for testing
mock_renv <- function() {
  # Define mock functions
  mock_restore <- function(...) {
    message("Mock renv::restore() called")
    return(TRUE)
  }
  
  mock_snapshot <- function(...) {
    message("Mock renv::snapshot() called")
    return(TRUE)
  }
  
  # Return mocked functions
  return(list(
    restore = mock_restore,
    snapshot = mock_snapshot
  ))
}

#' Setup mock packages for testing
setup_mock_packages <- function() {
  # Define mock functions for packages that might not be installed
  if (!requireNamespace("mockery", quietly = TRUE)) {
    stub <- function(where, what, how) {
      warning("mockery not available, using simple stub")
      # Simple stub implementation
      body(where) <- substitute({
        ORIGINAL_FUN <- what
        what <- how
        RESULT <- body_of_function
        what <- ORIGINAL_FUN
        return(RESULT)
      }, list(body_of_function = body(where)))
      return(where)
    }
    assign("stub", stub, envir = globalenv())
  }
}
