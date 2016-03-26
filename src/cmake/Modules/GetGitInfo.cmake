# This script finds a description of the last git commit (GIT_COMMIT_ID) and the current branch (GIT_BRANCH_NAME).
# It then calls configure_file to convert to provided template file (INFILE) into the desired output (OUTFILE).
# The script can be run at build time (e.g., through add_custom_target) to ensure the git info is always up to date.

find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --always --dirty
                  OUTPUT_VARIABLE GIT_COMMIT_ID
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET
                 )
  execute_process(COMMAND ${GIT_EXECUTABLE} name-rev --name-only HEAD
                  OUTPUT_VARIABLE GIT_BRANCH_NAME
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET
                 )
endif()
configure_file("${INFILE}" "${OUTFILE}")
