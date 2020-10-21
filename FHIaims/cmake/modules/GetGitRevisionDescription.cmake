# - Returns a version string from Git
#
# These functions force a re-configure on each git commit so that you can
# trust the values of the variables in your build system.
#
#  get_git_head_revision(<refspecvar> <hashvar> [<additional arguments to git describe> ...])
#
# Returns the refspec and sha hash of the current head revision
#
# Requires CMake 2.6 or newer (uses the 'function' command)
#
# Original Author:
# 2009-2010 Ryan Pavlik <rpavlik@iastate.edu> <abiryan@ryand.net>
# http://academic.cleardefinition.com
# Iowa State University HCI Graduate Program/VRAC
#
# Copyright Iowa State University 2009-2010.
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

if(__get_git_revision_description)
  return()
endif()
set(__get_git_revision_description YES)

# We must run the following at "include" time, not at function call time,
# to find the path to this module rather than the path to a calling list file
get_filename_component(_gitdescmoddir ${CMAKE_CURRENT_LIST_FILE} PATH)

function(get_git_head_revision _hashvar)
  set(GIT_PARENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  set(GIT_DIR "${GIT_PARENT_DIR}/.git")
  while(NOT EXISTS "${GIT_DIR}")	# .git dir not found, search parent directories
    set(GIT_PREVIOUS_PARENT "${GIT_PARENT_DIR}")
    get_filename_component(GIT_PARENT_DIR ${GIT_PARENT_DIR} PATH)
    if(GIT_PARENT_DIR STREQUAL GIT_PREVIOUS_PARENT)
      # We have reached the root directory, we are not in git
      set(${_refspecvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
      set(${_hashvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
      return()
    endif()
    set(GIT_DIR "${GIT_PARENT_DIR}/.git")
  endwhile()
  # check if this is a submodule
  if(NOT IS_DIRECTORY ${GIT_DIR})
    file(READ ${GIT_DIR} submodule)
    string(REGEX REPLACE "gitdir: (.*)\n$" "\\1" GIT_DIR_RELATIVE ${submodule})
    get_filename_component(SUBMODULE_DIR ${GIT_DIR} PATH)
    get_filename_component(GIT_DIR ${SUBMODULE_DIR}/${GIT_DIR_RELATIVE} ABSOLUTE)
  endif()
  set(GIT_DATA "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/git-data")
  if(NOT EXISTS "${GIT_DATA}")
    file(MAKE_DIRECTORY "${GIT_DATA}")
  endif()

  if(NOT EXISTS "${GIT_DIR}/HEAD")
    return()
  endif()
  set(HEAD_FILE "${GIT_DATA}/HEAD")
  configure_file("${GIT_DIR}/HEAD" "${HEAD_FILE}" COPYONLY)

  configure_file("${_gitdescmoddir}/GetGitRevisionDescription.cmake.in"
    "${GIT_DATA}/grabRef.cmake"
    @ONLY)
  include("${GIT_DATA}/grabRef.cmake")

  string(SUBSTRING "${HEAD_HASH}" 0 9 HEAD_HASH)
  set(${_hashvar} "${HEAD_HASH}" PARENT_SCOPE)
endfunction()
