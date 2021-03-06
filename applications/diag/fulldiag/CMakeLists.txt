#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LAPACK_FOUND)
  add_definitions(${LAPACK_DEFINITIONS})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_library(fulldiag_impl STATIC factory.C)
  target_link_libraries(fulldiag_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

  add_executable(fulldiag fulldiag.C)
  add_executable(fulldiag_evaluate fulldiag_evaluate.C)
  target_link_libraries(fulldiag fulldiag_impl)
  target_link_libraries(fulldiag_evaluate fulldiag_impl)
  install(TARGETS fulldiag fulldiag_evaluate RUNTIME DESTINATION bin COMPONENT applications)
else(LAPACK_FOUND)
  message(STATUS "fulldiag will not be built since lapack library has not been found")
endif(LAPACK_FOUND)
