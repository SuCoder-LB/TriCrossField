cmake_minimum_required(VERSION 3.25)
project(TriCrossField)

set(CMAKE_CXX_STANDARD 14)

include_directories(extern extern/vcglib)

set(DIRS ${DIRS};extern)
set(DIRS ${DIRS};extern/vcglib)

file(GLOB_RECURSE file_h *.h  extern/*.h)

foreach(DIR ${DIRS})
    string(REGEX REPLACE "\\+" "\\\\+" DIR ${DIR})
    source_group(${DIR} REGULAR_EXPRESSION ${DIR}/.*)
endforeach()

add_executable(${CMAKE_PROJECT_NAME}
        tcf_internal.cpp
        tcf_linear_system.cpp
        demo.cpp
        extern/vcglib/wrap/ply/plylib.cpp
        ${file_h}
        )


if(CYGWIN OR MSYS)
    # detect if we use the MinGW compilers on Cygwin - if we do, handle the build
    # as a pure Windows build and make cmake find pure Windows import libraries
    # (.lib)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
            CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpmachine
                OUTPUT_VARIABLE CXX_COMPILER_MACHINE
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(CXX_COMPILER_MACHINE MATCHES "mingw")
            set(WIN32 1)
            add_definitions(-DWIN32)
        endif()
    endif()
endif()


if(MSVC)
    add_definitions(/bigobj -DNOMINMAX -D_CRT_SECURE_NO_DEPRECATE -D_CRT_SECURE_NO_WARNINGS
            -D_SCL_SECURE_NO_DEPRECATE)
    add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/source-charset:utf-8>")
    #    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2")
    #    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MD")
    #    string(PREPEND "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif()


if(WIN32)
    add_definitions(-D_USE_MATH_DEFINES)
endif()