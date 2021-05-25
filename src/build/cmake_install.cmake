# Install script for directory: /home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages" TYPE DIRECTORY FILES "/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/magnum" REGEX "/magnum\\/magneto\\_cuda\\.py\\ magnum\\/\\_magneto\\_cuda\\.so\\ magnum\\/magneto\\_cpu\\.py\\ magnum\\/\\_magneto\\_cpu\\.so\\ CMakeLists\\.txt$" EXCLUDE)
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/magnum" TYPE FILE FILES
    "/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto_cuda.py"
    "/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/_magneto_cuda.so"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto/cmake_install.cmake")
  INCLUDE("/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto/bindings/cmake_install.cmake")
  INCLUDE("/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto/matrix/cmake_install.cmake")
  INCLUDE("/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto/math/cmake_install.cmake")
  INCLUDE("/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto/math/conv/cmake_install.cmake")
  INCLUDE("/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto/math/conv/kernels/cmake_install.cmake")
  INCLUDE("/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto/mmm/cmake_install.cmake")
  INCLUDE("/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto/evolver/cmake_install.cmake")
  INCLUDE("/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/magneto/mesh/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/twinkler/my_mycromagnumGPU_switch/MycroMagnum-master/src/build/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
