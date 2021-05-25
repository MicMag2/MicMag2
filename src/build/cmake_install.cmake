# Install script for directory: /localscratch/twinkler/MicMag2/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages" TYPE DIRECTORY FILES "/localscratch/twinkler/MicMag2/src/magnum" REGEX "/magnum\\/magneto\\_cuda\\.py\\ magnum\\/\\_magneto\\_cuda\\.so\\ magnum\\/magneto\\_cpu\\.py\\ magnum\\/\\_magneto\\_cpu\\.so\\ CMakeLists\\.txt$" EXCLUDE)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages/magnum" TYPE FILE FILES
    "/localscratch/twinkler/MicMag2/src/build/magneto_cuda.py"
    "/localscratch/twinkler/MicMag2/src/build/_magneto_cuda.so"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/localscratch/twinkler/MicMag2/src/build/magneto/cmake_install.cmake")
  include("/localscratch/twinkler/MicMag2/src/build/magneto/bindings/cmake_install.cmake")
  include("/localscratch/twinkler/MicMag2/src/build/magneto/matrix/cmake_install.cmake")
  include("/localscratch/twinkler/MicMag2/src/build/magneto/math/cmake_install.cmake")
  include("/localscratch/twinkler/MicMag2/src/build/magneto/math/conv/cmake_install.cmake")
  include("/localscratch/twinkler/MicMag2/src/build/magneto/math/conv/kernels/cmake_install.cmake")
  include("/localscratch/twinkler/MicMag2/src/build/magneto/mmm/cmake_install.cmake")
  include("/localscratch/twinkler/MicMag2/src/build/magneto/evolver/cmake_install.cmake")
  include("/localscratch/twinkler/MicMag2/src/build/magneto/mesh/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/localscratch/twinkler/MicMag2/src/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
