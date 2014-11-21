project(diameter)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated -Wno-sign-conversion -Wno-conversion-null")
add_graphlab_executable(diameter diameter.cpp)
