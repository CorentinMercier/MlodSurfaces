cmake_minimum_required(VERSION 3.9)
project(glad C)

add_library(glad STATIC src/glad.c include/glad/glad.h include/KHR/khrplatform.h)
target_include_directories(glad PUBLIC include)
if(NOT WIN32)
	target_link_libraries(glad PUBLIC ${CMAKE_DL_LIBS})
endif()
