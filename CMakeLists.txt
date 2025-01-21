cmake_minimum_required(VERSION 3.13.4)
project(satellite-tracker)

# Настройки C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_AUTOMOC ON)  # Включаем автоматическую поддержку MOC

# Настройки компилятора
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -pedantic")
endif()

# Находим Qt
find_package(Qt6 REQUIRED COMPONENTS Core)

# Скачиваем SGP4 если его нет
include(FetchContent)
FetchContent_Declare(
    sgp4
    GIT_REPOSITORY https://github.com/AcNasDev/sgp4_.git
    GIT_TAG        master
)

# Делаем SGP4 доступным
FetchContent_GetProperties(sgp4)
if(NOT sgp4_POPULATED)
    FetchContent_Populate(sgp4)
    add_subdirectory(${sgp4_SOURCE_DIR}/libsgp4 ${sgp4_BINARY_DIR}/libsgp4)
endif()

# Исходные файлы
set(SOURCES
    src/main.cpp
)

set(HEADERS
    src/SatelliteTracker.h
)

# Добавляем исполняемый файл
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_link_libraries(${PROJECT_NAME} PRIVATE
    Qt6::Core
    sgp4
)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${sgp4_SOURCE_DIR}/libsgp4
)
