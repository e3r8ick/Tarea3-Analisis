# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build"

# Include any dependencies generated for this target.
include test/CMakeFiles/tester.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/tester.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/tester.dir/flags.make

test/CMakeFiles/tester.dir/testAllocator.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testAllocator.cpp.o: ../test/testAllocator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/tester.dir/testAllocator.cpp.o"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testAllocator.cpp.o -c "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testAllocator.cpp"

test/CMakeFiles/tester.dir/testAllocator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testAllocator.cpp.i"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testAllocator.cpp" > CMakeFiles/tester.dir/testAllocator.cpp.i

test/CMakeFiles/tester.dir/testAllocator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testAllocator.cpp.s"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testAllocator.cpp" -o CMakeFiles/tester.dir/testAllocator.cpp.s

test/CMakeFiles/tester.dir/testAllocator.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testAllocator.cpp.o.requires

test/CMakeFiles/tester.dir/testAllocator.cpp.o.provides: test/CMakeFiles/tester.dir/testAllocator.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testAllocator.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testAllocator.cpp.o.provides

test/CMakeFiles/tester.dir/testAllocator.cpp.o.provides.build: test/CMakeFiles/tester.dir/testAllocator.cpp.o


test/CMakeFiles/tester.dir/testMain.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testMain.cpp.o: ../test/testMain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object test/CMakeFiles/tester.dir/testMain.cpp.o"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testMain.cpp.o -c "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testMain.cpp"

test/CMakeFiles/tester.dir/testMain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testMain.cpp.i"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testMain.cpp" > CMakeFiles/tester.dir/testMain.cpp.i

test/CMakeFiles/tester.dir/testMain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testMain.cpp.s"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testMain.cpp" -o CMakeFiles/tester.dir/testMain.cpp.s

test/CMakeFiles/tester.dir/testMain.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testMain.cpp.o.requires

test/CMakeFiles/tester.dir/testMain.cpp.o.provides: test/CMakeFiles/tester.dir/testMain.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testMain.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testMain.cpp.o.provides

test/CMakeFiles/tester.dir/testMain.cpp.o.provides.build: test/CMakeFiles/tester.dir/testMain.cpp.o


test/CMakeFiles/tester.dir/testMatrix.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testMatrix.cpp.o: ../test/testMatrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object test/CMakeFiles/tester.dir/testMatrix.cpp.o"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testMatrix.cpp.o -c "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testMatrix.cpp"

test/CMakeFiles/tester.dir/testMatrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testMatrix.cpp.i"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testMatrix.cpp" > CMakeFiles/tester.dir/testMatrix.cpp.i

test/CMakeFiles/tester.dir/testMatrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testMatrix.cpp.s"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testMatrix.cpp" -o CMakeFiles/tester.dir/testMatrix.cpp.s

test/CMakeFiles/tester.dir/testMatrix.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testMatrix.cpp.o.requires

test/CMakeFiles/tester.dir/testMatrix.cpp.o.provides: test/CMakeFiles/tester.dir/testMatrix.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testMatrix.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testMatrix.cpp.o.provides

test/CMakeFiles/tester.dir/testMatrix.cpp.o.provides.build: test/CMakeFiles/tester.dir/testMatrix.cpp.o


test/CMakeFiles/tester.dir/testRootFinders.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testRootFinders.cpp.o: ../test/testRootFinders.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object test/CMakeFiles/tester.dir/testRootFinders.cpp.o"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testRootFinders.cpp.o -c "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testRootFinders.cpp"

test/CMakeFiles/tester.dir/testRootFinders.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testRootFinders.cpp.i"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testRootFinders.cpp" > CMakeFiles/tester.dir/testRootFinders.cpp.i

test/CMakeFiles/tester.dir/testRootFinders.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testRootFinders.cpp.s"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test/testRootFinders.cpp" -o CMakeFiles/tester.dir/testRootFinders.cpp.s

test/CMakeFiles/tester.dir/testRootFinders.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testRootFinders.cpp.o.requires

test/CMakeFiles/tester.dir/testRootFinders.cpp.o.provides: test/CMakeFiles/tester.dir/testRootFinders.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testRootFinders.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testRootFinders.cpp.o.provides

test/CMakeFiles/tester.dir/testRootFinders.cpp.o.provides.build: test/CMakeFiles/tester.dir/testRootFinders.cpp.o


# Object files for target tester
tester_OBJECTS = \
"CMakeFiles/tester.dir/testAllocator.cpp.o" \
"CMakeFiles/tester.dir/testMain.cpp.o" \
"CMakeFiles/tester.dir/testMatrix.cpp.o" \
"CMakeFiles/tester.dir/testRootFinders.cpp.o"

# External object files for target tester
tester_EXTERNAL_OBJECTS =

test/tester: test/CMakeFiles/tester.dir/testAllocator.cpp.o
test/tester: test/CMakeFiles/tester.dir/testMain.cpp.o
test/tester: test/CMakeFiles/tester.dir/testMatrix.cpp.o
test/tester: test/CMakeFiles/tester.dir/testRootFinders.cpp.o
test/tester: test/CMakeFiles/tester.dir/build.make
test/tester: src/libanpi.a
test/tester: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
test/tester: /usr/lib/x86_64-linux-gnu/libboost_system.so
test/tester: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so
test/tester: test/CMakeFiles/tester.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable tester"
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tester.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/tester.dir/build: test/tester

.PHONY : test/CMakeFiles/tester.dir/build

test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testAllocator.cpp.o.requires
test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testMain.cpp.o.requires
test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testMatrix.cpp.o.requires
test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testRootFinders.cpp.o.requires

.PHONY : test/CMakeFiles/tester.dir/requires

test/CMakeFiles/tester.dir/clean:
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" && $(CMAKE_COMMAND) -P CMakeFiles/tester.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/tester.dir/clean

test/CMakeFiles/tester.dir/depend:
	cd "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code" "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/test" "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build" "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test" "/home/vmontero/Desktop/U/Anpi/Tarea 3/Tarea3-Analisis/code/build/test/CMakeFiles/tester.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : test/CMakeFiles/tester.dir/depend

