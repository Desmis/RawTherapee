include(FetchContent)
find_package(PkgConfig REQUIRED)

include(CheckCSourceCompiles)

# Use a macro as find_package() or pkg_check_modules() may not set variables
# at the global/parent scope.
macro(rt_setup_dependencies)
    rt_fetch_content()

    find_package(ATOMIC)
    find_package(JPEG REQUIRED)
    find_package(PNG REQUIRED)
    find_package(TIFF 4.0.4 REQUIRED)
    find_package(ZLIB REQUIRED)

    # Gtk version shall be greater than 3.24.3 for fixed Hi-DPI support
    pkg_check_modules(GTK REQUIRED IMPORTED_TARGET gtk+-3.0>=3.24.3)
    pkg_check_modules(GTKMM REQUIRED IMPORTED_TARGET gtkmm-3.0>=3.24)

    if(GTK_VERSION VERSION_GREATER "3.24.1" AND GTK_VERSION VERSION_LESS "3.24.7")
        if(GTK_VERSION VERSION_EQUAL "3.24.5")
            set(CERTAINTY "known to")
        else()
            set(CERTAINTY "likely to")
        endif()
        message(
            WARNING
                "\nWarning! You are using GTK+ version "
                ${GTK_VERSION}
                " which is "
                ${CERTAINTY}
                " have an issue where combobox menu scroll-arrows are missing when a Gtk::ComboBox list does not fit vertically on the screen. As a result, users of your build will not be able to select items in the following comboboxes: Processing Profiles, Film Simulation, and the camera and lens profiles in Profiled Lens Correction.\nIt is recommended that you either downgrade GTK+ to <= 3.24.1 or upgrade to >= 3.24.7."
        )
    endif()

    # These should be transitively included from PkgConfig::GTKMM.
    pkg_check_modules(CAIROMM REQUIRED IMPORTED_TARGET cairomm-1.0)
    pkg_check_modules(GIO REQUIRED IMPORTED_TARGET gio-2.0>=2.48)
    pkg_check_modules(GIOMM REQUIRED IMPORTED_TARGET giomm-2.4>=2.48)
    pkg_check_modules(GLIB2 REQUIRED IMPORTED_TARGET glib-2.0>=2.48)
    pkg_check_modules(GLIBMM REQUIRED IMPORTED_TARGET glibmm-2.4>=2.48)
    pkg_check_modules(GOBJECT REQUIRED IMPORTED_TARGET gobject-2.0>=2.48)
    pkg_check_modules(GTHREAD REQUIRED IMPORTED_TARGET gthread-2.0>=2.48)
    pkg_check_modules(SIGC REQUIRED IMPORTED_TARGET sigc++-2.0>=2.3.1)

    pkg_check_modules(EXIV2 REQUIRED IMPORTED_TARGET exiv2>=0.24)
    pkg_check_modules(EXPAT REQUIRED IMPORTED_TARGET expat>=2.1)
    pkg_check_modules(IPTCDATA REQUIRED IMPORTED_TARGET libiptcdata)
    pkg_check_modules(LENSFUN REQUIRED IMPORTED_TARGET lensfun>=0.2)
    pkg_check_modules(RSVG REQUIRED IMPORTED_TARGET librsvg-2.0>=2.52)

    pkg_check_modules(LCMS REQUIRED IMPORTED_TARGET lcms2>=2.6)
    # By default, little-cms2 uses the 'register' keyword which is deprecated
    # and removed in c++17. This definition forces it not to use the keyword.
    target_compile_definitions(PkgConfig::LCMS INTERFACE "CMS_NO_REGISTER_KEYWORD")

    # Check before OpenMP as rt_setup_openmp() handles fftw3f + omp integration
    pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
    if(OPTION_OMP)
        # Assumes FFTW3F is the package name for fftw3f
        rt_setup_openmp()
    endif()

    rt_setup_jxl()

    if(WITH_SYSTEM_KLT)
        find_package(KLT REQUIRED)
    endif()
    if(WITH_SYSTEM_LIBRAW)
        pkg_check_modules(LIBRAW REQUIRED IMPORTED_TARGET libraw_r>=0.21)
    endif()

    # Check for libcanberra-gtk3 (sound events on Linux):
    if(USE_LIBCANBERRA)
        pkg_check_modules(CANBERRA REQUIRED IMPORTED_TARGET libcanberra-gtk3)
        target_compile_definitions(PkgConfig::CANBERRA INTERFACE "USE_CANBERRA")
    endif()
endmacro()

function(rt_fetch_content)
    # fmt::fmt
    FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG 11.1.4
        GIT_SHALLOW ON
    )

    # Add all FetchContent-declared libraries here.
    # Don't use FetchContent_Declare after this.
    FetchContent_MakeAvailable(
        fmt
    )
endfunction()

macro(rt_setup_openmp)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        # Check for libfftw3f_omp
        if(NOT FFTW3F_FOUND)
            message(FATAL_ERROR "Missing PkgConfig::FFTW3F target")
        endif()

        # Prepare for check_c_source_compiles()
        set(CMAKE_REQUIRED_INCLUDES ${FFTW3F_INCLUDE_DIRS})
        set(CMAKE_REQUIRED_LIBRARIES)
        foreach(l ${FFTW3F_LIBRARIES})
            set(_f "_f-NOTFOUND")
            find_library(_f ${l} PATHS ${FFTW3F_LIBRARY_DIRS})
            list(APPEND CMAKE_REQUIRED_LIBRARIES ${_f})
        endforeach()

        check_c_source_compiles(
            "
            #include <fftw3.h>
            int main()
            {
                fftwf_init_threads();
                fftwf_plan_with_nthreads(1);
                return 0;
            }
            "
            FFTW3F_SUPPORTS_MULTITHREADING
        )
        if(FFTW3F_SUPPORTS_MULTITHREADING)
            target_compile_definitions(PkgConfig::FFTW3F INTERFACE "RT_FFTW3F_OMP")
            return()
        endif()

        find_library(fftw3f_omp fftw3f_omp PATHS ${FFTW3F_LIBRARY_DIRS})
        if(fftw3f_omp)
            target_link_libraries(PkgConfig::FFTW3F INTERFACE ${fftw3f_omp})
            target_compile_definitions(PkgConfig::FFTW3F INTERFACE "RT_FFTW3F_OMP")
        endif()
    endif()
endmacro()

macro(rt_setup_jxl)
    if(NOT (WITH_JXL OR WITH_JXL STREQUAL "AUTO"))
        set(JXL_VERSION "Disabled")
    else()
        if(NOT WITH_JXL STREQUAL "AUTO")
            set(JXL_REQUIRED "REQUIRED")
        else()
            set(JXL_REQUIRED "")
        endif()

        pkg_check_modules(JXL ${JXL_REQUIRED} IMPORTED_TARGET libjxl)
        pkg_check_modules(JXLTHREADS ${JXL_REQUIRED} IMPORTED_TARGET libjxl_threads)

        if(JXL_FOUND)
            target_compile_definitions(PkgConfig::JXL INTERFACE "LIBJXL")
        else()
            set(JXL_VERSION "Disabled, Auto")
        endif()
    endif()
endmacro()
