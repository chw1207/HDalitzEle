# https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/-/blob/master/cmake/Findxgboost.cmake

find_path(xgboost_INCLUDE_DIR xgboost/c_api.h)
find_library(xgboost_LIBRARY NAMES xgboost)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(xgboost
    REQUIRED_VARS xgboost_INCLUDE_DIR xgboost_LIBRARY)
mark_as_advanced(xgboost_FOUND xgboost_INCLUDE_DIR xgboost_LIBRARY)

if (xgboost_FOUND)
    set(xgboost_INCLUDE_DIRS ${xgboost_INCLUDE_DIR})
    set(xgboost_LIBRARIES ${xgboost_LIBRARY})
endif()

if (xgboost_FOUND AND NOT TARGET xgboost::xgboost)
    add_library(xgboost::xgboost SHARED IMPORTED ${xgboost_LIBRARY})
    set_target_properties(xgboost::xgboost PROPERTIES
        IMPORTED_LOCATION ${xgboost_LIBRARY}
        INTERFACE_INCLUDE_DIRECTORIES ${xgboost_INCLUDE_DIR}
    )
endif()

