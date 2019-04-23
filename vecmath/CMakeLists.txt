project(vecmath)

set(SRCS
    Matrix2f.cpp
    Matrix3f.cpp
    Matrix4f.cpp
    Quat4f.cpp
    Vector2f.cpp
    Vector3f.cpp
    Vector4f.cpp
)

set(HEADERS
    include/Matrix2f.h
    include/Matrix3f.h
    include/Matrix4f.h
    include/Quat4f.h
    include/Vector2f.h
    include/Vector3f.h
    include/Vector4f.h
    include/vecmath.h
)

add_library(vecmath ${SRCS} ${HEADERS})
target_include_directories(vecmath PUBLIC include)
