set(SOURCES ${sources} mgr_plugin.cpp)
set(HEADERS ${headers} mgr_plugin.h)

add_meshlab_plugin(mgr_projekt ${SOURCES} ${HEADERS})

