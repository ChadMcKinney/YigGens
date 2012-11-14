#-------------------------------------------------
#
# YigGens
#
#-------------------------------------------------

QT       -= core

TEMPLATE = lib
CONFIG -= qt plugin

mac {
    TARGET = YigGens.scx

    ##########################################################
    # MAKE SURE THIS POINTS TO YOUR MacOSX sdk!!!!!
    ##########################################################
    QMAKE_MAC_SDK = /Developer/SDKs/MacOSX10.6.sdk

    # This make it actually named what we want rather than libxxx.scx.dylib
    QMAKE_POST_LINK = mv -f $(TARGET) $${TARGET}
    QMAKE_CXXFLAGS += -DSC_DARWIN -DNOVA_SIMD
}

linux-g++ {
    TARGET = YigGens
    DEFINES += __LINUX__
    DEFINES += SC_AUDIO_API=SC_AUDIO_API_JACK
}

##########################################################
# CHANGE THIS PATH TO YOUR SUPERCOLLIDER SOURCE PATH!!!!!
##########################################################
SUPERCOLLIDER_SOURCE = /home/octopian/Documents/source/SuperCollider/SuperCollider-3.6-beta/supercollider

# Includes for sc sources
INCLUDEPATH += $${SUPERCOLLIDER_SOURCE}
INCLUDEPATH += $${SUPERCOLLIDER_SOURCE}/include
INCLUDEPATH += $${SUPERCOLLIDER_SOURCE}/external_libraries
INCLUDEPATH += $${SUPERCOLLIDER_SOURCE}/external_libraries/libsndfile
INCLUDEPATH += $${SUPERCOLLIDER_SOURCE}/include/common
INCLUDEPATH += $${SUPERCOLLIDER_SOURCE}/include/server
INCLUDEPATH += $${SUPERCOLLIDER_SOURCE}/include/plugin_interface
INCLUDEPATH += $${SUPERCOLLIDER_SOURCE}/include/lang


SOURCES += \
    source/YigGens/YigGens.cpp




