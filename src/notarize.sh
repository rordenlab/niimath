#!/bin/bash

#Before running, install a profile named "niimath"
# https://blog.smittytone.net/2022/06/09/wwdc-22-notarise-macos-command-line-apps-more-quickly/
# xcrun notarytool store-credentials niimath --apple-id robert@sc.edu --team-id 69BBQ234R --password abcd-efgh-bork-zork

COMPANY_NAME=mricro
APP_NAME=niimath
APP_DIR=macos
# security -v find-identity -p codesigning
# security find-identity
if [[ -z "$APPLE_ID_APP" ]]; then
    echo "APPLE_ID_APP environment variable required"
    echo "to find your ID, run:"
    echo " security find-identity"
    echo "then export this as a variable:"
    echo " export APPLE_ID_APP=\"Developer ID Application: John Doe (69BBQ123RR)\""
    exit 1
fi


cd "$(dirname "$0")"

mkdir ${APP_DIR}

DFLAGS="-DHAVE_ZLIB -DFSLSTYLE -DPIGZ -DREJECT_COMPLEX -DNII2MESH -DHAVE_64BITS -DHAVE_BUTTERWORTH -DHAVE_FORMATS -DHAVE_TENSOR -DHAVE_DTIFIT -DHAVE_CONFORM -DHAVE_BMP -DHAVE_ALLINEATE"
SRCS="niimath.c MarchingCubes.c meshify.c quadric.c base64.c radixsort.c fdr.c bwlabel.c bw.c core.c tensor.c dtifit.c core32.c core64.c conform.c unifize.c filter.c bmp.c spng.c nifti_io.c"
AL_SRCS="allineate.c powell_newuoa.c"

build_arch() {
    local target="$1" minver="$2" output="$3"
    # Compile allineate separately with -ffast-math (scoped)
    gcc -O3 -ffast-math -fno-finite-math-only ${DFLAGS} -target "$target" -mmacosx-version-min="$minver" -c ${AL_SRCS}
    # Compile everything else without -ffast-math, link allineate objects
    gcc -sectcreate TEXT info_plist Info.plist -O3 ${DFLAGS} ${SRCS} allineate.o powell_newuoa.o -lm -lz -target "$target" -mmacosx-version-min="$minver" -o "$output"
    rm -f allineate.o powell_newuoa.o
    strip "./$output"
}

# Compile x86
build_arch "x86_64-apple-macos10.12" "10.12" "niimathX86"

# Compile ARM
build_arch "arm64-apple-macos11" "11.0" "niimathARM"

# Create the universal binary
lipo -create -output ./${APP_DIR}/${APP_NAME} niimathX86 niimathARM
rm ./niimathX86; rm ./niimathARM

#code sign executable
codesign --timestamp --options=runtime -s "${APPLE_ID_APP}" -v ./${APP_DIR}/${APP_NAME}

#create a DMG
hdiutil create -volname ${APP_NAME} -srcfolder ./${APP_DIR} -ov -format UDZO -layout SPUD -fs HFS+J  ${APP_NAME}_macOS.dmg
xcrun notarytool submit ${APP_NAME}_macOS.dmg  --keychain-profile ${APP_NAME}  --wait

#create a PKG
if [[ -z "$APPLE_ID_INSTALL" ]]; then
    echo "APPLE_ID_INSTALL environment variable required"
    echo "to find your ID, run:"
    echo " security find-identity"
    echo "then export this as a variable:"
    echo " export APPLE_ID_INSTALL=\"Developer ID Installer: John Doe (69BBQ123RR)\""
    exit 1
fi
pkgbuild --identifier "com.${COMPANY_NAME}.${APP_NAME}.pkg" --sign "${APPLE_ID_INSTALL}" --timestamp --root $APP_DIR --install-location /usr/local/bin/ ${APP_NAME}_macOS.pkg
xcrun notarytool submit ${APP_NAME}_macOS.pkg  --keychain-profile ${APP_NAME}  --wait
