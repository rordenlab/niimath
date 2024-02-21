#!/bin/bash

#Before running, install a profile named "niimath"
# https://blog.smittytone.net/2022/06/09/wwdc-22-notarise-macos-command-line-apps-more-quickly/
# xcrun notarytool store-credentials niimath --apple-id robert@sc.edu --team-id 69BBQ234R --password abcd-efgh-bork-zork

COMPANY_NAME=mricro
APP_NAME=niimath
APP_DIR=macos
# security -v find-identity -p codesigning
# security find-identity

APPLE_ID_INSTALL="3rd Party Mac Developer Installer: Christopher Rorden (68BQDQS28R)"

cd "$(dirname "$0")"

mkdir ${APP_DIR}
# Compile x86
gcc -sectcreate TEXT info_plist Info.plist -O3 -DHAVE_ZLIB -DFSLSTYLE -DPIGZ -DREJECT_COMPLEX -DNII2MESH niimath.c MarchingCubes.c meshify.c quadric.c base64.c radixsort.c fdr.c bwlabel.c bw.c core.c tensor.c core32.c core64.c niftilib/nifti2_io.c znzlib/znzlib.c -I./niftilib -I./znzlib -lm -lz  -target x86_64-apple-macos10.12 -mmacosx-version-min=10.12 -o niimathX86
strip ./niimathX86

# Compile ARM
gcc -sectcreate TEXT info_plist Info.plist -O3 -DHAVE_ZLIB -DFSLSTYLE -DPIGZ -DREJECT_COMPLEX -DNII2MESH niimath.c MarchingCubes.c meshify.c quadric.c base64.c radixsort.c fdr.c bwlabel.c bw.c core.c tensor.c core32.c core64.c niftilib/nifti2_io.c znzlib/znzlib.c -I./niftilib -I./znzlib -lm -lz  -target arm64-apple-macos11 -mmacosx-version-min=11.0 -o niimathARM
strip ./niimathARM

# Create the universal binary
lipo -create -output ./${APP_DIR}/${APP_NAME} niimathX86 niimathARM
rm ./niimathX86; rm ./niimathARM

#pkgbuild --identifier "com.${COMPANY_NAME}.${APP_NAME}.pkg" --sign "${APPLE_ID_INSTALL}" --timestamp --root $APP_DIR --install-location /usr/local/bin/ ${APP_NAME}.pkg

hdiutil create -volname ${APP_NAME} -srcfolder ./${APP_DIR} -ov -format UDZO -layout SPUD -fs HFS+J  ${APP_NAME}_macOS.dmg

xcrun notarytool submit ${APP_NAME}_macOS.dmg  --keychain-profile ${APP_NAME}  --wait
