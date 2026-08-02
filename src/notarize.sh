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

# THIS SCRIPT BUILDS THE SHIPPED, BSD-2-BRANDED macOS BINARY. Never add anything from
# src/GPL/ to SRCS, and never add -DHAVE_GPL (or any other copyleft define) to DFLAGS: doing
# so links GPL-licensed code into an artifact distributed as BSD-2, and the binary's own
# version string would still end in " BSD". If you want the copyleft payload, build with
# GPL=1 (make) or -DENABLE_GPL=ON (CMake) instead -- that is a separate, GPL-2-or-later
# artifact and is not what gets notarized. This rule used to be stated only as a note about
# -bandpass/bw.c; that op was retired and the note went with it, so it is restated here in
# general form. release_smoke.py's copyleft check is the backstop, not the rule.
DFLAGS="-DHAVE_ZLIB -DFSLSTYLE -DPIGZ -DREJECT_COMPLEX -DNII2MESH -DHAVE_64BITS -DHAVE_FORMATS -DHAVE_TENSOR -DHAVE_DTIFIT -DHAVE_QC -DHAVE_CONFORM -DHAVE_BMP -DHAVE_ALLINEATE -DHAVE_ROMEO -DHAVE_MEDIC"
SRCS="niimath.c MarchingCubes.c meshify.c quadric.c base64.c radixsort.c fdr.c bwlabel.c core.c tensor.c dtifit.c qc.c core32.c core64.c conform.c unifize.c filter.c bmp.c spng.c nifti_io.c medic.c"
AL_SRCS="allineate.c powell_newuoa.c coreg_fast.c reface.c"
build_arch() {
    local target="$1" minver="$2" output="$3"
    # -moco and -stc build for BOTH slices: their Apple-Silicon-only gate was lifted once CI began
    # checking both numerically (release_smoke.py) on every shipped target.
    # -skullstrip is DELIBERATELY absent here: it is OFF by default pending the native-release
    # gate in skullstrip_plan.md, so the notarised release must not ship it.
    # Add "-DHAVE_SKULLSTRIP" and "skullstrip.c"
    # to the two strings below when it is promoted to default-on.
    local arch_dflags="${DFLAGS} -DHAVE_MOCO -DHAVE_STC" arch_srcs="${SRCS} moco.c stc.c"
    # Whole-program -ffast-math, matching the Makefile/CMake/WASM release contract so every
    # shipped artifact shares one FP behavior; -fno-finite-math-only preserves NaN/Inf.
    # (allineate no longer needs a separate scoped compile — everything is fast-math now.)
    # ROMEO is compiled SEPARATELY, strict-FP, for this architecture and linked as an object:
    # it must not join the whole-program fast-math source line (see src/romeo.c for the measured
    # consequence -- weight bins move and regions shift by 2*pi).
    gcc -O3 -fno-fast-math -ffp-contract=off ${arch_dflags} -c romeo.c -target "$target" -mmacosx-version-min="$minver" -o "romeo_$target.o"
    gcc -sectcreate TEXT info_plist Info.plist -O3 -ffast-math -fno-finite-math-only ${arch_dflags} ${arch_srcs} ${AL_SRCS} "romeo_$target.o" -lm -lz -target "$target" -mmacosx-version-min="$minver" -o "$output"
    rm -f "romeo_$target.o"
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
