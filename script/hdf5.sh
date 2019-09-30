#!/bin/sh
#  Copyright Synge Todo 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

VERSION=1.8.15-patch1

PREFIX="$1"
BUILD_DIR="$2"
SRC_DIR="$3"

if test -z "$BUILD_DIR"; then
  echo "$0 prefix build_dir [src_dir]"
  exit 127
fi

LOG="$0.log.$$"
echo "executing $0 $*" | tee "$LOG"

URL="http://ftp.hdfgroup.org/ftp/HDF5/releases/hdf5-$VERSION/src/hdf5-$VERSION.tar.gz"
SRC="$SRC_DIR/hdf5-$VERSION.tar.gz"

echo "cleaning up..." | tee -a "$LOG"
if test -d "$BUILD_DIR"; then
  rm -rf "$BUILD_DIR/hdf5-$VERSION"
else
  mkdir -p "$BUILD_DIR"
fi

echo "retrieving source files..." | tee -a "$LOG"
if test -n "$SRC_DIR" && test -f "$SRC"; then
  (cd "$BUILD_DIR" && tar zxf $SRC) 2>&1 | tee -a "$LOG"
else
  CURL=`which curl`
  if test -z "$CURL"; then
    echo "curl utility not found" | tee -a "$LOG"
    exit 127
  fi
  (cd "$BUILD_DIR" && "$CURL" "$URL" | tar zxf -) 2>&1 | tee -a "$LOG"
fi

( \
echo "configuring..." && \
(cd "$BUILD_DIR/hdf5-$VERSION" && ./configure --prefix="$PREFIX") && \
echo "building..." && \
(cd "$BUILD_DIR/hdf5-$VERSION" && make) && \
echo "installing..." && \
(cd "$BUILD_DIR/hdf5-$VERSION" && make install) && \
echo "cleaning up..." && \
(cd "$BUILD_DIR/hdf5-$VERSION" && make distclean) && \
echo "configuring with Fortran..." && \
(cd "$BUILD_DIR/hdf5-$VERSION" && ./configure --prefix="$PREFIX" --enable-fortran) && \
echo "building with Fortran..." && \
(cd "$BUILD_DIR/hdf5-$VERSION" && make) && \
echo "installing with Fortran..." && \
(cd "$BUILD_DIR/hdf5-$VERSION" && make install) && \
echo "cleaning up..." && \
rm -rf "$BUILD_DIR/hdf5-$VERSION" && \

echo "done" \
) 2>&1 | tee -a "$LOG"

echo "log file = $LOG" | tee -a "$LOG"
