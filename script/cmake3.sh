#!/bin/sh
#  Copyright Synge Todo 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

VERSION=3.3.1

PREFIX="$1"
BUILD_DIR="$2"
SRC_DIR="$3"

if test -z "$BUILD_DIR"; then
  echo "$0 prefix build_dir [src_dir]"
  exit 127
fi

LOG="$0.log.$$"
echo "executing $0 $*" | tee "$LOG"

URL="http://www.cmake.org/files/v3.3/cmake-$VERSION.tar.gz"
SRC="$SRC_DIR/cmake-$VERSION.tar.gz"

echo "cleaning up..." | tee -a "$LOG"
if test -d "$BUILD_DIR"; then
  rm -rf "$BUILD_DIR/cmake-$VERSION"
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
(cd "$BUILD_DIR/cmake-$VERSION" && ./bootstrap --prefix="$PREFIX") && \
echo "building..." && \
(cd "$BUILD_DIR/cmake-$VERSION" && make) && \
echo "installing..." && \
(cd "$BUILD_DIR/cmake-$VERSION" && make install) && \
echo "cleaning up..." && \
rm -rf "$BUILD_DIR/cmake-$VERSION" && \
echo "done" \
) 2>&1 | tee -a "$LOG"

echo "log file = $LOG" | tee -a "$LOG"
