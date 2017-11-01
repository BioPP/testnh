#! /bin/sh
arch=`uname -m`
version=1.2.0-1

strip TestNH/testnh
strip TestNH/mapnh
strip TestNH/partnh
strip TestNH/randnh
tar cvzf testnh-${arch}-bin-static-${version}.tar.gz TestNH/testnh TestNH/mapnh TestNH/partnh TestNH/randnh

