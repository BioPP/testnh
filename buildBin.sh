#! /bin/sh
arch=`uname -m`
version=2.3.2

strip TestNH/testnh
strip TestNH/mapnh
strip TestNH/partnh
strip TestNH/randnh
tar cvzf testnh-${arch}-bin-static-${version}.tar.gz TestNH/testnh TestNH/mapnh TestNH/partnh TestNH/randnh

