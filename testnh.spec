%define _basename testnh
%define _version 1.1.0
%define _release 1
%define _prefix /usr

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: Julien Dutheil
Source: http://biopp.univ-montp2.fr/repos/sources/testnh/%{_basename}-%{_version}.tar.gz
Summary: The TestNH package
Group: Productivity/Scientific/Other

Requires: libbpp-phyl9 = 2.2.0
Requires: libbpp-seq9 = 2.2.0
Requires: libbpp-core2 = 2.2.0

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.6.0
BuildRequires: gcc-c++ >= 4.0.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core2 = 2.2.0
BuildRequires: libbpp-core-devel = 2.2.0
BuildRequires: libbpp-seq9 = 2.2.0
BuildRequires: libbpp-seq-devel = 2.2.0
BuildRequires: libbpp-phyl9 = 2.2.0
BuildRequires: libbpp-phyl-devel = 2.2.0


AutoReq: yes
AutoProv: yes
%if 0%{?mdkversion} >= 201100 || %{?distribution} == "Mageia"
BuildRequires: xz
%define zipext xz
%else
%if 0%{?mdkversion}
BuildRequires: lzma
%define zipext lzma
%else
BuildRequires: gzip
%define zipext gz
%endif
%endif

%description
Includes programs:
 - TestNH for performing Bowker's test for non-stationary evolution
 - MapNH for mapping substitutions and cluster branches according to substitution process
 - PartNH for fitting and selecting non-homogeneous models of sequence evolution
 - RandNH for generating random non-homogeneous models of sequence evolution

%prep
%setup -q

%build
CFLAGS="-I%{_prefix}/include $RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix}"
if [ %{_lib} == 'lib64' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DLIB_SUFFIX=64"
fi
if [ %{zipext} == 'lzma' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=lzma -DDOC_COMPRESS_EXT=lzma"
fi
if [ %{zipext} == 'xz' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=xz -DDOC_COMPRESS_EXT=xz"
fi

cmake $CMAKE_FLAGS .
make
make info

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/bin/testnh
%{_prefix}/bin/mapnh
%{_prefix}/bin/partnh
%{_prefix}/bin/randnh
%{_prefix}/share/info/testnh.info.%{zipext}
%{_prefix}/share/man/man1/testnh.1.%{zipext}
%{_prefix}/share/man/man1/mapnh.1.%{zipext}
%{_prefix}/share/man/man1/partnh.1.%{zipext}
%{_prefix}/share/man/man1/randnh.1.%{zipext}

%changelog
* Mon Oct 06 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.1.0-1
- Bug fixes and output format update.
* Tue Feb 09 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.0.0-1
- Initial release.

