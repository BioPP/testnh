%define _basename testnh
%define _version 2.3.2
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
%define compress_program xz
%else
%if 0%{?mdkversion}
BuildRequires: lzma
%define compress_program lzma
%else
BuildRequires: gzip
%define compress_program gzip
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
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DCOMPRESS_PROGRAM=%{compress_program}"
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/bin/*
%{_prefix}/share/info/*.info*
%{_prefix}/share/man/man1/*.1*

%changelog
* Mon Oct 06 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.1.0-1
- Bug fixes and output format update.
* Tue Feb 09 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.0.0-1
- Initial release.

