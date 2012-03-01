%define name testnh
%define version 1.0.0
%define release 1
%define _prefix /usr

Summary: The TestNH package.
Name: %{name}
Version: %{version}
Release: %{release}
Vendor: Julien Dutheil
Source: http://download.gna.org/comap/%{name}-%{version}.tar.gz
License: CeCILL 2
Group: Productivity/Scientific/Other
BuildRoot: %{_builddir}/%{name}-root
Packager: Julien Dutheil
Prefix: %{_prefix}
AutoReq: yes
AutoProv: yes

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
cmake $CMAKE_FLAGS .
make
make info

%install
rm -rf $RPM_BUILD_ROOT
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
%{_prefix}/share/info/testnh.info.gz
%{_prefix}/share/man/man1/testnh.1.gz
%{_prefix}/share/man/man1/mapnh.1.gz
%{_prefix}/share/man/man1/partnh.1.gz
%{_prefix}/share/man/man1/randnh.1.gz

%changelog
* Tue Feb 09 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- TestNH 1.0.0 release

