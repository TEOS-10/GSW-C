

# SPDX-FileCopyrightText: 2024 Peter Urban
#
# SPDX-License-Identifier: BSD-3-Clause

# 'TEOS-10 V3.06.16-0 GSW Oceanographic Toolbox in C'
# 'http://teos-10.org'

# --- Project ---
# Define project meta data
project(
  'gswteos-10',
  ['c', 'cpp'], #language is c on linux/mac and cpp on windows
  license: 'BSD-3-Clause',

  version: 'v3.06.16-0',
  default_options: ['warning_level=2', 'buildtype=release'],
  meson_version: '>=1.5.1' #latest tested version,,,,
)

project_name = 'gswteos-10'
description = 'TEOS-10 V3.06.16-0 GSW Oceanographic Toolbox in C'
url = 'http://teos-10.org'

# --- dependencies ---
# libm (math library)
# Required is set to false because libm is only necessary and available if
# it is not already included in the compiler (e.g. gcc/linux)
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required: false)

# --- compile options ---
c_compile_args = []
cpp_compile_args = []

# Force compile language to be c++ on windows (cl and clang-cl) (or if force-cpp option is enabled)
link_language = 'c'
if cc.get_id() == 'cl' or cc.get_id() == 'msvc' or cc.get_id() == 'clang-cl'
  message('Setting /Tp flag to interpret c files as c++ files on windows')
  cpp_compile_args += ['/TP'] #interpret c files as c++ files on windows
  c_compile_args += ['/TP'] #interpret c files as c++ files on windows
  link_language = 'cpp'
elif get_option('force-cpp').enabled()
  if cc.get_id() == 'gcc' or cc.get_id() == 'clang'
    c_compile_args += ['-xc++']
  else
    error('Unhandled compiler id: ' + cc.get_id() + '! Cannot enable option force-cpp')
  endif
endif

# --- sources ---
program_sources = ['gsw_check_functions.c']

library_sources = [
  'gsw_oceanographic_toolbox.c',
  'gsw_saar.c',
]

export_headers = ['gswteos-10.h']

# --- define library and dependencies ---
gswteos_10_lib = library(
  project_name,
  library_sources,
  c_args: c_compile_args,
  cpp_args: cpp_compile_args,
  dependencies: [m_dep],
  link_language: link_language,
  install: true,
)

gswteos_10_dep = declare_dependency(
  dependencies: [m_dep],
  link_with: [gswteos_10_lib],
  include_directories: include_directories('.'),
)

# --- Installation ---
# install headers
install_headers(export_headers)

# pkgconfig file
pkg = import('pkgconfig')
pkg.generate(
  description: description,
  url: url,
  version: meson.project_version(),
  name: project_name,
  libraries: gswteos_10_lib,
)

# --- define tests ---
gswteos_10_exe = executable(
  'gsw_check',
  sources: [program_sources],
  dependencies: [gswteos_10_dep],
  c_args: c_compile_args,
  cpp_args: cpp_compile_args,
  link_language: link_language,
)

test('run_gsw_check', gswteos_10_exe)
