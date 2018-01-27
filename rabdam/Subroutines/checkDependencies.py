
# RABDAM
# Copyright (C) 2018 Garman Group, University of Oxford

# This file is part of RABDAM.

# RABDAM is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# RABDAM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General
# Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


def check_RABDAM_dependencies():
    # Checks whether the system has the required packages / programs installed
    # in order to be able to run RABDAM

    import imp
    import os
    import platform
    import pkg_resources  # Not in standard Python library, but currently I
    # can't find an alternative method to get the version of a Python package

    dependencies_met = True

    # Checks python packages called by RABDAM
    python_package_list = ['numpy', 'matplotlib', 'pandas', 'requests',
                           'seaborn']
    for package in python_package_list:
        try:
            imp.find_module(package)
            if package in ['pandas', 'seaborn']:
                version = pkg_resources.get_distribution('%s' % package).version
                version_list = list(version)
                count = version_list.count('.')
                if count < 2:
                    version = version + '.0'
                version = int(version.replace('.', ''))
                if package == 'pandas':
                    if version < 180:
                        dependencies_met = False
                        print('\nThe pandas Python data analysis library is '
                              'installed on the system,\n'
                              'but a more recent version is required to run '
                              'RABDAM.'
                              '\nUpdate via "pip install pandas --upgrade"')
                elif package == 'seaborn':
                    if version < 70:
                        dependencies_met = False
                        print('\nThe seaborn Python plotting library is '
                              'installed on the system,\n'
                              'but a more recent version is required to run '
                              'RABDAM.'
                              '\nUpdate via "pip install seaborn --upgrade"')

        except ImportError:
            dependencies_met = False
            print('Python package %s not installed' % package)
            print('Install via "pip install %s"' % package)

    # Checks that CCP4 has been downloaded and is accessible in the command
    # prompt / terminal being used
    operating_system = platform.system().strip()
    if operating_system.lower() == 'windows':
        os.system('pdbcur > test.txt')
    else:
        os.system('pdbcur > test.txt 2>&1')

    test_file = open('test.txt', 'r')
    test_file_lines = test_file.readlines()
    if test_file_lines == []:
        dependencies_met = False
        print('\nCCP4 program suite not found. If you have downloaded the\n'
              'CCP4 suite, ensure that you are running RABDAM in a\n'
              'terminal / command prompt that can run CCP4 programs')
    else:
        test_file_line = test_file_lines[0]
        if 'pdbcur:commandnotfound' in test_file_line.replace(' ', ''):
            dependencies_met = False
            print('\nCCP4 program suite not found. If you have downloaded the\n'
                  'CCP4 suite, ensure that you are running RABDAM in a\n'
                  'terminal / command prompt that can run CCP4 programs')
    test_file.close()
    os.remove('test.txt')

    if dependencies_met is True:
        print('\nSystem meets all RABDAM dependencies')
