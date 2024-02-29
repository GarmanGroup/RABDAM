
# RABDAM
# Copyright (C) 2024 Garman Group, University of Oxford

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
    """
    Checks whether the system has the required packages / programs installed
    in order to be able to run RABDAM
    """

    import imp
    import pkg_resources  # Not in standard Python library, but currently I
    # can't find an alternative method to get the version of a Python package

    dependencies_met = True

    # Checks python packages called by RABDAM
    python_package_dict = {}
    with open('requirements.txt', 'r') as f:
        for line in f.readlines():
            line = line.replace('\n', '').strip()
            if line == '':
                pass
            elif any(x in line for x in ['>=', '>', '==', '<', '<=']):
                package = line.split('>')[0].split('=')[0].split('<')[0].strip()
                version = line.split('>')[-1].split('=')[-1].split('<')[-1].strip()
                python_package_dict[package] = version
            else:
                python_package_dict[line] = None

    for package, req_version in python_package_dict.items():
        try:
            imp.find_module(package)
        except ImportError:
            dependencies_met = False
            print('Python package %s not installed' % package)
            print('Install via "pip install %s"' % package)

        if not req_version is None:
            act_version = pkg_resources.get_distribution('%s' % package).version.strip()
            act_version_list = act_version.split('.')
            req_version_list = req_version.split('.')

            if len(act_version_list) != 3:
                dependencies_met = False
                print(
                    'Expected {} package version to be in format X.Y.Z, instead'
                    ' it is formatted as {} - unable to recognise if the '
                    'package version meets RABDAM\'s requirements or '
                    'not'.format(package, act_version)
                )

            else:
                for i, req_val in enumerate(req_version_list):
                    req_val = int(req_val)
                    act_val = int(act_version_list[i])

                    if act_val > req_val:
                        break
                    elif act_val < req_val:
                        dependencies_met = False
                        print(
                            '\nThe {} package is installed on the system, but a'
                            ' more recent version is required to run RABDAM.'
                            '\nUpdate via "pip install {} --upgrade"'.format(package, package)
                        )
                    else:
                        if i == 2:
                            break
                        else:
                            continue

    if dependencies_met is True:
        print('\nSystem meets all RABDAM dependencies')
