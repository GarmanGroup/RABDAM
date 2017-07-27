
def check_RABDAM_dependencies():
    # Checks whether the system has the required packages / programs installed
    # in order to be able to run RABDAM

    import imp
    import os
    import pkg_resources

    dependencies_met = True

    # Checks python packages called by RABDAM
    python_package_list = ['numpy', 'matplotlib', 'pandas', 'requests',
                           'seaborn']
    for package in python_package_list:
        try:
            imp.find_module(package)
            if package in ['pandas', 'seaborn']:
                version = pkg_resources.get_distribution('%s' % package).version
                version.lstrip('0')
                version = int(version.replace('.', ''))
                if package == 'pandas':
                    if version < 180:
                        dependencies_met = False
                        print ('\nThe pandas Python data analysis library is '
                               'installed on the system,\n'
                               'but a more recent version is required to run '
                               'RABDAM.'
                               '\nUpdate via "pip install pandas --upgrade"')
                elif package == 'seaborn':
                    if version < 7:
                        dependencies_met = False
                        print ('\nThe seaborn Python plotting library is '
                               'installed on the system,\n'
                               'but a more recent version is required to run '
                               'RABDAM.'
                               '\nUpdate via "pip install seaborn --upgrade"')

        except ImportError:
            dependencies_met = False
            print ('Python package %s not installed' % package)
            print('Install via "pip install %s"' % package)

    # Checks that CCP4 has been downloaded and is accessible in the command
    # prompt / terminal being used
    os.system('pdbcur > test.txt')
    test_file = open('test.txt', 'r')
    test_file_lines = test_file.readlines()
    if test_file_lines == []:
        dependencies_met = False
        print ('\nCCP4 program suite not found. If you have downloaded the\n'
               'CCP4 suite, ensure that you are running RABDAM in a\n'
               'terminal / command prompt that can run CCP4 programs')
    else:
        test_file_line = test_file_lines[0]
        if 'pdbcur:command not found' in test_file_line.replace(' ', ''):
            dependencies_met = False
            print ('\nCCP4 program suite not found. If you have downloaded the\n'
                   'CCP4 suite, ensure that you are running RABDAM in a\n'
                   'terminal / command prompt that can run CCP4 programs')
    test_file.close()
    os.remove('test.txt')

    if dependencies_met is True:
        print '\nSystem meets all RABDAM dependencies'
