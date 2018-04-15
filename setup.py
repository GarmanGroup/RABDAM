
from setuptools import setup, find_packages

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except:
    long_description = open('README.md').read()

setup(
    name='rabdam',
    packages=find_packages(),
    package_data={'rabdam': ['Subroutines/HTML_stylesheet.css',
                             'Subroutines/HTML_stylesheet.js']},
    version='1.2.0',
    description=('RABDAM: A program to identify and quantify specific '
                 'radiation damage within individual protein crystal '
                 'structures'),
    long_description=long_description,
    author=('Kathryn Shelley, Thomas Dixon and Jonathan Brooks-Bartlett in '
            'the lab of Professor Elspeth Garman, University of Oxford'),
    author_email='kathryn.l.shelley@gmail.com',
    url='https://github.com/GarmanGroup/RABDAM',
    download_url='https://github.com/GarmanGroup/RABDAM/archive/1.1.tar.gz',
    license='LGPL v3',
    keywords=['radiation damage', 'specific damage', 'atomic Bfactors',
              'atomic displacement parameters', 'BDamage', 'Bnet'],
    install_requires=['numpy', 'matplotlib', 'seaborn', 'pandas', 'requests',
                      'setuptools'],
    classifiers=['Programming Language :: Python'],
    python_requires=('>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, '
                     '!=3.5.*, <4'),
    entry_points={'console_scripts': ['rabdam = rabdam.rabdam:main']}
)
