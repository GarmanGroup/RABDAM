
from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='rabdam',
    packages=find_packages(),
    package_data={'rabdam': ['Subroutines/HTML_stylesheet.css',
                             'Subroutines/HTML_stylesheet.js']},
    version='1.6.4',
    description=('RABDAM: A program to identify and quantify specific '
                 'radiation damage within individual protein crystal '
                 'structures'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=('Kathryn Shelley, Thomas Dixon and Jonathan Brooks-Bartlett in '
            'the lab of Professor Elspeth Garman, University of Oxford'),
    author_email='kathryn.l.shelley@gmail.com',
    url='https://github.com/GarmanGroup/RABDAM',
    license='LGPL v3',
    keywords=['radiation damage', 'specific damage', 'atomic Bfactors',
              'atomic displacement parameters', 'BDamage', 'Bnet'],
    install_requires=['numpy>=1.15.0', 'matplotlib>=2.2.0', 'scipy>=1.1.0',
                      'pandas>=0.24.1', 'requests', 'setuptools'],
    classifiers=['Programming Language :: Python'],
    python_requires=('>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, '
                     '!=3.5.*, <4'),
    entry_points={'console_scripts': ['rabdam = rabdam.rabdam:main']}
)
