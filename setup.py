
from setuptools import setup, find_packages

def readme():
    with open('README.md', 'r') as info:
        return info.read()

setup(
    name='rabdam',
    packages=find_packages(),
    package_data={'rabdam': ['Subroutines/HTML_stylesheet.css',
                             'Subroutines/HTML_stylesheet.js']},
    version='1.0.0',
    description=('RABDAM: A program to identify and quantify specific '
                 'radiation damage within individual protein crystal '
                 'structures'),
    long_description=readme(),
    author=('Kathryn Shelley, Thomas Dixon and Jonathan Brooks-Bartlett in '
            'the lab of Professor Elspeth Garman, University of Oxford'),
    author_email='kathryn.l.shelley@gmail.com',
    url='https://github.com/GarmanGroup/RABDAM',
    download_url='https://github.com/GarmanGroup/RABDAM/archive/1.0.tar.gz',
    license='LGPL v3',
    keywords=['radiation damage', 'specific damage', 'atomic Bfactors',
              'atomic displacement parameters', 'BDamage', 'Bnet'],
    install_requires=['numpy', 'matplotlib', 'seaborn', 'pandas', 'requests',
                        'setuptools'],
    classifiers=['Programming Language :: Python :: 3.6'],
    python_requires='>=3.6',
    entry_points={'console_scripts': ['rabdam = rabdam.rabdam:main']}
)
