#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read()

test_requirements = [
    'pytest', 'coverage', "flake8"
]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='cerebra',
    version='1.1.3',
    description="finds mutants in your scRNA-seq experiment",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Lincoln Harris",
    author_email='ljharris018@gmail.com',
    url='https://github.com/czbiohub/cerebra',
    packages=[
        'cerebra',
    ],
    package_dir={'cerebra':
                     'cerebra'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT",
    zip_safe=False,
    keywords='cerebra',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Operating System :: Unix",
        "Operating System :: MacOS",
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    entry_points={
        'console_scripts': [
            'cerebra = cerebra.commandline:cli'
        ]
    },
    test_suite='tests',
    tests_require=test_requirements
)
