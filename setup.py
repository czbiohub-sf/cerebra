#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.md') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read()

test_requirements = [
    'pytest', 'coverage', "flake8"
]

setup(
    name='cerebra',
    version='0.1.0',
    description="finds mutants in your scRNA-seq experiment",
    long_description=readme + '\n\n' + history,
    author="Lincoln Harris",
    author_email='lincoln.harris@czbiohub.org',
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
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    entry_points={
        'console_scripts': [
            'cerebra = cerebra.commandline:cli'
        ]
    },
    test_suite='tests',
    tests_require=test_requirements
)
