from setuptools import setup, find_packages

setup(
    name='HATKit',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[],
    entry_points={
        'console_scripts': [
            'hatkit = hatkit.cli:main',
        ],
    },
    author='Christopher M. Rogers',
    description='HemTools Auxiliary Transcription Factor Kit',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/HATKit',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)
