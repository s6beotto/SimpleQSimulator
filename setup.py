from setuptools import setup

setup(
    name='simpleqsimulator',
    version='0.1',
    description='A simple quantum simulator',
    author='Benedikt Otto',
    author_email='s6beotto@uni-bonn.de',
    packages=['simpleqsimulator'],
    entry_points={
        'console_scripts': [
            'sqs-shell = simpleqsimulator.shell:start',
        ]},
    install_requires=['IPython', 'numpy', 'scipy'],
)
