from setuptools import setup, find_packages
import versioneer

requirements = [
    'numpy',
    'pandas',
    'networkx',
    'scikit-learn',
    'markov_clustering',
    'pyteomics',
    'parmap',
    # 'tcrdist3'
    # package requirements go here
]

setup(
    name='clusTCR',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="tcr clustering",
    license="MIT",
    author="Sebastiaan Valkiers & Max Van Houcke",
    author_email='sebastiaan.valkiers@uantwerpen.be',
    url='https://github.com/svalkiers/clustcr',
    packages=find_packages(exclude=('analyses')),
    package_data={
        'clustcr': ['../data/*', '../data/*/*']
    },
    include_package_data=True,
    install_requires=requirements,
    keywords='https://github.com/svalkiers/clustcr',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
