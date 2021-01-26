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
    'faiss'
    # 'tcrdist3'
    # package requirements go here
]

setup(
    name='https://github.com/svalkiers/clustcr',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="tcr clustering",
    license="MIT",
    author="Max Van Houcke",
    author_email='max@sg.ua',
    url='https://github.com/maxvanhoucke/https://github.com/svalkiers/clustcr',
    packages=find_packages(),
    install_requires=requirements,
    keywords='https://github.com/svalkiers/clustcr',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
