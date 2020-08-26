''' setup for ocean_model_topog_generator '''
import setuptools

setuptools.setup(
    name="OMtopogen",
    version="0.0.1",
    author="",
    author_email="",
    description=(""),
    license="",
    keywords="",
    url="",
    packages=['OMtopogen'],
    scripts=['OMtopogen/create_topog_refinedSampling.py',
             'OMtopogen/ice9.py',
             'OMtopogen/merge_topog_tiles.py']
)
