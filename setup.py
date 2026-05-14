from setuptools import setup

setup(
   name='scalable_smr',
   version='1.0',
   description='A simple design of a scalable Small Modular Reactor',
   author='Gasper Srsen',
   author_email='gasper.srsen@ijs.si',
   packages=['scalable_smr'],  #same as name
   install_requires=['numpy', 'scipy', 'matplotlib'], #external packages as dependencies
#    scripts=[
#             'scripts/self.mats["cool"]',
#             'scripts/skype',
#            ]
)
